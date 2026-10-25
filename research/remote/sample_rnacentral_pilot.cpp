#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

uint64_t fnv1a(const std::string& text) {
    uint64_t value = 1469598103934665603ULL;
    for (unsigned char c : text) {
        value ^= static_cast<uint64_t>(c);
        value *= 1099511628211ULL;
    }
    return value;
}

uint64_t splitmix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

std::string record_id(const std::string& header) {
    const std::size_t begin = !header.empty() && header[0] == '>' ? 1 : 0;
    const std::size_t end = header.find_first_of(" \t", begin);
    return header.substr(begin, end == std::string::npos ? end : end - begin);
}

struct Options {
    std::string output_dir;
    std::size_t shards = 16;
    std::size_t min_length = 20;
    std::size_t max_length = 1024;
    double sample_probability = 0.05419;
    double max_ambiguous_fraction = 0.05;
    uint64_t seed = 20261025ULL;
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        auto require_value = [&]() -> std::string {
            if (++index >= argc) {
                throw std::runtime_error("Missing value for " + argument);
            }
            return argv[index];
        };
        if (argument == "--output-dir") {
            options.output_dir = require_value();
        } else if (argument == "--shards") {
            options.shards = std::stoull(require_value());
        } else if (argument == "--min-length") {
            options.min_length = std::stoull(require_value());
        } else if (argument == "--max-length") {
            options.max_length = std::stoull(require_value());
        } else if (argument == "--sample-probability") {
            options.sample_probability = std::stod(require_value());
        } else if (argument == "--max-ambiguous-fraction") {
            options.max_ambiguous_fraction = std::stod(require_value());
        } else if (argument == "--seed") {
            options.seed = std::stoull(require_value());
        } else {
            throw std::runtime_error("Unknown argument: " + argument);
        }
    }
    if (options.output_dir.empty()) {
        throw std::runtime_error("--output-dir is required");
    }
    if (options.shards == 0 || options.min_length > options.max_length) {
        throw std::runtime_error("Invalid shard count or length bounds");
    }
    if (!(options.sample_probability > 0.0 && options.sample_probability <= 1.0)) {
        throw std::runtime_error("sample probability must be in (0, 1]");
    }
    if (!(options.max_ambiguous_fraction >= 0.0 && options.max_ambiguous_fraction <= 1.0)) {
        throw std::runtime_error("max ambiguous fraction must be in [0, 1]");
    }
    return options;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        std::vector<std::ofstream> outputs(options.shards);
        std::vector<uint64_t> shard_records(options.shards, 0);
        std::vector<uint64_t> shard_nt(options.shards, 0);
        for (std::size_t shard = 0; shard < options.shards; ++shard) {
            std::ostringstream path;
            path << options.output_dir << "/pilot_shard_" << std::setw(2)
                 << std::setfill('0') << shard << ".fa";
            outputs[shard].open(path.str(), std::ios::binary);
            if (!outputs[shard]) {
                throw std::runtime_error("Cannot open " + path.str());
            }
        }

        const long double scaled = static_cast<long double>(options.sample_probability) *
                                   static_cast<long double>(std::numeric_limits<uint64_t>::max());
        const uint64_t threshold = static_cast<uint64_t>(scaled);
        uint64_t total_records = 0;
        uint64_t eligible_records = 0;
        uint64_t selected_records = 0;
        uint64_t selected_nt = 0;

        std::string header;
        std::string sequence;
        std::string line;
        auto emit = [&]() {
            if (header.empty()) {
                return;
            }
            ++total_records;
            const std::size_t length = sequence.size();
            if (length < options.min_length || length > options.max_length) {
                return;
            }
            std::size_t canonical = 0;
            for (unsigned char raw : sequence) {
                const char base = static_cast<char>(std::toupper(raw));
                canonical += (base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'U');
            }
            const double ambiguous_fraction = length == 0
                ? 1.0
                : 1.0 - static_cast<double>(canonical) / static_cast<double>(length);
            if (ambiguous_fraction > options.max_ambiguous_fraction) {
                return;
            }
            ++eligible_records;
            const std::string id = record_id(header);
            const uint64_t hash = splitmix64(fnv1a(id) ^ options.seed);
            if (hash > threshold) {
                return;
            }
            const std::size_t shard = splitmix64(hash) % options.shards;
            outputs[shard] << header << '\n';
            for (std::size_t offset = 0; offset < sequence.size(); offset += 80) {
                outputs[shard] << sequence.substr(offset, 80) << '\n';
            }
            ++selected_records;
            selected_nt += length;
            ++shard_records[shard];
            shard_nt[shard] += length;
            if (total_records % 5000000ULL == 0) {
                std::cerr << "records=" << total_records
                          << " eligible=" << eligible_records
                          << " selected=" << selected_records << '\n';
            }
        };

        while (std::getline(std::cin, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if (!line.empty() && line[0] == '>') {
                emit();
                header = line;
                sequence.clear();
            } else {
                for (unsigned char raw : line) {
                    if (!std::isspace(raw)) {
                        sequence.push_back(static_cast<char>(std::toupper(raw)));
                    }
                }
            }
        }
        emit();

        const std::string manifest_path = options.output_dir + "/pilot_manifest.json";
        std::ofstream manifest(manifest_path);
        manifest << "{\n"
                 << "  \"schema_version\": 1,\n"
                 << "  \"seed\": " << options.seed << ",\n"
                 << "  \"sample_probability\": " << std::setprecision(12)
                 << options.sample_probability << ",\n"
                 << "  \"min_length\": " << options.min_length << ",\n"
                 << "  \"max_length\": " << options.max_length << ",\n"
                 << "  \"max_ambiguous_fraction\": " << options.max_ambiguous_fraction << ",\n"
                 << "  \"input_records\": " << total_records << ",\n"
                 << "  \"eligible_records\": " << eligible_records << ",\n"
                 << "  \"selected_records\": " << selected_records << ",\n"
                 << "  \"selected_nt\": " << selected_nt << ",\n"
                 << "  \"shards\": [\n";
        for (std::size_t shard = 0; shard < options.shards; ++shard) {
            manifest << "    {\"shard\": " << shard
                     << ", \"records\": " << shard_records[shard]
                     << ", \"nt\": " << shard_nt[shard] << "}"
                     << (shard + 1 == options.shards ? "\n" : ",\n");
        }
        manifest << "  ]\n}\n";
        std::cerr << "complete records=" << total_records
                  << " eligible=" << eligible_records
                  << " selected=" << selected_records
                  << " selected_nt=" << selected_nt << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "sample_rnacentral_pilot: " << error.what() << '\n';
        return 2;
    }
}
