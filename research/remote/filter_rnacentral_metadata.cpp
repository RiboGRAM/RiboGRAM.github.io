#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Options {
    std::string ids;
    std::string rfam;
    std::string id_mapping;
    std::string taxonomy_nodes;
    std::string output_dir;
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) {
                throw std::runtime_error("Missing value for " + key);
            }
            return argv[i];
        };
        if (key == "--ids") options.ids = value();
        else if (key == "--rfam") options.rfam = value();
        else if (key == "--id-mapping") options.id_mapping = value();
        else if (key == "--taxonomy-nodes") options.taxonomy_nodes = value();
        else if (key == "--output-dir") options.output_dir = value();
        else throw std::runtime_error("Unknown option: " + key);
    }
    if (options.ids.empty() || options.rfam.empty() || options.id_mapping.empty() ||
        options.taxonomy_nodes.empty() || options.output_dir.empty()) {
        throw std::runtime_error(
            "Required: --ids --rfam --id-mapping --taxonomy-nodes --output-dir"
        );
    }
    return options;
}

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.emplace_back(line.substr(begin, end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

bool gz_getline(gzFile handle, std::string& line) {
    line.clear();
    constexpr int chunk_size = 1 << 16;
    char chunk[chunk_size];
    while (true) {
        char* result = gzgets(handle, chunk, chunk_size);
        if (result == nullptr) return !line.empty();
        line += chunk;
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
            if (!line.empty() && line.back() == '\r') line.pop_back();
            return true;
        }
        if (gzeof(handle)) return true;
    }
}

std::unordered_map<std::string, std::uint32_t> load_ids(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open IDs: " + path);
    std::unordered_map<std::string, std::uint32_t> result;
    result.reserve(2400000);
    std::string line;
    while (std::getline(input, line)) {
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        const std::uint64_t index = std::stoull(line.substr(0, tab));
        if (index > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("Global index exceeds uint32");
        }
        result.emplace(line.substr(tab + 1), static_cast<std::uint32_t>(index));
    }
    std::cerr << "loaded_ids=" << result.size() << '\n';
    return result;
}

class Taxonomy {
  public:
    explicit Taxonomy(const std::string& nodes_path) {
        std::ifstream input(nodes_path);
        if (!input) throw std::runtime_error("Cannot open taxonomy nodes: " + nodes_path);
        std::string line;
        while (std::getline(input, line)) {
            const std::size_t first = line.find('\t');
            if (first == std::string::npos) continue;
            const std::size_t second_begin = line.find_first_of("0123456789", first);
            if (second_begin == std::string::npos) continue;
            const std::size_t second_end = line.find('\t', second_begin);
            const std::uint32_t child = static_cast<std::uint32_t>(std::stoul(line.substr(0, first)));
            const std::uint32_t parent = static_cast<std::uint32_t>(
                std::stoul(line.substr(second_begin, second_end - second_begin))
            );
            const std::size_t need = static_cast<std::size_t>(std::max(child, parent)) + 1;
            if (parents_.size() < need) parents_.resize(need, 0);
            parents_[child] = parent;
        }
        domains_.resize(parents_.size(), unknown);
        set_domain(2, bacteria);
        set_domain(2157, archaea);
        set_domain(2759, eukaryota);
        set_domain(10239, viruses);
        std::cerr << "taxonomy_nodes_capacity=" << parents_.size() << '\n';
    }

    std::uint8_t domain(std::uint32_t taxid) {
        if (taxid >= parents_.size() || taxid == 0) return unknown;
        if (domains_[taxid] != unknown) return domains_[taxid];
        std::vector<std::uint32_t> path;
        path.reserve(64);
        std::uint32_t cursor = taxid;
        std::uint8_t answer = other;
        for (int depth = 0; depth < 256; ++depth) {
            if (cursor >= parents_.size() || cursor == 0) break;
            if (domains_[cursor] != unknown) {
                answer = domains_[cursor];
                break;
            }
            path.push_back(cursor);
            const std::uint32_t parent = parents_[cursor];
            if (parent == cursor) break;
            cursor = parent;
        }
        for (std::uint32_t node : path) domains_[node] = answer;
        return answer;
    }

    static std::uint8_t mask_for(std::uint8_t domain_code) {
        if (domain_code < bacteria || domain_code > viruses) return 0;
        return static_cast<std::uint8_t>(1u << (domain_code - 1));
    }

  private:
    static constexpr std::uint8_t unknown = 0;
    static constexpr std::uint8_t bacteria = 1;
    static constexpr std::uint8_t archaea = 2;
    static constexpr std::uint8_t eukaryota = 3;
    static constexpr std::uint8_t viruses = 4;
    static constexpr std::uint8_t other = 5;
    std::vector<std::uint32_t> parents_;
    std::vector<std::uint8_t> domains_;

    void set_domain(std::uint32_t taxid, std::uint8_t value) {
        if (taxid < domains_.size()) domains_[taxid] = value;
    }
};

struct MappingAggregate {
    bool selected = false;
    std::uint32_t index = 0;
    std::uint64_t rows = 0;
    std::uint8_t taxon_mask = 0;
    std::unordered_set<std::uint32_t> taxids;
    std::unordered_set<std::string> sources;
    std::unordered_map<std::string, std::uint32_t> type_counts;

    void reset(bool is_selected, std::uint32_t global_index) {
        selected = is_selected;
        index = global_index;
        rows = 0;
        taxon_mask = 0;
        taxids.clear();
        sources.clear();
        type_counts.clear();
    }
};

void write_mapping_row(
    std::ofstream& output,
    const std::string& urs,
    const MappingAggregate& aggregate
) {
    if (!aggregate.selected) return;
    std::vector<std::string> sources(aggregate.sources.begin(), aggregate.sources.end());
    std::sort(sources.begin(), sources.end());
    std::string primary_type;
    std::uint32_t primary_count = 0;
    for (const auto& item : aggregate.type_counts) {
        if (item.second > primary_count ||
            (item.second == primary_count && item.first < primary_type)) {
            primary_type = item.first;
            primary_count = item.second;
        }
    }
    output << aggregate.index << '\t' << urs << '\t' << aggregate.rows << '\t'
           << aggregate.taxids.size() << '\t' << static_cast<unsigned>(aggregate.taxon_mask)
           << '\t' << primary_type << '\t';
    for (std::size_t i = 0; i < sources.size(); ++i) {
        if (i) output << ';';
        output << sources[i];
    }
    output << '\n';
}

void filter_rfam(
    const Options& options,
    const std::unordered_map<std::string, std::uint32_t>& ids
) {
    gzFile input = gzopen(options.rfam.c_str(), "rb");
    if (!input) throw std::runtime_error("Cannot open Rfam annotations: " + options.rfam);
    std::ofstream output(fs::path(options.output_dir) / "rfam_filtered.tsv");
    if (!output) throw std::runtime_error("Cannot create filtered Rfam output");
    std::string line;
    std::uint64_t rows = 0;
    std::uint64_t kept = 0;
    while (gz_getline(input, line)) {
        ++rows;
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        const auto found = ids.find(line.substr(0, tab));
        if (found != ids.end()) {
            output << found->second << '\t' << line << '\n';
            ++kept;
        }
        if (rows % 10000000ULL == 0) {
            std::cerr << "rfam_rows=" << rows << " kept=" << kept << '\n';
        }
    }
    gzclose(input);
    std::cerr << "rfam_complete rows=" << rows << " kept=" << kept << '\n';
}

void aggregate_mapping(
    const Options& options,
    const std::unordered_map<std::string, std::uint32_t>& ids,
    Taxonomy& taxonomy
) {
    gzFile input = gzopen(options.id_mapping.c_str(), "rb");
    if (!input) throw std::runtime_error("Cannot open ID mapping: " + options.id_mapping);
    std::ofstream output(fs::path(options.output_dir) / "id_mapping_summary.tsv");
    if (!output) throw std::runtime_error("Cannot create mapping summary output");
    output << "global_index\turs\tmapping_rows\tunique_taxids\ttaxon_mask\tprimary_rna_type\tsources\n";

    std::string current_urs;
    MappingAggregate aggregate;
    std::string line;
    std::uint64_t rows = 0;
    std::uint64_t selected_rows = 0;
    std::uint64_t selected_ids = 0;
    while (gz_getline(input, line)) {
        ++rows;
        const std::size_t first_tab = line.find('\t');
        if (first_tab == std::string::npos) continue;
        const std::string urs = line.substr(0, first_tab);
        if (urs != current_urs) {
            if (!current_urs.empty()) {
                if (aggregate.selected) ++selected_ids;
                write_mapping_row(output, current_urs, aggregate);
            }
            current_urs = urs;
            const auto found = ids.find(urs);
            aggregate.reset(found != ids.end(), found == ids.end() ? 0 : found->second);
        }
        if (!aggregate.selected) continue;

        ++selected_rows;
        ++aggregate.rows;
        const std::vector<std::string> fields = split_tabs(line);
        if (fields.size() > 1 && !fields[1].empty()) aggregate.sources.insert(fields[1]);
        if (fields.size() > 4 && !fields[4].empty()) ++aggregate.type_counts[fields[4]];
        if (fields.size() > 3 && !fields[3].empty()) {
            try {
                const std::uint32_t taxid = static_cast<std::uint32_t>(std::stoul(fields[3]));
                aggregate.taxids.insert(taxid);
                aggregate.taxon_mask |= Taxonomy::mask_for(taxonomy.domain(taxid));
            } catch (const std::exception&) {
            }
        }
        if (rows % 25000000ULL == 0) {
            std::cerr << "mapping_rows=" << rows << " selected_rows=" << selected_rows
                      << " selected_ids=" << selected_ids << '\n';
        }
    }
    if (!current_urs.empty()) {
        if (aggregate.selected) ++selected_ids;
        write_mapping_row(output, current_urs, aggregate);
    }
    gzclose(input);
    std::cerr << "mapping_complete rows=" << rows << " selected_rows=" << selected_rows
              << " selected_ids=" << selected_ids << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        fs::create_directories(options.output_dir);
        const auto ids = load_ids(options.ids);
        Taxonomy taxonomy(options.taxonomy_nodes);
        filter_rfam(options, ids);
        aggregate_mapping(options, ids, taxonomy);
        std::ofstream success(fs::path(options.output_dir) / "METADATA_SUCCESS");
        success << "ok\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "filter_rnacentral_metadata: " << error.what() << '\n';
        return 2;
    }
}
