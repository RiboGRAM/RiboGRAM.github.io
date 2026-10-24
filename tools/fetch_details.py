"""Mirror the public RiboGRAM per-family interpretation bundles for GitHub Pages."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import time
import urllib.request
from pathlib import Path


DEFAULT_BASE = "https://www.ultrarnalab.com/ribogram_assets/explorer/data/"


def load_manifest(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def fetch_one(url: str, destination: Path, retries: int = 4) -> tuple[str, int]:
    if destination.exists() and destination.stat().st_size:
        return "cached", destination.stat().st_size

    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(retries):
        temporary = destination.with_suffix(destination.suffix + ".part")
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "RiboGRAM-Atlas/1.0"})
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read()
            with gzip.GzipFile(fileobj=__import__("io").BytesIO(payload)) as handle:
                json.load(handle)
            temporary.write_bytes(payload)
            temporary.replace(destination)
            return "downloaded", len(payload)
        except Exception as error:  # pragma: no cover - network recovery path
            last_error = error
            temporary.unlink(missing_ok=True)
            time.sleep(0.4 * (attempt + 1))
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("data/manifest.json.gz"))
    parser.add_argument("--output", type=Path, default=Path("data/details"))
    parser.add_argument("--base-url", default=DEFAULT_BASE)
    parser.add_argument("--workers", type=int, default=20)
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    jobs: list[tuple[str, Path]] = []
    for family in manifest["families"]:
        relative = family["detail_url"]
        jobs.append((args.base_url + relative, args.output / Path(relative).name))

    downloaded = cached = bytes_total = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(fetch_one, url, path): path for url, path in jobs}
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            state, size = future.result()
            downloaded += state == "downloaded"
            cached += state == "cached"
            bytes_total += size
            if index % 200 == 0 or index == len(jobs):
                print(
                    f"{index:,}/{len(jobs):,} bundles · {bytes_total / 1024 / 1024:.1f} MiB "
                    f"· {downloaded:,} downloaded · {cached:,} cached",
                    flush=True,
                )


if __name__ == "__main__":
    main()
