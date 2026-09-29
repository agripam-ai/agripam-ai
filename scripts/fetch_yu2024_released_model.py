#!/usr/bin/env python3
"""Fetch only the three required files from the versioned Yu et al. Zenodo ZIP.

HTTP byte ranges avoid downloading the complete 200 MB archive.  ZIP offsets,
uncompressed sizes and SHA-256 values pin the exact archived release.
"""
from __future__ import annotations

import argparse
import hashlib
import struct
import urllib.request
import zlib
from pathlib import Path


URL = "https://zenodo.org/api/records/10262867/files/BarquistLab/CRISPRi_guide_efficiency_bacteria-publish.zip/content"
ENTRIES = {
    "CRISPRi_headers.sav": (153112841, 417, 2077, "0e4cf56a486705ad7a140ab20f6a476edbdfb067358e5bbf6f0a7b157415abbe"),
    "CRISPRi_model.sav": (153113396, 20969545, 57029462, "04bf43713d1ed8f89b016d1de16e9be11b259e2af8b7beebf6e6e694a8ad453d"),
    "pur_gRNAs.tsv": (200431270, 105855, 316646, "31f7a6ca3039429ce092aea61c111b29be0553711ddc1fcf0c28d220028d427a"),
}


def request_range(start: int, end: int) -> bytes:
    request = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}"})
    with urllib.request.urlopen(request) as response:
        data = response.read()
    expected = end - start + 1
    if len(data) != expected:
        raise ValueError(f"Server returned {len(data)} bytes; expected {expected} for {start}-{end}.")
    return data


def fetch_entry(offset: int, compressed_size: int) -> bytes:
    header = request_range(offset, offset + 511)
    signature, _, _, method, _, _, _, archived_size, _, name_len, extra_len = struct.unpack(
        "<4s5H3L2H", header[:30]
    )
    if signature != b"PK\x03\x04" or method != 8 or archived_size != compressed_size:
        raise ValueError("Pinned ZIP entry metadata does not match the Zenodo archive.")
    data_start = offset + 30 + name_len + extra_len
    return zlib.decompress(request_range(data_start, data_start + compressed_size - 1), -15)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name, (offset, compressed, size, digest) in ENTRIES.items():
        content = fetch_entry(offset, compressed)
        if len(content) != size or hashlib.sha256(content).hexdigest() != digest:
            raise ValueError(f"Integrity check failed for {name}.")
        (args.output / name).write_bytes(content)
        print(f"verified {name}: {size} bytes, sha256={digest}")


if __name__ == "__main__":
    main()
