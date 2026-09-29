#!/usr/bin/env python3
"""Fetch the two pinned crisprHAL release tables required by the benchmark."""
from __future__ import annotations

import argparse
import hashlib
import urllib.request
from pathlib import Path


FILES = {
    "TevSpCas9_training_data.csv": (
        "https://raw.githubusercontent.com/tbrowne5/crisprHAL/main/data/TevSpCas9_training_data.csv",
        "383c404ed7305fb975bc0d1582eb3eb2ef0dae4fb88cc7b9a20dc019ff8f0aaa",
    ),
    "TevSpCas9_testing_data.csv": (
        "https://raw.githubusercontent.com/tbrowne5/crisprHAL/main/data/TevSpCas9_testing_data.csv",
        "971321d6174fc3981e4a0b48453416791a9a602b92be3eb9a3c96d033750d2ce",
    ),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("data/external_validation/raw/CRISPRHAL_TEV_RELEASED_HOLDOUT"))
    args = parser.parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, (url, expected) in FILES.items():
        payload = urllib.request.urlopen(url, timeout=120).read()
        observed = hashlib.sha256(payload).hexdigest()
        if observed != expected:
            raise RuntimeError(f"SHA-256 mismatch for {name}: expected {expected}, observed {observed}")
        (args.output_dir / name).write_bytes(payload)
        print(f"{name}\t{len(payload)} bytes\t{observed}")


if __name__ == "__main__":
    main()
