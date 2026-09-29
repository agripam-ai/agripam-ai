#!/bin/zsh
set -u
set -o pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
DATASETS="$PROJECT/.tools/bin/datasets"
PART_DIR="$PROJECT/config/caudoviricetes_refseq_chunks/chunk_04_parts"
PACKAGE_DIR="$PROJECT/data/ncbi_packages/caudoviricetes_refseq_chunks/chunk_04_parts"
LOG_DIR="$PROJECT/metadata/caudoviricetes_download_logs/chunk_04_parts"

mkdir -p "$PACKAGE_DIR" "$LOG_DIR"

for part_file in "$PART_DIR"/part_*; do
    part="${part_file:t}"
    zip_file="$PACKAGE_DIR/${part}.zip"
    log_file="$LOG_DIR/${part}.log"
    requested="$(wc -l < "$part_file" | tr -d ' ')"

    if [[ -s "$zip_file" ]] && unzip -tq "$zip_file" >/dev/null 2>&1; then
        observed="$(
            unzip -p "$zip_file" \
              'ncbi_dataset/data/genomic.fna' |
            grep -c '^>'
        )"

        if [[ "$observed" -eq "$requested" ]]; then
            print "Already complete: $part ($observed/$requested)"
            continue
        fi
    fi

    if [[ -e "$zip_file" ]]; then
        timestamp="$(date '+%Y%m%dT%H%M%S')"
        mv "$zip_file" "${zip_file}.incomplete_${timestamp}"
    fi

    print "Downloading $part ($requested accessions)"

    "$DATASETS" download virus genome accession \
        --inputfile "$part_file" \
        --include genome \
        --filename "$zip_file" \
        > "$log_file" 2>&1

    exit_code=$?

    if [[ "$exit_code" -eq 0 ]] && unzip -tq "$zip_file" >/dev/null 2>&1; then
        observed="$(
            unzip -p "$zip_file" \
              'ncbi_dataset/data/genomic.fna' |
            grep -c '^>'
        )"
    else
        observed=0
    fi

    if [[ "$exit_code" -eq 0 && "$observed" -eq "$requested" ]]; then
        print "Finished $part: complete ($observed/$requested)"
    else
        print "Finished $part: failed ($observed/$requested)"
    fi
done
