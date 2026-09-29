#!/bin/zsh
set -u
set -o pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
DATASETS="$PROJECT/.tools/bin/datasets"
CHUNK_DIR="$PROJECT/config/caudoviricetes_refseq_chunks"
PACKAGE_DIR="$PROJECT/data/ncbi_packages/caudoviricetes_refseq_chunks"
LOG_DIR="$PROJECT/metadata/caudoviricetes_download_logs"
STATUS_FILE="$PROJECT/results/caudoviricetes_download_status.tsv"

mkdir -p "$PACKAGE_DIR" "$LOG_DIR" "$PROJECT/results"

if [[ ! -e "$STATUS_FILE" ]]; then
    print "chunk\tstatus\texit_code\trequested\tobserved\tzip_file\tlog_file" > "$STATUS_FILE"
fi

for chunk_file in "$CHUNK_DIR"/chunk_*; do
    chunk="${chunk_file:t}"
    zip_file="$PACKAGE_DIR/${chunk}.zip"
    log_file="$LOG_DIR/${chunk}.log"

    if awk -F '\t' -v chunk="$chunk" \
        'NR > 1 && $1 == chunk && $2 == "complete" {
            found=1
        }
        END { exit !found }' \
        "$STATUS_FILE"; then
        print "Already recorded: $chunk"
        continue
    fi

    requested="$(wc -l < "$chunk_file" | tr -d ' ')"

    if [[ -s "$zip_file" ]] && unzip -tq "$zip_file" >/dev/null 2>&1; then
        observed="$(
            unzip -p "$zip_file" \
              'ncbi_dataset/data/genomic.fna' |
            grep -c '^>'
        )"

        if [[ "$requested" -eq "$observed" ]]; then
            print "$chunk\tcomplete\t0\t$requested\t$observed\t$zip_file\t$log_file" >> "$STATUS_FILE"
            print "Validated existing $chunk: $observed sequences"
            continue
        fi
    fi

    if [[ -e "$zip_file" ]]; then
        timestamp="$(date '+%Y%m%dT%H%M%S')"
        mv "$zip_file" "${zip_file}.incomplete_${timestamp}"
        print "Preserved incomplete archive for $chunk"
    fi

    print "Downloading $chunk ($requested accessions)"

    "$DATASETS" download virus genome accession \
        --inputfile "$chunk_file" \
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

    if [[ "$exit_code" -eq 0 && "$requested" -eq "$observed" ]]; then
        download_status="complete"
    else
        download_status="failed"
    fi

    print "$chunk\t$download_status\t$exit_code\t$requested\t$observed\t$zip_file\t$log_file" >> "$STATUS_FILE"
    print "Finished $chunk: $download_status ($observed/$requested)"
done

print
print "Download summary:"
cut -f2 "$STATUS_FILE" | tail -n +2 | sort | uniq -c
print "Status table: $STATUS_FILE"
