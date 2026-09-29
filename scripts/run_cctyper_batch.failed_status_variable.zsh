#!/bin/zsh
set -u
set -o pipefail

PROJECT="$(cd "$(dirname "$0")/.." && pwd)"
MAMBA="$PROJECT/.tools/bin/micromamba"
ENV_PREFIX="$PROJECT/.envs/cctyper-local"
DB="$PROJECT/data/databases/cctyper_1.8.0"
INPUT_DIR="$PROJECT/data/genomes_raw"
OUTPUT_ROOT="$PROJECT/crispr_detection/cctyper_1.8.0"
LOG_ROOT="$PROJECT/crispr_detection/logs"
STATUS="$PROJECT/results/cctyper_batch_status.tsv"

mkdir -p "$OUTPUT_ROOT" "$LOG_ROOT" "$PROJECT/results"

if [[ -e "$STATUS" ]]; then
    print -u2 "ERROR: Status file already exists: $STATUS"
    print -u2 "Move or rename it before starting a new batch."
    exit 1
fi

print "accession\tstatus\texit_code\tstarted\tfinished\toutput_directory\tlog_file" > "$STATUS"

for fasta in "$INPUT_DIR"/GCF_*.fna; do
    accession="${${fasta:t}%.fna}"
    output="$OUTPUT_ROOT/$accession"
    log="$LOG_ROOT/${accession}.log"
    started="$(date '+%Y-%m-%dT%H:%M:%S%z')"

    print "Starting $accession"

    if [[ -e "$output" ]]; then
        finished="$(date '+%Y-%m-%dT%H:%M:%S%z')"
        print "$accession\tnot_run_existing\tNA\t$started\t$finished\t$output\t$log" >> "$STATUS"
        print "Skipped $accession because its output directory already exists."
        continue
    fi

    CCTYPER_DB="$DB" "$MAMBA" run \
        --prefix "$ENV_PREFIX" \
        cctyper \
        "$fasta" \
        "$output" \
        --threads 4 \
        --prodigal single \
        --seed 42 \
        --keep_tmp \
        --no_plot \
        --simplelog \
        > "$log" 2>&1

    exit_code=$?
    finished="$(date '+%Y-%m-%dT%H:%M:%S%z')"

    if [[ "$exit_code" -eq 0 && -s "$output/crisprs_all.tab" && -s "$output/cas_operons.tab" ]]; then
        status="completed"
    else
        status="failed"
    fi

    print "$accession\t$status\t$exit_code\t$started\t$finished\t$output\t$log" >> "$STATUS"
    print "Finished $accession: $status"
done

print
print "Batch summary:"
cut -f2 "$STATUS" | tail -n +2 | sort | uniq -c
print "Status table: $STATUS"
