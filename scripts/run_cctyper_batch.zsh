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
STATUS_FILE="$PROJECT/results/cctyper_batch_status.tsv"

mkdir -p "$OUTPUT_ROOT" "$LOG_ROOT" "$PROJECT/results"

if [[ ! -e "$STATUS_FILE" ]]; then
    print "accession\tstatus\texit_code\tstarted\tfinished\toutput_directory\tlog_file" > "$STATUS_FILE"
fi

for fasta in "$INPUT_DIR"/GCF_*.fna; do
    accession="${${fasta:t}%.fna}"
    output="$OUTPUT_ROOT/$accession"
    log="$LOG_ROOT/${accession}.log"
    started="$(date '+%Y-%m-%dT%H:%M:%S%z')"

    if awk -F '\t' -v accession="$accession" \
        'NR > 1 && $1 == accession { found=1 } END { exit !found }' \
        "$STATUS_FILE"; then
        print "Already recorded: $accession"
        continue
    fi

    print "Starting $accession"

    if [[ -e "$output" ]]; then
        finished="$(date '+%Y-%m-%dT%H:%M:%S%z')"

        if [[ -s "$output/crisprs_all.tab" && -s "$output/cas_operons.tab" ]]; then
            run_status="completed_existing"
            exit_code=0
        else
            run_status="incomplete_existing"
            exit_code="NA"
        fi

        print "$accession\t$run_status\t$exit_code\t$started\t$finished\t$output\t$log" >> "$STATUS_FILE"
        print "Finished $accession: $run_status"
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

    if [[ "$exit_code" -eq 0 ]]; then
        run_status="completed"
    else
        run_status="failed"
    fi

    print "$accession\t$run_status\t$exit_code\t$started\t$finished\t$output\t$log" >> "$STATUS_FILE"
    print "Finished $accession: $run_status"
done

print
print "Batch summary:"
cut -f2 "$STATUS_FILE" | tail -n +2 | sort | uniq -c
print "Status table: $STATUS_FILE"
