#!/usr/bin/env bash
set -euo pipefail

# ---- user settings ----
SLURM_SCRIPT="sc_Allelome.PRO2.slurm"
FILELIST="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/9w/sc_filelist.txt"
CHUNK_SIZE=200
POLL_SECONDS=120

# If first 200 already finished successfully, start from 200.
# Use 0 to start from the beginning.
START_OFFSET=200

# Set to 1 if you want the wrapper to stop when any task in a chunk fails.
STOP_ON_FAILURE=1

# ---- sanity checks ----
if [[ ! -f "$SLURM_SCRIPT" ]]; then
    echo "ERROR: Cannot find SLURM script: $SLURM_SCRIPT" >&2
    exit 1
fi

if [[ ! -f "$FILELIST" ]]; then
    echo "ERROR: Cannot find file list: $FILELIST" >&2
    exit 1
fi

TOTAL=$(wc -l < "$FILELIST")
echo "Total samples in file list: $TOTAL"
echo "Chunk size: $CHUNK_SIZE"
echo "Starting offset: $START_OFFSET"
echo

submit_chunk () {
    local offset="$1"
    echo "Submitting chunk with OFFSET=$offset (lines $((offset+1))-$((offset+CHUNK_SIZE)))"

    # sbatch output is usually: "Submitted batch job 123456"
    # or "Submitted batch job 123456 on cluster serial"
    local sbatch_out
    sbatch_out=$(sbatch --export=OFFSET="$offset" "$SLURM_SCRIPT")
    echo "$sbatch_out"

    local jobid
    jobid=$(awk '/Submitted batch job/ {print $4}' <<< "$sbatch_out")

    if [[ -z "${jobid:-}" ]]; then
        echo "ERROR: Could not parse job ID from sbatch output" >&2
        exit 1
    fi

    echo "Chunk job ID: $jobid"
    echo "$jobid"
}

wait_for_job () {
    local jobid="$1"

    echo "Waiting for job $jobid to leave the queue..."
    while squeue --clusters=serial --partition=serial_std -h -j "$jobid" | grep -q .; do
        sleep "$POLL_SECONDS"
    done

    echo "Job $jobid has left the queue."
}

summarise_job () {
    local jobid="$1"

    echo
    echo "Accounting summary for job $jobid:"
    sacct -M serial -r serial_std -j "$jobid" --array \
        --format=JobID,JobName%24,State,ExitCode,Elapsed,MaxRSS
    echo
}

chunk_has_failures () {
    local jobid="$1"

    # Look for clearly bad end states in any task
    local bad_count
    bad_count=$(
        sacct -M serial -r serial_std -n -j "$jobid" --array --format=State \
        | awk '
            {
              gsub(/[[:space:]]+/, "", $0)
              if ($0 ~ /FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED/) c++
            }
            END {print c+0}
        '
    )

    [[ "$bad_count" -gt 0 ]]
}

for (( offset=START_OFFSET; offset<TOTAL; offset+=CHUNK_SIZE )); do
    jobid=$(submit_chunk "$offset" | tail -n1)
    wait_for_job "$jobid"
    summarise_job "$jobid"

    if [[ "$STOP_ON_FAILURE" -eq 1 ]]; then
        if chunk_has_failures "$jobid"; then
            echo "Detected failed/cancelled/timed-out/OOM tasks in chunk job $jobid."
            echo "Stopping here so you can inspect logs before submitting the next chunk."
            exit 1
        fi
    fi

    echo "Chunk OFFSET=$offset completed successfully."
    echo
done

echo "All chunks submitted and completed."
