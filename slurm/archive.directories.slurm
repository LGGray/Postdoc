#!/bin/bash
#SBATCH -J archive_array
#SBATCH -o archive_%A_%a.out
#SBATCH -t 06:00:00
#SBATCH --mem=4G
#SBATCH -c 8
#SBATCH -p serial
#SBATCH --array=1-3%2          # edit range to match line count; %2 = max 2 concurrent (be gentle on the filesystem)
set -euo pipefail

LIST=/path/to/dirs_to_archive.txt
STORE=/path/to/archive_store
mkdir -p "$STORE"

SRC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
TAG=$(echo "$SRC" | sed 's#^/##; s#/#__#g')   # sanitize full path -> unique flat name
ARCHIVE="$STORE/${TAG}.tar"

find "$SRC" -type f | wc -l  > "${ARCHIVE}.manifest_count"
du -sb "$SRC" | cut -f1      > "${ARCHIVE}.manifest_size"
echo "$SRC" > "${ARCHIVE}.src"

tar --create --verify --file="$ARCHIVE" "$SRC"

pigz -p 8 "$ARCHIVE"
gzip -t "${ARCHIVE}.gz"

echo "OK: $SRC -> ${ARCHIVE}.gz (verified)"
echo "To authorize deletion: touch ${ARCHIVE}.gz.CONFIRM_DELETE"
