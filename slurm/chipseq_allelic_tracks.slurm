#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J chip_AS_track
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8    # serial_std caps at 16, but 8 is the real need: haplotag is
                             # single-threaded per contig and 8 workers already cover the
                             # makespan, which is bound by chr1 alone
#SBATCH --mem=24G
#SBATCH --time=8:00:00       # serial_std hard cap is 24h; the run is ~2-4h, and a smaller
                             # ask backfills far sooner

# Allele-specific ChIP-seq tracks for IGV.
#
#   merged.clean.sorted.bam
#     -> whatshap haplotag  (B6/CAST SNPs as a pseudo-phased VCF, one chromosome per core)
#     -> whatshap split     -> <MARK>_B6.bam   (H1 = REF = B6   = maternal)
#                           -> <MARK>_CAST.bam (H2 = ALT = CAST = paternal)
#     -> bamCoverage        -> <MARK>_B6.bw    (plotted up)
#                           -> <MARK>_CAST.bw  (plotted down, same scale factor)
#
# Both bigWigs use ONE shared scale factor, 1e6 / (N_B6 + N_CAST). Normalising each
# track to its own library size would divide the allelic imbalance away, which is the
# only thing the track is meant to show. With a shared denominator the two bar heights
# are directly comparable and their sum is the assignable coverage in CPM.
#
# --export=NONE, so arguments are positional:
#   $1 MARK     H3K27ac | H3K4me3            (default H3K27ac -> heart_H3K27ac/)
#   $2 FRAGLEN  read extension in bp         (default 188; single-end ChIP, so bamCoverage
#                                             needs an explicit fragment length. 188 is the
#                                             mean of the MACS2 "alternative fragment length"
#                                             for the two H3K27ac replicates, 183 and 193.
#                                             H3K4me3: check its own *_peaks.xls.)
#   $3 BINSIZE  bigWig bin size in bp        (default 25)
#   $4 MAXV     y-axis limit in the IGV session XML (default 10)
#
#   sbatch slurm/chipseq_allelic_tracks.slurm H3K27ac

set -euo pipefail
date

MARK=${1:-H3K27ac}
FRAGLEN=${2:-188}
BINSIZE=${3:-25}
MAXV=${4:-10}

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

######------ paths ------######
ROOT=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
WORK=$ROOT/heart_${MARK}
BAM=$WORK/merged.clean.sorted.bam                       # replicate-merged, -F 1804 -q 30, blacklist-filtered
REF=$ROOT/bwa-index/GRCm39.primary_assembly.genome.fa
SNPBED=$ROOT/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed
VCF=$ROOT/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.vcf.gz
OUT=$WORK/igv
TMP=${SCRATCH:-$WORK}/whatshap_${MARK}_${SLURM_JOB_ID:-manual}

NPROC=${SLURM_CPUS_PER_TASK:-8}
JOBS=$(( NPROC < 8 ? NPROC : 8 ))                       # haplotag processes in flight; more
                                                        # than 8 buys little with 20 contigs
THRD_PER_JOB=$(( NPROC / JOBS )); THRD_PER_JOB=$(( THRD_PER_JOB > 0 ? THRD_PER_JOB : 1 ))

mkdir -p "$OUT" "$TMP"

for f in "$BAM" "$REF" "$SNPBED"; do
  [[ -f $f ]] || { echo "missing input: $f" >&2; exit 1; }
done

######------ environment ------######
# whatshap, samtools, htslib and deeptools all live in RNAseq, so no env switching.
conda activate RNAseq
for tool in whatshap samtools bgzip tabix bamCoverage; do
  command -v "$tool" >/dev/null || { echo "$tool not found in the RNAseq env" >&2; exit 1; }
done
echo -e "\nwhatshap version: $(whatshap --version)"

######------ SNP bed -> pseudo-phased VCF (built once, then reused) ------######
# SNPfile_C57BL6_NJ_CAST_EiJ.bed is chr / pos / pos+1 / <REF><ALT> / 0 / .
# Column 2 is ALREADY the 1-based position, not a 0-based BED start: checked against
# GRCm39.primary_assembly.genome.fa (chr1 3069198 "AC" -> reference base A at 3069198,
# T at 3069199). Do not "fix" this to $3.
# Column 4 is B6 base then CAST base, matching Allelome.PRO2's A1 = B6 convention, so
# GT 0|1 makes haplotype 1 = B6 and haplotype 2 = CAST.
if [[ ! -f ${VCF}.tbi ]]; then
  echo -e "\n######## Building $VCF ########\n"
  {
    echo '##fileformat=VCFv4.2'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    # single pass, unlike the column-by-column paste in AS_BW_whatshap.sh: rows whose
    # 4th field is not two bases are dropped whole, so chrom/pos can never slip out of
    # register against REF/ALT
    awk 'BEGIN{OFS="\t"} $4 ~ /^[ACGT][ACGT]$/ {
           print $1, $2, ".", substr($4,1,1), substr($4,2,1), "50.0", ".", ".", "GT", "0|1"
         }' "$SNPBED"
  } | bgzip -c > "$VCF"
  tabix -p vcf "$VCF"
fi

[[ -f ${REF}.fai ]] || samtools faidx "$REF"
[[ -f ${BAM}.bai ]] || samtools index -@ "$NPROC" "$BAM"

######------ haplotag, one chromosome per core ------######
# --regions makes haplotag fetch only that contig, so each process reads its own slice
# of the BAM and writes its own read->haplotype list. The tagged BAMs are throwaway:
# only the lists are used downstream.
mapfile -t CHROMS < <(tabix -l "$VCF")
echo -e "\nhaplotagging ${#CHROMS[@]} contigs, $JOBS at a time"

tag_chrom () {
  set -euo pipefail
  chrom="$1"
  whatshap haplotag \
    --reference "$REF" \
    --regions "$chrom" \
    --ignore-read-groups --ignore-linked-read --skip-missing-contigs \
    --output-threads "$THRD_PER_JOB" \
    --output-haplotag-list "$TMP/${chrom}.haplotag.tsv" \
    -o "$TMP/${chrom}.tagged.bam" \
    "$VCF" "$BAM"
  rm -f "$TMP/${chrom}.tagged.bam"
  echo "  done $chrom"
}
export -f tag_chrom
export REF VCF BAM TMP THRD_PER_JOB

printf '%s\n' "${CHROMS[@]}" | xargs -n1 -P "$JOBS" -I{} bash -c 'tag_chrom "$@"' _ {}

# one list for the whole genome: keep the first header, drop the rest
first="$TMP/${CHROMS[0]}.haplotag.tsv"
{ head -n1 "$first"; cat "$TMP"/*.haplotag.tsv | grep -v '^#'; } > "$TMP/haplotag_list.tsv"
echo -e "\ntagged reads: $(grep -vc '^#' "$TMP/haplotag_list.tsv") in list, \
H1 $(awk -F'\t' '$2=="H1"' "$TMP/haplotag_list.tsv" | wc -l), \
H2 $(awk -F'\t' '$2=="H2"' "$TMP/haplotag_list.tsv" | wc -l)"

######------ split into the two alleles ------######
# H1 = REF = B6 = maternal (the active X in this cross), H2 = ALT = CAST = paternal (Xi).
# Reads that never overlap a SNP are unassignable and are dropped here by design.
echo -e "\n######## Splitting ########\n"
whatshap split \
  --output-h1 "$TMP/${MARK}_B6.bam" \
  --output-h2 "$TMP/${MARK}_CAST.bam" \
  "$BAM" "$TMP/haplotag_list.tsv"

for allele in B6 CAST; do
  samtools sort -@ "$NPROC" -m 1G -o "$OUT/${MARK}_${allele}.bam" "$TMP/${MARK}_${allele}.bam"
  samtools index -@ "$NPROC" "$OUT/${MARK}_${allele}.bam"
  rm -f "$TMP/${MARK}_${allele}.bam"
done

N_B6=$(samtools view -c -@ "$NPROC" "$OUT/${MARK}_B6.bam")
N_CAST=$(samtools view -c -@ "$NPROC" "$OUT/${MARK}_CAST.bam")
SF=$(awk -v a="$N_B6" -v b="$N_CAST" 'BEGIN{ printf "%.10f", 1000000/(a+b) }')
echo -e "\nassigned reads: B6 $N_B6, CAST $N_CAST"
awk -v a="$N_B6" -v b="$N_CAST" 'BEGIN{ printf "genome-wide B6 fraction: %.4f (escape on chrX is 1 - this, per peak)\n", a/(a+b) }'
echo "shared scale factor: $SF"

######------ bigWigs ------######
echo -e "\n######## bigWig ########\n"

bamCoverage -b "$OUT/${MARK}_B6.bam"   -o "$OUT/${MARK}_B6.bw" \
  --binSize "$BINSIZE" --extendReads "$FRAGLEN" --scaleFactor "$SF"  -p "$NPROC"
bamCoverage -b "$OUT/${MARK}_CAST.bam" -o "$OUT/${MARK}_CAST.bw" \
  --binSize "$BINSIZE" --extendReads "$FRAGLEN" --scaleFactor "-$SF" -p "$NPROC"

######------ IGV session ------######
# Relative paths: copy the whole igv/ folder to the Mac and open the .xml in IGV.
cat > "$OUT/${MARK}_allelic_session.xml" <<XML
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" locus="chrX:102503978-102526839" version="8">
  <Resources>
    <Resource path="${MARK}_B6.bw"/>
    <Resource path="${MARK}_CAST.bw"/>
  </Resources>
  <Panel name="DataPanel">
    <Track id="${MARK}_B6.bw" name="${MARK} B6 (maternal, Xa)" color="0,90,181"
           renderer="BAR_CHART" windowFunction="mean" autoScale="false" height="70">
      <DataRange minimum="0.0" maximum="${MAXV}.0" baseline="0.0" type="LINEAR"/>
    </Track>
    <Track id="${MARK}_CAST.bw" name="${MARK} CAST (paternal, Xi)" color="220,50,32"
           renderer="BAR_CHART" windowFunction="mean" autoScale="false" height="70">
      <DataRange minimum="-${MAXV}.0" maximum="0.0" baseline="0.0" type="LINEAR"/>
    </Track>
  </Panel>
  <Panel name="FeaturePanel">
    <Track id="Reference sequence" name="Reference sequence"/>
  </Panel>
</Session>
XML

rmdir "$TMP" 2>/dev/null || echo "note: $TMP not empty, left in place"

echo -e "\nWrote:"
ls -lh "$OUT"
date
