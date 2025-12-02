
### H3K4me3 CHIP-seq ###

# Create bwa index
# bwa-mem2 index GRCm39.primary_assembly.genome.fa

# Merge replicates
# cat ENCFF771WXQ.fastq.gz ENCFF684VBM.fastq.gz ENCFF375HQY.fastq.gz ENCFF908NUT.fastq.gz ENCFF648ICY.fastq.gz ENCFF606CJL.fastq.gz ENCFF619KGH.fastq.gz ENCFF410YWZ.fastq.gz > rep_1.fastq.gz
# cat ENCFF669RHY.fastq.gz ENCFF480DXP.fastq.gz ENCFF636DAW.fastq.gz ENCFF624URI.fastq.gz ENCFF800WAS.fastq.gz ENCFF550MNG.fastq.gz ENCFF097ZHF.fastq.gz ENCFF437ENE.fastq.gz > rep_2.fastq.gz

# Map reads
bwa-mem2 mem -t 8 ../bwa-index/GRCm39.primary_assembly.genome.fa rep_1.fastq.gz | samtools view -bSh > rep_1.bam
# Sort reads
samtools sort -o rep_1.sorted.bam -@ 8 rep_1.bam
# Filter reads: remove unmapped, secondary, QC-fail, duplicates; MAPQ>=30
samtools view -h -b -F 1804 -q 30 rep_1.sorted.bam > rep_1.filtered.bam
# Remove blacklisted regions
bedtools intersect -v -abam rep_1.filtered.bam -b blacklist_GRCm39.bed > rep_1.clean.bam

macs2 callpeak \
  -t rep_1.clean.bam \
  -f BAM \
  -g mm \
  --broad \
  -n H3K4me3_rep1 \
  --outdir macs2

# Add score column to MACS2 summits bed file for Allelome.PRO
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,"."}' \
  macs2/H3K4me3_rep1_summits.bed \
  > macs2/H3K4me3_rep1_summits.bed6


# Map reads
bwa-mem2 mem -t 8 ../bwa-index/GRCm39.primary_assembly.genome.fa rep_2.fastq.gz | samtools view -bSh > rep_2.bam
# Sort reads
samtools sort -o rep_2.sorted.bam -@ 8 rep_2.bam
# Filter reads: remove unmapped, secondary, QC-fail, duplicates; MAPQ>=30
samtools view -h -b -F 1804 -q 30 rep_2.sorted.bam > rep_2.filtered.bam
# Remove blacklisted regions
bedtools intersect -v -abam rep_2.filtered.bam -b blacklist_GRCm39.bed > rep_2.clean.bam

macs2 callpeak \
  -t rep_2.clean.bam \
  -f BAM \
  -g mm \
  --broad \
  -n H3K4me3_rep2 \
  --outdir macs2


# Intersect peaks from replicates
bedtools intersect \
  -a macs2/H3K4me3_rep1_peaks.broadPeak \
  -b macs2/H3K4me3_rep2_peaks.broadPeak \
  > macs2/H3K4me3_consensus_peaks.bed

# Add score column to MACS2 summits bed file for Allelome.PRO
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,"."}' \
  macs2/H3K4me3_rep2_summits.bed \
  > macs2/H3K4me3_rep2_summits.bed6

# Allelome.PRO with MAC2 peaks as annotation
SNPfile=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK_mod/00_src/Allelome.PRO2.sh \
  -i rep_1.clean.bam \
  -a macs2/H3K4me3_consensus_peaks.bed \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o .

# Allelome.PRO with MAC2 peaks as annotation
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK_mod/00_src/Allelome.PRO2.sh \
  -i rep_2.clean.bam \
  -a macs2/H3K4me3_consensus_peaks.bed \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o .


### H3K27ac CHIP-seq ###

# merge replicates
cat ENCFF216FFO.fastq.gz ENCFF517MJZ.fastq.gz ENCFF905MLO.fastq.gz ENCFF903ERX.fastq.gz ENCFF582ECO.fastq.gz ENCFF046LTE.fastq.gz ENCFF958XKB.fastq.gz ENCFF117KYH.fastq.gz > rep_1.fastq.gz
cat ENCFF731YYJ.fastq.gz ENCFF993HAM.fastq.gz  ENCFF149IQD.fastq.gz ENCFF072FEA.fastq.gz ENCFF334RGO.fastq.gz ENCFF641KXJ.fastq.gz ENCFF686LHJ.fastq.gz ENCFF698BDJ.fastq.gz > rep_2.fastq.gz

# Map reads
bwa-mem2 mem -t 8 ../bwa-index/GRCm39.primary_assembly.genome.fa rep_1.fastq.gz | samtools view -bSh > rep_1.bam
# Sort reads
samtools sort -o rep_1.sorted.bam -@ 8 rep_1.bam
# Filter reads: remove unmapped, secondary, QC-fail, duplicates; MAPQ>=30
samtools view -h -b -F 1804 -q 30 rep_1.sorted.bam > rep_1.filtered.bam
# Remove blacklisted regions
bedtools intersect -v -abam rep_1.filtered.bam -b blacklist_GRCm39.bed > rep_1.clean.bam

macs2 callpeak \
  -t rep_1.clean.bam \
  -f BAM \
  -g mm \
  --broad \
  -n H3K27ac_rep1 \
  --outdir macs2


# Map reads
bwa-mem2 mem -t 8 ../bwa-index/GRCm39.primary_assembly.genome.fa rep_2.fastq.gz | samtools view -bSh > rep_2.bam
# Sort reads
samtools sort -o rep_2.sorted.bam -@ 8 rep_2.bam
# Filter reads: remove unmapped, secondary, QC-fail, duplicates; MAPQ>=30
samtools view -h -b -F 1804 -q 30 rep_2.sorted.bam > rep_2.filtered.bam
# Remove blacklisted regions
bedtools intersect -v -abam rep_2.filtered.bam -b blacklist_GRCm39.bed > rep_2.clean.bam

macs2 callpeak \
  -t rep_2.clean.bam \
  -f BAM \
  -g mm \
  --broad \
  -n H3K27ac_rep2 \
  --outdir macs2


# Intersect peaks from replicates
bedtools intersect \
  -a macs2/H3K27ac_rep1_peaks.broadPeak \
  -b macs2/H3K27ac_rep2_peaks.broadPeak \
  > macs2/H3K27ac_consensus_peaks.bed

# Allelome.PRO with MAC2 peaks as annotation
SNPfile=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK_mod/00_src/Allelome.PRO2.sh \
  -i rep_1.clean.bam \
  -a macs2/H3K27ac_consensus_peaks.bed \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o .

# Allelome.PRO with MAC2 peaks as annotation
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK_mod/00_src/Allelome.PRO2.sh \
  -i rep_2.clean.bam \
  -a macs2/H3K27ac_consensus_peaks.bed \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o .

