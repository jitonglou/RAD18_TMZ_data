#!/bin/sh
#SBATCH -J YangGBM_PD20_alah_alah # jobname
#SBATCH -p general # Partition to submit to (queue)
#SBATCH -N 1 # number of nodes, Ensure that all cores are on one machine
#SBATCH -n 1 # number of cores
#SBATCH --mem 4000 # memory pool for all cores
#SBATCH -t 0-04:00 # time (D-HH:MM)
#SBATCH -o /proj/yanglb/users/jitong/GBM_PD20/logs/all.alah_alah.%j.out # output file, %N name of nodelist, %j jobID
#SBATCH -e /proj/yanglb/users/jitong/GBM_PD20/logs/all.alah_alah.%j.err # output file, %N name of nodelist, %j jobID

### set directories
main_dir=/proj/yanglb
rawdata_dir=/proj/yanglb/users/jitong/GBM_PD20/raw_data

user_dir=${main_dir}/users/jitong
proj_dir=${user_dir}/GBM_PD20

ref_dir=${proj_dir}/ref
records_dir=${proj_dir}/records
results_dir=${proj_dir}/results
temp_dir=${proj_dir}/temp
logs_dir=${proj_dir}/logs


# Tools
THREADS=8
RAM=4 # in GB
FASTQC=/nas/longleaf/apps/fastqc/0.11.8/FastQC/fastqc
trim_dir=/proj/yanglb/users/jitong/Downloads/Trimmomatic-0.39
BWA=/nas/longleaf/apps/bwa/0.7.17/bin/bwa
PICARD_JAR=/nas/longleaf/apps/picard/2.23.4/picard-2.23.4/picard.jar
GATK=/nas/longleaf/apps/gatk/4.1.9.0/gatk-4.1.9.0/gatk
MULTIQC=/nas/longleaf/apps/multiqc/1.7/bin/multiqc
snpeff_dir=/nas/longleaf/apps/snpeff/4.3t/snpEff
samtools_dir=/nas/longleaf/apps/samtools/1.11
PERL=/nas/longleaf/apps/perl/5.18.2/bin/perl

# genome reference and exon cover region
REF=${user_dir}/Downloads/b37/Homo_sapiens_assembly19.fasta # broadinstitute b37
BED=${ref_dir}/S07604514_Covered.bed
KNOWN_VARIATION=${user_dir}/Downloads/GRCh37_dbSNP151/common_all_20180423.vcf.gz


## sample_name is like EV11_CKDN200001674-1A_HNLGYDSXX_L1_1.fq.gz
name=RT2;
fullname=RT2_CKDN210000105-1A_HVFN3DSXY_L4;
sample_dir=${rawdata_dir}/${name}


timelog=${logs_dir}/all_${name}.log

# output directories
fastqc_dir=${results_dir}/QC/fastqc
bams_dir=${results_dir}/bams
bqsr_dir=${results_dir}/QC/BQSR
picardmetric_dir=${results_dir}/QC/PicardMetric
multiqc_dir=${results_dir}/QC/multiqc
mutect2_dir=${results_dir}/Mutect2

#### workflow start ####
echo "Start at -> `date`; " > $timelog

#### step 1 - FastQC ####
## Quality control for raw reads 
[ ! -d ${results_dir}/QC ] && mkdir ${results_dir}/QC
[ ! -d ${fastqc_dir} ] && mkdir ${fastqc_dir}

echo "FastQC analysis start -> `date`; " >> $timelog

$FASTQC -t $THREADS \
-o ${fastqc_dir} \
${sample_dir}/${fullname}_1.fq.gz \
${sample_dir}/${fullname}_2.fq.gz \
> ${records_dir}/fastqc.raw.${name}.log \
2>&1

echo "FastQC analysis done -> `date`; " >> $timelog
# -t threads


#### step 2 - Trimmomatic ####
## Discard short reads and reads with insufficient base qualities
## http://www.usadellab.org/cms/?page=trimmomatic
## Truncation might cause error when sorting the BAM in step 5 !!!!! https://github.com/nanoporetech/medaka/issues/176
echo "Trim short reads start -> `date`; " >> $timelog

java -jar $trim_dir/trimmomatic-0.39.jar PE -threads $THREADS -phred33 \
${sample_dir}/${fullname}_1.fq.gz \
${sample_dir}/${fullname}_2.fq.gz \
${temp_dir}/${name}.R1.passed.fastq.gz \
${temp_dir}/${name}.R1.not_passed.fastq.gz \
${temp_dir}/${name}.R2.passed.fastq.gz \
${temp_dir}/${name}.R2.not_passed.fastq.gz \
LEADING:25 TRAILING:25 MINLEN:50 \
SLIDINGWINDOW:10:25 \
ILLUMINACLIP:$trim_dir/adapters/TruSeq3-PE-2.fa:2:30:10 \
> ${records_dir}/trimmomatic.${name}.log \
2>&1

echo "Trim short reads done -> `date`; " >> $timelog
# PE pair end
# Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
# Remove leading low quality or N bases (below quality 3) (LEADING:3)
# Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
# Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
# Drop reads below the 36 bases long (MINLEN:36)


#### step 3 - FastQC ####
## Qulity control for trimmed reads
echo "FastQC trimmed reads start -> `date`; " >> $timelog

$FASTQC -t $THREADS \
-o ${fastqc_dir} \
${temp_dir}/${name}.R1.passed.fastq.gz \
${temp_dir}/${name}.R2.passed.fastq.gz \
> ${records_dir}/fastqc.passed.${name}.log \
2>&1

echo "FastQC trimmed reads done -> `date`; " >> $timelog
# -t threads

#### step 4 - Align with bwa-mem ####
[ ! -d ${results_dir}/bams ] && mkdir ${results_dir}/bams

## Build index for reference genome
echo "bwa ref index start -> `date`; " >> $timelog
$BWA index $REF -t $THREADS
2> ${records_dir}/bwa.ref.index.log
echo "bwa ref index done -> `date`; " >> $timelog

## tumor ##
echo "bwa mem start -> `date`; " >> $timelog

# bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
$BWA mem $REF -t $THREADS -v 1 \
${temp_dir}/${name}.R1.passed.fastq.gz \
${temp_dir}/${name}.R2.passed.fastq.gz \
2> ${records_dir}/bwa.${name}.log \
> ${temp_dir}/${name}.aligned.sam

echo "bwa mem done -> `date`; " >> $timelog
# -v verbose level
# 2> : save standard error to file

#### step 5 - Postprocess of aligned reads ####
## Convert sam to bam using picard (1 hour)
echo "sam to bam start -> `date`; " >> $timelog

java -jar $PICARD_JAR CleanSam \
I=${temp_dir}/${name}.complete.aligned.sam \
2> ${records_dir}/SamToBam.${name}.log \
O=${temp_dir}/${name}.aligned.cleaned.bam

echo "sam to bam done -> `date`; " >> $timelog

## Sort bam files using samtools
echo "sort bam start -> `date`; " >> $timelog

$samtools_dir/bin/samtools sort -@ $THREADS \
${temp_dir}/${name}.aligned.cleaned.bam \
2> ${records_dir}/SortBam.${name}.log \
-o ${temp_dir}/${name}.aligned.cleaned.sorted.bam

echo "sort bam done -> `date`; " >> $timelog

## Add reading groups using picard (1 hour)
echo "Add reading group start -> `date`; " >> $timelog

java -jar $PICARD_JAR AddOrReplaceReadGroups \
I=${temp_dir}/${name}.aligned.cleaned.sorted.bam \
2> ${records_dir}/AddReadGroups.${name}.log \
O=${temp_dir}/${name}.aligned.cleaned.sorted.readgroups.bam \
RGID=${name} RGLB=Lib_${name} RGPL=Illumina RGPU=Run_${name} RGSM=Tumor

echo "Add reading group done -> `date`; " >> $timelog
# -ID 1 -LB Lib1 -PL ILLUMINA -PU Run1 -SM $type

## Mark duplicated reads using picard (need about 50GB memory, 30 minutes)
[ ! -d ${results_dir}/QC/MarkDuplicates ] && mkdir ${results_dir}/QC/MarkDuplicates
MAX_RECORDS_IN_RAM=$(expr $RAM \* 250000)

echo "Mark duplicates start -> `date`; " >> $timelog

java -jar $PICARD_JAR MarkDuplicates \
I=${temp_dir}/${name}.aligned.cleaned.sorted.readgroups.bam \
2> ${records_dir}/MarkDuplicates.${name}.log \
O=${temp_dir}/${name}.aligned.cleaned.sorted.readgroups.marked.bam \
METRICS_FILE=${results_dir}/QC/MarkDuplicates/${name}.duplicate_metrics.txt \
REMOVE_DUPLICATES=false ASSUME_SORTED=true \
MAX_RECORDS_IN_RAM=$MAX_RECORDS_IN_RAM

echo "Mark duplicates end -> `date`; " >> $timelog

## Validate bam again using Picard (10 minutes)
echo "Validate bam start -> `date`; " >> $timelog

java -jar $PICARD_JAR ValidateSamFile \
I=${temp_dir}/${name}.aligned.cleaned.sorted.readgroups.marked.bam \
MODE=SUMMARY \
2> ${records_dir}/ValidateBam.${name}.log \

echo "Validate bam done -> `date`; " >> $timelog

#### step 6 - Base recalibration #### (4 hours)
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
KNOWN_VARIATION=${user_dir}/Downloads/GRCh37_dbSNP151/common_all_20180423.vcf.gz

[ ! -d ${results_dir}/QC/BQSR ] && mkdir ${results_dir}/QC/BQSR

# Build index for vcf, which is needed for BQSR
$GATK IndexFeatureFile -F $KNOWN_VARIATION \
2> ${records_dir}/KnownVariationIndex.${name}.log
# -F feature file

## Analyze patterns of covariation in the sequence dataset
echo "BQSR pre start -> `date`; " >> $timelog

$GATK BaseRecalibrator \
-I ${temp_dir}/${name}.aligned.cleaned.sorted.readgroups.marked.bam \
-R $REF \
--known-sites $KNOWN_VARIATION \
--use-original-qualities \
2> ${records_dir}/BQSR.pre.${name}.log \
-O ${bqsr_dir}/${name}.GATK4.pre.recal.table # in r table format

echo "BQSR pre done -> `date`; " >> $timelog
# --use-original-qualities / -OQ: Use the base quality scores from the OQ tag
# This flag tells GATK to use the original base qualities (that were in the data before BQSR/recalibration) which are stored in the OQ tag, if they are present, rather than use the post-recalibration quality scores. If no OQ tag is present for a read, the standard qual score will be used.
# This flag is used if a samples has already been processed by BQSR and we are re-running BQSR, we want the tool to use the original qualities. 


## Apply the recalibration to your sequence data
echo "Apply recalibration start -> `date`; " >> $timelog

$GATK ApplyBQSR \
-I ${temp_dir}/${name}.aligned.cleaned.sorted.readgroups.marked.bam \
-R $REF \
-bqsr ${bqsr_dir}/${name}.GATK4.pre.recal.table \
2> ${records_dir}/BQSR.apply.${name}.log \
-O ${bams_dir}/${name}.ready.bam

echo "Apply recalibration done -> `date`; " >> $timelog

## Post BQSR table pattern
echo "BQSR post start -> `date`; " >> $timelog

$GATK BaseRecalibrator \
-I ${bams_dir}/${name}.ready.bam \
-R $REF \
--known-sites $KNOWN_VARIATION \
--use-original-qualities \
2> ${records_dir}/BQSR.post.${name}.log \
-O ${bqsr_dir}/${name}.GATK4.post.recal.table

echo "BQSR post done -> `date`; " >> $timelog

$samtools_dir/bin/samtools index -@ $THREADS ${bams_dir}/${name}.ready.bam
rm ${bams_dir}/${name}.ready.bai


#### step 7 - Quality control of the alignments ####
[ ! -d ${results_dir}/QC/PicardMetric ] && mkdir ${results_dir}/QC/PicardMetric
[ ! -d ${results_dir}/QC/multiqc ] && mkdir ${results_dir}/QC/multiqc

echo "picard CollectSequencingArtifactMetrics start -> `date`; " >> $timelog

java -jar $PICARD_JAR CollectSequencingArtifactMetrics \
R=$REF \
I=${bams_dir}/${name}.ready.bam \
2> ${records_dir}/picard.Artifact.Metrics.${name}.log \
O=${picardmetric_dir}/${name}.ready.bam.artifacts

echo "picard CollectSequencingArtifactMetrics end -> `date`; " >> $timelog


echo "picard CollectMultipleMetrics start -> `date`; " >> $timelog

java -jar $PICARD_JAR CollectMultipleMetrics \
R=$REF \
I=${bams_dir}/${name}.ready.bam \
2> ${records_dir}/picard.Multiple.Metrics.${name}.log \
O=${picardmetric_dir}/${name}.ready.bam.metrics

echo "picard CollectMultipleMetrics end -> `date`; " >> $timelog

$samtools_dir/bin/samtools idxstats ${bams_dir}/${name}.ready.bam \
> ${picardmetric_dir}/${name}.ready.bam.idxstats

## Calculate metrics for WES, including sequencing coverage
# # convert bed to interval list (add header lines)  
# java -jar $PICARD_JAR BedToIntervalList \
      # I=$BED \
      # O=${ref_dir}/list.interval_list \
      # SD=${user_dir}/Downloads/b37/Homo_sapiens_assembly19.dict
# # SD: The sequence dictionary, or BAM/VCF/IntervalList from which a dictionary can be extracted.
INTERVAL=${ref_dir}/list.interval_list

echo "picard CollectHsMetrics start -> `date`; " >> $timelog

java -jar $PICARD_JAR CollectHsMetrics \
SAMPLE_SIZE=100000 \
R=$REF \
I=${bams_dir}/${name}.ready.bam \
2> ${records_dir}/picard.Hs.Metrics.${name}.log \
O=${picardmetric_dir}/${name}.ready.bam.metrics \
BAIT_INTERVALS=$INTERVAL \
TARGET_INTERVALS=$INTERVAL

echo "picard CollectHsMetrics end -> `date`; " >> $timelog


#### step 8 - SNVs and small indels ####
[ ! -d ${results_dir}/Mutect2 ] && mkdir ${results_dir}/Mutect2

echo "Mutect2 caller start -> `date`; " >> $timelog

$GATK Mutect2 \
--native-pair-hmm-threads $THREADS \
-R $REF \
-I ${bams_dir}/${name}.ready.bam \
2> ${records_dir}/mutect2.call.${name}.log \
-O ${mutect2_dir}/${name}.m2.vcf \
-bamout ${mutect2_dir}/${name}.m2.bam

echo "Mutect2 caller done -> `date`; " >> $timelog

#### step 9 - Remove probable technical or germline artifacts ####
echo "GATK FilterMutectCalls start -> `date`; " >> $timelog

$GATK FilterMutectCalls \
--reference $REF \
--variant ${mutect2_dir}/${name}.m2.vcf \
2> ${records_dir}/mutect2.FilterMutectCalls.${name}.log \
--output ${mutect2_dir}/${name}.m2.filt.vcf

echo "GATK FilterMutectCalls done -> `date`; " >> $timelog

#### step 10 - Filter FFPR artifacts (C/T) and oxidative DNA damage artifacts (G/T) ####
## Execute this step only if there is sufficient evidence that these samples are affected by one of these technical artifacts
## need do step 7 CollectSequencingArtifactMetrics first

echo "No artifact filtering start -> `date`; " >> $timelog

cp ${mutect2_dir}/${name}.m2.filt.vcf ${mutect2_dir}/${name}.m2.filt.AM.vcf
cp ${mutect2_dir}/${name}.m2.filt.AM.vcf ${mutect2_dir}/${name}.m2.filt.AM.filtered.vcf


#### step 11 - filter out all indels >10 bp as follows ####
echo "GATK SelectVariants filter long indels start -> `date`; " >> $timelog
$GATK SelectVariants --max-indel-size 10 \
-V ${mutect2_dir}/${name}.m2.filt.AM.filtered.vcf \
2> ${records_dir}/gatk.FilterByOrientationBias.${name}.log \
-O ${mutect2_dir}/${name}.m2.filt.AM.filtered.selected.vcf

echo "GATK SelectVariants filter long indels done -> `date`; " >> $timelog

#### step 12 - Use additional filters to decrease the false-positive rate ####
# apply filters for mutant allele frequency (≥10%), coverage at particular positions in tumor and
# normal samples (≥10×) and supporting reads for the mutation in the tumor sample (at least three)

echo "SnpSift filter false-positive mutations start -> `date`; " >> $timelog

cat ${mutect2_dir}/${name}.m2.filt.AM.filtered.selected.vcf \
| java -jar $snpeff_dir/SnpSift.jar filter \
"((FILTER = 'PASS') & (GEN[Tumor].AF >= 0.1) & \
(GEN[Tumor].AD[0] + GEN[Tumor].AD[1] >= 10) & \
(GEN[Tumor].AD[1] >= 3))" \
2> ${records_dir}/SnpSift.filter.AFAD.${name}.log \
> ${mutect2_dir}/${name}.m2.postprocessed.vcf
echo "SnpSift filter false-positive mutations done -> `date`; " >> $timelog

#### step 13 - To further reduce false-positive callings, compare the SNVs and indels to known polymorphisms #### (1 minute)
## my filtering, MPOS=-2147483648 in the first vcf (calculated by GATK, indicates missing value) will cause bcftools crash
echo "SnpSift filter negative MPOS (median distance from end of read) start -> `date`; " >> $timelog

cat ${mutect2_dir}/${name}.m2.postprocessed.vcf \
| java -jar $snpeff_dir/SnpSift.jar filter \
"(MPOS >= 0)" \
2> ${records_dir}/SnpSift.filter.MPOS.${name}.log \
> ${mutect2_dir}/${name}.m2.postprocessed.mpos.vcf

$samtools_dir/htslib/bin/bgzip ${mutect2_dir}/${name}.m2.postprocessed.mpos.vcf
$samtools_dir/htslib/bin/tabix -p vcf ${mutect2_dir}/${name}.m2.postprocessed.mpos.vcf.gz

echo "SnpSift filter negative MPOS (median distance from end of read) done -> `date`; " >> $timelog

ALTERNATIVE_VARIATION=${results_dir}/Mutect2/RRef.m2.postprocessed.mpos.vcf.gz

echo "bcftools isec compare to known polymorphisms(RRef) start -> `date`; " >> $timelog

$samtools_dir/bin/bcftools isec -C -c none -O z -w 1 \
-o ${mutect2_dir}/${name}.m2.postprocessed.mpos.mutation_removed.vcf.gz \
${mutect2_dir}/${name}.m2.postprocessed.mpos.vcf.gz \
$ALTERNATIVE_VARIATION

echo "bcftools isec compare to known polymorphisms(RRef) done -> `date`; " >> $timelog
# http://samtools.github.io/bcftools/bcftools.html#common_options
# -C, --complement: output positions present only in the first file but missing in the others
# -c, collapse: none, only records with identical REF and ALT alleles are compatible
# -O, --output-type: z, compressed VCF
# -w, --write: list of input files to output given as 1-based indices.

$samtools_dir/bin/bcftools norm -m -any \
${mutect2_dir}/${name}.m2.postprocessed.mpos.mutation_removed.vcf.gz \
-O z -o ${mutect2_dir}/${name}.Mutect2.vcf.gz

gunzip -f ${mutect2_dir}/${name}.Mutect2.vcf.gz

echo "All is done -> `date`; " >> $timelog
