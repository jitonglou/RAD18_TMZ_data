VCFTOOLSISEC=/nas/longleaf/apps/vcftools/0.1.15/bin/vcf-isec
samtools_dir=/nas/longleaf/apps/samtools/1.11
GATK=/nas/longleaf/apps/gatk/4.1.9.0/gatk-4.1.9.0/gatk
bedtools_dir=/nas/longleaf/apps/bedtools/2.29.0/bin

rawdata_dir=/proj/yanglb/users/jitong/GBM_PD20/raw_data
mutect2_dir=/proj/yanglb/users/jitong/GBM_PD20/results/Mutect2
BED=/proj/yanglb/users/jitong/GBM_PD20/ref/Yang_WES_GBM_exon_b37.bed

[ ! -d ${mutect2_dir}/filtered_vcfs ] && mkdir ${mutect2_dir}/filtered_vcfs
[ ! -d ${mutect2_dir}/filtered_vcfs/exon ] && mkdir ${mutect2_dir}/filtered_vcfs/exon
[ ! -d ${mutect2_dir}/filtered_vcfs/SNVs ] && mkdir ${mutect2_dir}/filtered_vcfs/SNVs
[ ! -d ${mutect2_dir}/filtered_vcfs/INDELs ] && mkdir ${mutect2_dir}/filtered_vcfs/INDELs

for ii in `ls ${rawdata_dir} | sed -n '2,50p'`; do 
  # # $samtools_dir/htslib/bin/bgzip ${mutect2_dir}/${ii}.m2.postprocessed.mpos.vcf
  # # $samtools_dir/htslib/bin/tabix -p vcf ${mutect2_dir}/${ii}.m2.postprocessed.mpos.mutation_removed.vcf.gz
  gunzip -f ${mutect2_dir}/${ii}.m2.postprocessed.mpos.mutation_removed.vcf.gz
  cp ${mutect2_dir}/${ii}.m2.postprocessed.mpos.mutation_removed.vcf \
  ${mutect2_dir}/filtered_vcfs/${ii}.m2.postprocessed.mpos.mutation_removed.vcf
done

VCFLIST=`ls ${mutect2_dir} | grep "mpos.vcf.gz" | cut -d "." -f1-6 | uniq | sed -n '1,36p'`
echo $VCFLIST | cut -d " " -f36-37


cd ${mutect2_dir}/filtered_vcfs
for ii in `ls *.vcf | cut -d "." -f1-5 | sed -n '2,50p'`; do
  head -162 ${ii}.vcf > exon/${ii}.exon.vcf;
  $bedtools_dir/intersectBed \
  -a ${ii}.vcf \
  -b $BED >> exon/${ii}.exon.vcf;
  $GATK SelectVariants --select-type-to-include SNP \
  -V exon/${ii}.exon.vcf \
  -O SNVs/snvs.${ii}.exon.vcf;
  $GATK SelectVariants --select-type-to-include INDEL \
  -V exon/${ii}.exon.vcf \
  -O INDELs/indels.${ii}.exon.vcf; 
done