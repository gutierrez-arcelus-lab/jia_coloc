#!/usr/bin/sh

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=8:00:00
#SBATCH -p bch-compute
#SBATCH --array=1
#SBATCH --job-name=LD1KGP
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/LD1KGP-%A-%a

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export R_X=4.1

IDS=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/colocalization/finemap/data/europeans_samples.txt
REGIONSFILE=${SLURM_SUBMIT_DIR}/data/gwas_region.txt
REGION=$( awk -v ARR="$SLURM_ARRAY_TASK_ID" 'FNR == ARR' $REGIONSFILE )
CHR=$( echo $REGION | cut -d':' -f1 )
VCF1K=/reference_databases/1000G_VCF/GRCh38/Genotype_VCFs/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.recalibrated_variants.vcf.gz
PREFIX=${SLURM_SUBMIT_DIR}/data/${REGION}
VCFOUT=${PREFIX}.vcf.gz

# IMPORTANT:
## Selecting variants that are polymorphich (non-monomorphic for either REF or ALT);
## Excluding singletons;
## collapse variants after these filterings so we can recover more matches to the GWAS
## for not excluding variants that are multiallelics just because one occurrence of a 2nd ALT;
## Finally select biallelics

bcftools view --threads $SLURM_CPUS_PER_TASK -r $REGION $VCF1K |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $IDS --force-samples |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --min-ac 1:nref - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -e INFO/AC==INFO/AN - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -e INFO/AC==1 - |\
    bcftools norm --threads $SLURM_CPUS_PER_TASK -m +both - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x FORMAT --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --genotype ^miss -m2 -M2 -O z -o $VCFOUT -

Rscript ${SLURM_SUBMIT_DIR}/subset_vcf.R $VCFOUT ${PREFIX}_JIApositions.vcf.gz

plink --vcf ${PREFIX}_JIApositions.vcf.gz \
    --r square spaces \
    --threads $SLURM_CPUS_PER_TASK \
    --freqx \
    --out $PREFIX

rm ${PREFIX}.nosex
