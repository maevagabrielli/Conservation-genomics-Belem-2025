#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -e het.e
#SBATCH -o het.o

### subset set of snp to run the pca, take one SNP every 50kb, this was already run
#module load bioinfo/VCFtools/0.1.16
#vcftools --gzvcf /jarvis/cold/backup/usr/gabrielli/podarcis/raffonei/variant_calling/species/filtration/merge/Podarcis.all.allchr.allsites.miss0.75.mod.vcf.gz --thin 50000 --recode --recode-INFO-all --stdout | bgzip >Podarcis.all.allchr.miss0.75.thinned.vcf.gz

### get number of SNP
# in the vcf, it corresponds to the number of lines without the headers (starting by #)
zcat ../data/Podarcis.all.allchr.miss0.75.thinned.vcf.gz | grep -v ^# | wc -l
# for each group: need to create a subset of SNPs segregating in each of the group (use the option mac 1 of vcftools: at least one copy of the minor allele, which implies that this is a biallelic site)
# load vcftools
module load bioinfo/VCFtools/0.1.16
# loop over the 3 populations/species
for pop in LC ST WG
do
#create files with individual name for each group
cat ../data/Podarcis_inds_sp_pop.txt | awk -v pop=$pop '$3==pop' | cut -f 1 >$pop.list
#extract individuals from this group and with mac 1
vcftools --gzvcf ../data/Podarcis.all.allchr.miss0.75.thinned.vcf.gz --keep $pop.list --mac 1 --recode --out Podarcis.$pop.allchr.miss0.75.thinned.mac1
#nb of SNPs:
echo pop: $pop
echo nb of SNP: $(cat Podarcis.$pop.allchr.miss0.75.thinned.mac1.recode.vcf | grep -v ^# | wc -l)
done

### Compute individual count of heterozygotes
# load vcftools
module load bioinfo/VCFtools/0.1.16
# output statistics about homozygous and heterozygous counts with the --het option
vcftools --gzvcf ../data/Podarcis.all.allchr.miss0.75.thinned.vcf.gz --het --out Podarcis.all.allchr.miss0.75.thinned
# recapt het count
cat Podarcis.all.allchr.miss0.75.thinned.het | tail -n +2 | awk '{print $1 "\t" $4-$2}' >stats_het.txt #tail -n +2 removes the 1st line; awk print the first column (individual name) and the difference between the number of genotypes and the number of homozygous sites = number of heterozygous
