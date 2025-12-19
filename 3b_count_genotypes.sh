#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH -e 3b_count_genotypes.sh.e
#SBATCH -o 3b_count_genotypes.sh.o

# load modules
module load bioinfo/bedtools/2.30.0 
module load bioinfo/VCFtools/0.1.16

### Here, the ancestral allele is the one of the reference genome (0) and 1 is the derived allele

### 1) Restrict to CDS
# get CDS coordinates
cat ../snpEff/data/rPodRaf1/genes.gff | awk '$3=="CDS"' | awk '{print $1 "\t" $4 "\t" $5}' | bedtools sort -i - | bedtools merge -i - >CDS.bed
# length of CDS ? cat CDS.bed | awk '{sum+=$3-$2}END{print sum}' #34618025
# get snp coordinates in bed format
zcat Podarcis.all.allchr.miss0.75.thinned.ann.vcf.gz | grep -v ^# | awk '{print $1 "\t" $2-1 "\t" $2}' >snp.bed
# intersect with CDS to get coordinates of snp within CDS
bedtools intersect -a snp.bed -b CDS.bed >snpCDS.bed

### 2) Output vcf for CDS
bedtools intersect -header -a Podarcis.all.allchr.miss0.75.thinned.ann.vcf.gz -b snpCDS.bed -u >Podarcis.all.allchr.miss0.75.thinned.ann.CDS.vcf

### Create a vcf for each population/pop
for pop in LC ST WG
do 
vcftools --vcf Podarcis.all.allchr.miss0.75.thinned.ann.CDS.vcf --keep $pop.list --recode --recode-INFO-all --out Podarcis.$pop.allchr.miss0.75.thinned.ann.CDS
done

### 3) count genotypes
for pop in LC ST WG # loop over the different populations
do
vcf=Podarcis.$pop.allchr.miss0.75.thinned.ann.CDS.recode.vcf #define the population vcf to use
max=$(cat $vcf | grep -v ^# | head -n 1 | awk '{print NF}') #get the number of fields in the vcf to know when to stop the column loop (individual loop)
for effect in synonymous missense HIGH #loop over the different effects
do
# create empty files before the loop
>nb_der_hom_"$effect"_$pop.txt
>nb_anc_hom_"$effect"_$pop.txt
>nb_het_"$effect"_$pop.txt
>nb_miss_"$effect"_$pop.txt
>nb_tot_"$effect"_$pop.txt
for ((id=10; id<=$max; id++)) #loop over individuals: in the vcf, individual columns are between column 10 and the max number of columns
do
ind=$(head -n 10000 $vcf | grep "#CHR" | cut -f $id)
nbanc=$(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/' | grep -v ^"#" | cut -f $id | grep -c -e ^"0/0" -e ^"0|0") #ind nb of ancestral homozygous
nbder=$(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/' | grep -v ^"#" | cut -f $id | grep -c -e ^"1/1" -e ^"1|1") #ind nb of derived homozygous
nbhet=$(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/' | grep -v ^"#" | cut -f $id | grep -c -e ^"0/1" -e ^"0|1" -e ^"1|0") #ind nb of heterozygous genotypes
nbmiss=$(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/' | grep -v ^"#" | cut -f $id | grep -c -e ^"\./\." -e ^"\.|\.") #ind nb of missing genotypes
nbtot=$(($nbder+$nbanc+$nbhet+$nbmiss)) #ind nb of genotypes, check that it is the same for all individuals of a given population
echo $nbder >>nb_der_hom_"$effect"_$pop.txt #write nb of derived homozygous
echo $nbhet >>nb_het_"$effect"_$pop.txt #write nb of heterozygous genotypes
echo $nbanc >>nb_anc_hom_"$effect"_$pop.txt #write nb of ancestral homozygous
echo $nbmiss >>nb_miss_"$effect"_$pop.txt #write nb of missing genotypes
echo $nbtot >>nb_tot_"$effect"_$pop.txt #write nb of genotypes
done
done
done

### 4) recap counts
#paste together the different files to do a recap file for each effect and population
for effect in synonymous missense HIGH
do
for pop in LC ST WG
do
>temp1
>temp2
vcf=Podarcis.$pop.allchr.miss0.75.thinned.ann.CDS.recode.vcf
max=$(cat $vcf | grep -v ^# | head -n 1 | awk '{print NF}')
for ((id=10; id<=$max; id++))
do
ind=$(head -n 10000 $vcf | grep "#CHR" | cut -f $id)
echo $ind >>temp1
done
paste nb_der_hom_"$effect"_$pop.txt nb_het_"$effect"_$pop.txt nb_miss_"$effect"_$pop.txt nb_tot_"$effect"_$pop.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $1*2+$2 "\t" $4}' >>temp2
echo -e "ind \t der.hom \t het \t miss \t nb.all \t tot.sites" >recap_count_"$effect"_$pop.txt
paste temp1 temp2 >>recap_count_"$effect"_$pop.txt
rm temp1 temp2
done
done
