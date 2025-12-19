#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH -e 4b_count_genotypes.sh.e
#SBATCH -o 4b_count_genotypes.sh.o

#set working directory
#WD=~/work/podarcis/raffonei/workshop2025/genetic_load/SnpEff/polarization
#cd $WD

# load modules
module load bioinfo/bedtools/2.30.0 
module load bioinfo/VCFtools/0.1.16

# 1) extract positions OF but not FB
cat tagmiss | grep "OF" >tagANC #OF=Outgroup Fixed, so site have been polarized
cat tagANC | grep -v "FBXXXX" | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $4}' >OFsnp.bed #FB=Fixated in Both populations (SMALL and LARGE), -v = do not show, so output polymorphic site

#2) restrict to CDS
# the file CDS.bed has already been created in the previous count script
bedtools intersect -a OFsnp.bed -b ../CDS.bed >OFsnpCDS.bed

#3) Redo the vcf
bedtools intersect -header -a ../Podarcis.all.allchr.miss0.75.thinned.ann.vcf.gz -b OFsnpCDS.bed -u >Podarcis.all.allchr.miss0.75.CDS.thinned.ann.OF.snp.vcf

#4) Separate the vcf for each pop/sp
for pop in LC ST WG
do
vcftools --vcf Podarcis.all.allchr.miss0.75.CDS.thinned.ann.OF.snp.vcf --keep $pop.list --recode --recode-INFO-all --out Podarcis.$pop.allchr.miss0.75.CDS.thinned.ann.OF.snp
done

#5) count genotypes depending on the ancestral state
for pop in LC ST WG
do
vcf=Podarcis.$pop.allchr.miss0.75.CDS.thinned.ann.OF.snp.recode.vcf
max=$(cat $vcf | grep -v ^# | head -n 1 | awk '{print NF}') #get the number of fields in the vcf to know when to stop the column loop (individual loop)
for effect in synonymous missense HIGH # loop over the different effects
do
# create empty files
>nb_der_hom_"$effect"_$pop.txt
>nb_anc_hom_"$effect"_$pop.txt
>nb_het_"$effect"_$pop.txt
>nb_miss_"$effect"_$pop.txt
>nb_tot_"$effect"_$pop.txt
for ((id=10; id<=$max; id++))
do
#2 types of derived homozygous genotypes depending on the ancestral state (0 or 1), checked intersecting the vcf with the bed file with ancestral state
nbder00=$(cat OFsnpCDS.bed | awk '$4==1' | bedtools intersect -a <(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/') -b - | cut -f $id | grep -c -e ^"0/0" -e ^"0|0")
nbder11=$(cat OFsnpCDS.bed | awk '$4==0' | bedtools intersect -a <(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/') -b - | cut -f $id | grep -c -e ^"1/1" -e ^"1|1")
nbder=$(($nbder00+$nbder11))
#2 types of ancestral homozygous genotypes depending on the ancestral state (0 or 1), checked intersecting the vcf with the bed file with ancestral state
nbanc11=$(cat OFsnpCDS.bed | awk '$4==1' | bedtools intersect -a <(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/') -b - | cut -f $id | grep -c -e ^"1/1" -e ^"1|1")
nbanc00=$(cat OFsnpCDS.bed | awk '$4==0' | bedtools intersect -a <(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/') -b - | cut -f $id | grep -c -e ^"0/0" -e ^"0|0")
nbanc=$(($nbanc00+$nbanc11))
# then, heterozygous, missing, and total counts are not taking ancestral state into account
nbhet=$(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/' | grep -v ^"#" | cut -f $id | grep -c -e ^"0/1" -e ^"0|1" -e ^"1|0")
nbmiss=$(cat $vcf | awk -F'[,]' '/^#/ || $1~/'$effect'/' | grep -v ^"#" | cut -f $id | grep -c ^"\./\.")
nbtot=$(($nbder+$nbanc+$nbhet+$nbmiss))
# write the different files
echo $nbder >>nb_der_hom_"$effect"_$pop.txt
echo $nbhet >>nb_het_"$effect"_$pop.txt
echo $nbanc >>nb_anc_hom_"$effect"_$pop.txt
echo $nbmiss >>nb_miss_"$effect"_$pop.txt
echo $nbtot >>nb_tot_"$effect"_$pop.txt
done
done
done

# 6) recap counts
for effect in synonymous missense HIGH
do
>temp1
for pop in LC ST WG
do
vcf=Podarcis.$pop.allchr.miss0.75.CDS.thinned.ann.OF.snp.recode.vcf
max=$(cat $vcf | grep -v ^# | head -n 1 | awk '{print NF}') #get the number of fields in the vcf to know when to stop the column loop (individual loop)
for ((id=10; id<=$max; id++))
do
ind=$(head -n 10000 $vcf | grep "#CHR" | cut -f $id)
echo $ind >>temp1
done
echo -e "ind \t der.hom \t het \t miss \t nb.all \t tot.sites" >recap_count_"$effect"_$pop.txt
paste temp1 nb_der_hom_"$effect"_$pop.txt nb_het_"$effect"_$pop.txt nb_miss_"$effect"_$pop.txt nb_tot_"$effect"_$pop.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $2*2+$3 "\t" $5}' >>recap_count_"$effect"_$pop.txt
rm temp1
done
done
