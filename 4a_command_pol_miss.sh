#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH -e 4a_command_pol_miss.sh.e
#SBATCH -o 4a_command_pol_miss.sh.o

WD=~/work/podarcis/raffonei/workshop2025/genetic_load/polarization
cd $WD

# define the groups of interest for the polarization, the small and large populations/sp of study and the 2 outgroups
SMALL=raffonei
LARGE=waglerianus
OUT1=tiliguerta
OUT2=siculus
GZVCF=../../data/Podarcis.all.allchr.miss0.75.thinned.vcf.gz

zcat $GZVCF | head -n 10000 | grep "#CHR" >ind_vcf.list
for pop in $SMALL $LARGE $OUT1 $OUT2
do
>$pop.id
cat ind_species.txt | awk -v pop="$pop" '$2==pop' | cut -f 1 | while read ind
do
id=$(cat ind_vcf.list | tr -s '\t' '\n' | grep -n "$ind" | cut -f 1 -d ":") #individual column number
idpyt=$(($id-1)) #nb in python (0-based)
echo -n "$idpyt " >>$pop.id
done
echo >>$pop.id
done

python3 polarization_miss_ancestralallele.py -vcf $GZVCF -out tagmiss -OUT1 $(cat "$OUT1".id) -OUT2 $(cat "$OUT2".id) -SMALL $(cat "$SMALL".id) -LARGE $(cat "$LARGE".id)

