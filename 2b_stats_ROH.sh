#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -e stats_ROHan.e
#SBATCH -o stats_ROHan.o

###compute the percentage of ROH for each individual run in the window100000 folder
echo -e "ind \t PCT_ROH \t min_PCT_ROH \t max_PCT_ROH" >pct_roh_calc.txt
ls window100000/*.summary.txt | cut -f 2 -d "/" | cut -f 1 -d "." >temp1

cat window100000/*.summary.txt | grep "Segments in ROH" | grep -v "%" | cut -f 2 | cut -f 1 -d " " >temp2
cat window100000/*.summary.txt | grep "Segments in ROH" | grep -v "%" | cut -f 2 | cut -f 2 -d " " | cut -f 1 -d "," | sed 's/(//g' >temp3
cat window100000/*.summary.txt | grep "Segments in ROH" | grep -v "%" | cut -f 2 | cut -f 2 -d " " | cut -f 2 -d "," | sed 's/)//g' >temp4

for file in window100000/*.mid.hmmp.gz; do zcat $file  | tail -n +2 | wc -l; done >temp5

paste temp1 temp2 temp3 temp4 temp5 | awk '{print $1 "\t" $2/($5*100000)*100 "\t" $3/($5*100000)*100 "\t" $4/($5*100000)*100}' >>pct_roh_calc.txt
rm temp*

### calculate ROH length
# load bedtools, used to merge consecutive windows in ROH
module load bioinfo/bedtools/2.30.0
ls window100000/*summary.txt | cut -f 1 -d "." | cut -f 2 -d "/" | while read ind
do
	zcat window100000/$ind.mid.hmmp.gz | awk '$4 >= 0.99 && $5!="1"' | bedtools merge -i - | awk '{print $3-$2}' >length_ROH_$ind.txt 
	#awk: probability of ROH higher than 0.99 and probability of non-ROH different than 1 as gives a non-0 very small value for p(ROH) when p(non-ROH)=1 that awk cannot read
	# bedtools merge: merge windows if consecutive
	# awk '{print $3-$2}' gives the length of the window
done
