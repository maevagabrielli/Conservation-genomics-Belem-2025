#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH -e 3a_SnpEff.sh.e
#SBATCH -o 3a_SnpEff.sh.o

### 1 change the config file adding the information of your genome
# path: ~/work/podarcis/raffonei/workshop2025/snpEff/snpEff.config, this was already done but you can check it

### 2 build the database
# load java (version needs to be higher than java21)
module load devel/java/24.0.2
#java -Xmx4g -jar ~/work/podarcis/raffonei/workshop2025/snpEff/snpEff.jar build -gff3 -v rPodRaf1, this was already run

### 3 annotate the vcf
java -Xmx4g -jar ../snpEff/snpEff.jar rPodRaf1 ../data/Podarcis.all.allchr.miss0.75.thinned.vcf.gz | gzip -c >Podarcis.all.allchr.miss0.75.thinned.ann.vcf.gz
