#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -e PCA.e
#SBATCH -o PCA.o

### PCA
# load plink
module load bioinfo/PLINK/1.90b7
# define species and individuals to work with
species=raffonei-waglerianus
# creates a list with all raffonei and waglerianus individuals
grep raffonei ../data/ind_species.txt | cut -f 1 >temp
grep waglerianus ../data/ind_species.txt | cut -f 1 >>temp
# use plink format for individuals names: needs 2 columsnn (Family ID (FID); Individual ID (IID)),so we duplicate the individual names
paste temp temp >$species.list
# run the PCA
plink --vcf ../data/Podarcis.all.allchr.miss0.75.thinned.vcf.gz --pca 100 --keep $species.list --mac 1 --out Podarcis.$species.allchr.miss0.75.thinned --allow-extra-chr
# and then plot in R with plot_pca_plink.R
