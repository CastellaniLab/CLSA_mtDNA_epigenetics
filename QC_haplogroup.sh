cd /home/jnguy254/projects/def-ccastel/SharedResources/CLSA_private/data/genomics

module load plink
module load java

#Convert bed file to vcf, remove indels and autosomal variants
plink --bfile clsa_gen_v3 --chr mt --snps-only 'just-acgt' --recode vcf --out clsa_qc

#Classify haplogroups using haplogrep
./haplogrep classify --in clsa_snps.vcf --format vcf --chip --out haplogroups.txt

mv haplogroups.txt ../../results/genomics
