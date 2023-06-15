#!/bin/bash
sniffles_path='/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/myelong/results_concatenated_samples/sniffles/patient_sniffles.vcf'
cutesv_path='/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/myelong/results_concatenated_samples/cutesv/patient_cutesv.vcf'

printf "INFO \t SNIFFLES \t CUTE\n" > bcftools/analysis.tsv

total_variants_snif=$(bcftools query -f '%CHROM\n' $sniffles_path | grep -v -c '^#')
total_variants_cute=$(bcftools query -f '%CHROM\n' $cutesv_path | grep -v -c '^#')

chromlist=($(bcftools query -f '%CHROM\n' $sniffles_path| uniq|grep 'chr'))

snps_snif=$(bcftools view -v snps $sniffles_path |grep -v -c '^#')
snps_cute=$(bcftools view -v snps $cutesv_path |grep -v -c '^#')

indels_snif=$(bcftools view -v indels $sniffles_path |grep -v -c '^#')
indels_cute=$(bcftools view -v indels $cutesv_path |grep -v -c '^#')

bnd_snif=$(bcftools view -v bnd $sniffles_path |grep -v -c '^#')
bnd_cute=$(bcftools view -v bnd $cutesv_path |grep -v -c '^#')

del_snif=$(bcftools query -f '%ID\n' $sniffles_path |grep 'Sniffles2.DEL'| grep -v -c '^#')
del_cute=$(bcftools query -f '%ID\n' $cutesv_path |grep 'cuteSV.DEL'| grep -v -c '^#')

ins_snif=$(bcftools query -f '%ID\n' $sniffles_path |grep 'Sniffles2.INS'| grep -v -c '^#')
ins_cute=$(bcftools query -f '%ID\n' $cutesv_path |grep 'cuteSV.INS'| grep -v -c '^#')

inv_snif=$(bcftools query -f '%ID\n' $sniffles_path |grep 'Sniffles2.INV'| grep -v -c '^#')
inv_cute=$(bcftools query -f '%ID\n' $cutesv_path |grep 'cuteSV.INV'| grep -v -c '^#')

dup_snif=$(bcftools query -f '%ID\n' $sniffles_path |grep 'Sniffles2.DUP'| grep -v -c '^#')
dup_cute=$(bcftools query -f '%ID\n' $cutesv_path |grep 'cuteSV.DUP'| grep -v -c '^#')

del_snif_float=$del_snif.0
del_cute_float=$del_cute.0

ins_snif_float=$ins_snif.0
ins_cute_float=$ins_cute.0


indels_ratio_snif=$(awk "BEGIN {print $del_snif_float/$ins_snif_float}")
indels_ratio_cute=$(awk "BEGIN {print $del_cute_float/$ins_cute_float}")
#indels_ratio_cute=$(($del_cute_float / $ins_cute_float))


printf "total variants \t $total_variants_snif \t $total_variants_cute \n" >> results/bcftools/analysis.tsv

for chrom in ${chromlist[@]}
do
count_snif=$(bcftools view -t $chrom $sniffles_path | grep -v -c '^#')
count_cute=$(bcftools view -t $chrom $cutesv_path | grep -v -c '^#')
printf "$chrom\t$count_snif\t$count_cute\n" >> results/bcftools/analysis.tsv
done

printf "SNPs count\t $snps_snif \t $snps_cute\n" >> results/bcftools/analysis.tsv
printf "Indels count\t $indels_snif \t $indels_cute\n" >> results/bcftools/analysis.tsv
printf "Deletions count\t $del_snif \t $del_cute\n" >> results/bcftools/analysis.tsv
printf "Insertions count\t $ins_snif \t $ins_cute\n" >> results/bcftools/analysis.tsv
printf "Inversions count\t $inv_snif \t $inv_cute\n" >> results/bcftools/analysis.tsv
printf "Breakends count\t $bnd_snif \t $bnd_cute\n" >> results/bcftools/analysis.tsv
printf "Duplications count\t $dup_snif \t $dup_cute\n" >> results/bcftools/analysis.tsv
printf "DEL/INS ratio\t $indels_ratio_snif \t $indels_ratio_cute\n" >> results/bcftools/analysis.tsv
