#!/bin/bash
sniffles_path='/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/myelong/results_concatenated_samples/sniffles/patient_sniffles.vcf'
cutesv_path='/storage/scratch01/groups/bu/myeloma_altum/agostoli/board_alegosto01/myelong/results_concatenated_samples/cutesv/patient_cutesv.vcf'

printf "INFO \t SNIFFLES \t CUTE\n" > "${snakemake_output[0]}"

total_variants_snif=$(bcftools query -f '%CHROM\n' "${snakemake_input[sniffles]}" | grep -v -c '^#')
total_variants_cute=$(bcftools query -f '%CHROM\n' "${snakemake_input[cutesv]}" | grep -v -c '^#')

chromlist=($(bcftools query -f '%CHROM\n' "${snakemake_input[sniffles]}"| uniq|grep 'chr'))

snps_snif=$(bcftools view -v snps "${snakemake_input[sniffles]}" |grep -v -c '^#')
snps_cute=$(bcftools view -v snps "${snakemake_input[cutesv]}" |grep -v -c '^#')

indels_snif=$(bcftools view -v indels "${snakemake_input[sniffles]}" |grep -v -c '^#')
indels_cute=$(bcftools view -v indels "${snakemake_input[cutesv]}" |grep -v -c '^#')

bnd_snif=$(bcftools view -v bnd "${snakemake_input[sniffles]}" |grep -v -c '^#')
bnd_cute=$(bcftools view -v bnd "${snakemake_input[cutesv]}" |grep -v -c '^#')

del_snif=$(bcftools query -f '%ID\n' "${snakemake_input[sniffles]}" |grep 'Sniffles2.DEL'| grep -v -c '^#')
del_cute=$(bcftools query -f '%ID\n' "${snakemake_input[cutesv]}" |grep 'cuteSV.DEL'| grep -v -c '^#')

ins_snif=$(bcftools query -f '%ID\n' "${snakemake_input[sniffles]}" |grep 'Sniffles2.INS'| grep -v -c '^#')
ins_cute=$(bcftools query -f '%ID\n' "${snakemake_input[cutesv]}" |grep 'cuteSV.INS'| grep -v -c '^#')

inv_snif=$(bcftools query -f '%ID\n' "${snakemake_input[sniffles]}" |grep 'Sniffles2.INV'| grep -v -c '^#')
inv_cute=$(bcftools query -f '%ID\n' "${snakemake_input[cutesv]}" |grep 'cuteSV.INV'| grep -v -c '^#')

dup_snif=$(bcftools query -f '%ID\n' "${snakemake_input[sniffles]}" |grep 'Sniffles2.DUP'| grep -v -c '^#')
dup_cute=$(bcftools query -f '%ID\n' "${snakemake_input[cutesv]}" |grep 'cuteSV.DUP'| grep -v -c '^#')

del_snif_float=$del_snif.0
del_cute_float=$del_cute.0

ins_snif_float=$ins_snif.0
ins_cute_float=$ins_cute.0


indels_ratio_snif=$(awk "BEGIN {print $del_snif_float/$ins_snif_float}")
indels_ratio_cute=$(awk "BEGIN {print $del_cute_float/$ins_cute_float}")
#indels_ratio_cute=$(($del_cute_float / $ins_cute_float))


printf "total variants \t $total_variants_snif \t $total_variants_cute \n" >> "${snakemake_output[0]}"

for chrom in ${chromlist[@]}
do
count_snif=$(bcftools view -t $chrom "${snakemake_input[sniffles]}" | grep -v -c '^#')
count_cute=$(bcftools view -t $chrom "${snakemake_input[cutesv]}" | grep -v -c '^#')
printf "$chrom\t$count_snif\t$count_cute\n" >> "${snakemake_output[0]}"
done

printf "SNPs count\t $snps_snif \t $snps_cute\n" >> "${snakemake_output[0]}"
printf "Indels count\t $indels_snif \t $indels_cute\n" >> "${snakemake_output[0]}"
printf "Deletions count\t $del_snif \t $del_cute\n" >> "${snakemake_output[0]}"
printf "Insertions count\t $ins_snif \t $ins_cute\n" >> "${snakemake_output[0]}"
printf "Inversions count\t $inv_snif \t $inv_cute\n" >> "${snakemake_output[0]}"
printf "Breakends count\t $bnd_snif \t $bnd_cute\n" >> "${snakemake_output[0]}"
printf "Duplications count\t $dup_snif \t $dup_cute\n" >> "${snakemake_output[0]}"
printf "DEL/INS ratio\t $indels_ratio_snif \t $indels_ratio_cute\n" >> "${snakemake_output[0]}"
