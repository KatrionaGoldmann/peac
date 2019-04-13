for c in {1..9};
do
  echo $c
  #cp "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr${c}.vcf.gz" "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr${c}_temp.vcf.gz"
  #bcftools view "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr${c}_temp.vcf.gz" | awk '{gsub(/^0/,""); print;}' > "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr${c}_temp2.vcf.gz"
  mv "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr${c}_temp2.vcf.gz" "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr${c}_temp.vcf.gz"
done
