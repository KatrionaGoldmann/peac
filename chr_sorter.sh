for c in {1..9};
do
  echo $c
  cp "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr${c}.vcf.gz" "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_temp_chr${c}.vcf.gz"
  echo '  - tidying'
  bcftools view "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_temp_chr${c}.vcf.gz" | awk '{gsub(/^0/,""); print;}' | sed 's/##contig=<ID=0/##contig=<ID=/' > "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr${c}.vcf"
  #mv "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr${c}_temp2.vcf" "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_temp_chr${c}.vcf"
  #gunzip "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_temp_chr${c}.vcf"
  ##mv "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_temp_chr${c}.vcf.gz" "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_temp_chr${c}.vcf"
  echo '   - zipping'
  bgzip -f "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr${c}.vcf"
  tabix -p vcf "/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr${c}.vcf.gz"
done
