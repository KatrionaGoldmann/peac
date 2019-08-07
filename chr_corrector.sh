#cd /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_reworked/imputed_vcf_kg

# run as sudo bash -c "./chr_corrector.sh"

for i in {1..2}
do
  echo "Welcome $i times"
  gunzip -k -f /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_reworked/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz
  sed -n 's/^[a-z]*0*//p ' /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_reworked/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr$i.vcf > /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf
  gzip /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf
  tabix /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz
done
