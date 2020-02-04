# Katriona kgoldmann
#
# This script fixed the format of chromosome numbers to be consistent between
# files (reference and PEAC). Switches between 1 and 01 etc. 

#cd /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_reworked/imputed_vcf_kg

# run as sudo bash -c "./chr_corrector.sh"

for i in {1..9}
do
  echo "Chromosome $i - fixing the chroms"
  gunzip -k -f /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz
  sed -n 's/^[a-z]*0*//p ' /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr$i.vcf > /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf
  bgzip -c  /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf > /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz
  rm /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf
  rm /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr$i.vcf
  tabix -p vcf /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz
done

for i in {10..22}
do
  echo "Chromosome $i - just copying"
  cp /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz
  cp /home/kgoldmann/Documents/gcpeac/Kevin/PEAC-eQTL/ManchesterGenotyping/PEAC/imputed_vcf_kg/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz.tbi /home/kgoldmann/tmp/PEAC_1000G_Phase3_Oct14_chr$i.vcf.gz.tbi
done
