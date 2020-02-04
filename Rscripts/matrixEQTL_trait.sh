for i in {4..4}
   do
     echo "PEAC_matqtl_trait.R Synovium $i" #check call to r script is as expected
     Rscript PEAC_matqtl_trait.R "Synovium" $i
   done


for i in {4..4}
    do
      echo "PEAC_matqtl_trait.R Blood $i" #check call to r script is as expected
      Rscript PEAC_matqtl_trait.R "Blood" $i
    done
