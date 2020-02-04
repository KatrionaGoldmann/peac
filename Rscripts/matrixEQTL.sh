for i in {1..5}
   do
     echo "PEAC_matqtl.R $i Synovium" #check call to r script is as expected
     Rscript PEAC_matqtl.R $i "Synovium"
   done


for i in {1..4}
    do
      echo "PEAC_matqtl.R $i Blood" #check call to r script is as expected
      Rscript PEAC_matqtl.R $i "Blood"
    done
