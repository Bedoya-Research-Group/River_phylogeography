#!/bash/bin

for K in {1..10}

do

admixture -C 0.0001 --cv=10 marathrum_chiriqui_filtered_final_bin_plink.bed $K | tee marathrum_chiriqui_filtered_final_bin_plink.${K}.out
grep -h CV marathrum_chiriqui_filtered_final_bin_plink*.out > marathrum_chiriqui_filtered_final_bin_plink.bed_chooseK.txt
done
