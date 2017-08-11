while read -u11 line ; 
do 
    printf "$line \n\n" ; 
    ./sgeSubvolReconFlexible.sh ADNI_T1_baselines mar12_2016 $line wholevol itersLS_K5_d1_3_15_dec12 latentSubspaceR &
    printf "\n\n\n" 
done 11< /data/vision/polina/projects/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/selSids_fiveStars_from813_do50.txt 
