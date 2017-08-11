datatype=$1
datapath="/data/vision/polina/projects/patchSynthesis/data"

if [ $datatype = "stroke" ] ; then
     
    for i in 10591 10557 ; do
        dr="${datapath}/stroke/subvols/wholevol/mar12_2016//recons_hierLS_K5_dec06_ds_hier5_2_9_latentSubspaceR/$i/" 
        LSR=`ls $dr | wc -l` 
        dr="${datapath}/stroke/subvols/wholevol/mar12_2016//recons_hierLS_K5_dec06_ds_hier5_2_9_latentMissingR/$i/"
        LMR=`ls $dr | wc -l`
        echo $i LSR: $LSR LMR:$LMR 
     done
else
    while read -u11 i ; do
        dr="${datapath}/ADNI_T1_baselines/subvols/wholevol/mar12_2016//recons_itersLS_K5_d1_3_15_dec12_latentSubspaceR/$i/" 
        LSR=`ls $dr | wc -l` 
        dr="${datapath}/ADNI_T1_baselines/subvols/wholevol/mar12_2016//recons_itersLS_K5_d1_3_15_dec12_latentMissingR/$i/"
        LMR=`ls $dr | wc -l`
        echo $i LSR: $LSR LMR:$LMR 
    done 11< /data/vision/polina/projects/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/selSids_fiveStars_from813_do50.txt #_first50
fi
