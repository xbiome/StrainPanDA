#!/bin/bash
if [ $# -ne 1 ]
then
    echo ""
    echo "StrSIM_path must be provided, such as ../StrainPanDA/SimStr/"
   
    exit 1
fi

StrSIM_path=$(readlink -f $1)

mkdir -p eval_results/FigS2A_sJSD/

rm eval_results/FigS2A_sJSD/significant_pT_merge.txt

#cross group eval for specific datasets
for strn in $(echo "2 4 6 8")
do
    mkdir eval_results/FigS2A_sJSD/${strn}str_cross_eval
    rm -rf eval_results/FigS2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
    rm -rf eval_results/FigS2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt
    
    for LV in $(echo "ErrFree ARTError 1x_pWGS")
    do 
        echo $LV
        echo -e "strainpanda_eval/strainpanda_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_StrainPanDA_sorted_JSD_all.csv\t${LV}-StrainPanDA" >> eval_results/FigS2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
        echo -e "strainpanda_eval/strainpanda_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_StrainPanDA_sorted_RA_all.csv\t${LV}-StrainPanDA" >> eval_results/FigS2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt
        
        echo -e "strainest_eval/strainest_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_EST_sorted_JSD_all.csv\t${LV}-StrainEst" >> eval_results/FigS2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
        echo -e "strainest_eval/strainest_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_EST_sorted_RA_all.csv\t${LV}-StrainEst" >> eval_results/FigS2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt

        echo -e "pstrain_eval/pstrain_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_PS_sorted_JSD_detected.csv\t${LV}-PStrain" >> eval_results/FigS2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
        echo -e "pstrain_eval/pstrain_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_PS_sorted_RA_all.csv\t${LV}-PStrain" >> eval_results/FigS2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt
        
    done

    sed -i "s/1x_pWGS-/WGSError-/g" eval_results/FigS2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
    sed -i "s/ErrFree-/AErrorFree-/g" eval_results/FigS2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
    
    sed -i "s/1x_pWGS-/WGSError-/g" eval_results/FigS2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt
    sed -i "s/ErrFree-/AErrorFree-/g" eval_results/FigS2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt
    
    Rscript ${StrSIM_path}/SIM_evaluate_across_methods_sJSD.R sJSD_all_list_f=eval_results/FigS2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt metadata_f=metadata.txt sAb_all_list_f=eval_results/FigS2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt out_name=eval_results/FigS2A_sJSD/${strn}str_cross_eval/FigS2_${strn}str

Rscript ${StrSIM_path}/stat_across_group.R in_matr_f=eval_results/FigS2A_sJSD/${strn}str_cross_eval/FigS2_${strn}str_sJSD_matr.txt
echo ${strn}str >> eval_results/FigS2A_sJSD/significant_pT_merge.txt
cat eval_results/FigS2A_sJSD/${strn}str_cross_eval/FigS2_${strn}str_sJSD_matr.txt_significantPair.txt >> eval_results/FigS2A_sJSD/significant_pT_merge.txt

done


#plot FigS2A_sJSD seqerr
cd eval_results/FigS2A_sJSD/
echo -e "2str_cross_eval//FigS2_2str_sJSD_matr.txt\t2str" > JSD_all_file_list.txt
echo -e "4str_cross_eval//FigS2_4str_sJSD_matr.txt\t4str" >> JSD_all_file_list.txt
echo -e "6str_cross_eval//FigS2_6str_sJSD_matr.txt\t6str" >> JSD_all_file_list.txt
echo -e "8str_cross_eval//FigS2_8str_sJSD_matr.txt\t8str" >> JSD_all_file_list.txt

Rscript ${StrSIM_path}/all_group_MCC_paired_Fig2B.R MCC_all_list_f=JSD_all_file_list.txt dist_type=sJSD
