#!/bin/bash
if [ $# -ne 1 ]
then
    echo ""
    echo "StrSIM_path must be provided, such as ../StrainPanDA/SimStr/"
   
    exit 1
fi

StrSIM_path=$(readlink -f $1)

#cross group eval for specific datasets
for strn in $(echo "2 4 6 8")
do
mkdir -p eval_results/Fig2A_sJSD/${strn}str_cross_eval

echo -e "strainpanda_eval/strainpanda_${strn}str_1x_pWGS_eval/Ecoli99_WGS_${strn}str_1x_pWGS_StrainPanDA_sorted_JSD_all.csv\tStrainPanDA" > eval_results/Fig2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
echo -e "strainpanda_eval/strainpanda_${strn}str_1x_pWGS_eval/Ecoli99_WGS_${strn}str_1x_pWGS_StrainPanDA_sorted_RA_all.csv\tStrainPanDA" > eval_results/Fig2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt

echo -e "strainest_eval/strainest_${strn}str_1x_pWGS_eval/Ecoli99_WGS_${strn}str_1x_pWGS_EST_sorted_JSD_all.csv\tStrainEst" >> eval_results/Fig2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
echo -e "strainest_eval/strainest_${strn}str_1x_pWGS_eval/Ecoli99_WGS_${strn}str_1x_pWGS_EST_sorted_RA_all.csv\tStrainEst" >> eval_results/Fig2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt

echo -e "pstrain_eval/pstrain_${strn}str_1x_pWGS_eval/Ecoli99_WGS_${strn}str_1x_pWGS_PS_sorted_JSD_all.csv\tPStrain" >> eval_results/Fig2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt
echo -e "pstrain_eval/pstrain_${strn}str_1x_pWGS_eval/Ecoli99_WGS_${strn}str_1x_pWGS_PS_sorted_RA_all.csv\tPStrain" >> eval_results/Fig2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt

Rscript $StrSIM_path/SIM_evaluate_across_methods_sJSD.R sJSD_all_list_f=eval_results/Fig2A_sJSD/${strn}str_cross_eval/sJSD_all_file_list.txt metadata_f=metadata.txt sAb_all_list_f=eval_results/Fig2A_sJSD/${strn}str_cross_eval/sAb_all_file_list.txt out_name=eval_results/Fig2A_sJSD/${strn}str_cross_eval/Fig2A_sJSD_${strn}str
done

#plot Fig2A_sJSD
echo -e "2str_cross_eval/Fig2A_sJSD_2str_sAB_melt.txt\t2str" > eval_results/Fig2A_sJSD/all_abun_list_f.txt
echo -e "4str_cross_eval/Fig2A_sJSD_4str_sAB_melt.txt\t4str" >> eval_results/Fig2A_sJSD/all_abun_list_f.txt
echo -e "6str_cross_eval/Fig2A_sJSD_6str_sAB_melt.txt\t6str" >> eval_results/Fig2A_sJSD/all_abun_list_f.txt
echo -e "8str_cross_eval/Fig2A_sJSD_8str_sAB_melt.txt\t8str" >> eval_results/Fig2A_sJSD/all_abun_list_f.txt

echo -e "bench" > eval_results/Fig2A_sJSD/sample_order.txt
echo -e "StrainPanDA" >>eval_results/Fig2A_sJSD/sample_order.txt
echo -e "StrainEst" >>eval_results/Fig2A_sJSD/sample_order.txt
echo -e "PStrain" >>eval_results/Fig2A_sJSD/sample_order.txt

cd eval_results/Fig2A_sJSD/
Rscript $StrSIM_path/constrain_type_stackplot-Fig2A_sJSD.R all_abun_list_f=all_abun_list_f.txt sample_order_f=sample_order.txt

