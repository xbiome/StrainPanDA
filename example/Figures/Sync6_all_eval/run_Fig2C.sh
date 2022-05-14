#!/bin/bash
if [ $# -ne 1 ]
then
    echo ""
    echo "StrSIM_path must be provided, such as ../StrainPanDA/SimStr/"
   
    exit 1
fi


mkdir -p eval_results/Fig2C_sJSD/

StrSIM_path=$1

rm -rf eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
rm -rf eval_results/Fig2C_sJSD/sAb_all_file_list.txt

for MT in $(echo "Bifidobacterium-longum Clostridium-difficile Escherichia-coli Enterococcus-faecalis Faecalibacterium-prausnitzii Prevotella-copri")
do
echo -e "strainpanda_eval/${MT}_eval/${MT}_StrainPanDA_sorted_RA_all.csv\t${MT}-StrainPanDA" >> eval_results/Fig2C_sJSD/sAb_all_file_list.txt
echo -e "strainpanda_eval/${MT}_eval/${MT}_StrainPanDA_sorted_JSD_all.csv\t${MT}-StrainPanDA" >> eval_results/Fig2C_sJSD/sJSD_all_file_list.txt

echo -e "strainest_eval/${MT}_eval/${MT}_EST_sorted_JSD_all.csv\t${MT}-StrainEST" >> eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
echo -e "strainest_eval/${MT}_eval/${MT}_EST_sorted_RA_all.csv\t${MT}-StrainEST" >> eval_results/Fig2C_sJSD/sAb_all_file_list.txt

echo -e "pstrain_eval/${MT}_eval/${MT}_PS_sorted_JSD_all.csv\t${MT}-PStrain" >> eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
echo -e "pstrain_eval/${MT}_eval/${MT}_PS_sorted_RA_all.csv\t${MT}-PStrain" >> eval_results/Fig2C_sJSD/sAb_all_file_list.txt

done

sed -i "s/Bifidobacterium-longum-/Bifidobacterium_longum-/g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
sed -i "s/Clostridium-difficile-/Clostridium_difficile-/g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
sed -i "s/Escherichia-coli-/Escherichia_coli-/g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
sed -i "s/Enterococcus-faecalis-/Enterococcus_faecalis-/g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
sed -i "s/Faecalibacterium-prausnitzii-/Faecalibacterium_prausnitzii-/g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
sed -i "s/Prevotella-copri-/Prevotella_copri-/g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt

sed -i "s/Bifidobacterium-longum-/Bifidobacterium_longum-/g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt
sed -i "s/Clostridium-difficile-/Clostridium_difficile-/g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt
sed -i "s/Escherichia-coli-/Escherichia_coli-/g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt
sed -i "s/Enterococcus-faecalis-/Enterococcus_faecalis-/g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt
sed -i "s/Faecalibacterium-prausnitzii-/Faecalibacterium_prausnitzii-/g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt
sed -i "s/Prevotella-copri-/Prevotella_copri-/g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt

#for Clostridium-difficile Faecalibacterium-prausnitzii Prevotella-copri, EST results were using StrainPanDA results for generating figure in order to avoid error. But in fact they should be NA.
sed -i "s|strainest_eval/Clostridium-difficile_eval/Clostridium-difficile_EST_sorted|strainpanda_eval/Clostridium-difficile_eval/Clostridium-difficile_StrainPanDA_sorted|g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
sed -i "s|strainest_eval/Faecalibacterium-prausnitzii_eval/Faecalibacterium-prausnitzii_EST_sorted|strainpanda_eval/Faecalibacterium-prausnitzii_eval/Faecalibacterium-prausnitzii_StrainPanDA_sorted|g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt
sed -i "s|strainest_eval/Prevotella-copri_eval/Prevotella-copri_EST_sorted|strainpanda_eval/Prevotella-copri_eval/Prevotella-copri_StrainPanDA_sorted|g" eval_results/Fig2C_sJSD/sJSD_all_file_list.txt

sed -i "s|strainest_eval/Clostridium-difficile_eval/Clostridium-difficile_EST_sorted|strainpanda_eval/Clostridium-difficile_eval/Clostridium-difficile_StrainPanDA_sorted|g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt
sed -i "s|strainest_eval/Faecalibacterium-prausnitzii_eval/Faecalibacterium-prausnitzii_EST_sorted|strainpanda_eval/Faecalibacterium-prausnitzii_eval/Faecalibacterium-prausnitzii_StrainPanDA_sorted|g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt
sed -i "s|strainest_eval/Prevotella-copri_eval/Prevotella-copri_EST_sorted|strainpanda_eval/Prevotella-copri_eval/Prevotella-copri_StrainPanDA_sorted|g" eval_results/Fig2C_sJSD/sAb_all_file_list.txt

Rscript $StrSIM_path/SIM_evaluate_across_methods_sJSD.R sJSD_all_list_f=eval_results/Fig2C_sJSD/sJSD_all_file_list.txt metadata_f=metadata.txt sAb_all_list_f=eval_results/Fig2C_sJSD/sAb_all_file_list.txt out_name=eval_results/Fig2C_sJSD/Fig2C_sJSD_4str
Rscript ${StrSIM_path}/stat_across_group.R in_matr_f=eval_results/Fig2C_sJSD/Fig2C_sJSD_4str_sJSD_matr.txt

cd eval_results/Fig2C_sJSD/

echo -e "Fig2C_sJSD_4str_sJSD_matr.txt\t4str" > JSD_all_file_list.txt

Rscript $StrSIM_path/all_group_MCC_paired_Fig2B.R MCC_all_list_f=JSD_all_file_list.txt dist_type=sJSD