#!/bin/bash
if [ $# -ne 1 ]
then
    echo ""
    echo "StrSIM_path must be provided, such as ../StrainPanDA/SimStr/"
   
    exit 1
fi

StrSIM_path=$(readlink -f $1)

mkdir -p eval_results/FigS4CR_MCC/

rm eval_results/FigS4CR_MCC/significant_pT_merge.txt

#cross group eval for specific datasets
for strn in $(echo "2 4 6 8")
do
    mkdir eval_results/FigS4CR_MCC/${strn}str_cross_eval
    rm -rf eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
    rm -rf eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
    rm -rf eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
    rm -rf eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
    rm -rf eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
    rm -rf eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt
    
    for LV in $(echo "1x_MI 5x_MI 10x_MI 25x_MI 100x_MI 1x_pWGS")
    do 
        echo $LV
        echo -e "strainpanda_eval/strainpanda_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_StrainPanDA_stat_all.csv\t${LV}-StrainPanDA" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
        echo -e "strainpanda_eval/strainpanda_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_StrainPanDA_dis_all.csv\t${LV}-StrainPanDA" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
        echo -e "strainpanda_eval/strainpanda_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_StrainPanDA_rank_abun_diff_all.csv\t${LV}-StrainPanDA" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
        echo -e "strainpanda_eval/strainpanda_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_StrainPanDA_rank_abun_all.csv\t${LV}-StrainPanDA" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
        echo -e "strainpanda_eval/strainpanda_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_StrainPanDA_abun_all.csv\t${LV}-StrainPanDA" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
        echo -e "strainpanda_eval/strainpanda_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_StrainPanDA_abun_diff_all.csv\t${LV}-StrainPanDA" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt
        
        echo -e "strainest_eval/strainest_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_EST_stat_all.csv\t${LV}-StrainEst" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
        echo -e "strainest_eval/strainest_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_EST_dis_all.csv\t${LV}-StrainEst" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
        echo -e "strainest_eval/strainest_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_EST_rank_abun_diff_all.csv\t${LV}-StrainEst" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
        echo -e "strainest_eval/strainest_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_EST_rank_abun_all.csv\t${LV}-StrainEst" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
        echo -e "strainest_eval/strainest_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_EST_abun_all.csv\t${LV}-StrainEst" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
        echo -e "strainest_eval/strainest_${strn}str_${LV}_eval/Ecoli99_WGS_${strn}str_${LV}_EST_abun_diff_all.csv\t${LV}-StrainEst" >> eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt
        
    done

    sed -i "s/1x_pWGS-/BG0-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
    sed -i "s/1x_pWGS-/BG0-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
    sed -i "s/1x_pWGS-/BG0-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
    sed -i "s/1x_pWGS-/BG0-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
    sed -i "s/1x_pWGS-/BG0-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
    sed -i "s/1x_pWGS-/BG0-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt
    
    sed -i "s/1x_MI-/BG01F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
    sed -i "s/1x_MI-/BG01F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
    sed -i "s/1x_MI-/BG01F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
    sed -i "s/1x_MI-/BG01F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
    sed -i "s/1x_MI-/BG01F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
    sed -i "s/1x_MI-/BG01F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt

    sed -i "s/25x_MI-/BG25F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
    sed -i "s/25x_MI-/BG25F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
    sed -i "s/25x_MI-/BG25F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
    sed -i "s/25x_MI-/BG25F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
    sed -i "s/25x_MI-/BG25F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
    sed -i "s/25x_MI-/BG25F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt

    sed -i "s/5x_MI-/BG05F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
    sed -i "s/5x_MI-/BG05F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
    sed -i "s/5x_MI-/BG05F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
    sed -i "s/5x_MI-/BG05F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
    sed -i "s/5x_MI-/BG05F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
    sed -i "s/5x_MI-/BG05F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt

    sed -i "s/10x_MI-/BG10F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
    sed -i "s/10x_MI-/BG10F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
    sed -i "s/10x_MI-/BG10F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
    sed -i "s/10x_MI-/BG10F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
    sed -i "s/10x_MI-/BG10F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
    sed -i "s/10x_MI-/BG10F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt

    sed -i "s/100x_MI-/BGx100F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt
    sed -i "s/100x_MI-/BGx100F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt
    sed -i "s/100x_MI-/BGx100F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt
    sed -i "s/100x_MI-/BGx100F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt
    sed -i "s/100x_MI-/BGx100F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt
    sed -i "s/100x_MI-/BGx100F-/g" eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt
    
    Rscript ${StrSIM_path}/SIM_evaluate_across_methods.R stat_all_list_f=eval_results/FigS4CR_MCC/${strn}str_cross_eval/stat_all_file_list.txt metadata_f=metadata.txt dis_all_list_f=eval_results/FigS4CR_MCC/${strn}str_cross_eval/dis_all_file_list.txt RAD_all_list_f=eval_results/FigS4CR_MCC/${strn}str_cross_eval/rank_abun_diff_all_file_list.txt RA_all_list_f=eval_results/FigS4CR_MCC/${strn}str_cross_eval/prof_all_file_list.txt Ab_all_list_f=eval_results/FigS4CR_MCC/${strn}str_cross_eval/Abun_all_file_list.txt AD_all_list_f=eval_results/FigS4CR_MCC/${strn}str_cross_eval/abun_diff_all_file_list.txt out_name=eval_results/FigS4CR_MCC/${strn}str_cross_eval/FigS4CR_MCC_${strn}str

Rscript ${StrSIM_path}/stat_across_group.R in_matr_f=eval_results/FigS4CR_MCC/${strn}str_cross_eval/FigS4CR_MCC_${strn}str_MCC_matr.txt
echo ${strn}str >> eval_results/FigS4CR_MCC/significant_pT_merge.txt
cat eval_results/FigS4CR_MCC/${strn}str_cross_eval/FigS4CR_MCC_${strn}str_MCC_matr.txt_significantPair.txt >> eval_results/FigS4CR_MCC/significant_pT_merge.txt


done


#plot Fig2C_sJSD-s1_seq_err_3M
cd eval_results/FigS4CR_MCC/

rm MCC_all_file_list.txt
for StrN in $(echo "2str 4str 6str 8str")
do
    echo -e "${StrN}_cross_eval//FigS4CR_MCC_${StrN}_MCC_matr.txt\t${StrN}" >> MCC_all_file_list.txt

done

Rscript ${StrSIM_path}/all_group_MCC_paired_Fig2C_2M.R MCC_all_list_f=MCC_all_file_list.txt dist_type=MCC
