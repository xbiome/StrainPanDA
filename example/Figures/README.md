# List of codes

## NOTE！All the shell scripts or R codes (in htmls), should be run under the SimStr environment



## Supporting evaluation files:

Simulated datasets:

 - [Ecoli_sim_all_eval](Ecoli_sim_all_eval): evaluation files, codes and results from synthetic data of *E. coli* strains
 - [Sync6_all_eval](Sync6_all_eval): evaluations from synthetic data of 6 gut bacterial species

Under each folder:

 - `eval_results`: results for different figures in seperated folders
 - `simulation_matrix`: simulation matrix used to generate basic simulation
 - `strainpanda_eval`: primary evaluation result of StrainPanDA outputs
 - `strainest_eval`: primary evaluation result of StrainEst outputs
 - `pstrain_eval`: primary evaluation result of PStrain outputs

Example to generate strainpanda evaluation

```sh
time Rscript ../../../SimStr/SIM_evaluate_unifrac_denovo.R \
  bench_tb=simulation_matrix/benchmark_comp_2str_all_str_prof.csv \
  method_anno_tb=../Ecoli_sim_all_csv/strainpanda_csv/panphlan_2str_1x_pWGS_Escherichia-coli-202009.xstrain_str_anno_prof.csv \
  denovo_ref=StrainPanDA_ref_str_list.txt tree_nwk=Ecoli99_parsnp.tree \
  name_out=strainpanda_eval/strainpanda_2str_1x_pWGS_eval/Ecoli99_WGS_2str_1x_pWGS_StrainPanDA
```

Example to generate strainest evaluation

```sh
time Rscript ../../../SimStr/SIM_evaluate_unifrac_denovo.R \
  bench_tb=simulation_matrix/benchmark_comp_2str_all_str_prof.csv \
  method_anno_tb=../Ecoli_sim_all_csv/strainest_tsv/strainest_2str_1x_pWGS.tsv \
  denovo_ref=StrainEst_ref_str_list.txt \
  tree_nwk=Ecoli99_parsnp.tree \
  name_out=strainest_eval/strainest_2str_1x_pWGS_eval/Ecoli99_WGS_2str_1x_pWGS_EST
```

Example to generate pstrain evaluation

```sh
time Rscript ../../../SimStr/SIM_evaluate_unifrac_denovo.R \
  bench_tb=simulation_matrix/benchmark_comp_2str_all_str_prof.csv \
  method_anno_tb=../Ecoli_sim_all_csv/pstrain_csv/pstrain_2str_1x_pWGS.csv denovo_ref=StrainEst_ref_str_list.txt \
  tree_nwk=Ecoli99_parsnp.tree \
  name_out=pstrain_eval/pstrain_2str_1x_pWGS_eval/Ecoli99_WGS_2str_1x_pWGS_PS \
  str_match=FALSE
```


## Main figures：

### Figure 2:

Fig2A: 
 - Working directory: [Ecoli_sim_all_eval](Ecoli_sim_all_eval)
 - Command: `sh run_Fig2A.sh ../../../SimStr/ 2>run_Fig2A.log`
 - Results: [eval_results/Fig2A_sJSD](Ecoli_sim_all_eval/eval_results/Fig2A_sJSD)
 - Main output: [All_groups_abun_stackplot_by_method_ordered.pdf](Ecoli_sim_all_eval/eval_results/Fig2A_sJSD/All_groups_abun_stackplot_by_method_ordered.pdf)

Fig2B:
 - Working directory: [Ecoli_sim_all_eval](Ecoli_sim_all_eval)
 - Command: `sh run_Fig2B.sh ../../../SimStr/ 2>run_Fig2B.log`
 - Results: [eval_results/Fig2B_sJSD](Ecoli_sim_all_eval/eval_results/Fig2B_sJSD)
 - Main output: [All_groups_sJSD_boxplot.pdf](Ecoli_sim_all_eval/eval_results/Fig2B_sJSD/All_groups_sJSD_boxplot.pdf), [All_groups_sJSD_y01_boxplot.pdf](Ecoli_sim_all_eval/eval_results/Fig2B_sJSD/All_groups_sJSD_y01_boxplot.pdf) and [All_groups_sJSD_yauto_boxplot.pdf](Ecoli_sim_all_eval/eval_results/Fig2B_sJSD/All_groups_sJSD_yauto_boxplot.pdf)

Fig2C:
 - Working directory: [Sync6_all_eval](Sync6_all_eval)
 - Command: `sh run_Fig2C.sh ../../../SimStr/ 2>run_Fig2C.log`
 - Results: [eval_results/Fig2C_sJSD](Sync6_all_eval/eval_results/Fig2C_sJSD)
 - Main output: [All_groups_sJSD_boxplot.pdf](Sync6_all_eval/eval_results/Fig2C_sJSD/All_groups_sJSD_boxplot.pdf), [All_groups_sJSD_y01_boxplot.pdf](Sync6_all_eval/eval_results/Fig2C_sJSD/All_groups_sJSD_y01_boxplot.pdf) and [All_groups_sJSD_yauto_boxplot.pdf](Sync6_all_eval/eval_results/Fig2C_sJSD/All_groups_sJSD_yauto_boxplot.pdf)
 - Statistical tests: [Fig2C_sJSD_4str_sJSD_matr.txt_significantPair.txt](Sync6_all_eval/eval_results/Fig2C_sJSD/Fig2C_sJSD_4str_sJSD_matr.txt_significantPair.txt)

Fig2D: 
 - Working directory: [Ecoli_sim_all_eval](Ecoli_sim_all_eval)
 - Command: refer to [run_Fig2D_FigS5_FigS6.html](https://github.com/xbiome/StrainPanDA/blob/main/example/Figures/Ecoli_sim_all_eval/run_Fig2D_FigS5_FigS6.html)
 - Results: [eval_results/Fig2D_heatmap](Ecoli_sim_all_eval/eval_results/Fig2D_heatmap)
 - Main output: [Fig2D_strain_nb_SIM_allgene_pangenome_heatmap.pdf](Ecoli_sim_all_eval/eval_results/Fig2D_heatmap/Fig2D_strain_nb_SIM_allgene_pangenome_heatmap.pdf)

Fig2E: 
 - Working directory: [Sync6_all_eval](Sync6_all_eval)
 - Command: refer to [run_Fig2E.html](https://htmlpreview.github.io/?https://github.com/xbiome/StrainPanDA/blob/main/example/Figures/Sync6_all_eval/run_Fig2E.html)
 - Results:[eval_results/Fig2E_AUPRC](Sync6_all_eval/eval_results/Fig2E_AUPRC)
 - Main output: [Fig2E_AUC_wRC_y01_boxplot.pdf](Sync6_all_eval/eval_results/Fig2E_AUPRC/Fig2E_AUC_wRC_y01_boxplot.pdf)
 - Statistical tests: [Fig2E_AUC_wRC_stat_p.csv](Sync6_all_eval/eval_results/Fig2E_AUPRC/Fig2E_AUC_wRC_stat_p.csv)
 - Inputs: [eval_results/Fig2E_AUPRC/input_files_6species](Sync6_all_eval/eval_results/Fig2E_AUPRC/input_files_6species)
 - Supporting files to generate input files: [eval_results/Fig2E_AUPRC/supporting_files_6species](Sync6_all_eval/eval_results/Fig2E_AUPRC/supporting_files_6species), see [run_Fig2E_support_generate.html](https://htmlpreview.github.io/?https://github.com/xbiome/StrainPanDA/blob/main/example/Figures/Sync6_all_eval/run_Fig2E_support_generate.html) for an example of Bifidobacterium-longum.

### Figure 3:

FigA: 
 - Working directory: [MI_study](MI_study)
 - Command: refer to [Fig3AB_S10_S11_mother_infantt.html](https://htmlpreview.github.io/?https://github.com/xbiome/StrainPanDA/blob/main/example/Figures/MI_study/Fig3AB_S10_S11_mother_infant.html)
 - Results: [Bifidobacterium-longum_out/simplex_plot](MI_study/Bifidobacterium-longum_out/simplex_plot)
   - The `*_only.pdf` were used in Fig3A (M: mother, B: newborn, 4M: 4 month, 12M: 12 months)

FigA: 
 - Working directory: [MI_study](MI_study)
 - Command: refer to [Fig3AB_S10_S11_mother_infantt.html](https://htmlpreview.github.io/?https://github.com/xbiome/StrainPanDA/blob/main/example/Figures/MI_study/Fig3AB_S10_S11_mother_infant.html)
 - Results: [MI_study/Bifidobacterium-longum_out/simplex_plot](MI_study/Bifidobacterium-longum_out/simplex_plot)
   - The `*_only.pdf` were used in Fig3A (M: mother, B: newborn, 4M: 4 month, 12M: 12 months)

FigB: 
 - Working directory: [MI_study](MI_study)
 - Command: refer to [Fig3AB_S10_S11_mother_infantt.html](https://htmlpreview.github.io/?https://github.com/xbiome/StrainPanDA/blob/main/example/Figures/MI_study/Fig3AB_S10_S11_mother_infant.html)
 - Results: [MI_study/Bifidobacterium-longum_out/strain_pair_dist_breastStop](Bifidobacterium-longum_out/strain_pair_dist_breastStop)
   - `BreastStop_both_strain1_StopvsCont_yfull_boxplot.pdf`, `BreastStop_both_strain2_StopvsCont_yfull_boxplot.pdf`, `BreastStop_both_strain3_StopvsCont_yfull_boxplot.pdf` were used in Fig3B.

FigC:
 - Working directory: [MI_study](MI_study)
 - Command: refer to [Fig3AB_S10_S11_mother_infantt.html](https://htmlpreview.github.io/?https://github.com/xbiome/StrainPanDA/blob/main/example/Figures/MI_study/run_Fig3C_S12_S13_S14_heatmap.html)
 - Results: [MI_study/Bifidobacterium-longum_out/Fig3C_heatmap](Bifidobacterium-longum_out/Fig3C_heatmap)

### Figure 4：
 - Working directory: [kong_et_al_IBD_FMT](kong_et_al_IBD_FMT)
 - Command: refer to [analysis.Rmd](https://github.com/xbiome/StrainPanDA/blob/main/example/Figures/kong_et_al_IBD_FMT/analysis.Rmd)

## Supplementary Figures：

### Supplementary Figure 2

Supplementary Figure 2A: in Ecoli_sim_all_eval, use the following command: sh run_FigS2A.sh ../../../SimStr/ 2>run_FigS2A.log

​	Results are in eval_results/FigS2A_JSD

​	The main output is: All_groups_sJSD_boxplot.pdf; All_groups_sJSD_y01_boxplot.pdf and All_groups_sJSD_yauto_boxplot.pdf showed different way to set the y-axis range.

​	The statistical significant values are listed in significant_pT_merge.txt 

Supplementary Figure 2B: in Ecoli_sim_all_eval, use the following command: sh run_FigS2B.sh ../../../SimStr/ 2>run_FigS2B.log

​	Results are in eval_results/FigS2B_MCC

​	The main output is: All_groups_MCC_boxplot.pdf; All_groups_MCC_y01_boxplot.pdf and All_groups_MCC_yauto_boxplot.pdf showed different way to set the y-axis range.

​	The statistical significant values are listed in significant_pT_merge.txt 

### Supplementary Figure 3:

Supplementary Figure 3A: in Ecoli_sim_all_eval, use the following command: sh run_FigS3A.sh ../../../SimStr/ 2>run_FigS3A.log

​	Results in eval_results/FigS3A_JSD

​	The main output is: All_groups_sJSD_boxplot.pdf; All_groups_sJSD_y01_boxplot.pdf and All_groups_sJSD_yauto_boxplot.pdf showed different way to set the y-axis range.

​	The statistical significant values are listed in significant_pT_merge.txt 

Supplementary Figure 3B: in Ecoli_sim_all_eval, use the following command: sh run_FigS3B.sh ../../../SimStr/ 2>run_FigS3B.log

​	Results are in eval_results/FigS3B_MCC

​	The main output is: All_groups_MCC_boxplot.pdf; All_groups_MCC_y01_boxplot.pdf and All_groups_MCC_yauto_boxplot.pdf showed different way to set the y-axis range.

​	The statistical significant values are listed in significant_pT_merge.txt 

### Supplementary Figure 4: 	

Left panels: in Ecoli_sim_all_eval, use the following command:sh run_FigS4AL.sh ../../../SimStr/ 2>run_FigS4AL.log ; sh run_FigS4BL.sh ../../../SimStr/ 2>run_FigS4BL.log ; sh run_FigS4CL.sh ../../../SimStr/ 2>run_FigS4CL.log 

​	Results in eval_results/FigS4AL_JSD , FigS4BL_JSD and FigS4CL_JSD

​	The main output is: All_groups_sJSD_boxplot.pdf; All_groups_sJSD_y01_boxplot.pdf and All_groups_sJSD_yauto_boxplot.pdf showed different way to set the y-axis range.

​	The statistical significant values are listed in significant_pT_merge.txt 

Right panels: in Ecoli_sim_all_eval, use the following command:sh run_FigS4AR.sh ../../../SimStr/ 2>run_FigS4AR.log ; sh run_FigS4BR.sh ../../../SimStr/ 2>run_FigS4BR.log ; sh run_FigS4CR.sh ../../../SimStr/ 2>run_FigS4CR.log 

​	Results  are in eval_results/FigS4AR_MCC , FigS4BR_MCC , FigS4CR_MCC

​	The main output is: All_groups_MCC_boxplot.pdf; All_groups_MCC_y01_boxplot.pdf and All_groups_MCC_yauto_boxplot.pdf showed different way to set the y-axis range.

​	The statistical significant values are listed in significant_pT_merge.txt 

### Supplementary Figure 5:

In Ecoli_sim_all_eval, see run_Fig2D_FigS5_FigS6.html for detail.

​	Results are in eval_results/FigS5_boxplot

​	All other dependes for 2, 6 and 8 strains were provided in the folder as well. panphlan_single_errfree_reference_gene_presence_absence.tsv is universal for all 8 strains.

### Supplementary Figure 6: 

In Ecoli_sim_all_eval, see run_Fig2D_FigS5_FigS6.html for detail.

​	Result is in eval_results/FigS6_PRplot/

### Supplementary Figure 7:

In Ecoli_sim_all_eval, see FigS7_O104_genes_evaluate_StrainPanDA.html for detail.

​	Result is in eval_results/FigS7_heatmap/

​	The main output is: StrainPanDA_O104_functional_gene_pangenome_heatmap.pdf

### Supplementary Figure 8:

FigS8 AB: In Ecoli_sim_all_eval, see FigS8_strainpanda_vs_panphlan_geneprofile.html for details
	Results and the folder of inputs are in eval_results/FigS8AB_panphlan/,

​	The main outputs are FigS8_strain_nb_SIM_vs_panphlan_dom_heatmap.pdf and XT_Pan_dom_paired_JD_boxplot.pdf	

FigS8 CD: In Ecoli_sim_all_eval, see FigS8_strainpanda_vs_panphlan3_geneprofile.html for details

​	Results and the folder of inputs are in eval_results/FigS8CD_panphlan3/

​	The main outputs are FigS8C_strain_to_pan3_nb_SIM_vs_pan3dom_heatmap.pdf and XT_Pan3_dom_paired_JD_boxplot.pdf

### Supplementary Figure 9:

In MI_study folder, see run_FigS9.html for details, 

​	Results are in Bifidobacterium-longum_out/FigS9

### Supplementary Figure 10:

In MI_study folder, see Fig3AB_S10_S11_mother_infant.html for details, in the generate Fig3B and FigS10 section

​	Results are in Bifidobacterium-longum_out/strain_pair_dist_breastStop

​	BreastStop_4M_strain1_BFFvsOthers_boxplot.pdf, BreastStop_4M_strain2_BFFvsOthers_boxplot.pdf, BreastStop_12M_strain1_BBFvsBBM_boxplot.pdf and BreastStop_12M_strain2_BBFvsBBM_boxplot.pdf were used in FigS10.

### Supplementary Figure 11:

In MI_study folder, see Fig3AB_S10_S11_mother_infant.html for details, in the generate FigS11 section

​	Results are in Bifidobacterium-longum_out/str_relative_abun_breastStop

​	Only the results of subspecies1 (strain1) were used in the FigS11

### Supplementary Figure 12:

In MI_study folder, see run_Fig3C_S12_S13_S14_heatmap.html for details, in the generate FigS12 section

​	Results are in Bifidobacterium-longum_out/FigS12_heatmap

### Supplementary Figure 13:

In MI_study folder, see run_Fig3C_S12_S13_S14_heatmap.html for details, in the generate FigS13 section

​	Results are in Bifidobacterium-longum_out/FigS13_heatmap

### Supplementary Figure 14:

In MI_study folder, see run_Fig3C_S12_S13_S14_heatmap.html for details, in the generate FigS14 section

​	Results are in Bifidobacterium-longum_out/FigS14_heatmap

### Supplementary Figure 15:

In kong_et_al_IBD_FMT folder, see analysis.html for details,

​	Inputs from StrainPanDA are in strain_results folder

### Supplementary Figure 16:

In kong_et_al_IBD_FMT folder, see analysis.html for details,

​	Inputs from StrainPanDA are in strain_results folder

### Supplementary Figure 17:

In kong_et_al_IBD_FMT folder, see analysis.html for details,

​	Inputs from StrainPanDA are in strain_results folder

### Supplementary Figure 18:

In sample_size_effect_eval folder, see run_FigS18_sample_size.html for details.

### Supplementary Figure 19:

In running_time_eval folder, see run_FigS19_running_time.html for details.
