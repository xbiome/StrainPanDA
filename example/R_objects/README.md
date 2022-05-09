# Examples in the StrainPanDA manuscript

## rds files of synthetic datasets were submitted to Zenodo

The Zenodo link is: .

## The list of rds files

### The rds filename of the synthetic data of *E. coli* strains is in the following format:

panphlan_${strn}str_${type}_Escherichia-coli-202009.rds

"panphlan" showing it is using the coverage matrix from PanPhlAn

${strn} is the number of strains. For example: 2, 4, 6 and 8.

${type} strands for the simulation type: 

​	ErrFree and ARTError are the error-free simulation and the simulation with sequencing errors.

​	02x_pWGS 05x_pWGS 1x_pWGS 2x_pWGS 5x_pWGS are pWGS data at 0.2x, 0.5x,1x, 2x, and 5x sequencing depth separately.

​	1x_IBD 5x_IBD 10x_IBD 25x_IBD 100x_IBD are WGSBG dataset having 1, 5, 10, 25 and 100 fold of IBD backgrounds with the pWGS data at 1x sequencing depth

​	1x_FMT 5x_FMT 10x_FMT 25x_FMT 100x_FMT are WGSBG dataset having 1, 5, 10, 25 and 100 fold of FMT backgrounds with the pWGS data at 1x sequencing depth

​	1x_ERR 5x_ERR 10x_ERR 25x_ERR 100x_ERR are WGSBG dataset having 1, 5, 10, 25 and 100 fold of MI backgrounds with the pWGS data at 1x sequencing depth

​	

### The rds file of O104

### The rds file of real dataset

