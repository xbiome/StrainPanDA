# StrainPanDA - A strain analysis pipeline based on pangenome
# FAQs (Frequently asked questions)


## Errors pulling docker images

### Reason: the authorization of docker is not set correctly to the user.
![image](docker_pull_error.png)

Solution1: sudo docker pull

```
sudo docker pull yuxiangtan/strainpanda-mapping:dev
sudo docker tag yuxiangtan/strainpanda-mapping:dev strainpanda-mapping:dev
```

Solution2: add docker group for users by root. See [here](https://docs.docker.com/engine/install/linux-postinstall/) for details.

## Failing to run PanPhlAN

### Main reason: the name of reference database is not correct.

The name of referece database must following the correct rule.

To double check, you can `ls -l <the path of error log>`

For example as the following:
![image](panphlan_error.png)


## Setup the number of processors.

You can set it by changing the cpus value in the conf/base.config file, following [this](https://www.nextflow.io/docs/latest/process.html#cpus )

## Get the functional annotation of gene families 

You can get the annotation of gene families in the prebuild databases by using the "anno.csv" files. In the real practice, you can use "merge" function in R to get the annotation of the gene family profile.

## How to analysis your own strains (strains used in your own experiments that are not in the database) using StrainPanDA?
First, you can compare the pangenome of the references in the database with your own strain, and calculate AJI (Average Jacard Index). If AJI >=99, your self-owned strain will be merged with the existing strains in the database during the reconstruction process, therefore database reconstruction is meanningless. If AJI < 99, you may consider rebuilding the database. Or without regeneration, you can compare the predicted gene profile of strains with you own strain to see their similarity. 

## Why some of the input sample was lost in the final result?
StrainPanDA has two filter steps: first, if the sample has insufficient coverage of the species, the sample will be filtered out. The second filtering rule is as following: Samples were filtered out if the number of gene families detected was below 0.9 × gmin (gmin is the minimum number of gene families found in all reference genomes). Therefore, the .counts.csv table will generally have all samples (but may be filtered if the sample does not have the species), and some samples may not be present in the strain composition result.

## In StrainPanDA, how to calculate the lower threshold of relative abundance of target species according to the amount of sequencing data?
According to the discussion in the paper, for a typical metagenomic data (6GB), assuming that the average genome size of the species is 6MB, then the target species should have at least 1% relative abundance, that is, 60MB data volume. This amount of data corresponds to a sequencing depth of approximately 10x, where the lower limit of expected detection of the strain is 10%.

The corresponding formula is: The relative abundance threshold of species Lspecies (%) = 100%/ (Sdata / Sgenome / (100%/Lstrain (%)) ), where the metagenomic data is Sdata (1GB=1000MB), the average genome size of strains is Sgenome (MB), and the expected lower limit of strain detection is Lstrain (%).Users can adjust the variables of this fomular to obtain the target value.
