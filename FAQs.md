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
