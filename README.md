# StrainPanDA -- A strain analysis pipeline based on pangenome

TODO: Intoduction/Abstract

## Dependencies

### Using the strainpanda container (recommended)
1. docker
1. **or** singularity (TODO: this should work with the docker container)

### Local installation
1. [Nextflow](https://www.nextflow.io/): full pipeline
1. [PanPhlAn](https://github.com/segatalab/panphlan): mapping to pangenome databases
1. R and [strainpandar](src/strainpandar): decomposing pangenome into strains

### Pre-built pangenome databases

TODO: url for downloading

## Installation

### Docker

```sh
docker build -t strainpanda:dev .
```

### Local installation

1. install nextflow

```sh
curl -s https://get.nextflow.io | bash
```

2. install PanPhlAn

https://github.com/SegataLab/panphlan/wiki#download-the-panphlan-software

3. install strainpandar

TODO: required packages `dplyr`, `nnlm`, `foreach`, `MASS`

```sh
# Setup ssh connection auth files (private repository).
creds = git2r::cred_ssh_key(publickey="~/.ssh/gitlab_key.pub",
    privatekey="~/.ssh/gitlab_key")

# Install 
devtools::install_git("git@gitlab.com:xbiome/system-development/xviz.git",
    #branch="develop",
    credentials=creds)
```

```sh
tar -czf strainpandar.tar.gz src/strainpandar
R CMD INSTALL strainpandar.tar.gz
```

## Run analysis

Download test data
```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/007/SRR5813297/SRR5813297_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/007/SRR5813297/SRR5813297_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/008/SRR5813298/SRR5813298_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/008/SRR5813298/SRR5813298_2.fastq.gz
```

### With docker

```sh
PATH_TO_REPO/main.nf -profile docker --ref_path PATH_TO_REFERENCE --path PATH_TO_READS --ref_list SPECIES_FILE 
```

### Local

```sh
PATH_TO_REPO/main.nf --ref_path PATH_TO_REFERENCE --path PATH_TO_READS --ref_list SPECIES_FILE
```
