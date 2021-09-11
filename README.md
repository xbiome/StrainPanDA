# StrainPanDA -- A strain analysis pipeline based on pangenome

StrainPanDA is a tool that deconvolute pangenome coverage into strain composition and strain gene profile.

## Dependencies

### Using the strainpanda container (recommended)
1. docker

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

Install required packages `dplyr`, `nnlm`, `foreach`, `MASS`

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
