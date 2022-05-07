# StrainPanDA -- A strain analysis pipeline based on pangenome

StrainPanDA is a tool that deconvolute pangenome coverage into strain composition and strain gene profile.

## Installation

### Local installation
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

2. install PanPhlAn (for mapping short reads to pan-genome databases)

https://github.com/SegataLab/panphlan/wiki#download-the-panphlan-software

3. install strainpandar

Install required packages `dplyr`, `foreach`, `MASS`, `NMF`

```sh
tar -czf strainpandar.tar.gz src/strainpandar
R CMD INSTALL strainpandar.tar.gz
```

### Using the strainpanda container (recommended)

We provide two docker files to build the two images required by the nextflow pipeline.

```
docker build -t strainpanda-mapping:dev . -f docker/Dockerfile_mapping
docker build -t strainpanda-strainpandar:dev . -f docker/Dockerfile_strainpandar
```

### Pre-built pangenome databases

TODO: url for downloading

## Run analysis

### Run full analysis
Download test data to a foler e.g. `reads`

```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/005/SRR5813295/SRR5813295_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/005/SRR5813295/SRR5813295_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/006/SRR5813296/SRR5813296_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/006/SRR5813296/SRR5813296_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/007/SRR5813297/SRR5813297_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/007/SRR5813297/SRR5813297_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/008/SRR5813298/SRR5813298_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/008/SRR5813298/SRR5813298_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/009/SRR5813299/SRR5813299_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/009/SRR5813299/SRR5813299_2.fastq.gz
```

Create a species list

```sh
echo "Faecalibacterium prausnitzii" > species_list.txt
```

#### With docker

```sh
PATH_TO_REPO/main.nf -profile docker --ref_path PATH_TO_REFERENCE --path reads/ --ref_list species_list.txt
```

#### Local

```sh
PATH_TO_REPO/main.nf --ref_path PATH_TO_REFERENCE --path reads/ --ref_list species_list.txt
```

### Run strainpandar


## Outputs

Ouput files from the above run:

```sh
strainpanda_out/
├── Faecalibacterium-prausnitzii-202009.counts.csv  ## Merged count matrix (gene family by sample)
├── Faecalibacterium-prausnitzii-202009_mapping     ## Sample specific count files
│   ├── SRR5813295_Faecalibacterium-prausnitzii-202009.csv.bz2
│   ├── SRR5813296_Faecalibacterium-prausnitzii-202009.csv.bz2
│   ├── SRR5813297_Faecalibacterium-prausnitzii-202009.csv.bz2
│   ├── SRR5813298_Faecalibacterium-prausnitzii-202009.csv.bz2
│   └── SRR5813299_Faecalibacterium-prausnitzii-202009.csv.bz2
├── Faecalibacterium-prausnitzii-202009_strainpandar_out ## strainpandar outputs
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.anno_strain_sample.pdf ## annotation to the closest reference
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.genefamily_strain.csv ## gene family-strain matrix
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.genefamily_strain.pdf ## heatmap visualization
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.rds ## R object contains strainpandar results
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.strain_sample.csv ## strain-sample matrix
│   └── Faecalibacterium-prausnitzii-202009.strainpanda.strain_sample.pdf ## barplot visualization
└── pipeline_info  ## pipeline run statistics
    ├── strainpanda_DAG.svg
    ├── strainpanda_report.html
    ├── strainpanda_timeline.html
    └── strainpanda_trace.txt
```