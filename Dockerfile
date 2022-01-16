############################################################
# Dockerfile to build strainpanda pipeline
# Based on ubuntu
# To build:
# docker build -t strainpanda .
############################################################

FROM conda/miniconda3

################# BEGIN INSTALLATION ######################
## install system pacakges
RUN apt-get update && apt-get -y install procps

## install panphlan
RUN conda update -n base -c defaults conda
RUN conda install samtools=0.1.19 -c bioconda
RUN conda install -y panphlan -c bioconda

## install r dependencies
RUN conda install r-base -c r
RUN conda install r-dplyr r-MASS r-foreach r-NNLM -c r
RUN conda install r-getopt -c r
RUN conda install r-pracma -c r
RUN conda install r-pheatmap r-reshape2 r-ggplot2 -c r
RUN conda install r-permute r-cluster -c r
RUN curl -O https://cran.r-project.org/src/contrib/Archive/vegan/vegan_2.5-6.tar.gz 
RUN R CMD INSTALL vegan_2.5-6.tar.gz
RUN conda install r-data.table r-r.utils -c r
RUN conda install r-nmf -c r

## install StrainPanDAr
COPY ./src/strainpandar	 /strainpandar
RUN tar -czf strainpandar.tar.gz strainpandar
RUN R CMD INSTALL strainpandar.tar.gz

## install minpath
COPY ./src/MinPath /MinPath
ENV MinPath /MinPath

RUN rm -r strainpandar && rm -r strainpandar.tar.gz && rm -r vegan_2.5-6.tar.gz

##################### INSTALLATION END #####################

