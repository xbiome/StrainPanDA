#!/bin/bash

REF=$1
REF_LIST=$2
READS=$3

ADD_ARGS=" -profile docker -resume"

nextflow main.nf $ADD_ARGS --ref_path $REF --path $READS --ref_list $REF_LIST
