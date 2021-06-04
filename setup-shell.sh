#!/bin/bash

conda create -y -n bootcampr -c r -c bioconda -c conda-forge r-irkernel hmmer r-biocmanager salmon

conda init bash
source ~/.bashrc
conda activate bootcampr

Rscript R-installs.R > R-script-output.txt 2>&1
