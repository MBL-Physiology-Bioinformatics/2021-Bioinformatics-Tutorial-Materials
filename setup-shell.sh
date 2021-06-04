#!/bin/bash

conda create -y -n bootcampr -c r -c bioconda -c conda-forge r-irkernel hmmer r-biocmanager salmon

conda activate bootcampr

Rscript R-installs.R