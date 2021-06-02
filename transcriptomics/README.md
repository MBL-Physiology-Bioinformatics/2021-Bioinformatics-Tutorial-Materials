# MBL 2021 Physiology Course: Transcriptomics

## Welcome to the `Salmon` and `DESeq2` quantification + differential expression module!

To get started with this dataset, first clone the repository, which will give you access to everything you need to complete this part of the course.

```
git clone https://github.com/AlexanderLabWHOI/2021-physiology-subsampling
```

Inside of this repository, you'll find a few relevant subfolders:

```
- 2021-physiology-subsampling
    |- complete-sequence-data
    |- salmon-output
    |- subsampled-data
    |- code
    |- guided-tutorial
```

Inside of the `complete-sequence-data` folder is the `fasta` file that you will use to generate an index with `Salmon`, and inside of the `salmon-output` folder, you'll find everything you need from the result of the final quantification step.

## Introduction to the dataset

![Zebrafish](/images/lax_30860_elife-30860-fig1-v1.jpg)

This is a time series of mRNA sequencing data from zebrafish embryonic development. The really neat thing about this dataset is that, in addition to a number of mRNA transcriptome sequencing files, the authors have taken the time to set their dataset up on an interactive web viewer. You can find that [here](https://www.ebi.ac.uk/gxa/experiments/E-ERAD-475/Results) - _after_ we do some initial exploration of the transcriptome dataset. Among other things, you can look at how the expression of a particular putative gene changes over time, and compare that gene to similarly expressed genes. This is a cool way to be able to visualize particular subsets of the data, once we use some scripts to get a more global feel for the data and to identify some of those important genes ourselves. 

- 5 biological replicates were sequenced
- Each of these 5 biological replicates consisted of a pool of 12 embryos that were sequenced

## Gene identification

The authors have used `Ensembl`, which is what the above viewer shows you, to identify gene models present in the RNA-Seq data for each detailed developmental stage during zebrafish rearing.

## Expression quantification

`Salmon`, a quantification tool for RNA-Seq data, was used to estimate the abundance of the identified genes in the raw data. This tool works faster than direct estimates because it uses _k_-mers, which are _k_ (where _k_ is some variable) sized fragments of the sequences which can be quickly chunked out from the sequence of interest and aligned + summed for quantification. This particular tool has been shown to be highly effective because of the care it takes to eliminate bias related to the sequencing process and the overall distribution of bases in the sequence files. You can read more about that [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/).

Because of the unique 2-step approach used by `Salmon`, where query sequences are first split into these _k_-mers and subsequently used to get the count/abundance data we're after, there are also two steps to actually running `Salmon`. For the differential expression analysis we'll be talking about in this course, we'll use the `Salmon` output provided by the authors of this paper for two developmental stages (128-cell blastula, early in development, and 4 days post-fertilization (larval), towards the end of the experiment). However, it's important to know how you would run `Salmon`, after identifying the gene models.

First, we create an index, which is where the _k_-mers are created. The authors mention that they use the Ensembl 90 dataset for this reference creation, using the flags `--threads 4 --seqBias --gcBias --libType ISR`. 

```
salmon index -t ensembl/ensembldatabase.fa.gz -i ensembl_90 --threads 4
```

Where `ensembl/ensembldatabase.fa.gz` would be a FASTA file containing the reference sequences, `-i ensembl_90` tells us the name of the index we wish to be created, and the other flags have to do with the steps `Salmon` takes to quality control the data and to pre-determine the type of library we're using to quantify the abundance of the transcripts against. The reason that this works is because we're quantifying against the same dataset that "we" (the authors) used to identify genes in the data.

### Practice `Salmon` index generation

In the base folder of this repository, there is a folder on path `complete-sequence-data/Danio_rerio.GRCz10.cds.all.fa.gz` that is the Ensembl transcripts identified from the zebrafish (_Danio rerio_) genome. We'll use this to practice generating an index with `Salmon`. This should be fairly easy to do in an interactive session. From the base of this `GitHub` repository (`2021-physiology-subsampling`), run the following command:

```
salmon index -t complete-sequence-data/Danio_rerio.GRCz10.cds.all.fa.gz -i ensembl_90_zebrafish --threads 4
```

After the command completes, you should find an index with the name `ensembl_90_zebrafish` in the current directory. 

### `Salmon` quantification step

Next, we quantify each of our sets of raw data files against the reference to get the counts + the "transcripts per million" (an abundance metric that is normalized by the length of each sequence and the overall size of the sample). 

```
salmon quant -i ensembl_90 -l ISR \
         -1 sample_file_forward.fastq.gz \
         -2 sample_file_reverse.fastq.gz \
         --seqBias --gcBias --libType ISR \
         -p 8 --validateMappings -o quants/sample_file_quant
```

Where `-i ensembl_90` references the name of the database we created, and the `-1` and `-2` flags show us where to find the forward and reverse read output of our sequencing project, which in our case we would separate out by zebrafish developmental phase.

## Differential expression analysis

In order to work on differential expression, we'll be working out of the sample notebook on path `2021-physiology-subsampling/guided-tutorial/deseq2-tutorial.ipynb`. Go ahead and spin up an R-based `Jupyter` notebook with this file loaded!