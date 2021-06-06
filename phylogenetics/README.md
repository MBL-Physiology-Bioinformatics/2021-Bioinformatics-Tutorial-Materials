# Bioinformatics Bootcamp - Phylogenetics

## Downloading data sets and installing software

### Data sets

#### Course materials

First we will make sure that our [course materials from GitHub](https://github.com/MBL-Physiology-Bioinformatics/2021-Bioinformatics-Tutorial-Materials) are up to date:

```
cd ~/2021-Bioinformatics-Tutorial-Materials
git pull
```

#### EukProt

[EukProt](https://doi.org/10.1101/2020.06.30.180687) is a database of genome-scale predicted proteins from diverse eukaryotes. We will use it in our exercises because we can download the proteins for about 750 species from [a FigShare repository](https://doi.org/10.6084/m9.figshare.12417881.v2) with a single `wget` command, which makes things convenient:

```
mkdir ~/data
cd ~/data
wget https://ndownloader.figshare.com/files/23580944
tar -xzf 23580944 &
```

Note that we've added an `&` following the `tar` command. This will allow us to continue on to other steps, while the proteins are uncompressed in the background. At some point, when the command finishes, you'll get a message sent to your terminal that looks something like this: `[1]+  Done tar -xzf 23580944`.

Next we'll fetch the metadata for the protein sequences, which contains information for each species' data file: its taxonomic lineage, the URL from which it was originally downloaded, the DOI of the paper that produced the data set, and so forth.

```
wget https://ndownloader.figshare.com/files/23580767
mv 23580767 EukProt_included_data_sets.v02.2020_06_30.txt
```

Finally we'll fetch some accessory files for EukProt that we'll use to annotate the phylogenetic trees that we plan to produce. These accessory files are stored in the [EukProt GitHub repository](https://github.com/beaplab/EukProt):

```
git clone https://github.com/beaplab/EukProt.git
```

The files we'll be using are in the `EukProt/iTOL` directory.

### Software

Most of the phylogenetics software we will be using is already installed, either because you created your XSEDE Jetstream instance using the 'Genomics Toolkit', or because you used `conda` to install them during the previous session. For `conda`, you'll just need to activate the environment into which the programs were installed:

```
conda activate bootcampr
```

There are still a few more programs that we will need, so we will install them now into our `~/bin` directory. We choose this directory because, when we try to run a command, the directory will automatically be searched for executable files:

```
mkdir ~/bin
echo $PATH
```

#### trimAl

[trimAl](http://trimal.cgenomics.org/) is a program to remove non-homologous sites from a sequence alignment. It is one of many approaches. We will use it in this tutorial because it is convenient and frequently used.

trimAl is written in a programming language that needs to be compiled from source code. The instance you are using already has the necessary tools to compile, so we only need to download the source and run the compile commands. We will make use of the `/tmp` directory to store the source code. The `/tmp` directory is convenient because its contents are automatically removed whenever the computer is rebooted. That means that after we compile `trimAl` and move the resulting executable file to our `~/bin` directory, we can leave the source code where it was, and it will be cleaned up for us.

```
cd /tmp
wget http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz
tar -xzf trimal.v1.2rev59.tar.gz
rm trimal.v1.2rev59.tar.gz
cd trimAl/source
make
mv trimal ~/bin
```

#### FastTree

[FastTree](http://www.microbesonline.org/fasttree/) is a program to reconstruct phylogenetic trees using an approximately maximum likelihood approach. As the name implies, it is optimized for speed, which makes it useful for this tutorial. It is distributed as a pre-compiled binary for Linux, so we don't even need to compile it:

```
cd ~/bin
wget http://www.microbesonline.org/fasttree/FastTreeMP
```

But wait, there is one more thing that we need to do. Even though the file we just downloaded is an executable, we need to tell our operating system about that. Run the `ls -l` command. You will see that `FastTreeMP` is not executable (it lacks the `x` flags). So let's make it executable:

```
chmod +x FastTreeMP
```

#### Phylogenetic Diversity Analyzer

[Phylogenetic Diversity Analyzer](http://www.cibiv.at/software/pda/) is a program that we will use to select a subset of proteins in a phylogenetic tree that represent the maximal phylogenetic diversity. It is already compiled, as was FastTree, but it comes with a bunch of files we don't need, so we'll download and unpack them in the `/tmp` directory. Here we'll use a trick to get ourselves to the `/tmp` directory that will "remember" the directory that we started from. We will "push" the `/tmp` directory onto a stack, then when we are done, we will "pop" the directory back off of the stack, which will return us to where we started (the `~/bin` directory):

```
pushd /tmp
wget http://www.cibiv.at/software/pda/download/pda-1.0.3/pda-1.0.3-Linux.tar.gz
tar -xzf pda-1.0.3-Linux.tar.gz
mv pda-1.0.3-Linux/bin/pda ~/bin
popd
```

#### InterProScan

[InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/) is a tool to annotate domains in protein sequences using a library of HMMs. We will be using the results from InterProScan in this tutorial, but we won't be downloading the software, because it is a rather large file and it has some dependencies (Python 3, Java 11) that are difficult to install.

Instead, when we need the results, we will find them (pre-computed) in the course GitHub repository.

## Exercises

### Exercise 1: Rad51

The first protein we'll investigate is Rad51. We chose this protein for our first example because it is defined by a [Pfam domain](http://pfam.xfam.org/family/rad51#tabview=tab4) and most Rad51 proteins have a simple domain structure (containing only the Rad51 domain). The first thing we need to do is download the [HMM for the domain](http://pfam.xfam.org/family/rad51#tabview=tab6):

```
mkdir ~/data/Rad51
cd ~/data/Rad51
wget http://pfam.xfam.org/family/PF08423/hmm
mv hmm Rad51.hmm
```

Next we will search the EukProt database with the Rad51 Pfam domain using `hmmsearch`. In order for convenience as part of this tutorial, we wrote a wrapper script to search each individual protein sequence file (there is one file per species) and combine the results. We will use that script here (it should finish the search in about 5 minutes):

```
~/2021-Bioinformatics-Tutorial-Materials/phylogenetics/wrap_hmmsearch.pl -hmm Rad51.hmm -output Rad51.fasta -target ~/data/proteins
```

We used the default options, which will retrieve the best hit (if any) for each species. Take a look at the output file `Rad51.fasta` to see the results.

Next, we need to align the sequences to determine which sites in each species are homologous. For this, we will use `hmmalign`, part of the HMMer package, which will use the Rad51 to create an alignment. The benefit of using `hmmalign` is that it runs very quickly and will very accurately determine which sites are homologous (because they match to the HMM), but the drawback is that it will not align sites outside the HMM. In this instance, it is not a big deal, since the Rad51 domain covers most of the Rad51 protein in most eukaryotic species. Note that we will use the `--trim` option, which removes sites outside the HMM from the alignment.

```
hmmalign --outformat afa --trim Rad51.hmm Rad51.fasta > Rad51.hmmalign.fasta
```

Next we will reconstruct a phylogenetic tree using FastTree approximate maximum likelihood:

```
FastTreeMP Rad51.hmmalign.fasta > Rad51.hmmalign.FastTree.tree
```

Finally, we will use a helper script to produce an annotation file that we will use when visualizing the tree with iTOL. In this case, the annotation file will color each species by its taxonomic lineage:

```
~/2021-Bioinformatics-Tutorial-Materials/phylogenetics/fasta_iTOL_taxonomy_annotation.pl -fasta Rad51.fasta -taxonomy ~/data/EukProt_included_data_sets.v02.2020_06_30.txt -colors ~/data/EukProt/iTOL/clade_colors.v02.2020_03_27.txt -output Rad51.iTOL_taxonomy_annotation.txt
```

Now, we need to transfer the tree file `Rad51.hmmalign.FastTree.tree` and the annotation file `Rad51.iTOL_taxonomy_annotation.txt` to our local computer.

To visualize the tree, let's go to [iTOL](https://itol.embl.de/), click on the "Upload" button, and browse to where we saved the tree on our local computer.

By default, the tree displays in a circular format. To add the color annotations by taxonomic lineage, we simply drag and drop the annotation file onto the tree.

### Exercise 2: Caspases and metacaspases

Let's start at the [InterPro web site entry for the Caspase-like domain superfamily](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR029030/). If we click on "Proteins", we can see that there are 56,000, which is clearly too many for our tree. Conveniently, we can filter these proteins by clicking on "Reviewed", to get the subset that have been expert reviewed. Let's click "Export" and then "FASTA" to download these to our local computer. Then we will transfer them into a new directory on our remote instance called `~/data/caspase`.

Instead of aligning with `hmmalign`, as we did in the first exercise, we will use MAFFT, since we do not have an HMM available, and because we want to align the full protein sequences:
```
mafft protein-matching-IPR029030.fasta > protein-matching-IPR029030.mafft.fasta
```

With the alignment, we can now make our own HMM that we will use to search the EukProt database. Since we have built the HMM ourselves, it won't have the gathering threshold that comes from Pfam. So we'll need to specify an E-value. Let's choose `1e-20` (this should take about 7-8 minutes):

```
hmmbuild -n IPR029030 IPR029030.hmm protein-matching-IPR029030.mafft.fasta
~/2021-Bioinformatics-Tutorial-Materials/phylogenetics/wrap_hmmsearch.pl -hmm IPR029030.hmm -output IPR029030.fasta -target ~/data/proteins -evalue 1e-20
```

Next let's align the sequences that we found with MAFFT. Although we have an HMM available, we want to align the full sequences, just in case there are some interesting parts of the proteins we found that lie outside the HMM that we built:
```
mafft IPR029030.fasta > IPR029030.mafft.fasta
```

And then we can trim the alignment with trimAl to remove sites where less than 20% of species are present:

```
trimal -gt 0.2 -in IPR029030.mafft.fasta -out IPR029030.mafft.trimal_gt0.2.fasta
```

With the trimmed alignment, we can then build a tree:

```
FastTreeMP IPR029030.mafft.trimal_gt0.2.fasta > IPR029030.mafft.trimal_gt0.2.tree
```

Finally, we run the helper script to produce colors by taxonomic lineage:

```
~/2021-Bioinformatics-Tutorial-Materials/phylogenetics/fasta_iTOL_taxonomy_annotation.pl -fasta IPR029030.fasta -taxonomy ~/data/EukProt_included_data_sets.v02.2020_06_30.txt -colors ~/data/EukProt/iTOL/clade_colors.v02.2020_03_27.txt -output IPR029030.iTOL_taxonomy_annotation.txt
```

In the course GitHub repository, there is also a Pfam domain annotation file available, which we pre-computed for you. The annotation file is formatted for import into iTOL. You can download this to your local machine:

```
~/2021-Bioinformatics-Tutorial-Materials/phylogenetics/IPR029030.iTOL_Pfam_annotation.txt
```

Visualizing the tree in iTOL is a lot of information all at once, and lots of the sequences seem pretty similar to one another. Let's try to reduce the number of sequences to a more managable number, by choosing a subset of a fixed size with maximal phylogenetic diversity with the Phylogenetic Diversity Analyzer tool. We can choose a subset of whichever size we like. Here let's choose 100:

```
pda -k 100 IPR029030.mafft.trimal_gt0.2.tree IPR029030.mafft.trimal_gt0.2.pda_k100
```

This prints various output information, but what we really want is the tree. So let's use the `grep` command to extract the tree and redirect it into a file. Notice that the tree line begins with an open parenthesis, followed by an EP number. So let's grep for that:

```
grep "^(EP" IPR029030.mafft.trimal_gt0.2.pda_k100 > IPR029030.mafft.trimal_gt0.2.pda_k100.tree
```

Now we can transfer this file to our local computer and open it in iTOL. We can use the same annotation files as for the full tree. iTOL will produce error messages that most of the annotated sequences are missing, but we can dismiss these messages, since this is what we expect.

### Exercise 3: RFX transcription factors

This is a group exercise. The idea is to test the hypothesis presented in this paper from 2010 that RFX transcription factors originated in the last common ancestor of Amoebozoa and Opisthokonta:

![RFX Evolution](https://www.pnas.org/content/pnas/107/29/12969/F1.large.jpg)

You should start with the [RFX DNA binding Pfam domain](http://pfam.xfam.org/family/PF02257).

Break into small groups, each of which may make different choices in the tree reconstruction process. Once we each have trees, we can compare and contrast how they are different depending on the choices that you made. We expect that you will encounter problems (this is the essence of bioinformatics!), but we will be here to help you.

Here are the steps you will want to follow:

1. Start with the RFX DNA-binding domain HMM
1. Search EukProt (how many hits? 1? 3? 5?)
1. Align (with hmmalign? with MAFFT?)
1. Trim with trimAl? Using which settings?
1. Build a tree with FastTree? Or RAxML (which requires a model of sequence evolution)?
1. Prune the tree with PDA? Or not?
1. Annotate the tree with InterProScan and by taxonomy
1. Visualize the tree with iTOL

The InterProScan annotation file (for the top 5 hits from every species) is available here:

```
~/2021-Bioinformatics-Tutorial-Materials/phylogenetics/RFX_DNA_binding.iTOL_Pfam_annotation.txt
```

For the groups that choose to use RAxML, here is a command line to get you started:

```
raxmlHPC-PTHREADS-AVX2 -m PROTGAMMALGF -s <alignment FASTA file> -T 6 -n RFX_DNA_binding -p <random number> -x <random number> -f a -N 100
```
