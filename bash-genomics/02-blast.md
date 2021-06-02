

## **BLAST**
The Basic Local Alignment Search Tool (BLAST) finds regions of local similarity between sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance of matches. BLAST can be used to infer functional and evolutionary relationships between sequences as well as help identify members of gene families.

#### BLAST Types

The programs in the BLAST+ suite can search for and against sequences in protein format  and in nucleotide format (A’s, C’s, T’s, and G’s). Depending on what type the query and subject sets are, different BLAST programs are used.
![Alt text](images/blast_types.png)

Two sequences (either nucleotide or protein) may be compared directly, but when comparing a nucleotide sequence to a protein sequence, we need to consider which reading frame of the nucleotide sequence corresponds to a protein. The **blastx** and **tblastn** programs do this by converting nucleotide sequences into protein sequences in all six reading frames and comparing against all of them. Generally such programs result in six times as much work to be done. The tblastx program compares nucleotide queries against nucleotide subjects, but it does so in protein space with all six conversions compared to all six on both sides.

Other  BLAST tools include psiblast, which produces an initial search and tweaks the scoring rules on the basis of the results; these tweaked scoring rules are used in a second search, generally finding even more matches. This process is repeated as many times as the user wishes, with more dissimilar matches being revealed in later iterations. The deltablast program considers a precomputed database of scoring rules for different types of commonly found (conserved) sequences. Finally, rpsblast searches for sequence matches against sets of profiles, each representing a collection of sequences.


----------
#### Running web-based blast searches


Go to [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
https://blast.ncbi.nlm.nih.gov/Blast.cgi

>A summer student recovered putative nitrogenase genes from marine metagenomic assemblies (Tara Oceans). Using web-based blast, can you help her/him to evaluate the results?
>
>\>Putative_nifH_seq1
VLGDVVCGGFSVPLRGGGREVVIIISGELMAIYAGNNILRAVKAIGGNPKILGYVVNMRG
VPNETSIVEAFATKTGVPILASIKRDPVTFKEAERLGGPIVRTLPDSEIADTFRQMAKRV
QLGAGDTVPRPIDSYDELFDMFMSFQDLEGETGKEGQASYLVKIPNDPVKRDKPIRISVY
GAGGIGKSTTSSNVSAALVLMGEKVFQIGCDPKHDSIANLCGGLKPTILDEIQKRGGARA
MTREILKELIFPGMNHNNHLYGSECGGPVPGKGCGGKGVDLAIKLLEDHGIMDEYQFSLI
LYDVLGDTVCGGFARPLKYTPQAYIIANGETSSLTQAMKIAQSIQVASSRGADVGLAGVI
NNMRGVPNEEAIVEEAFETVGIPVIHHIPRSILVQTAENLRTTVVEAFLESAQAKQYLEL
AEKMQQNTQRHALTREVLSSREIVAIVNKHA
\>Putative_nifH_seq2
MAKEIGIYGKGGIGKSTTVSNVSAATSELGFKVMQIGCDPKADSTNSLLGGKYIPTILDT
VLEADTVREYTEVDVSKVIFEGYNGIVCAECGGPDPGIGCAGRGVITAIELMKEQGAFDS
INPDFIFYDVLGDVVCGGFAMPLRQGICQQVYVIVSPNFASMYAANNLFKAIVRFGERGG
GQLAGLIANHINTPEEKDIVDDFAAITKTEVVEYIPHSMELARADLQGRTILESSPKSEL
AKTYRGLAAEIAKREQQIIPVPAEAPALRKWAQGWLERLHASALKE


#### Stand-alone blast
For large datasets, using the web-based blast to get functional (and/or taxonomic annotation) is not option. But, it can be installed and run locally.

We will use the same conda environment we created earlier *downloading* to install blast. Make sure that *downloading* is active and install blast:
```conda install -c bioconda blast```

Each of the various programs in the BLAST suite accepts a large number of options; run blastn -help to see them for the blastn program. Here is a summary of a few parameters that are most commonly used for blastn et al.:

> **-query**
	the name (or path) of the FASTA-formatted file to search for as query sequences.
>
>**-subject**
	the name (or path) of the FASTA-formatted file to search in as subject sequences.
>
> **-db** <database name>
    the name of the database to search against (as opposed to using -subject).
>
>**-evalue**
	only hits with E values smaller than this should be reported. For example: -evalue 0.001 or -evalue 1e-6 *evaluation of the evalue cut-off should be performed*
>
> **-num_threads** <integer>
    number of CPU cores to use.
>
>**-outfmt**
>	output format: The default, 0, provides a human-readable (but not programmatically parseable) text file. The values 6 and 7 produce tab-separated rows and columns in a text file, with 7 providing explanatory comment lines. Similarly, a value of 10 produces comma-separated output; 11 produces a format that can later be quickly turned into any other with another program called blast_formatter. Options 6, 7, and 10 can be highly configured in terms of what columns are shown.
>
>**-max_target_seqs**
	see Shah et al, 2018: "Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows"
>
>**-out**
Write the output to as opposed to the default of standard output.

#### BLAST Databases
Additionally, blast provides a tool called makeblastdb that converts a subject FASTA file into an indexed and quickly searchable (but not human-readable) version of the same information, stored in a set of similarly named files (often at least three ending in .pin, .psq, and .phr for protein sequences, and .nin, .nsq, and .nhr for nucleotide sequences). This set of files represents the “database,” and the database name is the shared file name prefix of these files.

We will use
```makeblastdb -in *fasta file* -out  *database name* -dbtype *type* -title *title*```  where *type* is one of prot or nucl, and *title* is a human-readable title

When using the -db option, the BLAST tools will search for the database files in three locations: (1) the present working directory, (2) your home directory, and (3) the paths specified in the $BLASTDB environment variable. If you created your database in other location, you have to provide the full path to this location.

>The file ""/proteins/For_database/dataset.faa"  contains protein sequences for all genes retrieved from several marine single cell genomes. Make a blast formatted database and scan this database for our three genes of interest (in the file proteins/query.faa). Before you begin make sure that you use SLURM (srun for this example). These kind of search are computationally heavy and cannot be handled by the login nodes.
For your queries use
evalue: 1e-10
number of threads: 2
and as outfmt: "6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore"

```srun -p scavenger --time=00:30:00 --ntasks-per-node=2 --pty bash```

> Do not forget to activate the conda environment where blast is installed

```makeblastdb -in dataset.faa -title dataset -dbtype prot -out dataset```


```blastp -query query.faa -db dataset -evalue 1e-10 -num_threads 4 -outfmt "6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out query_blastp.txt```

>Check the output file. What do you see?

## Other tools
There are a lot of other tools that someone can use to search sequences of interests against other sequences that are designed to fill in various "searching" needs such as HMMER (designed to detect remote homologs as sensitively as possible), LAST (handles big sequencing data i.e. can compare two vertebrate genomes), etc. We do not have time to cover them all but we will be happy to discuss them with you.



For more information on BLAST, see: Altschul et al. 1990 https://www.sciencedirect.com/science/article/pii/S0022283605803602

For this section material from Open Oregon State - Open Textbooks Sites (http://library.open.oregonstate.edu/computationalbiology/chapter/command-line-blast/) was used

