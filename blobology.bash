#!/bin/bash

## Blobology pipeline. 2013-09-12 Sujai Kumar 
## sujai.kumar@gmail.com github.com/sujaikumar/assemblage

## Needs:
## Read files (unless you are providing assembly as fasta and
## read alignments as BAM).
## ABySS in path (unless you are providing your own assembly)
## NCBI blast+ suite (2.2.28 or above, which provide taxids in hit)
## NCBI formatted nt blast database (because it includes taxon information)
## NCBI taxonomy dump (to calculate higher tax levels from hit taxids)
## samtools in path (to manipulate bam files)
## github assemblage scripts in path
## fastq-mcf from the ea-utils suite (version 1.1.2-537) if you want to
## quality- and adapter- trim your raw reads

####==========================================================================
#### Download read files from SRA
####==========================================================================

## Uncomment the lines below if you are replicating the results in the
## Blobology paper. Otherwise provide your own reads/assembly in the assembly/
## alignment steps.
## ERR138445 is a 300 bp paired-end library
## ERR138446 is a 600 bp paired-end library

# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR138/ERR138445/ERR138445_1.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR138/ERR138445/ERR138445_2.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR138/ERR138446/ERR138446_1.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR138/ERR138446/ERR138446_2.fastq.gz

## Uncomment if you want to do quality- and adapter- trimming using fastq-mcf

#fastq-mcf -o >(gzip -c >ERR138445_1.mcf.fastq.gz) -o >(gzip -c >ERR138445_2.mcf.fastq.gz) --max-ns 0 -l 50 -q 20 --qual-mean 20 -R adapters.fa \
#    ERR138445_1.fastq.gz ERR138445_2.fastq.gz &>ERR138445.mcf.err
#fastq-mcf -o >(gzip -c >ERR138446_1.mcf.fastq.gz) -o >(gzip -c >ERR138446_2.mcf.fastq.gz) --max-ns 0 -l 50 -q 20 --qual-mean 20 -R adapters.fa \
#    ERR138446_1.fastq.gz ERR138446_2.fastq.gz &>ERR138446.mcf.err

####==========================================================================
#### GENERAL CONFIG VARIABLES
####==========================================================================

## Num of processors. Or, uncomment the line below and enter manually.

NUMPROC=`grep "^processor" /proc/cpuinfo | tail -n 1 | awk '{print $3}'`
#NUMPROC=16

## K-mer size for prelim assembly. Not needed if you are providing your own
## assembly fasta file.
## You can crudely estimate this for your data using
## http://dna.med.monash.edu.au/~torsten/velvet_advisor/

KMER=61

## Location of local installation of nt blast database
## (not needed if using blast remotely, which is slower)

## The NCBI nt databases can be downloaded from
## ftp://ftp.ncbi.nlm.nih.gov/blast/db/ using the following command:

# wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"

BLASTDB=~/scratch/blastdb

## Location of NCBI tar gunzipped directory downloaded from
## ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
## Default: current directory

TAXDUMP=.

####==========================================================================
#### STEP 1, run ABySS assembler, or provide assembly fasta file
####==========================================================================

## If providing own assembly fasta file, uncomment and insert filename,
## and comment out the remaining commands in this section:

#ASSEMBLY=assembly.fasta

## If using ABySS, uncomment the commands below and insert your ownread files
## in the section 'lib="..."'.
## See http://www.bcgsc.ca/downloads/abyss/doc/#assemblingapaired-endlibrary
## for details. We also recommend running the
## command below on a cluster if possible (with np= the number of cores
## that you can request and j= number of processors on a single machine.
## e.g, if you have access to 6 nodes each with 8 processors, use np=48 j=8
## Even if you have only one machine with many cores, compile ABySS with
## mpi support so that it can run many processes in parallel.

abyss-pe name=pa${KMER} k=${KMER} np=60 j=12 lib="lib1 lib2" lib1="ERR138445_1.mcf.fastq.gz ERR138445_2.mcf.fastq.gz" lib2="ERR138446_1.mcf.fastq.gz ERR138446_2.mcf.fastq.gz"

## At the end of the abyss step, the assembled sequence will be in
## ${NAME}-scaffolds.fa so set the ASSEMBLY environment variable to that filename

ASSEMBLY=${NAME}-scaffolds.fa

####==========================================================================
#### STEP 2, find best blast hits of a random sample of contigs
####==========================================================================

RND_ASSEMBLY=minlen1k.random10k.$ASSEMBLY
fastaqual_select.pl -f $ASSEMBLY -s r -n 10000 -l 1000 >$RND_ASSEMBLY

## Change $BLASTDB to match your own path location of the nt blast database
## If running blast remotely on NCBI's servers, use -remote -db nt
## instead of -db $BLASTDB/nt
## -outfmt 6 qseqid staxids gives only a two-column output table with the
## query id and the taxid of the best hit

blastn -task megablast -query $RND_ASSEMBLY -db $BLASTDB/nt -evalue 1e-5 -num_threads $NUMPROC -max_target_seqs 1 -outfmt '6 qseqid staxids' -out $RND_ASSEMBLY.nt.1e-5.megablast

####==========================================================================
#### STEP 3, map reads back to assembly using bowtie2
####==========================================================================

## This step is not needed if you already have BAM files with the read alignments

bowtie2-build -o 3 -t 12 $ASSEMBLY $ASSEMBLY
bowtie2 -x $ASSEMBLY --very-fast-local -k 1 -t -p 12 -U ERR138445_1.mcf.fastq.gz,ERR138445_2.mcf.fastq.gz  | samtools view -S -b -T $ASSEMBLY - >$ASSEMBLY.ERR138445.bowtie2.bam
bowtie2 -x $ASSEMBLY --very-fast-local -k 1 -t -p 12 -U ERR138446_1.mcf.fastq.gz,ERR138446_2.mcf.fastq.gz  | samtools view -S -b -T $ASSEMBLY - >$ASSEMBLY.ERR138446.bowtie2.bam

####==========================================================================
#### STEP 4, combine BAM files and blast hits to create data file for plotting
####==========================================================================

## ABySS will create one or more BAM files. Replace *.bam with own BAM files
## if you have aligned them separately without using ABySS

gc_cov_annotate.pl --blasttaxid $RND_ASSEMBLY.nt.1e-5.megablast --assembly $ASSEMBLY --bam *.bam --out blobplot.txt

####==========================================================================
#### STEP 4, create blobplots using R, or visualise using blobsplorer
####==========================================================================

## Tested with R 2.15.2 and ggplot2 0.9.3.1
## If you used a different --out option for naming the output data file in
## gc_cov_annotate.pl, use that in place of blobplot.txt below
## 0.01 is the threshold - if fewer than this proportion of annotated contigs
## have a particular annotation, they are not shown in the legend

makeblobplot.R blobplot.txt 0.01 taxlevel_order
