Blobology
=========

Tools for making blobplots or Taxon-Annotated-GC-Coverage plots (TAGC plots) to visualise the contents of genome assembly data sets as a QC step

Blaxter Lab, Institute of Evolutionary Biology, University of Edinburgh

Bash/perl scripts for creating and collating GC, Coverage and taxon annotation for a preliminary assembly  

Scripts to accompany "Blobology: exploring raw genome data for contaminants, symbionts and parasites using taxon-annotated GC-coverage plots."
Sujai Kumar, Martin Jones, Georgios Koutsovoulos, Michael Clarke, Mark Blaxter
submitted to Frontiers in Bioinformatics and Computational Biology special issue : Quality assessment and control of high-throughput sequencing data

This is an update to the code at github.com/sujaikumar/assemblage which was used in my thesis. I could have updated the code in that repository, but enough things have changed (the basic file formats as well) that I thought it made sense to create a new repo.

Installation
============

    git clone git://github.com/blaxterlab/blobology.git

You also need the following software:

1. samtools (tested with version 0.1.19) http://sourceforge.net/projects/samtools/files/samtools/0.1.19/
2. ABySS (tested with version 1.3.6, compiled with mpi) - http://www.bcgsc.ca/platform/bioinfo/software/abyss/releases/1.3.6
3. R (tested with version 2.15.2)
4. ggplot2, an R graphics package (tested with ggplot2_0.9.3.1)

