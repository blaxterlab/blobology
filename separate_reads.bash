# Create proteobacteria database:

curl "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&db=nuccore&dopt=gilist&qty=2000000&filter=all&term=txid1224\[Organism:exp\]" >txid1224.gids
blastdb_aliastool -dbtype nucl -gilist txid1224.gids -db nt -out nt_Proteobacteria

# Create nematoda database:

curl "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&db=nuccore&dopt=gilist&qty=2000000&filter=all&term=txid6231\[Organism:exp\]" >txid6231.gids
blastdb_aliastool -dbtype nucl -gilist txid6231.gids -db nt -out nt_Nematoda

# Create actinobacteria database:

curl "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&db=nuccore&dopt=gilist&qty=10000000&filter=all&term=txid201174\[Organism:exp\]" >txid201174.gids
blastdb_aliastool -dbtype nucl -gilist txid201174.gids -db nt -out nt_Actinobacteria

# Do the blasts:

blastn -task megablast -query pa61-scaffolds.fa.l200 -db ~/scratch/blastdb/nt_Proteobacteria -outfmt 6 -evalue 1e-10 -max_target_seqs 1 -out pa61-scaffolds.fa.l200.nt_Proteobacteria.1e-10
    
blastn -task megablast -query pa61-scaffolds.fa.l200 -db ~/scratch/blastdb/nt_Nematoda -outfmt 6 -evalue 1e-10 -max_target_seqs 1 -out pa61-scaffolds.fa.l200.nt_Nematoda.1e-10

blastn -task megablast -query <(fastaqual_select.pl -f pa61-scaffolds.fa.l200 -s r) -db ~/scratch/blastdb/nt_Actinobacteria -outfmt 6 -evalue 1e-10 -max_target_seqs 1 -out pa61-scaffolds.fa.l200.nt_Actinobacteria.1e-10

# Separate by best blast hit:

cat pa61-scaffolds.fa.l200.nt_Actinobacteria.1e-10 pa61-scaffolds.fa.l200.nt_Proteobacteria.1e-10 >pa61-scaffolds.fa.l200.nt_contaminants.1e-10

blast_separate_taxa.pl -b1 pa61-scaffolds.fa.l200.nt_Nematoda.1e-10 -b2 pa61-scaffolds.fa.l200.nt_contaminants.1e-10

# Make list of contigids to be removed:

# based on megablast 1e-10 hit:

cut -f1 pa61-scaffolds.fa.l200.nt_contaminants.1e-10.only > toremove.contigids

# based on gc (col 3) and total cov (sum of cols 4 and 5):

awk '$3>=.55 && ($4+$5)<=50'  pa61-scaffolds.fa.bowtie2.txt | cut -f1 >>toremove.contigids
awk '$3>=.60 && ($4+$5)<=100' pa61-scaffolds.fa.bowtie2.txt | cut -f1 >>toremove.contigids
awk '$3>=.65'                 pa61-scaffolds.fa.bowtie2.txt | cut -f1 >>toremove.contigids
sort toremove.contigids | uniq | grep -v seqid >toremove.contigids.uniq; mv toremove.contigids.uniq toremove.contigids

awk '$2>=200' pa61-scaffolds.fa.bowtie2.txt | grep -v seqid | cut -f1 | ~/scripts/scripts/tools/fgrep.pl - -v -f toremove.contigids >tokeep.contigids

# extract reads mapping to these contigs

ASSEMBLY=pa61-scaffolds.fa.l200

for LIBNAME in ERR138446 ERR138445
do
    bowtie2_extract_reads_mapped_to_specific_contigs.pl -s <(samtools view $ASSEMBLY.$LIBNAME.bowtie2.bam) -id tokeep.contigids -out $LIBNAME.
done

