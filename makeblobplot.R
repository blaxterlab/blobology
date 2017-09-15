#!/usr/bin/env Rscript

# 2015-07-15
# works with R 3.0.2 and ggplot2_1.0.1
# Check ggplot2 help forums or contact sujai.kumar@gmail.com if something doesn't run
# because of updated programs/packages

# linda.bakker@wur.nl -- changes
# - additional script argument: sequence length cutoff & plot title
# - legend now includes the number of assigned sequences in a bin.
# - added colours to the plot
# - added a new bin "below threshold"
# - apply sequence cutoff steps before nr-of-colours cutoff 
#########################################################################


library(ggplot2)
library(reshape2)

subplotwidth=1200;
subplotheight=1000;

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --args1=file		- blobplot.txt file
      --args2=bin_cutoff	- nr of sequences that must be assigned to a specific classification bin, this is a percentage of the total nr of annotated sequences (i.e 0.01=1%)
      --args3=length_cutoff   - sequence length cutoff 
      --args4=taxlevel	- taxlevel to plot, taxlevel_superkingdom | taxlevel_phylum | taxlevel_order | taxlevel_species
      --args5=plottitle	- plot title
 
      Example:
      ./test.R blobplot.txt 0.01 200 taxlevel_order \"name of my plot\" \n\n")
 
  q(save="no")
}

#Load arguments
arg_input_file = args[1]
arg_ignore_below_prop = as.numeric(args[2])
arg_length_cutoff = as.numeric(args[3])
arg_taxlevel = args[4]
arg_plottitle = args[5]

orig <- read.delim(arg_input_file,header=TRUE,sep="\t")

# if cov_colnames >1, then create a new column cov_total = sum of all cov columns
if (length(grep("^cov_",colnames(orig))) >1) orig$cov_total = rowSums(orig[,grep("^cov_",colnames(orig))])

cov_colnames=colnames(orig)[grep("^cov_",colnames(orig))]
tax_colnames=colnames(orig)[grep("^taxlevel_",colnames(orig))]

numcols=length(cov_colnames)
taxlevel=arg_taxlevel;

m<-melt(orig,id.vars=c("seqid","len","gc",taxlevel),measure.vars=cov_colnames, variable.name="read_set", value.name="cov")

# a new level "Below threshold" is created:
levels(m[,taxlevel]) = c(levels(m[,taxlevel]), "Below threshold")

# filter for sequence length on annotated sequences:
m[,taxlevel][which((m$len<=arg_length_cutoff) & (m$taxlevel != "Not annotated"))]<-"Below threshold" 


# get the "annotated" sequences
annotated<-m[m[,taxlevel]!="Not annotated",]
total<-dim(annotated)[1]

# filtering: if the nr of sequences assigned to a classification bin is below the threshold (at least N% of annotated sequences should be assigned to the same classification), re-classify
levels(m[,taxlevel])[which(table(m[,taxlevel])<arg_ignore_below_prop*total)]<-"Below threshold" 

# restrict the plot to the number of colours (18 + grey for not annotated)
# (thanks to https://github.com/hobrien for the fix)
if (length(levels(m[,taxlevel])) > 19) {
  levels(m[,taxlevel])[which(table(m[,taxlevel])<=sort(as.numeric(table(m[,taxlevel])), decreasing=T)[18])]<-"Other"
}

taxnames=names(sort(table(m[,taxlevel]),decreasing=TRUE))
taxnames=c("Not annotated", taxnames[taxnames != "Not annotated"])
m[,taxlevel] <- factor(m[,taxlevel], levels = taxnames)

png(paste(arg_input_file,taxlevel,"png",sep="."), (numcols * subplotwidth), (1 * subplotheight) + 300, units="px",res=100)

theme_set(theme_bw())
colourtable=list(
c("#DDDDDD"),
c("#DDDDDD","#4477AA"),
c("#DDDDDD","#4477AA","#CC6677"),
c("#DDDDDD","#4477AA","#DDCC77","#CC6677"),
c("#DDDDDD","#4477AA","#117733","#DDCC77","#CC6677"),
c("#DDDDDD","#332288","#88CCEE","#117733","#DDCC77","#CC6677"),
c("#DDDDDD","#332288","#88CCEE","#117733","#DDCC77","#CC6677","#AA4499"),
c("#DDDDDD","#332288","#88CCEE","#44AA99","#117733","#DDCC77","#CC6677","#AA4499"),
c("#DDDDDD","#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#AA4499"),
c("#DDDDDD","#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499"),
c("#DDDDDD","#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#882255","#AA4499"),
c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#882255","#AA4499"),
c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499"),
c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499","#F9BD7E"),
c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499","#F9BD7E","#D24D3E"),
c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499","#F9BD7E","#D24D3E","#FB9A29"),
c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499","#F9BD7E","#D24D3E","#FB9A29","#CC4C02"),
c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499","#F9BD7E","#D24D3E","#FB9A29","#CC4C02","#662506"),
c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499","#F9BD7E","#D24D3E","#FB9A29","#CC4C02","#662506","#777777")
)

labels_info<- apply(as.data.frame(ftable(m[,taxlevel])),1,paste,collapse=" ")
g<-ggplot() + scale_colour_manual(values=colourtable[[length(levels(m[,taxlevel]))]], name="", limits=levels(m[,taxlevel]),labels=labels_info)


for (t in levels(m[,taxlevel])) {
  g <- g + geom_point(data=m[m[,taxlevel]==t,],aes_string(x="gc", y="cov", colour=taxlevel), size=2, alpha=I(1/3))
}
y_axis_breaks = c(10,15,25,50,100,200,500,1000);
g<-g +
  ggtitle(arg_plottitle) +
  #facet_wrap(~read_set, ncol=numcols) + 
  scale_y_log10(breaks = y_axis_breaks) + scale_x_continuous(limits=c(0, 1),breaks = seq(0,1,.1)) +
  labs(x="GC content", y="Read coverage") + 
  guides(colour = guide_legend(nrow=10, override.aes = list(alpha = 1,size=9))) + 
  theme (
    plot.title = element_text(colour = "black", size = 20, vjust = 0.5),
    strip.text.x = element_text(colour = "black", size = 20, vjust = 0.5),
    axis.text.x  = element_text(colour = "black", size = 20, vjust = 1),
    axis.text.y  = element_text(colour = "black", size = 15, vjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 15, vjust = 0),
    axis.title.y = element_text(colour = "black", size = 20, hjust = 0.5, vjust = 0.5, angle=90),
    legend.text  = element_text(colour = "black", size = 15, vjust = 0),
    legend.title = element_text(colour = "black", size = 20, vjust = 0, hjust = 0, lineheight=1),
    legend.justification=c(1,1), legend.position="bottom", legend.direction="horizontal"
  )
g
dev.off()
