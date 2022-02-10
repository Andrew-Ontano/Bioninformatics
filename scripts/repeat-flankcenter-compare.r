#!/usr/bin/env Rscript

#### Script for taking repeat-window outputs and plotting them
### requires: flank fasta file (fasta_flanks.py), center fasta file (fasta_extract_center.py), and lengt file (seq_length.py)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-f","--flanks"), type="character", default=NULL,
              help="input flank file' [default %default]",
              dest="flankFile"),
  make_option(c("-g","--center"), type="character", default=NULL,
              help="input center file' [default %default]",
              dest="centerFile"),
  make_option(c("-l","--length"), type="character", default=NULL,
              help="input length file' [default %default]",
              dest="lengthFile"),
  make_option(c("-c","--columns"), type="character", default=c(4),
              help="columns' [default %default]",
              dest="colsString"),
  make_option(c("-o","--output"), action="store_true", default=TRUE,
              help="write graph to output [default %default]",
              dest="save"),
  make_option(c("-s","--species"), type="character", default="none",
              help="species for output filename' [default %default]",
              dest="speciesName"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

options(error=traceback)
parser <- OptionParser(usage = "%prog -f flanks.tsv -g center.tsv -l lengths.tsv -c 4:(length(df)-1) -o -s species [options]",option_list=option_list)
opt = parse_args(parser)


if(opt$v){
  cat(paste0("PARAMETERS:\ninput flank file (-f): ", opt$flankFile,"\n"))
  cat(paste0("input center file (-g): ", opt$centerFile,"\n"))
  cat(paste0("input length file (-l): ", opt$lengthFile,"\n"))
  cat(paste0("collated columns (-c): ", opt$colsString,"\n"))
  cat(paste0("write output graph (-o): ", opt$save,"\n"))
  cat(paste0("species prefix (-s): ", opt$speciesName,"\n"))
}

#### function to make flank df more uniform
repeats_collate_flanks <- function(df){
  for(row in 1:nrow(df)){
    if(grepl("trailing",df[row,1])){
      df[row,1]=strsplit(df[row,1],"_")[[1]][1]
    }
    else{
      df[row,1]=strsplit(df[row,1],"_")[[1]][1]
    }
  }
  return(df %>% group_by(Scaffold) %>% summarise_all(.funs=sum))
}

#### function for binding length to df
repeats_bind_lengths <- function(df,lengths){
  df$Length = lengths$Length[match(df$Scaffold,lengths$Sequence)]
  #df$Length=rep(lengths$Length,each=2)
  return(df)
}

#### function for pivoting the columns you want to keep
# default is 4, which is the first tabulates column minimum
repeats_pivot <- function(df,cols=4:(ncol(df)-1)){
  return(df %>% tidyr::pivot_longer(cols))
}

#### function for full production
#inputs required: flank file, center file, length file, columns to collate
repeats_produce <- function(species,fileFlank,fileCenter,fileLength,colsString,save=FALSE){
  
  cols = eval(parse(text=colsString))
  flankValues = read.csv(fileFlank,sep="\t",stringsAsFactors = FALSE)
  lengthValues = read.csv(fileLength,sep="\t",stringsAsFactors = FALSE)
  centerValues = read.csv(fileCenter,sep="\t",stringsAsFactors = FALSE)

  #processing flanks
  flankValues = repeats_collate_flanks(flankValues)
  flankMax = max(rowSums(flankValues[,cols]))
  flankValues = repeats_bind_lengths(flankValues,lengthValues)
  flankValues = repeats_pivot(flankValues,cols)

  
  centerValues = repeats_bind_lengths(centerValues,lengthValues)
  centerMax = max(rowSums(centerValues[,cols]))
  centerValues = repeats_pivot(centerValues,cols)
  centerValues$value=-centerValues$value
  centerValues$Side="Center"
  flankValues$Side="Flank"
  df <- bind_rows(flankValues,centerValues)
  
  ggplot(df,aes(fill=name,y=value,x=reorder(Scaffold,desc(Length)))) + geom_bar(stat="identity",position="stack")
#  ylim=c(0,round_any(max(flankMax,centerMax),50,f=ceiling))
  
#  flankPlot<- ggplot(flankValues,aes(fill=name,y=value,x=reorder(Scaffold, desc(Length)))) +
#    geom_bar(position="stack",stat="identity") +
#    theme(axis.text.x=element_blank()) +
#    ggtitle("Flank sequences") +
#    ylim(ylim)
  
#  centerPlot<- ggplot(centerValues,aes(fill=name,y=value,x=reorder(Scaffold, desc(Length)))) +
#    geom_bar(position="stack",stat="identity") +
#    theme(axis.text.x=element_blank()) +
#    ggtitle("Central sequence") +
#    ylim(ylim)
  
  #ggarrange(flankPlot,centerPlot,ncol=1,nrow=2)
  if(save){
      ggsave(filename = paste0(species, "_innerouter.png"), width = 10, height = 6, units = "in", dpi = 600, limitsize = F)
  }
  
}

repeats_produce(opt$speciesName,opt$flankFile,opt$centerFile,opt$lengthFile,opt$colsString,opt$save)