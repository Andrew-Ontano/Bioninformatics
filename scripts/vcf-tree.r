#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

#Argument parser
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL, help="tsv file'", dest="input_filename"),
  make_option(c("-s","--sample"), type="character", default=NULL, help="Sample to plot relations", dest="sample"),
  make_option(c("-o","--output-file"), type="character", default=NULL, help="Output prefix", dest="output_filename"),
  make_option(c("-b","--baseline"), action="store_true", default=FALSE, help="Incorporate background baseline distances", dest="baseline"),
  make_option(c("-r","--resolution"), type="numeric", default=300, help="Plotting dpi resolution [default %default]", dest="resolution")
)
options(error=traceback)
parser = OptionParser(usage = "%prog -i data.tsv [options]", option_list=option_list)
opt = parse_args(parser)

if(is.null(opt$output_filename)){
  split_file = strsplit(opt$input_filename,"/")[[1]]
  output_file = c(split_file[length(split_file)],"png")
}else{
  split_file = strsplit(opt$output_filename,"\\.")[[1]]
  output_type = split_file[length(split_file)]
  if(output_type=="svg"){
    suppressPackageStartupMessages(library(svglite))
  }
  output_file = c(split_file[1], output_type)
}

if(is.null(opt$input_filename)){
  cat("Error: No tsv specified. See usage 'with vcf-tree.r -h'\n")
  quit()
}else{
  data = read.csv(opt$input_filename, sep="\t", header=TRUE, stringsAsFactors = FALSE)
}

data_long = data %>%
  pivot_longer(cols = names(data)[7:length(names(data))], names_to = "Pair", values_to = "Distance")

if(is.null(opt$sample)){
    p = ggplot(data_long%>%filter(Type=="Window"), aes(x = (Start+End)/2, y = Distance, color=Pair)) +
      geom_line(linewidth=.1) + scale_x_continuous(labels = scales::scientific_format(digits=2)) +
      geom_smooth(se = FALSE, linewidth = 0.5, method="loess") +
      theme_minimal()+theme(text = element_text(size=6), axis.text.x = element_text(size = 4)) + facet_wrap(~ Chromosome ,scales = "free_x")
}else{
    p = ggplot(data_long%>%filter(Type=="Window")%>%filter(stringr::str_detect(Pair, opt$sample)), aes(x = (Start+End)/2, y = Distance, color=Pair)) +
      geom_line(linewidth=.1) + scale_x_continuous(labels = scales::scientific_format(digits=2)) +
      geom_smooth(se = FALSE, linewidth = 0.5, method="loess") +
      theme_minimal()+theme(text = element_text(size=6), axis.text.x = element_text(size = 4)) + facet_wrap(~ Chromosome ,scales = "free_x")
}
if(opt$baseline){
  chromosomes = names(table(data$Chromosome))
  chromosomes = chromosomes[chromosomes!="All"]
  chromosomeDF = data.frame(matrix(ncol=2,nrow=0))
  names(chromosomeDF) = c("Chromosome","Mean")
  class(chromosomeDF$Chromosome) = "character"
  class(chromosomeDF$Mean) = "numeric"


  chromosomeDistances = data%>%filter(Type=='Chromosome')


  for(chromosome in chromosomes){
    mean = rowMeans((chromosomeDistances%>%filter(Chromosome==chromosome))[,7:length(names(data))])
    chromosomeDF%>%add_row(Chromosome=chromosome, Mean=mean[1])
  }

  p = p + stat_mean_line(linetype="dashed", color = "red")

}

ggsave(filename=paste0(output_file[1],'.',output_file[2]), device=output_file[2], width=10, height=6, units="in", dpi=opt$resolution, limitsize=F)
