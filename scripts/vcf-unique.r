#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

#Argument parser
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL, help="tsv file'", dest="input_filename"),
  make_option(c("-x","--index"), type="character", default=NULL, help="index file'", dest="index_filename"),
  make_option(c("-o","--output-file"), type="character", default=NULL, help="Output prefix", dest="output_filename"),
  make_option(c("-r","--resolution"), type="numeric", default=300, help="Plotting dpi resolution [default %default]", dest="resolution")
)
options(error=traceback)
parser = OptionParser(usage = "%prog -i data.tsv -x index.tsv -o output.png [options]", option_list=option_list)
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
}else if(is.null(opt$index_filename)){
  cat("Error: No index specified. See usage 'with vcf-tree.r -h'\n")
  quit()
}else{
  data = read.csv(opt$input_filename, sep="\t", header=TRUE, stringsAsFactors = FALSE)
}

# Read the TSV file (update the path to your file)
data <- read.table(opt$input_filename, header = TRUE, sep = "\t")
# Read denominator data
data.denominator <- read.table(opt$index_filename, header = TRUE, sep = '\t')
data["Baseline"] = data.denominator["Baseline"]
# Check the structure of the data
data_long <- data %>%
  pivot_longer(cols = colnames(data)[3]:colnames(data)[length(colnames(data))-1], names_to = "Sample", values_to = "Sample_Variant_Count")

p = ggplot(data_long, aes(x = Window, y = Sample, fill = Sample_Variant_Count/Baseline)) +
  geom_tile() +
  facet_wrap(~ Chrom, scales = "free_x") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Heatmap of Variant Counts Across Samples and Genomic Windows",
       x = "Genomic Window Start Position",
       y = "Sample",
       fill = "Variant Count")

ggsave(filename=paste0(output_file[1],'.',output_file[2]), device=output_file[2], width=10, height=6, units="in", dpi=opt$resolution, limitsize=F)
