#!/usr/bin/env Rscript
####
##### Repeats along (catermerize) sequence for multiple libraries

#load dependencies
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))


spectraPlot = function(values, tripletColors, legend=FALSE, facet=FALSE, frequencies=FALSE, ylims=TRUE){
	if(frequencies){
		p = ggplot() + geom_area(data=values, aes(fill=name, x=(End+Start)/2, y=value), stat="identity", position="stack")
	}else{
		p = ggplot() + geom_area(data=values, aes(fill=name, x=(End+Start)/2,y=value/(End-Start+1)), stat="identity", position="stack")
	}
	if(ylims){
		p = p + scale_y_continuous(limits=c(0,1), expand=c(0, 0))
	}else{
		p = p + scale_y_continuous(expand=c(0, 0))
	}
	
	p = p +
		scale_fill_manual(values=tripletColors) +
		scale_x_continuous(
			limits=c(min(values$Start), max(values$End)),
			breaks=waiver(),
			minor_breaks=waiver(),
			n.breaks=10,
			expand=c(0, 0)
		) +
		xlab("Window Position (nucleotide)") +
		ylab("Proportion") +
		theme_bw() +
		theme(
			axis.title.x = element_text(size=6),
			axis.title.y = element_text(size=6),
			legend.text = element_text(size=5),
			legend.key.size = unit(5,"pt"),
			legend.spacing = unit(-3,"pt"),
			axis.text.x = element_text(angle=0,vjust=0.5),
			line=element_blank(),
			axis.ticks = element_line(),
			legend.title= element_blank()
		)

	if(facet){
		p = p + facet_grid(rows = vars(values$Library))
	}
	if(!legend){
		p = p + theme(legend.position = "none")
	}
	return(p)
}

triplerColor = function(triplet){
	bases = c("A","C","G","T")
	colors = c("#DD0000","#990099","#00CCCC","#009966")
	return(colors[which(triplet == bases)])
}

#Argument parser
option_list <- list(
	make_option(c("-i","--input"), type="character", default=NULL, help="triplet tsv file' [default %default]", dest="input_filename"),
	make_option(c("-w", "--window-size"), type="character", default=NULL, help="if specified, generate plot only from the closest intervals of \"N,M\". Otherwise generate a plot for all positions [default %default]", dest="window_size"),
	make_option(c("-s", "--sequence"), type="character", default=NULL, help="if specified, generate plot only from the sequences with names in \"A,B,C\". Otherwise generate a plot for each sequence id [default %default]", dest="sequences"),
	make_option(c("-n", "--libraries"), type="character", default=NULL, help="if specified, generate plot only from the libraries with names in \"A,B,C\". Otherwise generate a plot for each library id [default %default]", dest="libraries"),
	make_option(c("-g", "--gff-file"), type="character", default=NULL, help="if specified, generate plot of overlapping gene annotations from supplied gff. [default %default]", dest="gffFile"),
	make_option(c("-t", "--gff-tracks"), type="character", default=NULL, help="curate which gff types to use in types \"A,B,C\". [default %default]", dest="gffTracks"),
	make_option(c("-l","--show-legend"), action="store_true", default=FALSE, help="display triplet color legend [default %default]", dest="legend"),
	make_option(c("-o","--output-file"), type="character", default=NULL, help="output filename [default %default]", dest="output_filename"),
	make_option(c("-r","--resolution"), type="numeric", default=300, help="plotting dpi resolution [default %default]", dest="resolution"),
	make_option(c("-f","--freq"), action="store_true", default=FALSE, help="data already transformed to frequencies [default %default]", dest="frequencies"),
	make_option(c("-y","--ylims"), action="store_false", default=TRUE, help="Limit results to y-axes between 0,1 [default %default]", dest="ylims")
)
options(error=traceback)
parser = OptionParser(usage = "%prog -i triplet.ysv [options]",option_list=option_list)
opt = parse_args(parser)

# prepare outputfile space
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

# gff preparation
gff = NULL
if(!is.null(opt$gffFile)){
	suppressPackageStartupMessages(library(ape))
	gff = read.gff(opt$gffFile)
	# filter
	if(!is.null(opt$gffTracks)){
		gff = gff %>% filter(type%in%as.vector(opt$gffTracks))
	}
	if(!is.null(opt$window_size)){
		coords = strsplit(opt$window_size,",")
		gff = gff %>% filter(start >= as.numeric(coords[[1]][1]))
		gff = gff %>% filter(end <= as.numeric(coords[[1]][2]))
	}
}

# read tsv values or exit with error
if(is.null(opt$input_filename)){
  cat("Error: No tsv specified. See usage 'with spectra-plot.r -h'\n")
  quit()
}else{
	values = read.csv(opt$input_filename,sep="\t",stringsAsFactors = FALSE)
}

# filter out wasteful data
if(!is.null(opt$libraries)){
	values = values %>% filter(Library%in%as.vector(opt$libraries))
}

if(!is.null(opt$sequences)){
	values = values %>% filter(Sequence%in%as.vector(opt$sequences))
}

if(!is.null(opt$window_size)){
	coords = strsplit(opt$window_size,",")
	values = values %>% filter(Start >= as.numeric(coords[[1]][1]))
	values = values %>% filter(End <= as.numeric(coords[[1]][2]))
}

lib.names = names(table(values$Library))
seq.names = names(table(values$Sequence))

if(length(lib.names)>1 && !is.null(opt$gffFile)){
	print("Warning: the gff track will be omitted because multiple libraries are being plotted.")
	gff = NULL
}

# Pivot table for stacking of columns
values = values %>% tidyr::pivot_longer(cols=starts_with(c("A","C","G","T")))

tripletNames=names(table(values$name))
tripletColors=sapply(tripletNames,triplerColor)

# write multiple plots in a single frame if sequence names are the same
for(seq in seq.names){
	seq.filename = paste0(output_file[1], '_', seq, '.', output_file[2])
	temp.values = values %>% filter(Sequence==seq)
	if(length(lib.names)>1){
		faceted=TRUE
		height.factor = length(names(table(temp.values$Library)))

	}else{
		faceted=FALSE
		height.factor = 1
	}
	
	p = spectraPlot(temp.values, tripletColors, legend=opt$legend, facet=faceted, frequencies=opt$frequencies, ylims=opt$ylims)
	
	if(!is.null(gff)){
		p = p + geom_segment(data=gff%>%filter(seqid==seq), aes(x=start, xend=end, y=-.03, yend=-.03,color=strand), size=4) + scale_y_continuous(limits=c(-.06, 1), expand=c(0, 0))
	}
	
	ggsave(filename=seq.filename,device=output_file[2], width=10, height=1+height.factor*2, units="in", dpi=opt$resolution, limitsize=F)
}

