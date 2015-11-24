#!/usr/bin/Rscript

package = require(pheatmap, quiet=T)

if (package==F) {
    cat("Package pheatmap is not installed. Exiting.\n")
    q()
}

require(RColorBrewer, quiet=T)
require(getopt, quiet=T)
opt1 <- c("help","h",0,"logical")
opt2 <- c("infile", "i", 1, "character")
opt3 <- c("outfile", "o", 1, "character")
opt4 <- c("sampledata", "s", 1, "character")
opt5 <- c("includecols", "I", 1, "character")
opt6 <- c("excludecols", "E", 1, "character")
opt7 <- c("clustercols", "c", 0, "logical")
opt8 <- c("clusterrows", "r", 0, "logical")
opt9 <- c("distmethod", "d", 1, "character")
opt10 <- c("genomedata", "G", 1, "character")
opt11 <- c("salinitysort", "S", 0, "logical")
opt12 <- c("groupsort", "O", 0, "logical")
opt13 <- c("height", "H", 1, "double")
opt14 <- c("png", "p", 0, "logical")
opt15 <- c("title", "t", 1, "character")
opt16 <- c("verbose", "v", 0, "logical")
opt17 <- c("zeroremove", "z", 0, "logical")
opt18 <- c("multiply", "m", 1, "double")

spec = matrix(c(opt1,opt2,opt3,opt4,opt5,opt6, opt7, opt8, opt9, opt10, opt11, opt12,opt13, opt14,opt15,opt16,opt17, opt18), byrow=T, ncol=4)

opt = getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help) || is.null(opt$infile)) {
        cat(getopt(spec, usage=TRUE))
        cat("Example: plot.FR.heatmap.R -i blastoutput.tab -o FRplot.pdf -s sampledata.tsv\n")
        q(status=1)
}

if (!is.null(opt$verbose)) {cat(paste(opt$infile),"\n", sep="")}

if (is.null(opt$outfile)) {opt$outfile = "FRplot"}
if (is.null(opt$clustercols)) {opt$clustercols = F}
if (is.null(opt$clusterrows)) {opt$clusterrows = F}
if (is.null(opt$distmethod)) {opt$distmethod = "euclidean"}
if (is.null(opt$height)) {opt$height=NA}
if (is.null(opt$title)) {opt$title=""}
title = gsub(".+\\/","",opt$outfile)
title = gsub("\\.pdf", "", title)
if (is.null(opt$png)) {pdf(paste(opt$outfile,".pdf",sep=""), useDingbats=F, title=title)
}else {png(paste(opt$outfile,".png",sep=""), res=300, width=2000, height=2000)}
fontsize = 6
fontsize_row = fontsize
fontsize_col = fontsize
## Set color gradient and make the low values more easy to distinguish
color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(180)
colorlow = colorRampPalette(c("darkgrey",color[1]))(20) ## Setting the lower scale to 10% of the total
color = append(colorlow,color)

## Read fragment recruitment matrix
frdata <- read.delim(opt$infile, header=T, row.names=1)

## Multiply the matrix if specified (e.g. by 1000 to get RPKB (reads per kb per 10000 reads of metagenome)
if (!is.null(opt$multiply)) {frdata <- frdata*opt$multiply}
if (!is.null(opt$includecols)) {
    include <- readLines(opt$includecols)
    cols <- which(colnames(frdata)%in%include)
    if (length(cols) >0) {frdata <- frdata[,cols]}
}

if (!is.null(opt$excludecols)) {
    exclude <- readLines(opt$excludecols)
    cols <- which(colnames(frdata)%in%exclude)
    if (length(cols) >0) {frdata <- frdata[,-cols]}
}

distance <- function(x, distmethod, dim) {
    if (distmethod=="spearman") {dist <- spearman(x,dim)
    }else {dist <- distmethod}
    dist
}

spearman <- function(x, dim="row") {
    if (dim=="row") {x.cor <- cor(t(x), method="spearman")
    }else {x.cor <- cor(x, method="spearman")}
    x.dist <- (-1*x.cor + 1)/2
    x.dist <- as.dist(x.dist)
    x.dist
}


## Read genome info if it exists, add names from column 2 in the file to the rownames in frdata
if (!is.null(opt$genomedata)) {
    genomedata <- read.delim(opt$genomedata, header=F, row.names=1, sep="\t")
    
    ## Sort genomedata by frdata
    genomedata <- genomedata[rownames(frdata),]
    genomedata <- matrix(genomedata)
    rownames(genomedata) <- rownames(frdata)
    ## Sort frdata by gemomedata
    genomedata <- genomedata[sort.list(as.vector(genomedata[,1])),]
    frdata <- frdata[names(genomedata),]
    rownames(frdata) <- paste(rownames(frdata),as.vector(genomedata), sep=" ")
}

if (ncol(frdata)>50) {fontsize_col=4}
if (nrow(frdata)>20) {fontsize_row=4}
if (nrow(frdata)>31) {fontsize_row=3}

## Read sample info if it exists, this part can be customized
## Now annotating samples by SampleGroup
## Minimum info should be sampleID SampleName SampleGroup
if (!is.null(opt$sampledata)) {
    ## Read sampledata
    sampledata <- read.delim(opt$sampledata, header=T, row.names=1)
    ## Make sample data conform
    sampledata <- sampledata[intersect(rownames(sampledata),colnames(frdata)),]
    ## Make frdata conform to sampledata
    frdata <- frdata[,intersect(colnames(frdata),rownames(sampledata))]
    ## Change colnames of frdata to SampleName
    colnames(frdata) <- as.vector(sampledata[match(colnames(frdata),rownames(sampledata)),1])
    ## Change rownames of sampledata
    rownames(sampledata) <- sampledata[,1]
    ## Delete zero-sum rows if specified
    if (!is.null(opt$zeroremove)) { frdata <- frdata[which(rowSums(frdata)!=0),]}
    ## Sort samples by SampleGroup or Salinity
    if (!is.null(opt$groupsort)) {sampledata <- sampledata[sort.list(sampledata[,2]),]}
    if (!is.null(opt$salinitysort)) {sampledata <- sampledata[sort.list(sampledata[,3]),]}
    frdata <- frdata[,rownames(sampledata)]
    
    ## Create annotations
    annotation <- cbind(as.vector(sampledata[,2]), round(sampledata[,3]))
    rownames(annotation) <- colnames(frdata)
    colnames(annotation) <- c("SampleGroup","Salinity")
    annotation[which(as.numeric(annotation[,2])>40),2] = ">40"
    annotation[which(is.na(annotation[,2])),2] <- "Unknown"
    annotation <- as.data.frame(annotation)
    
    ## Create salinity color ramp and mapping matrix
    colorramp <- colorRampPalette(rev(brewer.pal("Spectral",n=11)))
    salinitymat <- matrix(ncol=2, nrow=43) 
    salinitymat[,2] <- colorramp(43)
    salinitymat[2:42,1] <- seq(0,40)
    salinitymat[1,1] <- "Unknown"
    salinitymat[1,2] <- "white"
    salinitymat[43,1] <- ">40"
    Salinity = salinitymat[match(intersect(salinitymat[,1], as.vector(unique(annotation[,2]))),salinitymat[,1]),2]
    names(Salinity) = salinitymat[match(intersect(salinitymat[,1], as.vector(unique(annotation[,2]))),salinitymat[,1]),1]

    ## Create sample group color map
    ngroups = length(as.vector(unique(sampledata[,2])))
    if (ngroups<3) {n=3
    }else {n=ngroups}
    SampleGroup <- colorRampPalette(brewer.pal(name="Dark2", n=8))(n)
    if (ngroups<3) {SampleGroup <- SampleGroup[1:ngroups]}
    names(SampleGroup) <- as.vector(unique(sampledata[,2]))
    ann_colors  <- list(Salinity = Salinity, SampleGroup = SampleGroup)
    ## Plot
    pheatmap(frdata,main=opt$title,cellheight=opt$height, color=color,fontsize=fontsize, fontsize_row=fontsize_row, fontsize_col=fontsize_col,cluster_cols=opt$clustercols, 
             cluster_rows=opt$clusterrows,clustering_distance_rows=distance(frdata,opt$distmethod,"row"), clustering_distance_columns = distance(frdata, opt$distmethod, "col"), 
             annotation=annotation, annotation_colors=ann_colors)
    dev.off()
    colmat <- as.matrix(SampleGroup)
    rownames(colmat) <- names(SampleGroup)
    colnames(colmat) <- c("Color")

    write.table(colmat, "", sep="\t", quote=F)
    ## Calculate zero-sum samples
#    zerosum <- colnames(frdata)[which(colSums(frdata)==0)]
#    cat("Samples with no hits:\n")
#    cat(paste(paste(zerosum, collapse="\n"),"\n", sep=""))
    q()
}

## Delete zero-sum rows if specified
if (!is.null(opt$zeroremove)) { frdata <- frdata[which(rowSums(frdata)!=0),]}



pheatmap(frdata, main=opt$title, color=color,fontsize=fontsize,fontsize_row=fontsize_row, fontsize_col=fontsize_col, cluster_cols = opt$clustercols, 
         cluster_rows=opt$clusterrows, clustering_distance_rows = distance(frdata,opt$distmethod,"row"), clustering_distance_columns = distance(frdata,opt$distmethod,"col"))
dev.off()

## Calculate zero-sum samples
#zerosum <- colnames(frdata)[which(colSums(frdata)==0)]
#cat(paste(zerosum, collapse="\n"))
