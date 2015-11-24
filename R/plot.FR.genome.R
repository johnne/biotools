#!/usr/bin/env Rscript
require(RColorBrewer, quiet=T)
require(getopt, quiet=T)
require(ggplot2, quiet=T)

opt1 <- c("help","h",0,"logical")
opt2 <- c("blastfile", "b", 1, "character")
opt3 <- c("outfile", "o", 1, "character")
opt4 <- c("sampledata", "s", 1, "character")
opt5 <- c("groupcols", "C", 1, "character")
opt6 <- c("idmin", "i", 1, "double")
opt7 <- c("idmax", "I", 1, "double")
opt8 <- c("alignlen", "a", 1, "double")
opt9 <- c("eval", "e", 1, "double")
opt10 <- c("genomeinfo", "G", 1, "character")
opt11 <- c("width", "w", 1, "double")

spec = matrix(c(opt1,opt2,opt3,opt4,opt5,opt6, opt7, opt8, opt9, opt10, opt11), byrow=T, ncol=4)
opt = getopt(spec)
# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help) || is.null(opt$blastfile)) {
        cat(getopt(spec, usage=TRUE))
        q(status=1)
}

if (is.null(opt$idmin)) {opt$idmin <- 50}
if (is.null(opt$alignlen)) {opt$alignlen <- 200}
if (is.null(opt$eval)) {opt$eval <- 1e-5}
if (is.null(opt$width)) {opt$width <- 15}

## Assume genomeinfo has sizes in first column and eventual renaming labels in second
genomeinfo <- read.delim(opt$genomeinfo, header=F, row.names=1)
Genomes <- rownames(genomeinfo)

if (!is.null(opt$groupcols)) {
    Groupcols <- read.delim(opt$groupcols, header=T, row.names=1, sep="\t")
    cols <- as.vector(Groupcols[,1])
    names(cols) <- rownames(Groupcols)
    Groupcols <- cols
    groupcol <- Groupcols
}else{groupcol <- c()}

sampledata <- read.delim(opt$sampledata, header=T, row.names=1)

## Read and filter the input
dta <- read.delim(opt$blastfile, header=F, sep="\t")
colnames(dta) <- c("query", "subject", "id", "alignlen", "mismatches", "gaps", "qstart", "qend", "sstart", "ssend", "eval", "score")
goodrows <- intersect(intersect(which(dta$id>=opt$idmin), which(dta$alignlen>=opt$alignlen)),which(dta$eval<opt$eval))
genomes <- gsub(".+=","",dta$subject)
goodrows <- intersect(goodrows, which(genomes%in%Genomes))
samples <- gsub(".+=","",dta$query)
samples.u <- unique(samples)
samples.u <- intersect(samples.u,rownames(sampledata))
goodrows <- intersect(goodrows, which(samples%in%samples.u))
dta <- dta[goodrows,]
samples.u <- intersect(samples.u, unique(gsub(".+=","",dta$query)))


genomes <- unique(gsub(".+=","",dta$subject))

df <- data.frame(Genome=NA, Sample=NA, Group=NA, id=NA, Gene=NA)

for (sample in samples.u) {
    group <- as.vector(sampledata$SampleGroup[which(rownames(sampledata)==sample)])
    sampledta <- dta[grep(sample,dta$query),]
    genomes <- gsub(".+=","",sampledta$subject)
    positions <- as.numeric(gsub("PROKKA_MOD_([0-9]+).+","\\1",sampledta$subject))
    df.tmp <- matrix(data=NA, ncol=5, nrow=nrow(sampledta))
    colnames(df.tmp) <- colnames(df)
    df.tmp <- as.data.frame(df.tmp)
    df.tmp$Gene <- positions
    df.tmp$Genome <- genomes
    df.tmp$id <- sampledta$id
    df.tmp$Sample <- sample
    df.tmp$Group <- group
    df <- rbind(df,df.tmp)
}

df <- df[-1,]
df[,4] <- sapply(df[,4], as.numeric)
df[,5] <- sapply(df[,5], as.numeric)

if (is.null(opt$outfile)) {opt$outfile <- "FRclassic"}
pdf(paste(opt$outfile, "pdf", sep="."), useDingbats=F, width=opt$width)
s=16
for (genome in Genomes) {
    if (ncol(genomeinfo)>1) {
        info <- as.vector(genomeinfo[genome,2])
        info <- paste(c("(",info,")"), collapse="")
    }else{info <- ""}
    xlim <- c(1,as.vector(genomeinfo[genome,1]))
    ylim <- c(opt$idmin,100)
    df.part <- df[which(df$Genome==genome),]
    p = ggplot(df.part, aes(Gene,id))+geom_point(aes(colour=Group), alpha = 0.5) + labs(title=paste(c(genome,info,sep=" "), collapse=" ")) + ylab("% nucleotide identity")
    p = p + theme(plot.title=element_text(size=s+6), legend.text=element_text(size=s+2), legend.title=element_text(size=s+4), axis.text.x=element_text(size=s), axis.text.y=element_text(size=s), axis.title.x=element_text(size=s+4), axis.title.y=element_text(size=s+4))
    if (length(groupcol)>0) {
        groupcol.part <- groupcol[unique(df.part$Group)]
        p = p + scale_colour_manual('Group', values = groupcol.part)
    }
    p = p + guides(colour = guide_legend(override.aes = list(size=8)))
    print(p)
}

dev.off()

