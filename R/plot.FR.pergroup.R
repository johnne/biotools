#!/usr/bin/env Rscript

require(RColorBrewer, quiet=T)
require(getopt, quiet=T)
require(ggplot2, quiet=T)
opt1 <- c("help","h",0,"logical")
opt2 <- c("filepattern", "f", 1, "character")
opt3 <- c("outfile", "o", 1, "character")
opt4 <- c("sampledata", "s", 1, "character")
opt5 <- c("png", "p", 0, "logical")
opt6 <- c("title", "t", 1, "character")
opt7 <- c("idmax", "I", 1, "double")
opt8 <- c("idmin", "i", 1, "double")
opt9 <- c("idby", "b", 1, "double")
opt10 <- c("groups","G",1,"character")
opt11 <- c("groupcols", "C", 1, "character")
opt12 <- c("width", "w", 1, "double")
opt13 <- c("height", "H", 1, "double")
spec = matrix(c(opt1,opt2,opt3,opt4,opt5,opt6, opt7, opt8, opt9, opt10, opt11, opt12, opt13), byrow=T, ncol=4)
opt = getopt(spec)
# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help) || is.null(opt$filepattern)) {
        cat(getopt(spec, usage=TRUE))
        q(status=1)
}

if (is.null(opt$outfile)) {opt$outfile="FRgroupplot"}
if (is.null(opt$title)) {opt$title = opt$outfile}
if (is.null(opt$idmax)) {opt$idmax = 100}
if (is.null(opt$idmin)) {opt$idmin = 70}
if (is.null(opt$idby)) {opt$idby = 5}
if (is.null(opt$width)) {opt$width = 7}
if (is.null(opt$height)) {opt$height = 7}
if (!is.null(opt$groups)) {
    Groups <- readLines(opt$groups)
}else{Groups <- c()}
if (!is.null(opt$groupcols)) {
    Groupcols <- read.delim(opt$groupcols, header=T, row.names=1, sep="\t")
    cols <- as.vector(Groupcols[,1])
    names(cols) <- rownames(Groupcols)
    Groups <- names(cols)
    Groupcols <- cols
}else{Groupcols <- c()}

readData <- function(idseq, sampledata, filepattern) {
    df <- data.frame(id=NA,Station=NA,Genome=NA,Group=NA, NA)
    df <- as.data.frame(df)
    colnames(df) <- c("id","Station","Genome","Group","RPKB")
    for (id in idseq) {
        infile = gsub("\\[ID\\]",id,filepattern)
        dta <- read.delim(infile, header=T, row.names=1)
        dta <- dta*1000
        samples <- intersect(colnames(dta),rownames(sampledata))
        for (sample in samples) {
            ids <- rep(id,nrow(dta))
            s <- rep(sample, nrow(dta))
            genomes <- rownames(dta)
            groups <- rep(as.vector(sampledata[sample,2]), nrow(dta))
            rpkb <- dta[,sample]
            df.tmp <- cbind(ids, s, genomes, groups, rpkb)
            colnames(df.tmp) <- colnames(df)
            df <- rbind(df,df.tmp)
        }
    }
    df
}

sampledata <- read.delim(opt$sampledata, header=T, row.names=1)
groups <- as.vector(unique(sampledata[,2]))
if (length(Groups)>0) {groups <- groups[which(groups%in%Groups)]}
if (length(Groupcols)>0) {
    groupcol <- Groupcols[intersect(groups,names(Groupcols))]
}else {
    groupcol <- colorRampPalette(brewer.pal(name="Dark2", n=8))(length(groups))
    names(groupcol) <- groups
}

idseq = seq(opt$idmax, opt$idmin, by=-opt$idby)
idseq = c(seq(70,90, by=5),seq(91,100, by=1))
filepattern = opt$filepattern

df.detailed <- readData(idseq, sampledata, filepattern)
df.detailed <- df.detailed[-1,]
df.detailed <- df.detailed[which(df.detailed$Group%in%groups),]
df.detailed <- as.data.frame(df.detailed)
df.detailed[,1] <- sapply(df.detailed[,1], as.numeric)
df.detailed[,5] <- sapply(df.detailed[,5], as.numeric)

groups <- intersect(unique(as.vector(df.detailed$Group)),groups)
groupcol <- groupcol[groups]


## Summarize the detailed data into mean recruitment across genomes in each group
df <- data.frame(id=NA,Station=NA,Group=NA, NA)
df <- as.data.frame(df)
colnames(df) <- c("id","Station","Group","RPKB")

for (group in groups) {
    df.group <- df.detailed[which(df.detailed$Group==group),]
    stations <- unique(as.vector(df.group$Station))
    df.group.agg <- aggregate(df.group, by=list(df.group$id,df.group$Station), FUN=mean, na.rm=T)
    df.group.agg$Station <- df.group.agg$Group.2
    df.group.agg$Group <- group
    df.group.agg <- df.group.agg[,-c(1,2,5)]
    df <- rbind(df,df.group.agg)
}
df <- df[-1,]
df[,1] <- sapply(df[,1], as.numeric)
df[,4] <- sapply(df[,4], as.numeric)
outfile = paste(opt$outfile,"pdf",sep=".")
pdf(outfile, useDingbats=F, width=opt$width, height=opt$height)
p = ggplot(df, aes(factor(id), RPKB, fill=Group))+geom_boxplot(outlier.size=0.5)+xlab("%ID")+ylab("Recruitment per kb per 10K reads")
if (!is.null(opt$groupcols)) { p = p + scale_fill_manual('Group',values = groupcol)}
p
dev.off()
#g <- ggplot_build(p)
#print(unique(g$data[[1]]["fill"]))

## PLOT Estuary groups
df.part <- df.detailed[which(df.detailed$Group%in%c("Estuary Baltic","Estuary USA")),]
df.part <- df.part[which(df.part$Group%in%c("Estuary Baltic","Estuary USA")),]
df.part$Station <- as.vector(sampledata[df.part$Station,5])
clusters1 <- paste("CLUSTER",seq(1,10), sep="")
clusters2 <- paste("CLUSTER",seq(11,21), sep="")
clusters3 <- paste("CLUSTER",seq(22,30), sep="")
clusters3 <- append(clusters3, "BalticASM")
pdf(paste(opt$outfile,"percluster.pdf", sep="."), width=20, useDingbats=F)
ggplot(df.part[which(df.part$Genome%in%clusters1),], aes(factor(id), RPKB, fill=Station)) + geom_point(aes(colour=Station, shape=Group))+xlab("%ID")+ylab("Recruitment per kb per 10K reads") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(.~Genome)
ggplot(df.part[which(df.part$Genome%in%clusters2),], aes(factor(id), RPKB, fill=Station)) + geom_point(aes(colour=Station, shape=Group))+xlab("%ID")+ylab("Recruitment per kb per 10K reads") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(.~Genome)
ggplot(df.part[which(df.part$Genome%in%clusters3),], aes(factor(id), RPKB, fill=Station)) + geom_point(aes(colour=Station, shape=Group))+xlab("%ID")+ylab("Recruitment per kb per 10K reads") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(.~Genome)
dev.off()
