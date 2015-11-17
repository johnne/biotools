spearman <- function(x, ...) {
    x <- t(x)
    x.cor <- cor(x, method="spearman")
    x.dist <- (-1*x.cor + 1)/2
    x.dist <- as.dist(x.dist)
    attr(x.dist, "method") <- "spearman"
    x.dist
}
