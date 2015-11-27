barplot_lines <- function(data, bars, lwd=1, horiz=F) {
    ## EXAMPLE
    ##########
    #data <- cbind(c(10,30,40,20),c(40,30,20,10), c(10,10,50,30))
    #bars <- barplot(data)
    #barplot_lines(data,bars)

    #TODO: Add functionality for horizontal graph display
    # To insert lines between segments, first iterate over the columns in the data
    for (col in seq(1,ncol(data)-1)) {
        ## Get the x-coordinates for the end of this column and start of the next
        x <- c(bars[col]+0.5, bars[col+1]-0.5)
        ## Store the sum so far for each of the bars
        colsum_1 <- 0
        colsum_2 <- 0
        ## Iterate over the rows and connect lines
        for (row in seq(1,nrow(data))) {
            y <- c(data[row,col]+colsum_1, data[row, col+1]+colsum_2)
            colsum_1 <- colsum_1+data[row,col]
            colsum_2 <- colsum_2+data[row, col+1]
        }
    }
}
