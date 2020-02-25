library(data.table)
library(RColorBrewer) # for brewer.pal
library(ggplot2)
library(grDevices)
library(GenomicRanges)
#library(grid)
library(ggthemes)
library(scales)
source("lib/theme_publish.R") # defines theme_Publication()
source("lib/graphic.R")

# Suppl information from: https://academic.oup.com/bib/article/18/2/205/2562739
# Function require a vector with expression of one gene in different tissues.
# If expression for one tissue is not known, gene specificity for this gene is NA
# Minimum 2 tissues
fTau <- function(x){
    if(all(!is.na(x))){
        if(min(x, na.rm=TRUE) >= 0){
            if(max(x)!=0){
                x <- (1-(x/max(x)))
                res <- sum(x, na.rm=TRUE)
                res <- res/(length(x)-1)
            }else{
                res <- 0
            }
 		}else{
            res <- NA
            #print("Expression values have to be positive!")
 		} 
 	}else{
        res <- NA
        #print("No data for this gene avalable.")
    } 
    return(res)
}

