#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("draw sturcture figure for admixture result")

# Add command line arguments
p <- add_argument(p, "dir", help="input: Q matirx dir", type="character")
p <- add_argument(p, "sample", help="input: samplefile, ie  .nosex", type="character")
p <- add_argument(p, "sample_order", help="input: sample order file, define order in output figure ", type="character")
p <- add_argument(p, "output", help="output prefix", type="character")

# Parse the command line arguments
argv <- parse_args(p)

dir <- argv$dir
samp <- argv$sample
ord  <- argv$sample_order
outpre <- argv$output


#dir <- "result/"
#samp <- "./all.fam"
#ord  <- "./all.fam"
#outpre <- "test"



library(pophelper)
library(ggplot2)


Qfiles <- list.files( dir, pattern = "Q", full.names = T )
qlist <- readQ(Qfiles,  indlabfromfile=F)


label <- read.table( samp, header = F)
ordinfo <- read.table(ord, header = F, stringsAsFactors=F )

for( i in 1:length(qlist) ){
    rownames(qlist[[i]]) <- label$V1
    qlist[[i]] <- qlist[[i]][as.character(ordinfo[,1]),]
}

# head(qlist[[3]])

p1_sort  <- plotQ(sortQ(qlist)[1:length(qlist)], 
      sortind = "all",  ## ind排序 
      imgoutput = "join",  ## 一张图还是多张
      returnplot=T,  
      exportplot=F,
      clustercol=c("#ff0000","#2b92df", "#6a329f","#99ccff"),
      useindlab = T ,    ## 显示ind名称
      sharedindlab= F ,   ## ind出现一次
      showindlab=T  )

ggsave(filename = paste(outpre, "sorted.pdf",sep=".") , 
       p1_sort$plot[[1]],  
       width = 20,
       height = 8  
       )


p1_order  <- plotQ(sortQ(qlist)[1:length(qlist)], 
                  sortind = NA,  ## ind排序 
                  imgoutput = "join",  ## 一张图还是多张
                  returnplot=T,  
                  exportplot=F,
                  clustercol=c("#ff0000","#2b92df", "#6a329f","#99ccff"),
                  useindlab = T ,    ## 显示ind名称
                  sharedindlab= T ,   ## ind出现一次
                  showindlab=T  )

ggsave(filename = paste(outpre, "ordered.pdf",sep=".") , 
       p1_order$plot[[1]],  
       width = 20,
       height = 8  
)





