#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)# Create a parser
p <- arg_parser("draw PCA figure for plink pca  result")

# Add command line arguments
p <- add_argument(p, "eigenvec", help="input: eigenvec file", type="character" )
p <- add_argument(p, "Xaxis", help="input: which pc for x-axis  ", type="character")
p <- add_argument(p, "Yaxis", help="input: which pc for y-axis ", type="character")
p <- add_argument(p, "pop",  help="input: population table  ", type="character")
p <- add_argument(p, "output", help="output prefix", type="character")


# Parse the command line arguments
argv <- parse_args(p)

input <- argv$eigenvec
x <- as.numeric(argv$Xaxis)
y <- as.numeric(argv$Yaxis)
pop <- argv$pop
outpre <- argv$output

library(ggplot2)

vec <- read.table(  input, header = F , row.names = 1, sep = " ")
pop <- read.table(  pop, header = F, row.names = 1, sep = "\t")

population <- pop[rownames(vec),1]



p1<- ggplot(mapping = aes(x=vec[,x+1],y=vec[,y+1] ,colour = population))+ geom_point(size = 5) 
#p1 <- p1 +  stat_ellipse() + theme_bw()
p1 <- p1 + theme(panel.grid = element_blank(),
panel.background = element_blank(),
axis.line = element_line(),
axis.title.x = element_text(size = 14),
axis.title.y = element_text(size = 14),
axis.text.x = element_text(size = 12),
axis.text.y = element_text(size = 12))
p1 <- p1 +  xlab( paste("PC", x ,"(38.1%)",sep = "") )
p1 <- p1 +  ylab( paste("PC", y ,"(9.1%)",sep = "") )

pdf(file=paste(outpre, "pc", x, y, "pdf", sep = "."), height = 8, width = 10)
p1
dev.off()
png(file=paste(outpre, "pc", x, y, "png", sep = "."), height = 400, width = 500 )
p1
dev.off()
