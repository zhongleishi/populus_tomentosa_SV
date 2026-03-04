#!/usr/bin/env Rscript

args=commandArgs(TRUE)
if (length(args) != 2) {
          print ("usage: <gwas file(4 column)>   <out.prefix> ")
  q()
}

input <- args[1]
#max_y <- as.numeric(args[3])
prefix <- args[2]

library(CMplot)

gwas<- read.table(input, header = T);

cutoff <- 1/nrow(gwas) ###原始0.05

sigSNP <- gwas[gwas[,4] < cutoff,]

write.table( sigSNP, file = paste( prefix, "sigSite.out", sep =  "."), row.names = F, quote = F)


png(paste(prefix, "_manhattan_threshold.png", sep = ""), width=960, height=480)
CMplot(gwas,#[gwas$CHR==1,],只画第一条染色体
       plot.type = "m",
       LOG10 = T,
       col = c("#92C6A0", "#53A6D9"),
       cex = 0.3,
       ylab.pos = 2,
       #cex.axis = 1, 
       #chr.den.col=c("darkgreen", "yellow", "red"),##snp图
       threshold = c(1)/nrow(gwas), ## 显著性阈值，原始threshold = c(0.01,0.05)/nrow(gwas),
       threshold.col=c('grey','black'),  ## 阈值线颜色
       threshold.lty = c(1,2), ## 阈值线线型
       threshold.lwd = c(1,1), ## 阈值线粗细
       amplify = T,  ## 放大显著SNP
       signal.cex = c(1,1), ## 点大小
       signal.pch = c(20,20), ## 点形状
       signal.col = c("red","blue"), ## 点颜色
       file.output = F )
dev.off()

png(paste(prefix, "_manhattan.png", sep = ""), width=960, height=480)
CMplot(gwas,plot.type = "m",
       LOG10 = T,
       col = c("#92C6A0", "#53A6D9"),
       cex = 0.3,
       ylab.pos = 2, 
       #cex.axis = 1,
       threshold=NULL,
       file.output = F )

dev.off()

png(paste(prefix, "_qqplot.png", sep = ""), width=480, height=480)
CMplot(gwas,
       plot.type = "q", ## 绘制QQplot
       box=T, ## 是否加边框
       conf.int=T, ## 是否绘制置信区间
       conf.int.col=NULL, ## 置信区间颜色
       threshold.col="red", ## 对角线颜色
       threshold.lty=2,  ## 线型
       cex = 0.8,
       ylab.pos = 2, 
       #cex.axis = 1,
       main = "QQ-plot",
       file.output = F )
dev.off()

#gwas_chr14_region <- gwas[gwas$CHR == 14 & gwas$BP >= 13330001 & gwas$BP <= 13430000, ]

## 曼哈顿图1
pdf(paste(prefix, "_manhattan_threshold.pdf", sep = ""), width=18, height=5)
CMplot(gwas,
       plot.type = "m",
       LOG10 = T,
       col = c("#92C6A0", "#53A6D9"),
       cex = 0.3,
       ylab.pos = 2, 
       #cex.axis = 1,
       threshold = c(1)/nrow(gwas), ## 显著性阈值
       threshold.col=c('grey','black'),  ## 阈值线颜色
       threshold.lty = c(1,2), ## 阈值线线型
       threshold.lwd = c(1,1), ## 阈值线粗细
       amplify = T,  ## 放大显著SNP
       signal.cex = c(1,1), ## 点大小
       signal.pch = c(20,20), ## 点形状
       signal.col = c("red","blue"), ## 点颜色
       file.output = F )
dev.off()

## 曼哈顿图2

pdf(paste(prefix, "_manhattan.pdf", sep = ""), width=10, height=5)
CMplot(gwas,plot.type = "m",
       LOG10 = T,
       col = c("#92C6A0", "#53A6D9"),
       cex = 0.3,
       ylab.pos = 2, 
       #cex.axis = 1,
       threshold=NULL,
       file.output = F )

dev.off()

## QQplot 
pdf(paste(prefix, "_qqplot.pdf", sep = ""), width=10, height=10)
CMplot(gwas,
       plot.type = "q", ## 绘制QQplot
       box = T, ## 是否加边框
       conf.int=T, ## 是否绘制置信区间
       conf.int.col=NULL, ## 置信区间颜色
       threshold.col="red", ## 对角线颜色
       threshold.lty=2,  ## 线型
       cex = 0.8,
       ylab.pos = 2, 
       #cex.axis = 1,
       main = "QQ-plot",
       file.output = F )
dev.off()




