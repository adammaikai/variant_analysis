# RNA Editing Functions
library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(ggpubr)
library(maftools)
options(shiny.maxRequestSize=1000*1024^2)

ggViolin <- function(data){
  ggviolin(data, x="Condition", y="VAF", fill = "Condition") + 
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)+ theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggBox <- function(data){
  ggplot(data, aes(x=Variant_Classification, y=VAF, fill=Condition)) + geom_boxplot() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

topGenes <- function(data){
  maf.in <- read.maf(data)
  oncoplot(maf.in, showTumorSampleBarcodes = TRUE)
}

VariantsPerGene <- function(data){
  VariantsPerGene <- do.call(rbind, lapply(unique(data$Condition), function(i)
    cbind(Condition=i, data.frame(table(subset(data, Condition==i)$Hugo_Symbol)))))
  VariantsPerGene$Condition <- factor(VariantsPerGene$Condition, levels=levels(data$Condition))
  ggplot(VariantsPerGene, aes(x=log(Freq), fill=Condition)) + geom_density(alpha=0.7) + theme_bw() + ylab("Log Variants Per Gene")
}

MutatedGenesPerSampleBox <- function(data){
  MutatedGenesPerSample <- do.call(rbind, lapply(unique(data$Condition), function(i)
    cbind(Condition=i, aggregate(Hugo_Symbol ~ Tumor_Sample_Barcode, data=subset(data, Condition==i), length))))
  MutatedGenesPerSample$Condition <- factor(MutatedGenesPerSample$Condition, levels=levels(data$Condition))
  ggplot(MutatedGenesPerSample, aes(x=Condition, y=Hugo_Symbol, fill=Condition)) + geom_boxplot() + theme_bw() + ylab("Mutated Genes Per Sample")
}


