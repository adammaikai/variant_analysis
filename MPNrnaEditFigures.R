library(data.table)
library(maftools)
library(ggplot2)
library(ggpubr)
library(plyr)

# Combined Meta Data
comb.meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/metadata/DNA_RNA_combined_meta12052018.txt", stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")
comb.meta$Condition <- factor(as.character(comb.meta$Condition),
                              levels=c("Normal", "ET", "PV", "LR_MF_PostET", "Int_1_MF", "Int_2_MF",
                                       "Int_2_MF_PostET", "Int_2_MF_PostPV", "HR_MF", "HR_MF_PostPV",
                                       "sAML", "AML", "AP_CML", "BC_CML", "CP_CML", "CML", "Diseased"))
comb.meta$Treatment_type <- factor(as.character(comb.meta$Treatment_type),
                                   levels=c("None", "Untreated", "Jak2 Inhibitor", "SHH treated", "hydroxyurea", "TKI", "vidaza"))
# DNA Meta
meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/metadata/somatic_metadata_clean_12062018.txt", sep = "\t", stringsAsFactors = FALSE)
# RNA Meta 
rna.meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/metadata/RNAseq_with_controls_meta_20181130.txt", stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)

# ADAR p150 quantile data
load("~/mnt/workstation/mnt/data1/tomw/Holm_Jamieson_Analysis/ADARp150stemexp.rdata")
load("~/mnt/workstation/mnt/data1/tomw/Holm_Jamieson_Analysis/ADARp150progexp.rdata")
adar <- rbind(ADAR.p150.Stem.quartile.df, ADAR.p150.Prog.quartile.df)
adar <- subset(adar, !(grepl(".BM.", adar$Sample)))
adar$Sample <- gsub("\\.", "-", gsub("X|_ACAGTG|_GTGAAA|_S[0-9]|.PB", "", adar$Sample))
names(adar) <- c("Sample", "p150")

#mpn.maf <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/latest_mpn_combined/MPN_Normal_known_novel_merged_rnaedit_sites_04162019.maf", stringsAsFactors = FALSE, sep = "\t")
mpn.maf <- fread("~/mnt/workstation/mnt/data1/adam/jamieson/aws/rnaedit/mpn_combined/MPN_Normal_known_novel_merged_rnaedit_sites_04162019.maf", stringsAsFactors = FALSE, sep = "\t")
rnaedit.known.maf <- subset(mpn.maf, grepl("Known", Novelty) & Condition_code2 %in% c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "AML"))
rnaedit.known.maf[grepl("Mutation|Start|Splice", Variant_Classification)]$Variant_Classification <- "Missense_Mutation"
rnaedit.known.maf$Condition_code <- factor(rnaedit.known.maf$Condition_code, levels=c("Normal", "PV", "ET", 
                                                              "myelofibrosis", "high_risk_myelofibrosis", "int2_myelofibrosis", "int2_myelofibrosis_postET", 
                                                              "int2_myelofibrosis_postPV", "myelofibrosis_fromMDS", "myelofibrosis_postET", "myelofibrosis_postPV", 
                                                              "systemic_mastocytosis_myelofibrosis", "CML", "AML"))
rnaedit.known.maf$Condition_code2 <- factor(rnaedit.known.maf$Condition_code2, levels=c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "AML"))
#rnaedit.known.maf <- merge(rnaedit.known.maf, unique(comb.meta[c("Sample", "Tumor_Sample_Barcode", "Treatment_type", "Condition")]), by=c("Sample", "Tumor_Sample_Barcode"))
rnaedit.known.maf$VarID <- paste0(rnaedit.known.maf$Chromosome, ":", rnaedit.known.maf$Start_position)
# Figures
prog.all <- subset(rnaedit.known.maf, Cell_type=="Progenitor")
ggplot(prog.all, aes(x=Variant_Classification, y=VAF, fill=Condition_code2)) + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/rnaedit.prog.vaf.bycondition.region.pdf")
grant_subset <- subset(prog.all, Condition_code2 %in% c("Aged_Normal", "Young_Normal", "AML"))
ggviolin(grant_subset, x="Condition_code2", y="VAF", fill = "Condition_code2") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/grant_subset.rnaedit.violin.bycondition_10212019.pdf")
ggplot(grant_subset, aes(color=Condition_code2, x=VAF)) + geom_density(lwd=1.7) + theme_bw()
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/grant_subset.rnaedit.density.bycondition.pdf")
ggviolin(prog.all, x="Condition_code2", y="VAF", fill = "Condition_code2") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/rnaedit.violin.bycondition.pdf")
maf.p150 <- merge(rnaedit.known.maf, adar, by = "Sample")
ggplot(maf.p150, aes(x=Variant_Classification, y=VAF, fill=p150)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/rnaedit.p150byregion.pdf")

# Statistical tests for differences in VAF between phenotypes: Kolmogorov Smirnov Test and T Test
ks.phenotype <- sapply(levels(prog.all$Condition_code2)[2:7], function(i) ks.test(subset(prog.all, Condition_code2=="Aged_Normal")$VAF, subset(prog.all, Condition_code2==i)$VAF, alternative = "greater")$p.value) 
ttest.phenotype <- sapply(levels(prog.all$Condition_code2)[2:7], function(i) t.test(subset(prog.all, Condition_code2=="Aged_Normal")$VAF, subset(prog.all, Condition_code2==i)$VAF, alternative = "less")$p.value)
pv.phenotype.df <- data.frame(Kolmogorov_Smirnov=ks.phenotype, T_test=ttest.phenotype)
write.table(pv.phenotype.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/vaf_phenotype_statistical_signif.tsv", sep="\t", quote=FALSE, row.names=TRUE)

ks.risk <- sapply(levels(prog.all$Condition)[2:15], function(i) ks.test(subset(prog.all, Condition=="Normal")$VAF, subset(prog.all, Condition==i)$VAF, alternative = "greater")$p.value) 
ttest.risk <- sapply(levels(prog.all$Condition)[2:15], function(i) t.test(subset(prog.all, Condition=="Normal")$VAF, subset(prog.all, Condition==i)$VAF, alternative = "less")$p.value)
pv.risk.df <- data.frame(Kolmogorov_Smirnov=ks.risk, T_test=ttest.risk)
write.table(pv.risk.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/vaf_risk_statistical_signif.tsv", sep="\t", quote=FALSE, row.names=TRUE)

ks.treatment <- sapply(levels(prog.all$Treatment_type)[2:7], function(i) ks.test(subset(prog.all, Treatment_type=="None")$VAF, subset(prog.all, Treatment_type==i)$VAF, alternative = "greater")$p.value)
ttest.treatment <- sapply(levels(prog.all$Treatment_type)[2:7], function(i) t.test(subset(prog.all, Treatment_type=="None")$VAF, subset(prog.all, Treatment_type==i)$VAF, alternative = "less")$p.value)
pv.treatment.df <- data.frame(Kolmogorov_Smirnov=ks.treatment, T_test=ttest.treatment)
write.table(pv.treatment.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/vaf_treatment_statistical_signif.tsv", sep="\t", quote=FALSE, row.names=TRUE)

# Stratified by Variant Classification
ks.phenotype.by.vc <- sapply(levels(prog.all$Condition_code2)[2:7], function(i) lapply(unique(prog.all$Variant_Classification), function(j) ks.test(subset(prog.all, Variant_Classification==j & Condition_code2=="Aged_Normal")$VAF, subset(prog.all, Variant_Classification==j & Condition_code2==i)$VAF, alternative = "greater")$p.value)) 
row.names(ks.phenotype.by.vc) <- unique(prog.all$Variant_Classification)
write.table(ks.phenotype.by.vc, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/ks_phenotype_vc_statistical_signif.tsv", sep="\t", quote=FALSE, row.names=TRUE)

# ADAR expression correlation
load("~/mnt/workstation/mnt/data1/adam/jamieson/holm/rdata/ADAR.tx.lcpm.exp.rdata")
mean.vaf <- do.call(rbind, lapply(unique(rnaedit.known.maf$Sample), function(i) data.frame(Sample=i, mean.vaf=mean(subset(rnaedit.known.maf, Sample==i)$VAF))))
adar.exp <- do.call(cbind, ADAR.tx.lcpm.exp)
adar.exp <- as.data.frame(t(as.matrix(adar.exp)))
adar.exp$Sample <- gsub("\\.", "-", gsub("X|_ACAGTG|_GTGAAA|_S[0-9]|.PB|.BM", "", row.names(adar.exp)))
adar.exp <- merge(adar.exp, mean.vaf, by="Sample")
adar.exp <- merge(adar.exp, rna.meta, by="Sample")
adar.exp <- merge(adar.exp, adar, by="Sample")
# Stem=1, circle Progenitor=2, triangle
adar.exp$Cell_type <- as.numeric(factor(adar.exp$Cell_type, levels=c("Stem", "Progenitor")))
adar.exp$Condition_code2 <- factor(adar.exp$Condition_code2, levels=c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "AML"))
# both stem and prog
ggscatter(adar.exp, x="ENST00000368474", y="mean.vaf", color="Condition_code2", size=c(4, 4)[adar.exp$Cell_type], shape=c(15, 17)[adar.exp$Cell_type]) + ylab("Mean RNA editing VAF") + xlab("ADAR p150 log CPM expression") + stat_cor(label.x = .5) + guides(colour = guide_legend(override.aes = list(size=4)))
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/rnaedit.p150.meanvaf.corr.pdf")
ggscatter(adar.exp, x="ENST00000368471", y="mean.vaf", color="Condition_code2", size=c(4, 4)[adar.exp$Cell_type], shape=c(15, 17)[adar.exp$Cell_type]) + ylab("Mean RNA editing VAF") + xlab("ADAR p110 log CPM expression") + stat_cor(aes(color = cyl), label.x = .5) + guides(colour = guide_legend(override.aes = list(size=4)))
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/rnaedit.p110.meanvaf.corr.pdf")
ggplot(adar.exp, aes(p150, mean.vaf)) + geom_boxplot() + theme_classic() + ylab("Mean RNA editing VAF") + xlab("ADAR p150 expression quartile")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/rnaedit.p150.meanvaf.boxplot.pdf")

# CDK13 plots
cdk13.maf <- subset(rnaedit.known.maf, Hugo_Symbol=="CDK13")
cdk13.maf$Hugo_Symbol <- cdk13.maf$Protein_Change
cdk13 <- read.maf(cdk13.maf, clinicalData = rna.meta)
risk.colors <- RColorBrewer::brewer.pal(n = 7,name = "Paired")
names(risk.colors) <- c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "AML", "CML")
risk.order <- as.numeric(factor(names(risk.colors), levels=c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "AML", "CML")))
risk.colors = list(Condition_code2 = risk.colors)
dev.off(dev.list()["RStudioGD"])
oncoplot(maf =cdk13,
         clinicalFeatures = "Condition_code2", sortByAnnotation = TRUE, 
         annotationOrder = c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "AML", "CML"),
         annotationColor = risk.colors, drawColBar = FALSE)

# STAT plots
stat <- fread("~/mnt/workstation/mnt/data1/adam/jamieson/holm/expression/ADAR.tx.lcpm.exp.rdataMPN_STAT3_STAT5B_editingSites_AML_AN.csv", stringsAsFactors = FALSE)
stat$Hugo_Symbol <- paste(stat$Hugo_Symbol, stat$Location)
stat.in <- read.maf(stat, vc_nonSyn = unique(stat$Variant_Classification), clinicalData = rna.meta)
oncoplot(stat.in, clinicalFeatures = "Condition_code2")
stat <- fread("~/mnt/workstation/mnt/data1/adam/jamieson/holm/vcf/rnaedit/STAT3_STAT5B_rnaedit_sites_hg38pos_05292019.maf", stringsAsFactors = FALSE, sep = "\t")
stat <- subset(stat, Condition_code2 %in% c("Aged_Normal", "AML"))
stat$Hugo_Symbol <- paste0("chr", stat$Chromosome, ":", stat$Start_position)
#stat$published <- ifelse(stat$hg38_position %in% c(42318026, 42318064, 42318069, 42318091), yes="published", no="novel")
#stat.known <- subset(stat, published=="published" & Condition_code2 != "MDS")
stat.in <- read.maf(stat, vc_nonSyn = unique(stat$Variant_Classification))
clinicalData <- unique(stat[,c("Tumor_Sample_Barcode", "Condition_code2")])
tryCatch({dev.off(dev.list()["RStudioGD"])}, error=function(e){})
oncoplot(stat.in, annotationDat = clinicalData, clinicalFeatures = "Condition_code2", 
         genes=c("chr17:40470044", "chr17:40470109"), annotationOrder = c("AML", "Aged_Normal"), drawColBar = FALSE,
         sortByAnnotation = TRUE, showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE)

# lollipop
rkm <- read.maf(rnaedit.known.maf, clinicalData = rna.meta)
lollipopPlot(rkm, gene = "CDK13")


