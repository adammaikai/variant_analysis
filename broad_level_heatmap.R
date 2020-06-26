# Gistic broad-level events
# Hierarchical clustering of GISTIC output
library(reshape2)
library(data.table)
library(plyr)
library(ComplexHeatmap)
library(gplots)
library(circlize)
library(RColorBrewer)

gistic_broad <- read.csv("~/mnt/workstation/mnt/data1/adam/cleveland/20181106_cleveland_WGS_mouse_tumorigenesis_SV/gistic_1mb/broad_significance_results.txt", stringsAsFactors = FALSE, sep = "\t")
gistic_broad <- subset(gistic_broad, Amp.q.value < 0.01 | Del.q.value < 0.01)
gistic_broad <- data.frame(GISTIC2.0=ifelse(gistic_broad$Amp.q.value < 0.2, yes=gistic_broad$Amp.z.score, no=-gistic_broad$Del.z.score), "Chromosome Arm" = gistic_broad$Arm, check.names=FALSE)

sample_broad <- read.csv("~/mnt/workstation/mnt/data1/adam/cleveland/20181106_cleveland_WGS_mouse_tumorigenesis_SV/gistic_1mb/broad_values_by_arm.txt", stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)

data <- merge(sample_broad, gistic_broad, by="Chromosome Arm", all.x=TRUE)
data[is.na(data)] <- 0
chrOrder <- c(paste0(1:19, "q"))
data$`Chromosome Arm` <- factor(data$`Chromosome Arm`, levels=chrOrder)
data <- data[order(data$`Chromosome Arm`),]
row.names(data) <- data$`Chromosome Arm`
data$`Chromosome Arm` <- NULL
rotate <- function(x) t(apply(x, 2, rev))
rot.data <- rotate(rotate(rotate(data)))
rot.data.dist <- dist(rot.data)
df <- melt(as.matrix(rot.data.dist), varnames = c("row", "col"))
similarity <- subset(df, row == "GISTIC2.0" & col != "GISTIC2.0")
similarity$value.transformed <- 1 - (similarity$value - min(similarity$value))/(max(similarity$value) - min(similarity$value))

# Phenotype labels
meta <- read.csv("~/mnt/workstation/mnt/data1/adam/cleveland/20181106_cleveland_WGS_mouse_tumorigenesis_SV/metadata/cleveland_meta_07172019.txt", sep = "\t")
meta <- merge(meta, similarity, by.x="Sample.Name", by.y="col")
meta <- arrange(meta, split)
meta <- rbind.fill(meta, data.frame(Sample.Name="GISTIC2.0", condition="GISTIC2.0", split="GISTIC2.0", value=0))
condition.df <- data.frame(Condition=meta$condition, Similarity=meta$value.transformed, row.names=meta$Sample.Name)
condition.color <- list(Condition = c("-/-"="black", "-/-Doxycycline"="darkgray", "+/-Doxycycline"="lightgray", "GISTIC2.0"="white"))
colOrder <- as.character(meta$Sample.Name)

# Heatmap
set.seed(123)
conditionAnno <- rowAnnotation(df = condition.df, col = condition.color, show_annotation_name = FALSE)
hm <- Heatmap(rot.data[colOrder,],
              right_annotation = conditionAnno,
              row_split = meta$split,
              cluster_rows = function(x) hclust(dist(x, method="euclidean"), method="ward.D2"),
              cluster_row_slices = FALSE,
              cluster_columns = FALSE, show_column_names = TRUE,
              show_row_names = TRUE, row_names_side = "right",
              col = colorRamp2(c(-1.5, 0, 1.5), c("blue3", "white", "red3")),
              name = "")
draw(hm)


