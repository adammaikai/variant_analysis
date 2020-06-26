# RNA Edit Known and Novel Combined Figure
library(data.table)
library(maftools)

rna.meta <- read.csv("~/mnt/jam_db/projects/holm/rnameta_04122019.txt", stringsAsFactors = FALSE, sep = "\t")
comb.meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/metadata/DNA_RNA_combined_meta12052018.txt", stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")
comb.meta$Condition <- factor(as.character(comb.meta$Condition),
                              levels=c("Normal", "ET", "PV", "LR_MF_PostET", "Int_1_MF", "Int_2_MF",
                                       "Int_2_MF_PostET", "Int_2_MF_PostPV", "HR_MF", "HR_MF_PostPV",
                                       "sAML", "AML", "AP_CML", "BC_CML", "CP_CML", "CML", "Diseased"))
comb.meta$Treatment_type <- factor(as.character(comb.meta$Treatment_type),
                                   levels=c("None", "Untreated", "Jak2 Inhibitor", "SHH treated", "hydroxyurea", "TKI", "vidaza"))
meta <- merge(rna.meta, comb.meta, by=c("Sample", "Tumor_Sample_Barcode"))
meta <- subset(meta, Cell_type.x %in% c("Stem", "Progenitor"))

maf <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/latest_mpn_combined/MPN_Normal_known_novel_missense_rnaedit_sites_04122019.maf", stringsAsFactors = FALSE, sep = "\t")
maf <- unique(maf)
maf <- subset(maf, !(grepl("Novel", Novelty) & Hugo_Symbol=="MT-ND5"))
maf[grepl("Known", maf$Novelty) & !grepl("Non", maf$Novelty)]$Variant_Classification <- "Reported Alu"
maf[!grepl("Known", maf$Novelty) & !grepl("Non", maf$Novelty)]$Variant_Classification <- "Unreported Alu"
maf[grepl("Known", maf$Novelty) & grepl("Non", maf$Novelty)]$Variant_Classification <- "Reported Nonalu"
maf[!grepl("Known", maf$Novelty) & grepl("Non", maf$Novelty)]$Variant_Classification <- "Unreported Nonalu"

# oncoplot
maf <- subset(maf, Cell_type %in% c("Progenitor", "Stem"))
maf.in <- read.maf(maf, clinicalData = meta, vc_nonSyn = c("Reported Alu", "Unreported Alu", "Reported Nonalu", "Unreported Nonalu"))
Sample_type.colors <- RColorBrewer::brewer.pal(n = 7, name = 'Paired')
names(Sample_type.colors) <- c("Aged_normal_bone_marrow", "Young_normal_bone_marrow", "MF", "PV", "ET", "CML", "AML")
Treatment_type.colors <- RColorBrewer::brewer.pal(n = 7, name = 'Accent')
names(Treatment_type.colors) <- c("None", "Untreated", "hydroxyurea", "Jak2_Inhibitor", "TKI", "vidaza", "SHH_treated")
cell_type.colors <- gray.colors(2)
names(cell_type.colors) <- c("Stem", "Progenitor")
colors <- list(Sample_type = Sample_type.colors, Treatment_type=Treatment_type.colors, Cell_type.x=cell_type.colors)
vc.colors <- RColorBrewer::brewer.pal(n = 4, name = 'Set1')
names(vc.colors) <- c("Reported Alu", "Unreported Alu", "Reported Nonalu", "Unreported Nonalu")
dev.off(dev.list()["RStudioGD"])
oncoplot(maf = maf.in, clinicalFeatures = c("Sample_type", "Treatment_type", "Cell_type.x"), sortByAnnotation = TRUE, 
         showTumorSampleBarcodes = TRUE, annotationColor = colors, top = 25, fontSize = 6, 
         colors = vc.colors, annotationOrder = c("Aged_normal_bone_marrow", "Young_normal_bone_marrow", "MF", "PV", "ET", "CML", "AML"),
        bgCol="white", removeNonMutated = FALSE
         )
pdf("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures_update/science_figures_04082019/RNAEdit.known.novel.missense.combined.pdf")
