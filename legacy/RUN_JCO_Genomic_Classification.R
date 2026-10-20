######################################
##########                    ########
#### . . . - Timing paper - . . . ####
##########                    ########
######################################
#### . . . - Packages
library(dndscv)

#### . . . - Functions
source("~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/JCO_GenomicClassification_hg38.R")

######################################
##########                    ########
#### . . . - MMcart samples - . . . ####
##########                    ########
######################################

#### . . . - Metadata creation
cnv=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/CART_all_CNV.txt",stringsAsFactors = F)
cnv$cnv=1
colnames(cnv)[grep("ample|id",colnames(cnv),ignore.case = TRUE)]="sample"
cnv.1=unique(cnv[,c("sample","cnv")])

snv=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/all_snv_vaf.txt",
               stringsAsFactors = F)  

# snv$id=str_split_fixed(snv$sample,"_S_",2)[,2]
# snv$SAMPLE.tumor=paste(str_split_fixed(snv$id,"_",7)[,1],
#                        str_split_fixed(snv$id,"_",7)[,2],
#                        str_split_fixed(snv$id,"_",7)[,3],sep="_")

snv$snv=1
colnames(snv)[grep("ample|id",colnames(snv),ignore.case = TRUE)]="sample"
snv.1=unique(snv[,c("sample","snv")])
#snv.1$sample=gsub("IID_","",snv.1$sample)

clinic=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/all_clinical.txt",stringsAsFactors = F)
#clinic=clinic[complete.cases(clinic$ISABL_ID),]
clinic$clinic=1
colnames(clinic)[grep("ample|id",colnames(clinic),ignore.case = TRUE)]="sample"
clinic.1=unique(clinic[,c("sample","clinic")])
clinic.1$sample=gsub("IID_","",clinic.1$sample)

metadata=merge(cnv.1,snv.1,by="sample",all=T)
metadata=merge(metadata,clinic.1,by="sample",all=T)

for(i in 2:length(colnames(metadata))){
  metadata[,i]=ifelse(is.na(metadata[,i]),0,metadata[,i])
}

metadata$Genomic=rowSums(metadata[,2:3])
metadata$Genomic=ifelse(metadata$Genomic==2,"Yes","No")


write.table(metadata,
            "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/TEST_metadata.txt",
            col.names = T,
            row.names = F,
            quote = F,
            sep="\t")


#### . . . - Importing the CNV data - . . . ####
mm.cnv=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/CART_all_CNV.txt",stringsAsFactors = F)


#### . . . - Creating the input file for the JCO Genomic Classification - . . . ####
mm.cnv.1=inputCNV(mm.cnv,
                    SAMPLE.COL = "sample",
                    CHR.COL = "chr",
                    START.COL = "startpos",
                    END.COL = "endpos",
                    TOTAL.COL = "tot",
                    MINOR.COL = "min")

#dkrd.cnv.1=dkrd.cnv.1[,c("sample","chrom","start","end","total","minor")]

#### . . . - Creating the focal & large matrix - . . . ####
FocalLargeMatrix(CNV=mm.cnv.1,
                 CYTOBAND.FILE = "~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/zReferences/hg38/cytoBand_hg38.txt",
                 GISTIC.FILE = "~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/zReferences/hg38/Genomic_Coordinates_GISTICpeaks_hg38.txt",
                 DATABASE = "MMcart",
                 METADATA="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/TEST_metadata.txt",
                 FOCAL.LENGTH = 5000000,
                 OUT.DIR = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/")

#### . . . - HY identification - . . . ####
HYidentification(CNV=mm.cnv.1,
                 CYTOBAND.FILE = "~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/zReferences/hg38/cytoBand_hg38.txt",
                 ARMS = "~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/zReferences/hg38/arm_hg38.txt",
                 CHROM.LENGHT = "~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/zReferences/hg38/chrom_sizes_hg38.txt",
                 GISTIC.ARMS = "~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/zReferences/hg38/arm_level_hg38.txt",
                 DATABASE = "MMcart",
                 OUT.DIR = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/",
                 PERCENTAGE=0.6)

#### . . . - Cytogenetic matrix generation - . . . ####
CytogeneticMatrixGenerator(FocalLargeMATRIX="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/CNV_MMcartmatrix_focal_large_5Mb.txt",
                           HY.FILE="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/60%_ARMS_MMcarthyperdiploid.txt",
                           DATABASE="MMcart",
                           OUT.DIR="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/")



#### . . . - Mutations - . . . ####
snv=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/all_snv_vaf.txt",
               stringsAsFactors = F)
mutations=snv[,c("sample","chrom","pos","ref","alt")]
colnames(mutations)=c("sampleID","chr","pos","ref","mut")
dndsout = dndscv(mutations,
                 refdb = "~/Google Drive/My Drive/Bachisio/Desktop_OCT2023/Tools/R_functions/dnds/RefCDS_human_GRCh38_GencodeV18_recommended.rda",
                 cv=NULL)

driver=dndsout$sel_cv
driver.pos=driver[driver$qglobal_cv<0.1,]
driver.pos$gene_nam

annot.mut.all=dndsout$annotmuts
#annot.mut.pos=dn$annotmuts[dn$annotmuts$gene %in% driver.pos$gene_nam,]
annot.mut.all.ns=dndsout$annotmuts[dndsout$annotmuts$impact!="Synonymous",]

annot.mut.all.ns$sampleID=gsub("IID_","",annot.mut.all.ns$sampleID)

write.table(annot.mut.all.ns,
            "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/MMcart_dnds_ns_mutations.txt",
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")

MutationsMatrixGenerator(CODING.MUT = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/MMcart_dnds_ns_mutations.txt",
                         METADATA = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/TEST_metadata.txt",
                         UCSC.gene = "~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/zReferences/hg38/UCSC_hg38_ref_genes.dms",
                         DRIVERs = "~/Google Drive/My Drive/Bachisio/Documents/Script/2024-2025/R.FUNCTIONS/SHARED/MMGenomicClassification/zReferences/hg38/Driver_summary_hg38.txt",
                         GisticOtherAllPeaks = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/gistic.other_All_peaksMMcart.txt",
                         CNV = mm.cnv.1,
                         DATABASE = "MMcart",
                         OUT.DIR = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/")

snv$sample=gsub("IID_","",snv$sample)
snv=unique(snv)
APOBECanalyzer(SNV.OBJECT = snv,
               ITER = 10,
               SAMPLE = "sample",
               CHROM = "chrom",
               POS = "pos",
               REF = "ref",
               ALT = "alt",
               DATABASE = "MMcart",
               OUT.DIR = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/")


#### . . . - Final Matrix generator - . . . ####
FinalMatrixGenerator(INPUT.FILE.FOLDER="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/",
                     DATABASE="MMcart",
                     OUTPUT="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/MM_CART/2024_final_calls/JCOClassification/FinalMatrix_MMcart.txt")


translocations=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/PMMM.dkrd.txt",
                          stringsAsFactors = F)
translocations=translocations[,c("sample","t_CCND1","t_MAF","t_MMSET")]

chromothripsis=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/sv_to_annotation_dkrd.051921.txt",
                          stringsAsFactors = F)
chromothripsis=chromothripsis[,c("sample","PCAWG_class")]
chromothripsis=chromothripsis[chromothripsis$PCAWG_class=="chromothripsis",]
chromothripsis=unique(chromothripsis)
chromothripsis$sample=gsub("IID_","",chromothripsis$sample)

tra.chrom=merge(translocations,chromothripsis,by="sample",all.x=T)
tra.chrom$PCAWG_class=ifelse(is.na(tra.chrom$PCAWG_class),0,1)

colnames(tra.chrom)[5]="Chromothripsis"

matrix=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/FinalMatrix_DKRd.txt",
                  stringsAsFactors = F)
matrix.dkrd=merge(matrix,tra.chrom,by="sample",all.x=T)


#### . . . - DIS3 hotspot - . . . ####
matrix.dkrd.1=dis3_hs(MATRIX = matrix.dkrd,
                      CODING.MUT = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/DKRd_dnds_ns_mutations.txt")

matrix.dkrd.def=merge(matrix.dkrd,matrix.dkrd.1,by="sample")

colnames(matrix.dkrd.def)[131]="CNV.Sig"

for(i in 2:ncol(matrix.dkrd.def)){
  matrix.dkrd.def[,i]=ifelse(is.na(matrix.dkrd.def[,i]),0,matrix.dkrd.def[,i])
}

#### . . . - Classification & Heatmap - . . . ####
dkrd.class=MMGenoCla(input=matrix.dkrd.def)
y
y

dkrd.class$annotation$clusters.new=factor(dkrd.class$annotation$clusters.new,levels = c("HRD_RAS","HRD_Gains","HRD_Complex_Cytogenetic",
                                                                                        "CCND1_Simple","CCND1_Complex_Cytogenetic",
                                                                                        "MAF_and_or_HyperAPOBEC",
                                                                                        "NSD2_GainAmp1q_Del13q","NSD2_HRD","NSD2_Del13q",
                                                                                        "Multiple_Losses","Simple"))
space_heat=as.numeric(table(dkrd.class$annotation))

space_heat2=c(space_heat[1],
              space_heat[1] + space_heat[2],
              space_heat[1] + space_heat[2] + space_heat[3],
              space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4],
              space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5],
              space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6],
              space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7],
              space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8],
              space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8] + space_heat[9],
              space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8] + space_heat[9] + space_heat[10],
              space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8] + space_heat[9] + space_heat[10] + space_heat[11]
)

space.row=c(8,
            8 + 1,
            8 + 1 + 7,
            8 + 1 + 7 + 50,
            8 + 1 + 7 + 50 + 3,
            8 + 1 + 7 + 50 + 3 + 22,
            8 + 1 + 7 + 50 + 3 + 22 + 3,
            8 + 1 + 7 + 50 + 3 + 22 + 3 + 1,
            8 + 1 + 7 + 50 + 3 + 22 + 3 + 1 + 1)

pheatmap(as.matrix(t(dkrd.class$genomic)), annotation_col= dkrd.class$annotation,
         annotation_colors = dkrd.class$colors,
         cluster_cols = FALSE, show_colnames = F, fontsize_row = 5,
         cluster_rows =  FALSE, border_color = NA, legend = F,
         col=c("grey85","salmon1","darkred"),
         gaps_col = space_heat2,
         gaps_row = space.row
)

write.table(dkrd.class$genomic,
            "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/output/DKRd_matrix.txt",
            col.names = T,
            row.names = T,
            quote = F,
            sep="\t")
write.table(dkrd.class$annotation,
            "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/output/DKRd_GenomicClassification.txt",
            col.names = T,
            row.names = T,
            quote = F,
            sep="\t")






# ########################################
# ##########                    ##########
# #### . . . - German samples - . . . ####
# ##########                    ##########
# ########################################
# 
# #### . . . - Metadata creation
# cnv=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/Pre-Post-gain Project/german/cnv_german_hd6_all.txt",stringsAsFactors = F)
# cnv$cnv=1
# colnames(cnv)[grep("ample|id",colnames(cnv),ignore.case = TRUE)]="sample"
# cnv.1=unique(cnv[,c("sample","cnv")])
# 
# snv=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/Pre-Post-gain Project/MolTime/german_data_MolTime/mutindel.all.final.corrected.txt",
#                stringsAsFactors = F)  
# snv$snv=1
# colnames(snv)[grep("ample|id",colnames(snv),ignore.case = TRUE)]="sample"
# snv.1=unique(snv[,c("sample","snv")])
# snv.1$sample=gsub("IID_","",snv.1$sample)
# 
# clinic=read.delim("~/Desktop/Pre-Post-gain Project/master.file.new.txt",stringsAsFactors = F)
# clinic=clinic[clinic$cohort=="German",]
# clinic$clinic=1
# colnames(clinic)[grep("ample|id",colnames(clinic),ignore.case = TRUE)]="sample"
# clinic.1=unique(clinic[,c("sample","clinic")])
# 
# metadata=merge(cnv.1,snv.1,by="sample",all=T)
# metadata=merge(metadata,clinic.1,by="sample",all=T)
# 
# for(i in 2:length(colnames(metadata))){
#   metadata[,i]=ifelse(is.na(metadata[,i]),0,metadata[,i])
# }
# 
# metadata$Genomic=rowSums(metadata[,2:3])
# metadata$Genomic=ifelse(metadata$Genomic==2,"Yes","No")
# 
# metadata[metadata$Genomic=="No",]
# 
# metadata=metadata[!metadata$sample %in% c("K08K-1XZHPL_T1","K08K-MPT52W_T1","P044_BL","P346_BL","P505_BL","K08K-6KLSUS_BL","K08K-KXRUXD_BL",
#                                           "K08K-VHPB5B_BL","K08K-XWDDRB_BL"),]
# 
# write.table(metadata,
#             "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/German_metadata.txt",
#             col.names = T,
#             row.names = F,
#             quote = F,
#             sep="\t")
# 
# #### . . . - Importing the CNV data
# german.cnv=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/Pre-Post-gain Project/german/cnv_german_hd6_all.txt",
#                       stringsAsFactors = F)
# german.cnv$total=german.cnv$nMajor+german.cnv$nMinor
# length(unique(german.cnv$sample)) # 36
# 
# #### . . . - Creating the input file for the JCO Genomic Classification
# german.cnv.1=inputCNV(german.cnv,
#                       SAMPLE.COL = "sample",
#                       CHR.COL = "chr",
#                       START.COL = "start",
#                       END.COL = "end",
#                       TOTAL.COL = "total",
#                       MINOR.COL = "nMinor")
# 
# german.cnv.1=german.cnv.1[,c("sample","chrom","start","end","total","minor")]
# 
# #### . . . - Creating the focal & large matrix
# FocalLargeMatrix(CNV=german.cnv.1,
#                  CYTOBAND.FILE = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/Reference/cytoBand.txt",
#                  GISTIC.FILE = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/Reference/Genomic_Coordinates_GISTICpeaks.txt",
#                  DATABASE = "German_NDMM",
#                  FOCAL.LENGTH = 5000000,
#                  OUT.DIR = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/")
# 
# #### . . . - HY identification
# HYidentification(CNV=german.cnv.1,
#                  CYTOBAND.FILE = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/Reference/cytoBand.txt",
#                  ARMS = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/Reference/arm_hg19.txt",
#                  CHROM.LENGHT = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/Reference/hg19.chrom_sizes.txt",
#                  GISTIC.ARMS = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/Reference/arm_level.txt",
#                  DATABASE = "German_NDMM",
#                  OUT.DIR = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/",
#                  PERCENTAGE=0.6)
# 
# #### . . . - Cytogenetic matrix generation
# CytogeneticMatrixGenerator(FocalLargeMATRIX="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/CNV_German_NDMMmatrix_focal_large_5Mb.txt",
#                            HY.FILE="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/60%_ARMS_German_NDMMhyperdiploid.txt",
#                            DATABASE="German_NDMM",
#                            OUT.DIR="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/")
# 
# #### . . . - Mutations
# snv.german=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/Pre-Post-gain Project/MolTime/german_data_MolTime/mutindel.all.final.corrected.txt",
#                       stringsAsFactors = F) 
# length(unique(snv.german$sample)) # 365
# snv.german$sample[snv.german$sample %in% "K08K-1XZHPL_T1"]="K08K-1XZHPL_BL"
# snv.german$sample[snv.german$sample %in% "K08K-KXRUXD_T1"]="K08K-KXRUXD_BL"
# snv.german$sample[snv.german$sample %in% "K08K-6KLSUS_T1"]="K08K-6KLSUS_BL"
# snv.german$sample[snv.german$sample %in% "K08K-MPT52W_T1"]="K08K-MPT52W_BL"
# snv.german$sample[snv.german$sample %in% "K08K-VHPB5B_T1"]="K08K-VHPB5B_BL"
# snv.german$sample[snv.german$sample %in% "K08K-XWDDRB_T1"]="K08K-XWDDRB_BL"
# 
# 
# mutations.german=snv.german[,c("sample","chrom","pos","ref","alt")]
# colnames(mutations.german)=c("sampleID","chr","pos","ref","mut")
# dndsout.german = dndscv(mutations.german,
#                         refdb = "hg19")
# 
# driver.german=dndsout.german$sel_cv
# driver.pos.german=driver.german[driver.german$qglobal_cv<0.1,]
# driver.pos.german$gene_nam
# 
# annot.mut.all.german=dndsout.german$annotmuts
# #annot.mut.pos=dn$annotmuts[dn$annotmuts$gene %in% driver.pos$gene_nam,]
# annot.mut.all.ns.german=dndsout.german$annotmuts[dndsout.german$annotmuts$impact!="Synonymous",]
# 
# write.table(annot.mut.all.ns.german,
#             "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/German_NDMM_dnds_ns_mutations.txt",
#             col.names = T,
#             row.names = F,
#             quote = F,
#             sep = "\t")
# 
# MutationsMatrixGenerator(CODING.MUT = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/German_NDMM_dnds_ns_mutations.txt",
#                          METADATA = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/German_metadata.txt",
#                          UCSC.HG19 = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/Reference/UCSC_hg19_ref_genes.dms",
#                          DRIVERs = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/Reference/Driver_summary.txt",
#                          GisticOtherAllPeaks = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/gistic.other_All_peaksGerman_NDMM.txt",
#                          CNV = german.cnv.1,
#                          DATABASE = "German_NDMM",
#                          OUT.DIR = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/")
# 
# APOBECanalyzer(SNV.OBJECT = snv.german,
#                ITER = 10,
#                SAMPLE = "sample",
#                CHROM = "chrom",
#                POS = "pos",
#                REF = "ref",
#                ALT = "alt",
#                DATABASE = "German_NDMM",
#                OUT.DIR = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/")
# 
# 
# #
# FinalMatrixGenerator(INPUT.FILE.FOLDER="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/",
#                      DATABASE="German_NDMM",
#                      OUTPUT="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/FinalMatrix_German_NDMM.txt")
# 
# 
# translocations.german=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/bach_ID_final_clock.txt",
#                                  stringsAsFactors = F)
# translocations.german=translocations.german[translocations.german$group=="IGH_canon",]
# translocations.german=translocations.german[,c("sample","IGH_partner")]
# 
# translocations.german=as.data.frame.matrix(table(translocations.german$sample,translocations.german$IGH_partner))
# translocations.german$sample=rownames(translocations.german)
# translocations.german$MAF.new=rowSums(translocations.german[,c(3,4)])
# 
# translocations.german=translocations.german[,c("sample","CCND1","MAF.new","NSD2")]
# colnames(translocations.german)=c("sample","t_CCND1","t_MAF","t_MMSET")

# translocations.german.rr=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/rrmm.pm.matrix.txt",
#                                  stringsAsFactors = F)
# translocations.german.rr=translocations.german.rr[,c("sample","CCND1.IGH","MAF.IGH","NSD2.IGH")]
# colnames(translocations.german.rr)=c("sample","t_CCND1","t_MAF","t_MMSET")
# 
# translocations.german.1=rbind(translocations.german,translocations.german.rr)

# chromothripsis.german=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/Germans/NDMM_after_QC/UM/cnv_signature_prediction_value_GERMAN.txt",
#                                  stringsAsFactors = F)
# chromothripsis.german=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/missing_samples_chromothripsis.txt",
#                                  stringsAsFactors = F)
# #colnames(chromothripsis.german)[2]="Chromothripsis"
# chromothripsis.german=chromothripsis.german[,c("sample","Chromothripsis")]
# 
# tra.chrom.german=merge(translocations.german,chromothripsis.german,by="sample",all.y=T)
# 
# for(i in 2:4){
#   tra.chrom.german[,i]=ifelse(is.na(tra.chrom.german[,i]),0,tra.chrom.german[,i])
# }
# tra.chrom.german=tra.chrom.german[tra.chrom.german$Chromothripsis %in% c(0,1),]
# 
# matrix.german=read.delim("~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/FinalMatrix_German_NDMM.txt",
#                          stringsAsFactors = F)
# matrix.german=merge(matrix.german,tra.chrom.german,by="sample",all.x=T)
# dim(matrix.german) # 371

#table(matrix.german$Chromothripsis,useNA="ifany")
#miss.chromo.samp=matrix.german[,c("sample","Chromothripsis")]

#write.table(miss.chromo.samp,
#            "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/output/missing_samples_chromothripsis.txt",
#            col.names = T,
#            row.names = F,
#            quote = F,
#            sep="\t")

# 
# 
# 
# #### . . . - German
# matrix.german.1=dis3_hs(MATRIX = matrix.german,
#                         CODING.MUT = "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/matrix/German_NDMM_dnds_ns_mutations.txt")
# 
# matrix.german.def=merge(matrix.german,matrix.german.1,by="sample")
# 
# colnames(matrix.german.def)[131]="CNV.Sig"
# 
# for(i in 2:ncol(matrix.german.def)){
#   matrix.german.def[,i]=ifelse(is.na(matrix.german.def[,i]),0,matrix.german.def[,i])
#   matrix.german.def[,i]=as.numeric(as.character(matrix.german.def[,i]))
# }
# 
# german.class=MMGenoCla(input=matrix.german.def)
# y
# y
# 
# german.class$annotation$clusters.new=factor(german.class$annotation$clusters.new,levels = c("HRD_RAS","HRD_Gains","HRD_Complex_Cytogenetic",
#                                                                                             "CCND1_Simple","CCND1_Complex_Cytogenetic",
#                                                                                             "MAF_and_or_HyperAPOBEC",
#                                                                                             "NSD2_GainAmp1q_Del13q","NSD2_HRD","NSD2_Del13q","GainAmp1q_Del13q",
#                                                                                             "Multiple_Losses","Simple"))
# 
# space_heat=as.numeric(table(german.class$annotation))
# 
# 
# 
# space_heat2=c(space_heat[1],
#               space_heat[1] + space_heat[2],
#               space_heat[1] + space_heat[2] + space_heat[3],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8] + space_heat[9],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8] + space_heat[9] + space_heat[10],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8] + space_heat[9] + space_heat[10] + space_heat[11],
#               space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + space_heat[6] + space_heat[7] + space_heat[8] + space_heat[9] + space_heat[10] + space_heat[11] + space_heat[12]
# )
# 
# space.row=c(8,
#             8 + 1,
#             8 + 1 + 7,
#             8 + 1 + 7 + 50,
#             8 + 1 + 7 + 50 + 3,
#             8 + 1 + 7 + 50 + 3 + 22,
#             8 + 1 + 7 + 50 + 3 + 22 + 3,
#             8 + 1 + 7 + 50 + 3 + 22 + 3 + 1,
#             8 + 1 + 7 + 50 + 3 + 22 + 3 + 1 + 1)
# 
# pheatmap(as.matrix(t(german.class$genomic)), annotation_col= german.class$annotation,
#          annotation_colors = german.class$colors,
#          cluster_cols = FALSE, show_colnames = F, fontsize_row = 5,
#          cluster_rows =  FALSE, border_color = NA, legend = F,
#          col=c("grey85","salmon1","darkred"),
#          gaps_col = space_heat2,
#          gaps_row = space.row
# )
# 
# write.table(german.class$genomic,
#             "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/output/German_matrix.txt",
#             col.names = T,
#             row.names = T,
#             quote = F,
#             sep="\t")
# write.table(german.class$annotation,
#             "~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Pre_Post_gain/NatGen_NDMM_timing_Rebuttal/output/German_GenomicClassification.txt",
#             col.names = T,
#             row.names = T,
#             quote = F,
#             sep="\t")
# 
# 
