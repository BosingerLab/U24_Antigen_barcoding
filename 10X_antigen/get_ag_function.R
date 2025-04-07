library("Seurat")

setwd("/data/runs/Analysis_BLab/UM1/Analysis/antigen")

ids = c("RIr15","RLk15","RYs15","RHe16","275-12","RUp16")


vdj_path <-  "/data/runs/Analysis_BLab/UM1/Analysis/antigen/vdj_only/all/IgBLAST_AIRR/productive/filtered_paired"
ag_path <- "/data/runs/Analysis_BLab/UM1/Analysis/antigen/all_antigen"

list_ag <- list()
list_plots_norm <- list()
list_plots_raw <- list()

df_counts <- data.frame()

id <- "RIr15"

get_ag <- function(id){
  
  print(id)

  vdj_file <- paste0(id,".filtered_paired")
  antigen_file <- paste0(id,"_sample_feature_bc_matrix.h5")
  
  vdj <- read.table(file.path(vdj_path,vdj_file),sep="\t", header = T)
  head(vdj)
  vdj$Cell
  
  df_counts[id,"VDJ"] <<- nrow(vdj)

  
  counts <- Read10X_h5(file.path(ag_path,antigen_file))
  ag <- CreateSeuratObject(counts=(counts+1),assay = "HTO")
  
  df_counts[id,"Antigen"] <<- nrow(ag@meta.data)
  
  joint.bcs <- intersect(colnames(ag), vdj$Cell)
  
  ag.vdj <- ag[, joint.bcs]
  
  df_counts[id,"Antigen and VDJ"] <<- nrow(ag.vdj@meta.data)
  
  print(dim(ag.vdj))
  
  raw_counts <- t(as.data.frame(GetAssayData(ag.vdj,assay = "HTO", slot="counts")))
  head(raw_counts)
  raw_counts <- as.data.frame(raw_counts)
  head(raw_counts)
  colnames(raw_counts) <- paste0(colnames(raw_counts),"_raw_counts")
  ag.vdj <- AddMetaData(ag.vdj,raw_counts[as.character(colnames(ag.vdj)),])
  
  ag.vdj <- NormalizeData(ag.vdj, assay = "HTO", normalization.method = "CLR", margin = 2)
  norm_counts <- t(as.data.frame(GetAssayData(ag.vdj,assay = "HTO", slot="data")))
  norm_counts <- as.data.frame(norm_counts)
  colnames(norm_counts) <- paste0(colnames(norm_counts),"_norm_counts_2")
  ag.vdj <- AddMetaData(ag.vdj,norm_counts[as.character(colnames(ag.vdj)),])
  
  head(ag.vdj@meta.data)
  
  p1 <- FeatureScatter(ag.vdj,feature1 = "AF647.STV.BG505_raw_counts",feature2 = "AF546.STV.Biotin_raw_counts")
  p2 <- FeatureScatter(ag.vdj,feature1 = "AF488.STV.BG505_raw_counts",feature2 = "AF546.STV.Biotin_raw_counts")
  p3 <- FeatureScatter(ag.vdj,feature1 = "AF647.STV.BG505_raw_counts",feature2 = "AF488.STV.BG505_raw_counts") 
  
  out <- paste0(id,"_featurescatterraw.png")
  list_plots_raw[[id]] <<- p1+p2+p3 & NoLegend()
  ggsave(p1+p2+p3,file=out, width = 12, height= 4, units = "in")
  
  q_biotin <- quantile(ag.vdj$AF546.STV.Biotin_norm_counts_2, probs = 0.97)
  q_af647 <- quantile(ag.vdj$AF647.STV.BG505_norm_counts_2, probs = 0.03)
  q_af488 <- quantile(ag.vdj$AF488.STV.BG505_norm_counts_2, probs = 0.03)
  
  p1 <- FeatureScatter(ag.vdj,feature1 = "AF647-STV-BG505",feature2 = "AF546-STV-Biotin") + geom_vline(xintercept = q_af647, linetype='dashed', col = 'black') + geom_hline(yintercept = q_biotin, linetype='dashed', col = 'black') + labs(subtitle = id)
  p2 <- FeatureScatter(ag.vdj,feature1 = "AF488-STV-BG505",feature2 = "AF546-STV-Biotin") + geom_vline(xintercept = q_af488, linetype='dashed', col = 'black') + geom_hline(yintercept = q_biotin, linetype='dashed', col = 'black')
  p3 <- FeatureScatter(ag.vdj,feature1 = "AF647-STV-BG505",feature2 = "AF488-STV-BG505")  + geom_vline(xintercept = q_af647, linetype='dashed', col = 'black') + geom_hline(yintercept = q_af488, linetype='dashed', col = 'black')
  
  out <- paste0(id,"_featurescatter.png")
  list_plots_norm[[id]] <<- p1+p2+p3 & NoLegend() 
  ggsave(p1+p2+p3,file=out, width = 12, height= 4, units = "in")
  
  ag.vdj$Biotin_positive <- "No" 
  ag.vdj$AF647_positive <- "No" 
  ag.vdj$AF488_positive <- "No"
  
  ag.vdj$Biotin_positive[ag.vdj$AF546.STV.Biotin_norm_counts_2 > q_biotin] <- "Yes"
  ag.vdj$AF647_positive[ag.vdj$AF647.STV.BG505_norm_counts_2 > q_af647] <- "Yes"
  ag.vdj$AF488_positive[ag.vdj$AF488.STV.BG505_norm_counts_2 > q_af488] <- "Yes"
  
  print(table(ag.vdj$AF647_positive))
  print(table(ag.vdj$AF488_positive))
  print(table(ag.vdj$Biotin_positive))
  
  
  print(table(ag.vdj$AF647_positive,ag.vdj$AF488_positive))
  print(table(ag.vdj$AF647_positive,ag.vdj$Biotin_positive))
  print(table(ag.vdj$AF488_positive,ag.vdj$Biotin_positive))
  
  write.table(ag.vdj@meta.data,paste0(id,"ag_metadata.txt"),quote = F,sep="\t")
  
  row.names(vdj) <- vdj$Cell
  ag.vdj_meta <- merge(vdj,ag.vdj@meta.data, by = 0)
  ag.vdj_meta$Row.names <- NULL
  write.table(ag.vdj_meta,paste0(id,"_Ag_VDJ_combined.txt"),quote = F,sep="\t",row.names = F)
  
  df_counts[id,"Biotin positive"] <<- as.numeric(table(ag.vdj$Biotin_positive)[2])
  df_counts[id,"AF647-STV-BG505 positive"] <<- as.numeric(table(ag.vdj$AF647_positive)[2])
  df_counts[id,"AF488-STV-BG505 positive"] <<- as.numeric(table(ag.vdj$AF488_positive)[2])
  
  keep <- row.names(ag.vdj@meta.data[ag.vdj$Biotin_positive == "No" & ag.vdj$AF488_positive == "Yes"  & ag.vdj$AF647_positive == "Yes",])
  print(dim(keep))
  vdj2 <- vdj[vdj$Cell %in% keep,]
  print(dim(vdj2))
  
  df_counts[id,"BG505 double positive"] <<- nrow(vdj2)
  
  write.table(vdj2,paste0(id,"_VDJ_agspecific.txt"),quote = F,sep="\t",row.names = F)
  
  list_ag[[id]] <<- ag.vdj
}

lapply(ids,get_ag)

p <- wrap_plots(list_plots_norm[c("RUp16","RHe16","RIr15","275-12","RLk15","RYs15")], ncol = 2)
ggsave(plot = p, "FeatureScatter_ag_barcode.pdf", height = 8, width = 16, units = "in")

write.xlsx(df_counts,"MB_Counts.xlsx")
