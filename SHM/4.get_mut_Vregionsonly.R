library(alakazam)
library(shazam)
library(dplyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(ggbeeswarm)
library(tidyverse)

install.packages("tidyverse")

db <- read.table("All.DB_db-pass.tab", header = T, sep="\t")
dim(db)

db_mut <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT",
                                germlineColumn="GERMLINE_IMGT",
                                regionDefinition=IMGT_V_BY_SEGMENTS,
                                frequency = TRUE,
                                combine = TRUE,
                                nproc=8)

write.table(file = "SHM_V_region_beforeCDR3_allcolumns.txt",db_mut,sep="\t", quote = F, row.names = F)


mut <- db_mut[,c("SEQUENCE_ID","LOCUS","mu_freq")]

mut <- mut %>%
    separate("SEQUENCE_ID", into = c("AnimalID", "Contig_id"), sep = "_", extra = "merge", remove = F)

mut$CellID <- paste0(mut$AnimalID,"_",gsub("_contig.*","",mut$Contig_id))

head(mut)



cloneid <- read.table("Clone_CellID", header = F,sep="\t")
names(cloneid) <- c("CloneID","ContigID")
cloneid$CellID <- paste0(gsub("_\\d+","",cloneid$CloneID),"_",gsub("_contig.*","",cloneid$ContigID))
head(cloneid)

dim(mut)
dim(cloneid)
mut_cloneid <- merge(mut,cloneid, by = "CellID")
dim(mut_cloneid)
head(mut_cloneid)

mut_avg <- mut_cloneid %>%
  group_by(CloneID, LOCUS) %>%
  summarise(avg_mu_freq = mean(mu_freq, na.rm = TRUE), .groups = "drop")

mut_avg$AnimalID <- gsub("_\\d+","",mut_avg$CloneID)

head(mut_avg)

write.table(file = "SHM_Vregion_beforeCDR3_cells.txt",mut_cloneid,sep="\t", quote = F)

write.table(file = "SHM_V_region_beforeCDR3_avgclone.txt",mut_avg,sep="\t", quote = F)


# 275-25 - #705AA6
# RHe16 - #5182C3
# RIr15 - #5DBF89
# RLk15 - #243066
# RUp16 - #004243
# Rys15 - #F073AC

#mut_cloneid <- read.table("SHM_Vregion_beforeCDR3_cells.txt")
#mut_avg <- read.table("SHM_V_region_beforeCDR3_avgclone.txt")


p1 <- ggplot(mut_cloneid,aes(x=AnimalID,y=mu_freq)) + geom_violin(aes(fill=AnimalID)) + 
  scale_fill_manual(values = c("#705AA6","#5182C3","#5DBF89","#243066","#004243","#F073AC")) +
  geom_boxplot(width = 0.1) +
  facet_wrap(~LOCUS) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),  # Rotate x-axis labels
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),  # Increase facet label size
        strip.background = element_rect(fill = "white", color = "black"), 
        axis.title = element_text(size = 14)) +
  labs(
    x = "Animal ID", 
    y = "Mutation Frequency\n(V region - upto FWR3)", 
  )


p2 <- ggplot(mut_avg,aes(x=AnimalID,y=avg_mu_freq)) + geom_violin(aes(fill=AnimalID)) + 
  scale_fill_manual(values = c("#705AA6","#5182C3","#5DBF89","#243066","#004243","#F073AC")) +
  geom_boxplot(width = 0.1) +
  facet_wrap(~LOCUS) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),  # Rotate x-axis labels
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),  # Increase facet label size
        strip.background = element_rect(fill = "white", color = "black"), 
        axis.title = element_text(size = 14)) +
  labs(
    x = "Animal ID", 
    y = "Average mutation frequency per clone\n(V region - upto FWR3)", 
  )


p1/p2
