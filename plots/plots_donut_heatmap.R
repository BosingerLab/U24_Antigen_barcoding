library(stringr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)

setwd("/data/runs/Analysis_BLab/UM1/Analysis/Overlap/scPB_scMemory/Clones_MB/add_pb_cloneid/sorted")

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    plot.title=element_text(size=14, face="bold"),
    legend.position="none"
  )


files <- list.files(path = ".","*.txt")
files

#############################################################################################################################

# Clonal frequencies

lapply(files, function(x) {
  print(x)
  
  # x <- files[1]
  d <- read.table(x,header = T, sep ="\t")
  db <- as.data.frame(table(d$CloneID_MB))
  names(db) <- c("Clone","Freq")
  db$AnimalID <- str_split_fixed(db$Clone,"_",2)[,1]
  db$Total <- sum(db$Freq)
  db$Percentage <- (db$Freq/db$Total)*100
  db <- db[order(db$Percentage,decreasing = T),]
  head(db)
  db$Per <- round(db$Percentage,digits=1)
  db$Label <-  paste0(db$Per,"%")
  # db$Label[db$Per<1] <- ""
  db$Label[11:dim(db)[1]] <- ""
  
  out <- paste0(x,".png")
  out <- paste0("donut_charts/",out)
  
  write.table(db,paste0("donut_charts/",x,"_clonal_freq.txt"),sep="\t", quote = F, row.names = F)
  
  p <- ggplot(db,aes(x="",y=Percentage,fill=Percentage))+
    geom_bar(stat="identity", color = "black",width = 0.4) +
    blank_theme +
    coord_polar('y')+ theme(axis.ticks=element_blank(),
                            axis.text.y=element_blank(),
                            axis.text.x=element_text(colour='white', size = 0),
                            axis.title=element_blank()) +
    annotate("text",x=0.5,y=0.5,label=dim(db)[1],size=8) +
    scale_fill_gradient(low="white",high="steelblue4") +
    theme(plot.margin=unit(c(0,0,0,0),"mm"), plot.title = element_text(hjust=0.5, vjust=0.5, size=20),axis.text=element_text(size=9)) +
    geom_label_repel(aes(label= Label, y = cumsum(db$Percentage) - db$Percentage / 2), size = 5, direction = "y",fill = "white", hjust = 1, max.overlaps = 10, nudge_x = 0.3) +
    ggtitle(gsub("./","",strsplit(x,"_")[[1]][1]))
  p
  ggsave(plot=p, filename = out)
  
})



#############################################################################################################################

# V gene usage

combined_db <- do.call(rbind, lapply(files, function(x) {
  
  print(x)
  # x <- files[1]
  d <- read.table(x, header = T, sep = "\t")
  
  db <- as.data.frame(table(gsub("\\*.*", "", d$v_call)))
  names(db) <- c("V_gene", "Freq")
  db$AnimalID <- gsub("_clones.*", "", x)
  db$Total <- sum(db$Freq)
  db$Percentage <- (db$Freq / db$Total) * 100
  db <- db[order(db$Percentage, decreasing = T),]
  
  write.table(db, paste0("V_gene_usage/", x, "_vgeneusage.txt"), sep = "\t", quote = F, row.names = F)
  
  return(db)  # Return the db for each file
}))


library(reshape2)

# Reshape the combined_db so that each row is an animal and each column is a V gene
combined_db_wide <- dcast(combined_db, AnimalID ~ V_gene, value.var = "Percentage", fill = 0)

# Check the reshaped data
head(combined_db_wide)

# Calculate the total frequency of each V gene across all animals
v_gene_totals <- colSums(combined_db_wide[,-1])

# Order the V genes from most frequent to least
v_gene_order <- names(sort(v_gene_totals))

# Now, create a heatmap using ggplot2
heatmap_data <- melt(combined_db_wide, id.vars = "AnimalID", variable.name = "V_gene", value.name = "Percentage")

# Convert V_gene to a factor and order it by total frequency
heatmap_data$V_gene <- factor(heatmap_data$V_gene, levels = v_gene_order)

# Plot the heatmap
ggplot(heatmap_data, aes(x = AnimalID, y = V_gene, fill = Percentage)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of V Gene Frequency per Animal", x = "V Gene", y = "Animal ID", fill = "Frequency")



#############################################################################################################################

