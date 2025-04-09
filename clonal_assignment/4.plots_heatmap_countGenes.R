library(stringr)
library(ggplot2)
library(patchwork)
library(alakazam)
library(reshape2)
library(dplyr)
library(tidyr)


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    plot.title=element_text(size=14, face="bold"),
    legend.position="none"
  )


files <- list.files(path = ".","*.txt")
files


# Combined all changeO files

combined_data <- data.frame()

get_combined <- function(f){
  d <- read.table(f,header= T, sep="\t")
  d$Sample <- gsub("_clones.*","",f)
  combined_data <<- rbind(combined_data, d)
}

files

lapply(files,get_combined)

names(combined_data)

# Get counts for genes
count <- countGenes(combined_data,gene="v_call", mode="gene", groups= "Sample",clone = "CloneID_MB", remove_na = F)
write.table(count,"V_gene_counts.txt",quote = F, sep = "\t", row.names = F)


# Get order of V genes for plotting
gene_order <- count %>%
  group_by(gene) %>%
  summarise(total_clone_count = sum(clone_freq)) %>%
  arrange(desc(total_clone_count)) %>%
  pull(gene)

count$gene <- factor(count$gene, levels = rev(gene_order))

# Get a complete dataframe 
count_complete <- count %>%
  dplyr::ungroup() %>%
  tidyr::complete(Sample, gene, fill = list(clone_freq = NA))

# Plot the heatmap
ggplot(count_complete, aes(x = Sample, y = gene, fill = clone_freq)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue",na.value = "white") +
  theme_minimal() +
  theme(panel.grid.major=element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Heatmap of V Gene Frequency (clonal)", x = "V Gene", y = "Animal ID", fill = "Frequency")
