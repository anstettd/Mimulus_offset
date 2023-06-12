#setwd("/scratch/general/vast/u1123911/Mimulus")

library(dplyr)
library(rtracklayer)
library(tidyverse)
library(Biostrings)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")
BiocManager::install("limma")
BiocManager::install("GO.db")

library(limma)
library(GO.db)

#Read GFF3
Reference.curatedggf<- as.data.frame(readGFF("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/CE10.with_functional_annotation.sorted.gff3"))
Interest_SNPs <- read_csv("data/genomic_data/snp_set_env.csv")

genes_only <- filter(Reference.curatedggf,type == "gene")

GFF3.gr <- GRanges(seqnames = genes_only$seqid,
              ranges = IRanges(start = genes_only$start, end = genes_only$end),
              strand = "*",
              mcols = genes_only[,c("ID","IPR","GO","Pfam")])


# Separate 'chr_snp' column into 'chr' and 'snp' in the 'Interest_SNPs' data frame
Interest_SNPs_separated <- Interest_SNPs %>%
  separate(chr_snp, into = c("chr", "snp"), sep = "_(?=[^_]+$)", extra = "merge")
Interest_SNPs_separated$snp <- as.numeric(Interest_SNPs_separated$snp)

# Convert your SNP data frame to a GRanges object
snps <- GRanges(seqnames = Interest_SNPs_separated$chr, 
                ranges = IRanges(start = Interest_SNPs_separated$snp, end = Interest_SNPs_separated$snp))

# Find overlaps between genes and SNPs

overlaps <- findOverlaps(GFF3.gr, snps)

# Extract the overlapping elements from GFF3.gr and convert to a data frame
overlapping_elements <- as.data.frame(GFF3.gr[queryHits(overlaps)])

# Extract the overlapping SNPs and convert to a data frame
overlapping_snps <- as.data.frame(snps[subjectHits(overlaps)])

# Add an 'id' column to each data frame to enable merging
overlapping_elements$id <- seq_len(nrow(overlapping_elements))
overlapping_snps$id <- seq_len(nrow(overlapping_snps))

# Merge the data frames
Overlapping_SNPs_GFF3.df <- full_join(overlapping_elements, overlapping_snps, by = "id")

Overlapping_SNPs_GFF3.df$mcols.GO <- ifelse(Overlapping_SNPs_GFF3.df$mcols.GO=="GO:_",NA,
                                            Overlapping_SNPs_GFF3.df$mcols.GO)
Overlapping_SNPs_GFF3.df$mcols.GO <- ifelse(Overlapping_SNPs_GFF3.df$mcols.IPR=="_:_",NA,
                                            Overlapping_SNPs_GFF3.df$mcols.GO)


###Go terms Function and application to the dataframe

# Define the function
get_go_definitions <- function(go_terms) {
  # Use a regular expression to match the GO terms
  matches <- regmatches(go_terms, gregexpr("GO:\\d{7}", go_terms))
  
  # Retrieve the definitions for each GO term
  go_definitions <- sapply(unlist(matches), function(go_term) {
    # Retrieve the definition
    go_definition <- tryCatch(Term(go_term), error = function(e) return(NA))
    # Combine the GO term with its definition
    paste0(go_term, " (", go_definition, ")")
  })
  
  # Combine the GO terms and their definitions back into a single string
  go_definitions <- paste0(go_definitions, collapse = " | ")
  
  return(go_definitions)
}
# Apply the function to the 'mcols.GO' column of the dataframe and save the results in a new column
Overlapping_SNPs_GFF3.df$mcols.GO.definition <- sapply(Overlapping_SNPs_GFF3.df$mcols.GO, get_go_definitions)
Overlapping_SNPs_GFF3.df$interest_snp <- paste0(Overlapping_SNPs_GFF3.df$seqnames.y,"_",Overlapping_SNPs_GFF3.df$start.y)


Overlapping_SNPs_GFF3.df_curated <- 
  Overlapping_SNPs_GFF3.df[,c("interest_snp","mcols.ID","seqnames.x","start.x","end.x","mcols.GO.definition","mcols.IPR","mcols.Pfam")]

colnames(Overlapping_SNPs_GFF3.df_curated) <- c("SNP","Gene","Chr","Start","End","GO_definition","IPR","Pfam") 

write_csv(Overlapping_SNPs_GFF3.df_curated,"data/genomic_data/GO_overlaps.csv")

####614 snps list
colnames(Interest_SNPs) <- "SNP"
Interest_SNPs_full <- merge(Interest_SNPs,Overlapping_SNPs_GFF3.df_curated,by = "SNP", all.x=TRUE)

write_csv(Interest_SNPs_full,"data/genomic_data/Full_SNP_list_with_GO_overlaps.csv")
####Summary Tables

Genes_Table<- as.data.frame(sort(table(Overlapping_SNPs_GFF3.df_curated$Gene), decreasing = TRUE))
write_csv(Genes_Table,"data/genomic_data/Genes_Table.csv")

GO_Table<- as.data.frame(table(Overlapping_SNPs_GFF3.df_curated$Gene,
                                    Overlapping_SNPs_GFF3.df_curated$GO_definition))
GO_Table <- arrange(filter(GO_Table,Freq > 0), -Freq)

write_csv(Genes_Table,"data/genomic_data/GO_Genes_Table.csv")
