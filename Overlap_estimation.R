# ------------------------------------------------------- #
# ------------------------------------------------------- #
#  OTUs overlap between gut microbiota in vertebrates     #
# ------------------------------------------------------- #
# ------------------------------------------------------- #

# ----------------------#
# Script description    #
# ----------------------#

# 1. Objective
# 2. Data used
# 3. Load R packages
# 4. Input parameters
# 5. Select ASVs (to remove mitochondria and chloroplasts)
# 6. Load and merge OTU tables
# 7. Select samples : one (random) sample per host species 
# 8. Rarefy final table
# 9. Compute beta diversity and summarize results


# ------------- #
# 1. Objective  #
# ------------- #

# This script computes the overlap of OTUs between the gut microbiota of individuals 
# of different vertebrate genera 

# ------------- #
# 2. Data used  #
# --------------#

# Data set : 16S (V4) amplicon sequencing reads publicly available from the Song et al. study. 
# Associated paper title: Song, S. J., Sanders, J. G., Delsuc, F., Metcalf, J., Amato, K., Taylor, M. W., ... & Knight, R. (2020). Comparative analyses of vertebrate gut microbiomes reveal convergence between birds and bats. MBio, 11(1).
# Associated paper DOI: https://doi.org/10.1128/mBio.02901-19
# Description. This data set consists of the Earth Microbiome Project Release #1 with additional Vertebrate Microbiota samples
# The data consist of processed 16S (V4) amplicon sequencing reads. Reads were trimmed to 90pb and processed using Deblur as described in the SOng et al. methods. Several studies are pooled together and deblur has been run independently for each study. 

# Data (16S reads) location: Can be found on the Qiita website: https://qiita.ucsd.edu/study/description/11166.
# Metadata location: Can be found on the website of the publisher of the Song et al. paper as a Supplementary material: DATA SET S1: https://mbio.asm.org/content/mbio/11/1/e02901-19/DC1/embed/inline-supplementary-material-1.xlsx?download=true

# Downloading data 
# All biom tables can be downloaded directly from Qiita, the corresponding deblur table IDs are given below ("qiita_OTU_table_ID" vector) and in the method section. 


# ---------------------- #
#  3. Load R packages    #
# ---------------------- #

library(rbiom) # to load OTU table
library(tidyverse) # for tidy data manipulation 
library(vegan) # to compute Jaccard beta-diversity metric
library(janitor) # To produce summary tables
library(dada2) # to assign taxonomy (RDP classifier + SILVA 132 database)

# ------------------------ #
#  4. Input parameters     #
# ------------------------ #

# directories 
dir_path <- "/Users/fmazel/Data/Mammalian_microbiome_EMP/study_11166_030421-045713/" # path of the working directory where data has been dowloaded
directory_OTU_tables <- paste0(dir_path,"BIOM/") # path of the biom objects 

# Make resuls reproducible 
set.seed(99)

# Qiita IDs for the processed table used here
# We focus on all studies with the Song et al. dataset 
# and select the deblur OTU table using the 90pb reads
qiita_OTU_table_ID <- c("93862","93855","93819","93914","93900","93846","93851","94483") 

rarefaction_depth = 5000 # rarefaction depth

beta_metric = "jaccard" # beta-diversity metric used to quantify overlap between samples (used in the function vegan::vegdist)

# ----------------- #
#  5. Select ASVs   #
# ----------------- #

DNA_sequences <- list()
for (i in qiita_OTU_table_ID){
  DNA_sequences[[i]] <- read.fasta(paste0(directory_OTU_tables,i,"/reference-hit.seqs.fa"))
}
DNA_sequences <- unique(unlist(DNA_sequences))

OTU_Taxonomy <- assignTaxonomy(seqs=DNA_sequences, 
                               refFasta="/Users/fmazel/Data/SILVA/dada2_formated/silva_nr_v132_train_set.fa", 
                               multithread=3,
                               taxLevels =c("Kingdom", "Order", "Family"))


OTU_noMito_noChloro <- OTU_taxonomy[,c(1:3)] %>% 
  as_tibble() %>% 
  mutate(seqs = rownames(OTU_taxonomy)) %>% 
  subset(Kingdom == "Bacteria") %>% 
  subset(!Order=="Chloroplast") %>% 
  subset(!Family=="Mitochondria") %>% 
  pull(seqs)


# ------------------------------ #
#  6. Load and merge OTU tables  #
# ------------------------------ #

OTU_tables <- list()
for (i in qiita_OTU_table_ID){
  OTU_tables[[i]] <- read.biom(paste0(directory_OTU_tables,i,"/reference-hit.biom"))
}

# Retrieve sample names in OTU tables 
samplesID <- unlist(lapply(OTU_tables,sample.names))

# Merge tables
OTU_table <- as.matrix(OTU_tables[[1]]$counts)

for (i in qiita_OTU_table_ID[2:length(qiita_OTU_table_ID)]){
  OTU_table <- merge(OTU_table,as.matrix(OTU_tables[[i]]$counts),
                     by = "row.names",
                     all=TRUE,
                     suffixes = c("","duplicated")) # choose one sample from the one that have duplicated names between studies 
  row.names(OTU_table) <- OTU_table$Row.names
  OTU_table <- OTU_table[,-1]
} # Merge all OTUS table into one single table 

OTU_table[is.na(OTU_table)] = 0

# ------------------------------------------------ #
# 7. Select samples : one sample per host species  # 
# ------------------------------------------------ #
meta <- readxl::read_xlsx(paste0(dir_path,"inline-supplementary-material-1.xlsx"))

subset_meta <- meta %>% 
  subset(SampleID %in% samplesID) %>% # keep samples present in the OTU table
  subset(!is.na(Taxonomy_Species)) %>% # remove samples without data on host species 
  group_by(Taxonomy_Species) %>% 
  sample_n(1) # randomly sub-sample one sample per host species 

subset_samples <- subset_meta %>% 
  pull(SampleID)

# select OTUs with reads and that are not mitochondria nor chloroplasts
final_otus <- (rownames(OTU_table) %in% OTU_noMito_noChloro) & (rowSums(OTU_table)>0)

# select samples that were selected above 
final_samples <- intersect(subset_samples,colnames(OTU_table))

OTU_table <- OTU_table[final_otus, final_samples]

# ----------------------- #
#  8.  Rarefy final table
# ----------------------- #

rarefied_OTU_table <- vegan::rrarefy(x = t(OTU_table),
                                     sample=rarefaction_depth)

rarefied_OTU_table <- rarefied_OTU_table[rowSums(rarefied_OTU_table)==rarefaction_depth,] # remove sample with less that "rarefaction_depth" reads
dim(rarefied_OTU_table)
table(rowSums(rarefied_OTU_table))
sum(colSums(rarefied_OTU_table)>0)

# ------------------------- #
# 9. Compute beta diversity #
# ------------------------- #

Jaccard <- as.matrix(vegan::vegdist(rarefied_OTU_table,method = beta_metric,binary = T))


# Summrise overlap values
TaxonmyHost <- subset_meta %>% 
  subset(SampleID %in% rownames(rarefied_OTU_table)) %>% 
  select(SampleID, Taxonomy_Species, Taxonomy_Genus,Taxonomy_Order,Taxonomy_Class)

Beta_div_reformated <- Jaccard %>% 
  reshape2::melt(value.name = "Jaccard") %>% 
  left_join(TaxonmyHost, by = c("Var1"="SampleID")) %>%
  left_join(TaxonmyHost, by = c("Var2"="SampleID")) %>% 
  subset(!Var1==Var2) %>%  #remove self comparisons
  mutate(comparisonGenus=ifelse(Taxonomy_Genus.x==Taxonomy_Genus.y,"WithinGenus","BetweenGenus")) %>% 
  mutate(comparisonClass=ifelse(Taxonomy_Class.x==Taxonomy_Class.y,"WithinClass","BetweenClass"))

summaryBeta = Beta_div_reformated %>% 
  group_by(comparisonGenus,comparisonClass) %>% 
  summarise(Mean_Jaccard = 100*(1-round(mean(Jaccard),4)),
            Median_Jaccard = 100*(1-round(median(Jaccard),4)))


Order_summary <- TaxonmyHost %>% 
  group_by(Taxonomy_Order) %>% 
  summarise(n_species=n(),
            n_genus=n_distinct(Taxonomy_Genus),
            Class=unique(Taxonomy_Class)) %>% 
  adorn_totals("row")

write.csv(x=Order_summary, file = "Summary_Host_Taxonomy.csv")
write.csv(x=summaryBeta,file = "SummaryJaccard.csv")

