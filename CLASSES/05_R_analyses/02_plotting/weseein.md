---
title: "Class_exercise"
author: "Alexzanah Griffin"
date: "3/22/2023"
output: github_document
always_allow_html: true
---


# Load the libraries you need
# working directory = "../CLASSES/05_R_analyses/02_plotting
# Load functions you need "my_class_functions"
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(GenomicRanges)
source("../../../util/my_class_functions.R")
source("../../../util/plotting_functions.R")
source("../../../util/_setup.R")
```


# load in your peak files for each replicate of each protein
# Here I am starting to analyze my data for my proteins of interest:
# First I will read in each replicate file
```{r load in peak files}
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/atg/Analysis_of_genome_2023"
peak_path <- "CLASSES/03_Nextflow/00_my_chipseq/kurt/results/bwa/mergedLibrary/macs/broadPeak/"
broadpeakfilepath <- file.path(basepath, peak_path)

ARID3A_peak1 <- read_tsv("/scratch/Shares/rinnclass/CLASS_2023/atg/Analysis_of_genome_2023/CLASSES/03_Nextflow/00_my_chipseq/kurt/results/bwa/mergedLibrary/macs/broadPeak/ARID3A_R1_peaks.broadPeak", col_names = F)
ARID3A_peak2 <- read_tsv("/scratch/Shares/rinnclass/CLASS_2023/atg/Analysis_of_genome_2023/CLASSES/03_Nextflow/00_my_chipseq/kurt/results/bwa/mergedLibrary/macs/broadPeak/ARID3A_R2_peaks.broadPeak", col_names = F)

#renaming columns of peak files
names(ARID3A_peak1) <- c('chromosome', 'start', 'end', 'name', 'score', 'strand', 
                       'signalValue', 'pValue', 'qValue')
names(ARID3A_peak2) <- c('chromosome', 'start', 'end', 'name', 'score', 'strand', 
                       'signalValue', 'pValue', 'qValue')

# printing out a table of the number of peaks in each file (of All proteins):

peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

peak_num <- sapply(peak_list, length) %>% as.data.frame(row.names = T)

names(peak_num) <- c("num_peaks")

peak_num <- peak_num %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")

# Save the file 
write_csv(peak_num, "hw_results/num_peaks_df.csv")

# Peak widths 
ARID3A_peak1 <- ARID3A_peak1 %>%
   mutate(peak_width = end - start)
ARID3A_peak2 <- ARID3A_peak2 %>%
   mutate(peak_width = end - start)
```


# Now I am going to create consensus peaks for each protein
```{r consensus peaks}
consensouspeak_path <- "CLASSES/05_R_analyses/00_consensus_peaks/consensus_peaks"
consensusPeakPath <- file.path(basepath, peak_path)

# Unique protien names
dbps <- unique(sapply(names(peak_list), function(x) {
   unlist(strsplit(x, "_"))[1]
}))

# list of 
consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)

names(consensus_list) <- dbps

# export consensus peaks to results folder
for(i in 1:length(consensus_list)) {
rtracklayer::export(consensus_list[[i]], paste0("hw_results/consensus/", names(consensus_list)[i],"_consensus_peaks.bed")) }

```

# Now I am going to make my consensus peaks compatable with UCSC genome browser
```{r UCSC compatable files}
# list of consensus peak files in data set (from consensus_list)
consensus_file_list <- list.files("/scratch/Shares/rinnclass/CLASS_2023/atg/Analysis_of_genome_2023/CLASSES/05_R_analyses/02_plotting/hw_results/consensus", full.names = T, pattern = ".bed")

#view(consensus_file_list)

peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))
names(peaks) <- dbps

canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")
peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))

# name and export uscs peaks
formatted_peaks <- paste0("hw_results/consensus/uscs_", names(peaks), "_consensus_peaks.bed")
view(formatted_peaks)

for(i in 1:length(peaks)) {
  write.table(peaks[[i]], formatted_peaks[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}

headers <- paste0("track type=bed name=", names(peaks))
headers

formatted_peaks <- paste0("hw_results/consensus/formatted_peaks/", names(peaks), ".bed")
formatted_peaks
# print out consensus peak files in a results/UCSC directory


for(i in 1:length(peaks)) {
  writeLines(headers[[i]], formatted_peaks[[i]])
  write.table(peaks[[i]], formatted_peaks[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}

```

# I am curious if my proteins are transcription factors so I will use the annotations
# in a cell paper I found and see

```{r, error = TRUE}
# get genecode annotations
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
table(gencode_gr$type)

rtracklayer::export(gencode_genes, "hw_results/gene_annotations/gencode_genes.gtf")

# protien codind/ mRNA
mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"] 

rtracklayer::export(mrna_genes, "hw_results/gene_annotations/mrna_genes.gtf")
table(gencode_genes$gene_type)

# lncRNA
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 

rtracklayer::export(lncrna_genes, "hw_results/gene_annotations/lncrna_genes.gtf")

#  mRNA and lncRNA combines
mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
rtracklayer::export(mrna_lncrna_genes, "hw_results/gene_annotations/mrna_lncrna_genes.gtf")


# mRNA n lncRNA annotations as df
lncrna_mrna_genes <- rtracklayer::import("hw_results/gene_annotations/mrna_lncrna_genes.gtf")

lncrna_mrna_genes_df <- lncrna_mrna_genes %>% as.data.frame()


# Annotations from cell ppr
url <- "https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx"

destination_for_url <- "hw_results/TF_annotations.xlsx"

download.file(url, destination_for_url)

# Add annotations to num_df
human_tfs <- readxl::read_excel("hw_results/TF_annotations.xlsx",
                                sheet = 2, skip = 1)

names(human_tfs)[4] <- "is_tf"

length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$Name)))


human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(num_peaks_df$dbp), 1:4]

names(human_tfs) <- c("ensembl_id",
                      "dbp",
                      "dbd",
                      "tf")

num_peaks_df <- merge(num_peaks_df, human_tfs, all.x = T)

dim(num_peaks_df[is.na(num_peaks_df$tf),])

num_peaks_df <- num_peaks_df[,1:12]
write_csv(num_peaks_df, "hw_results/num_peaks_df.csv")


```

```{r set promotor regions}
lncrna_mrna_promoters <- promoters(lncrna_mrna_genes, upstream = 1000, downstream = 1000)

# try 500 next
#width(lncrna_mrna_promoters)

lncrna_mrna_promoters_df <-lncrna_mrna_promoters %>% as.data.frame()
rtracklayer::export(lncrna_mrna_promoters, "hw_results/gene_annotations/lncrna_mrna_promoters.gtf")

# RNA gene ids
lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
table(mrna_lncrna_genes$gene_type)

mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]

# Save me plzzzz... ok :)
write_csv(num_peaks_df, "hw_results/num_peaks_df.csv")

```


# Now I want to compare a protein with a previous analysis 
```{r}

# goto UCSC genome browser and load in a peak file for a given protein
# load in the data for the same protein from the previous analysis
# compare how your consensus peaks are similar or different to previous analyses

# how to include graphics with knit
knitr::include_graphics("/scratch/Shares/rinnclass/CLASS_2023/atg/Analysis_of_genome_2023/CLASSES/05_R_analyses/02_plotting/hw_shot.png")

#Consensus peaks align to the file I made

```


# Now I am going to determine how my peaks for each protein overlap annotations of the genome
# First I will find the overlaps between my consensus peaks with promoters of lncRNA and mRNA promoters

```{r}
# find overlaps of promoters for each protein

#consensus peaks
num_peaks_df <- data.frame("dbp" = names(consensus_list),
                           "num_peaks" = sapply(consensus_list, length))

num_peaks_df$total_peak_length <- sapply(consensus_list, function(x) sum(width(x)))


```

## results: 
#1) What can you determine from these overlaps?
#We can see infer that a protien is a transcripton factor based on if it is abundantly bound to promotor regions


# Now I want to compare the overlaps with lncRNA and mRNA promoters seperately 
```{r overlap rna with promoters}

promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_list, type = "counts")

num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])

num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])


```
## results:
# 1) What is the difference in overlaps between mRNA and lncRNA promoters
# overlaps in mRNS are coding for protiens, but lncRNAs done
# more than half of all promotor peaks 

# Now I am going to test if there is more binding over gene bodies than promoters
# I will seperate lncRNA and mRNA gene bodies to find the overlaps 

```{r overlapping gene bodies}

genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                consensus_list, 
                                                type = "counts")

num_peaks_df$peaks_overlapping_genebody <- 
  rowSums(genebody_peak_counts)

#peaks overlapping rna
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

num_peaks_df$peaks_overlapping_mrna_genebody <- 
  rowSums(genebody_peak_counts[,mrna_gene_ids])

write_csv(num_peaks_df, "hw_results/num_peaks_df.csv")

```
## results: 
# 1) Do my proteins have more overlaps with promoters or genebodies?
# Most of the peaks are bound to the genebody when compared to promoters


# It is nice and all to find overlaps, but I am interested in how many proteins
# bind a specific promoter. I will use my handy "occurence" parameter in 
# " count peaks per feature" 

```{r protiens to a single promoter}
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_list, 
                                               type = "occurrence")

stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

write.table(promoter_peak_occurence, "hw_results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# Make me a dataframe
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))

write_csv(peak_occurence_df, "hw_results/peak_occurence_dataframe.csv")

#dont forget to save the environment

```
## results: I find the max number of proteins on a promoter to be 1


# Now I want to start plotting my results
# First I will see if there is a realtionship between peak number and total DNA covered
```{r peak number to total genome covered}
#call for bigger data set:
big_num_peaks_df <- read_csv("/scratch/Shares/rinnclass/CLASS_2023/data/data/2022_num_peaks_df.csv")

#my data
ggplot(num_peaks_df, aes(x = num_peaks, 
                         y = total_peak_length)) +
  geom_point(color = 'purple') 

#save
ggsave("hw_results/Class_exercise_pics/peak_num_v_total_coverage.png")

![](/hw_results/Class_exercise_pics/are_they_tfs.png)

#big data
ggplot(big_num_peaks_df, aes(x = num_peaks, 
                         y = total_peak_length)) +
  geom_point(color = 'purple') 


ggsave("hw_results/Class_exercise_pics/big_peak_num_v_total_coverage.png")
```

# Now I want to color my plot by wether the protein is a TF or not.
```{r are they TFs?}
#my data
ggplot(num_peaks_df, aes(x = num_peaks, 
                 y = total_peak_length,
                )) +
  facet_wrap(tf ~ .) +
  geom_point() 

ggsave("hw_results/Class_exercise_pics/are_they_tfs.png")

#big data
ggplot(big_num_peaks_df, aes(x = num_peaks, 
                 y = total_peak_length,
                )) +
  facet_wrap(tf ~ .) +
  geom_point() 

ggsave("hw_results/Class_exercise_pics/big_are_they_tfs.png")
```

# I want to make a histogram of the number of peaks for each of my proteins

```{r}
# my data
ggplot(num_peaks_df, aes(x = num_peaks)) +
  geom_histogram(bins = 20)

ggsave("hw_results/Class_exercise_pics/hist_peak_for_one_protien.png")

#big data
ggplot(big_num_peaks_df, aes(x = num_peaks)) +
  geom_histogram(bins = 30)

ggsave("hw_results/Class_exercise_pics/big_hist_peak_for_one_protien.png")

```


# Now I want to facet this by the type of DNA binding domain my protein has.
```{r}
#my data
ggplot(num_peaks_df %>% filter(dbd %in% dbd),
       aes(x = num_peaks, y = total_peak_length )) +
  facet_wrap(tf ~ dbd) + 
  geom_point()

ggsave("hw_results/Class_exercise_pics/bind_domain_type.pdf")

#big data
ggplot(big_num_peaks_df %>% filter(dbd %in% dbd),
       aes(x = num_peaks, y = total_peak_length )) +
  facet_wrap(tf ~ dbd) + 
  geom_point()

ggsave("hw_results/Class_exercise_pics/big_bind_domain_type.pdf")

```


