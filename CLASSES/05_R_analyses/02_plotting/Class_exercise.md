Class\_exercise
================
Alexzanah Griffin
3/22/2023

# Load the libraries you need

# working directory = "../CLASSES/05\_R\_analyses/02\_plotting

# Load functions you need “my\_class\_functions”

`{r setup, include=FALSE} knitr::opts_chunk$set(echo = TRUE) library(tidyverse) library(GenomicRanges) library(ggplot2) library(GenomicRanges) source("../../../util/my_class_functions.R") source("../../../util/plotting_functions.R") source("../../../util/_setup.R")`

# load in your peak files for each replicate of each protein

# Here I am starting to analyze my data for my proteins of interest:

# proteinX, Y, Z …..

# First I will read in each replicate file

\`\`\`{r load in peak files} basepath &lt;-
“/scratch/Shares/rinnclass/CLASS\_2023/atg/Analysis\_of\_genome\_2023”
peak\_path &lt;-
“CLASSES/03\_Nextflow/00\_my\_chipseq/kurt/results/bwa/mergedLibrary/macs/broadPeak/”
broadpeakfilepath &lt;- file.path(basepath, peak\_path)

ARID3A\_peak1 &lt;-
read\_tsv(“/scratch/Shares/rinnclass/CLASS\_2023/atg/Analysis\_of\_genome\_2023/CLASSES/03\_Nextflow/00\_my\_chipseq/kurt/results/bwa/mergedLibrary/macs/broadPeak/ARID3A\_R1\_peaks.broadPeak”,
col\_names = F) ARID3A\_peak2 &lt;-
read\_tsv(“/scratch/Shares/rinnclass/CLASS\_2023/atg/Analysis\_of\_genome\_2023/CLASSES/03\_Nextflow/00\_my\_chipseq/kurt/results/bwa/mergedLibrary/macs/broadPeak/ARID3A\_R2\_peaks.broadPeak”,
col\_names = F)

\#renaming columns of peak files names(ARID3A\_peak1) &lt;-
c(‘chromosome’, ‘start’, ‘end’, ‘name’, ‘score’, ‘strand’,
‘signalValue’, ‘pValue’, ‘qValue’) names(ARID3A\_peak2) &lt;-
c(‘chromosome’, ‘start’, ‘end’, ‘name’, ‘score’, ‘strand’,
‘signalValue’, ‘pValue’, ‘qValue’)

# printing out a table of the number of peaks in each file (of All proteins):

peak\_list &lt;- import\_peaks(consensus\_file\_path =
broadpeakfilepath)

peak\_num &lt;- sapply(peak\_list, length) %&gt;%
as.data.frame(row.names = T)

names(peak\_num) &lt;- c(“num\_peaks”)

peak\_num &lt;- peak\_num %&gt;% rownames\_to\_column(var = “dbp”)
%&gt;% separate(col = dbp, into = c(‘dbp’, ‘replicate’), sep = "\_")

# Save the file

write\_csv(peak\_num, “hw\_results/num\_peaks\_df.csv”)

# Peak widths

ARID3A\_peak1 &lt;- ARID3A\_peak1 %&gt;% mutate(peak\_width = end -
start) ARID3A\_peak2 &lt;- ARID3A\_peak2 %&gt;% mutate(peak\_width = end
- start)



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

# Now I am going to make my consensus peaks compatable with UCSC genome browser

``` {r}
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

\`\`\`{r get annotations} \# get genecode annotations gencode\_gr &lt;-
rtracklayer::import(“/scratch/Shares/rinnclass/CLASS\_2023/data/data/genomes/gencode.v32.annotation.gtf”)

gencode\_genes &lt;- gencode\_gr\[gencode\_gr\$type == “gene”\]
table(gencode\_gr\$type)

rtracklayer::export(gencode\_genes,
“hw\_results/gene\_annotations/gencode\_genes.gtf”)

# protien codind/ mRNA

mrna\_genes &lt;- gencode\_genes\[gencode\_genes\$gene\_type %in%
“protein\_coding”\]

rtracklayer::export(mrna\_genes,
“hw\_results/gene\_annotations/mrna\_genes.gtf”)
table(gencode\_genes\$gene\_type)

# lncRNA

lncrna\_genes &lt;- gencode\_genes\[gencode\_genes\$gene\_type %in%
“lncRNA”\]

rtracklayer::export(lncrna\_genes,
“hw\_results/gene\_annotations/lncrna\_genes.gtf”)

# mRNA and lncRNA combines

mrna\_lncrna\_genes &lt;- gencode\_genes\[gencode\_genes\$gene\_type
%in% c(“protein\_coding”,“lncRNA”)\]
rtracklayer::export(mrna\_lncrna\_genes,
“hw\_results/gene\_annotations/mrna\_lncrna\_genes.gtf”)

# mRNA n lncRNA annotations as df

lncrna\_mrna\_genes &lt;-
rtracklayer::import(“hw\_results/gene\_annotations/mrna\_lncrna\_genes.gtf”)

lncrna\_mrna\_genes\_df &lt;- lncrna\_mrna\_genes %&gt;% as.data.frame()

# if you leave the object name you just created in the environment

# it will print out in the knit. For example :

# Annotations from cell ppr

url &lt;-
“<https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx>”

destination\_for\_url &lt;- “hw\_results/TF\_annotations.xlsx”

download.file(url, destination\_for\_url)

# Add annotations to num\_df

human\_tfs &lt;- readxl::read\_excel(“hw\_results/TF\_annotations.xlsx”,
sheet = 2, skip = 1)

names(human\_tfs)\[4\] &lt;- “is\_tf”

length(which(tolower(num\_peaks\_df$dbp) %in% tolower(human_tfs$Name)))

human\_tfs &lt;-
human\_tfs\[tolower(human\_tfs$Name) %in% tolower(num_peaks_df$dbp),
1:4\]

names(human\_tfs) &lt;- c(“ensembl\_id”, “dbp”, “dbd”, “tf”)

num\_peaks\_df &lt;- merge(num\_peaks\_df, human\_tfs, all.x = T)

dim(num\_peaks\_df\[is.na(num\_peaks\_df\$tf),\])

num\_peaks\_df &lt;- num\_peaks\_df\[,1:12\] write\_csv(num\_peaks\_df,
“hw\_results/num\_peaks\_df.csv”)


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
    write_csv(num_peaks_df, "results/num_peaks_df.csv")

# Now I want to compare a protein with a previous analysis

``` {r}
# goto UCSC genome browser and load in a peak file for a given protein
# load in the data for the same protein from the previous analysis
# compare how your consensus peaks are similar or different to previous analyses

# how to include graphics with knit
# knitr::include_graphics("picture name")
![UCSC comparison] ("hw_results/Class_exercise_pics/hw_shot.png")

```

# Now I am going to determine how my peaks for each protein overlap annotations of the genome

# First I will find the overlaps between my consensus peaks with promoters of lncRNA and mRNA promoters

``` {r}
# find overlaps of promoters for each protein

#consensus peaks
num_peaks_df <- data.frame("dbp" = names(consensus_list),
                           "num_peaks" = sapply(consensus_list, length))

num_peaks_df$total_peak_length <- sapply(consensus_list, function(x) sum(width(x)))

```

## results:

\#1) What can you determine from these overlaps? \#We can see infer that
a protien is a transcripton factor based on if it is abundantly bound to
promotor regions

# Now I want to compare the overlaps with lncRNA and mRNA promoters seperately

\`\`\`{r overlap rna with promoters}

promoter\_peak\_counts &lt;-
count\_peaks\_per\_feature(lncrna\_mrna\_promoters, consensus\_list,
type = “counts”)

num\_peaks\_df\$peaks\_overlapping\_promoters &lt;-
rowSums(promoter\_peak\_counts)

num\_peaks\_df\$peaks\_overlapping\_lncrna\_promoters &lt;-
rowSums(promoter\_peak\_counts\[,lncrna\_gene\_ids\])

num\_peaks\_df\$peaks\_overlapping\_mrna\_promoters &lt;-
rowSums(promoter\_peak\_counts\[,mrna\_gene\_ids\])

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

## results:

# 1) Do my proteins have more overlaps with promoters or genebodies?

# Most of the peaks are bound to the genebody when compared to promoters

# It is nice and all to find overlaps, but I am interested in how many proteins

# bind a specific promoter. I will use my handy “occurence” parameter in

# " count peaks per feature"

\`\`\`{r protiens to a single promoter} promoter\_peak\_occurence &lt;-
count\_peaks\_per\_feature(lncrna\_mrna\_promoters, consensus\_list,
type = “occurrence”)

stopifnot(all(colnames(promoter\_peak\_occurence) ==
lncrna\_mrna\_promoters\$gene\_id))

write.table(promoter\_peak\_occurence,
“hw\_results/lncrna\_mrna\_promoter\_peak\_occurence\_matrix.tsv”)

# Make me a dataframe

peak\_occurence\_df &lt;- data.frame(“gene\_id” =
colnames(promoter\_peak\_occurence), “gene\_name” =
lncrna\_mrna\_promoters$gene_name,  "gene_type" = lncrna_mrna_promoters$gene\_type,
“chr” = <lncrna_mrna_promoters@seqnames>,  
“1kb\_up\_tss\_start” = <lncrna_mrna_promoters@ranges@start>, “strand” =
<lncrna_mrna_promoters@strand>, “number\_of\_dbp” =
colSums(promoter\_peak\_occurence))

write\_csv(peak\_occurence\_df,
“hw\_results/peak\_occurence\_dataframe.csv”)

\#dont forget to save the environment

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
    #big data
    ggplot(big_num_peaks_df, aes(x = num_peaks, 
                             y = total_peak_length)) +
      geom_point(color = 'purple') 

# Now I want to color my plot by wether the protein is a TF or not.

\`\`\`{r are they TFs?} \#my data ggplot(num\_peaks\_df, aes(x =
num\_peaks, y = total\_peak\_length, )) + facet\_wrap(tf \~ .) +
geom\_point()

\#big data ggplot(big\_num\_peaks\_df, aes(x = num\_peaks, y =
total\_peak\_length, )) + facet\_wrap(tf \~ .) + geom\_point()


    # I want to make a histogram of the number of peaks for each of my proteins

    ```{r}
    # my data
    ggplot(num_peaks_df, aes(x = num_peaks)) +
      geom_histogram(bins = 20)

    #big data
    ggplot(big_num_peaks_df, aes(x = num_peaks)) +
      geom_histogram(bins = 30)

# Now I want to facet this by the type of DNA binding domain my protein has.

``` {r}
#my data
ggplot(num_peaks_df %>% filter(dbd %in% dbd),
       aes(x = num_peaks, y = total_peak_length )) +
  facet_wrap(tf ~ dbd) + 
  geom_point()

#big data
ggplot(big_num_peaks_df %>% filter(dbd %in% dbd),
       aes(x = num_peaks, y = total_peak_length )) +
  facet_wrap(tf ~ dbd) + 
  geom_point()
```

# Cool now I am ready to send my result to my collaborator as a

# Knitted document
