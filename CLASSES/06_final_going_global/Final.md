Final project
================
Alexzanah Griffin
4/21/2023

Analyze large data and load in data

``` r
#make paths

broadpeakfilepath <- "/scratch/Shares/rinnclass/CLASS_2023/data/data/peaks"

all_peaks_list <- import_peaks(consensus_file_path = broadpeakfilepath)

#peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)
```

``` r
# How does peak number and genome coverage compare (need gencode)

#Gene annotations
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
table(gencode_gr$type)
```

    ## 
    ##           gene     transcript           exon            CDS    start_codon 
    ##          60609         227462        1372308         761508          87662 
    ##     stop_codon            UTR Selenocysteine 
    ##          79913         310193            119

``` r
# exporting all genes file 
rtracklayer::export(gencode_genes, "results/gene_annotations/gencode_genes.gtf")

#mrna
mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"] 

rtracklayer::export(mrna_genes, "results/gene_annotations/mrna_genes.gtf")
table(gencode_genes$gene_type)
```

    ## 
    ##                          IG_C_gene                    IG_C_pseudogene 
    ##                                 14                                  9 
    ##                          IG_D_gene                          IG_J_gene 
    ##                                 37                                 18 
    ##                    IG_J_pseudogene                      IG_pseudogene 
    ##                                  3                                  1 
    ##                          IG_V_gene                    IG_V_pseudogene 
    ##                                144                                188 
    ##                             lncRNA                              miRNA 
    ##                              16849                               1881 
    ##                           misc_RNA                            Mt_rRNA 
    ##                               2212                                  2 
    ##                            Mt_tRNA             polymorphic_pseudogene 
    ##                                 22                                 42 
    ##               processed_pseudogene                     protein_coding 
    ##                              10171                              19965 
    ##                         pseudogene                           ribozyme 
    ##                                 18                                  8 
    ##                               rRNA                    rRNA_pseudogene 
    ##                                 52                                500 
    ##                             scaRNA                              scRNA 
    ##                                 49                                  1 
    ##                             snoRNA                              snRNA 
    ##                                942                               1901 
    ##                               sRNA                                TEC 
    ##                                  5                               1061 
    ##                          TR_C_gene                          TR_D_gene 
    ##                                  6                                  4 
    ##                          TR_J_gene                    TR_J_pseudogene 
    ##                                 79                                  4 
    ##                          TR_V_gene                    TR_V_pseudogene 
    ##                                106                                 33 
    ##   transcribed_processed_pseudogene     transcribed_unitary_pseudogene 
    ##                                495                                130 
    ## transcribed_unprocessed_pseudogene    translated_processed_pseudogene 
    ##                                923                                  2 
    ##  translated_unprocessed_pseudogene                 unitary_pseudogene 
    ##                                  2                                 98 
    ##             unprocessed_pseudogene                           vaultRNA 
    ##                               2631                                  1

``` r
#lncrna
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 
rtracklayer::export(lncrna_genes, "results/gene_annotations/lncrna_genes.gtf")

#both mrna and lncrna <3
mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
rtracklayer::export(mrna_lncrna_genes, "results/gene_annotations/mrna_lncrna_genes.gtf")

# mrna & lncrna annotation
lncrna_mrna_genes <- rtracklayer::import("results/gene_annotations/mrna_lncrna_genes.gtf")

#easier to see as df
lncrna_mrna_genes_df <- lncrna_mrna_genes %>% as.data.frame()

#promoter annotations set promoter size and save as df
lncrna_mrna_promoters <- promoters(lncrna_mrna_genes, upstream = 1000, downstream = 1000)

lncrna_mrna_promoters_df <-lncrna_mrna_promoters %>% as.data.frame()

rtracklayer::export(lncrna_mrna_promoters, "results/gene_annotations/lncrna_mrna_promoters.gtf")

#lncrna geneIDs
lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
table(mrna_lncrna_genes$gene_type)
```

    ## 
    ##         lncRNA protein_coding 
    ##          16849          19965

``` r
#mrna geneIDs
mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]
```

\#Make our peak data files

``` r
#how many peaks do i have in my files?
peak_num <- sapply(all_peaks_list, length) %>% as.data.frame(row.names = T)
```

    ## Warning in as.data.frame.integer(., row.names = T): 'row.names' is not a
    ## character vector of length 1064 -- omitting it. Will be an error!

``` r
names(peak_num) <- c("num_peaks")

peak_num <- peak_num %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")

#save this
#write_csv(peak_num, "results/num_peaks_df.csv")

#how many consensus peaks do I have (make consensus peak list)
dbps <- unique(sapply(names(all_peaks_list), function(x) {
   unlist(strsplit(x, "_"))[1]
}))

consensus_list <- lapply(dbps, consensus_from_reduced, all_peaks_list)

names(consensus_list) <- dbps

#how many consensus peaks do we have per gene
num_consensus_peaks <- sapply(consensus_list, length) %>% 
  as.data.frame() %>%
  rownames_to_column( var = "dbp") %>%
  dplyr::rename(consensus_peaks = ".")

#Add consensus into peak_num df
num_peaks <- sapply(consensus_list, length)

peak_num <- left_join(peak_num, num_consensus_peaks)
```

    ## Joining with `by = join_by(dbp)`

``` r
names(peak_num)[4] <- "total_consensus"
```

\#filter consensus peaks

``` r
#filter out data with less that 1000 peaks so make a filtered consensus peaks
filtered_consensus_list <- consensus_list[num_peaks >= 1000]


filtered_consensus_lower_list <- consensus_list[num_peaks < 1000]

#make a df
filtered_consensus_df <- sapply(filtered_consensus_list, length) %>% 
  as.data.frame() %>%
  rownames_to_column( var = "dbp") %>%
  dplyr::rename(consensus_peaks = ".")

filtered_consensus_lower_df <- sapply(filtered_consensus_lower_list, length) %>% 
  as.data.frame() %>%
  rownames_to_column( var = "dbp") %>%
  dplyr::rename(consensus_peaks = ".")


#save
#write_csv(filtered_consensus_df, "results/filtered_consensus_df.csv") 
#write_csv(filtered_consensus_lower_df,"results/consensus_lower_df.csv")
```

\#what is the peak genome coverage?

``` r
# Now we can do this: How does peak number and genome coverage compare (use consnesus_kist (grange obj) dont use dfs)
num_peaks_df <- data.frame("dbp" = names(filtered_consensus_list),
                           "num_peaks" = sapply(filtered_consensus_list, length))


num_peaks_df$total_peak_length <- sapply(filtered_consensus_list, function(x) sum(width(x)))


#Total genome coverage histogram

ggplot(num_peaks_df, aes(x = total_peak_length)) +
  geom_histogram(bins = 30) +
  ggtitle("peak genome coverage")
```

![](Final_files/figure-gfm/genome%20coverage-1.png)<!-- -->

``` r
#save
#ggsave("figures/peaks_v_genome_coverage_hist.pdf")

#peak distribution on dbps
 ggplot(num_peaks_df, aes(x = num_peaks)) + 
  geom_histogram(bins = 70, color = "black") +
   ggtitle("peaks_distribution")
```

![](Final_files/figure-gfm/genome%20coverage-2.png)<!-- -->

``` r
# saving
#ggsave("figures/peaks_distribution_hist.pdf")

#scatterplot
ggplot(num_peaks_df, aes(x = num_peaks, y = total_peak_length)) +
  geom_point() + 
  geom_smooth(method = "gam", se = TRUE, color = "black", lty = 2)+
         
  ylab("BP covered") +
  xlab("Number of peaks") +
  ggtitle("Peak count vs. total bases covered")
```

    ## `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'

![](Final_files/figure-gfm/genome%20coverage-3.png)<!-- -->

``` r
#ggsave("figures/peak_vs_coverage_scatter.pdf")
```

\#results \#\#The bulk of peaks are shorter and cover less of the
genome(&lt;50,000),But there are some (very few) peaks that are larger
and cover more of the genome. \#\# there is a linear relationship
between genome coverage and peak number.

``` r
#genebody peak counts
genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                filtered_consensus_list, 
                                                type = "counts")

#overlapping peaks to total genebodies
num_peaks_df$peaks_overlapping_genebody <- rowSums(genebody_peak_counts)

#lncrna to genebodies
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA gene bodies
num_peaks_df$peaks_overlapping_mrna_genebody <- rowSums(genebody_peak_counts[,mrna_gene_ids])

#write_csv(num_peaks_df, "results/num_peaks_df.csv")

#now lets look at it (peaks v. genebodies)
hist(num_peaks_df$peaks_overlapping_genebody)
```

![](Final_files/figure-gfm/peaks%20overlapping%20genebody-1.png)<!-- -->

``` r
#SAVE MEE
ggsave("figures/peak_vs_genebody.pdf")
```

    ## Saving 7 x 5 in image
    ## `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'

``` r
# all rna genebody overlaps
num_peaks_df$percent_total_genebody_overlaps <- (num_peaks_df$peaks_overlapping_genebody/36814)*100
percent_genebody_overlaps <- mean(num_peaks_df$percent_total_genebody_overlaps)

print(paste0(percent_genebody_overlaps, "%", " of peaks overlap total genebody"))
```

    ## [1] "62.2646781242222% of peaks overlap total genebody"

``` r
# % of lncrna genebody overlaps
num_peaks_df$percent_lncrna_genebody_overlaps <- (num_peaks_df$peaks_overlapping_lncrna_genebody/36814)*100
percent_lncrna_genebody_overlaps <- mean(num_peaks_df$percent_lncrna_genebody_overlaps)

print(paste0(percent_lncrna_genebody_overlaps, "%", " of peaks overlap lncrna genebody"))
```

    ## [1] "12.9164208257475% of peaks overlap lncrna genebody"

``` r
# % of mrna genebody overlaps
num_peaks_df$percent_mrna_genebody_overlaps <- (num_peaks_df$peaks_overlapping_mrna_genebody/36814)*100
percent_mrna_genebody_overlaps <- mean(num_peaks_df$percent_mrna_genebody_overlaps)

print(paste0(percent_mrna_genebody_overlaps, "%", " of peaks overlap mrna genebody"))
```

    ## [1] "49.3482572984747% of peaks overlap mrna genebody"

\#results \#\#62.2646781242222% of peaks overlap total genebody
\#\#12.9164208257475% of peaks overlap lncrna genebody
\#\#49.3482572984747% of peaks overlap mrna genebody \#\#genebody
ovelaps follow the same trend of promoter overlaps \#\#more overlaps are
on mRNA promoters that lncrna

\#peak promoter coverage?

``` r
# What is the distribution of promoter overlaps versus gene-bodies (hint hist)
#promoter and peak overlaps (matrix)
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "counts")

#matrix
promoter_peak_occurrence_matrix <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_list, type = "occurrence")

#make peak occurrnce df
promoter_peak_occurrence_df <-data.frame("gene_id" = colnames(promoter_peak_counts),
  "gene_name" = lncrna_mrna_promoters$gene_name,
  "gene_type" = lncrna_mrna_promoters$gene_type,
  "chr" = lncrna_mrna_promoters@seqnames,   
  "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
  "strand" = lncrna_mrna_promoters@strand,
  "number_of_dbp" = colSums(promoter_peak_counts))

#write_csv(promoter_peak_occurrence_df, "results/promoter_peak_occurrence_df.csv")
```

``` r
#peaks overlapping promoters
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

# lncrna promoter overlaps 
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])

# mrna promoter overlaps
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

#peaks vs promoters histogram
ggplot(num_peaks_df, aes(peaks_overlapping_promoters))+
  ggtitle("Number of peaks overlapping promoters")+
  geom_histogram(bins = 70, fill = "cyan4", color = "black")
```

![](Final_files/figure-gfm/peak%20overlaps%20with%20promoters-1.png)<!-- -->

``` r
#hist(num_peaks_df$peaks_overlapping_promoters)

#ggsave("figures/peaks_v_promoters.pdf")

#scatterplot
ggplot(num_peaks_df,
       aes(x = num_peaks, y = peaks_overlapping_promoters)) +
  xlab("Peaks per DBP") +
  ylab("Number of peaks overlapping promoters") +
  ggtitle(" Peaks on Promoters")+
  geom_point(color = "cyan4") +
  geom_abline(slope = 1, linetype="dashed") +
  geom_smooth(method = "lm", se=FALSE, formula = 'y ~ x',
              color = "firebrick4") 
```

![](Final_files/figure-gfm/peak%20overlaps%20with%20promoters-2.png)<!-- -->

``` r
#ggsave("figures/peaks_v_promoters_scatter.pdf")

#Density plot (may not need this)
#ggplot(promoter_peak_occurrence_df, aes(x = number_of_dbp)) +
#geom_density(alpha = 0.2, color = "black", fill = "cyan4") +
  #theme_paperwhite() +
#  xlab(expression("Number of DBPs")) +
 # ylab(expression("Density")) +
  #ggtitle("Promoter binding events",
   #       subtitle = "mRNA and lncRNA genes") 

# Avg promoter overlaps
#total_num_peaks <- sum(num_peaks_df$num_peaks)

# all promoter overlaps
num_peaks_df$percent_total_promoter_overlaps <- (num_peaks_df$peaks_overlapping_promoters/36814)*100
percent_promoter_overlaps <- mean(num_peaks_df$percent_total_promoter_overlaps)

print(paste0(percent_promoter_overlaps, "%", "of peaks overlap promoters"))
```

    ## [1] "27.2122081968311%of peaks overlap promoters"

``` r
# % of lncrna promoter overlaps
num_peaks_df$percent_lncrna_promoter_overlaps <- (num_peaks_df$peaks_overlapping_lncrna_promoters/36814)*100
percent_lncrna_promoter_overlaps <- mean(num_peaks_df$percent_lncrna_promoter_overlaps)

print(paste0(percent_lncrna_promoter_overlaps, "%", "of peaks overlap lncrna promoters"))
```

    ## [1] "6.62486844615484%of peaks overlap lncrna promoters"

``` r
# % of mrna promoter overlaps
num_peaks_df$percent_mrna_promoter_overlaps <- (num_peaks_df$peaks_overlapping_mrna_promoters/36814)*100
percent_mrna_promoter_overlaps <- mean(num_peaks_df$percent_mrna_promoter_overlaps)

print(paste0(percent_mrna_promoter_overlaps, "%", "of peaks overlap mrna promoters"))
```

    ## [1] "20.5873397506762%of peaks overlap mrna promoters"

\#results \#\#More peaks cover less DBP (linear relationship), most bind
to &lt;50,000 promoters, similar to the trend of genome coverage.
\#\#27.2122081968311% of peaks overlap promoters \#\#6.62486844615484%
of peaks overlap lncrna promoters \#\#20.5873397506762% of peaks overlap
mrna promoters

``` r
#lnc V mRNA 
lncrna_vs_mrna_overlaps <- num_peaks_df %>%
  dplyr::select(-peaks_overlapping_promoters) %>%
  pivot_longer(cols = peaks_overlapping_lncrna_promoters:peaks_overlapping_mrna_promoters,
               names_to = "gene_type",
               values_to = "peaks_overlapping_promoters") %>%
  mutate(gene_type = gsub("peaks_overlapping_", "", gene_type))

# plotting lncrna mrna
ggplot(lncrna_vs_mrna_overlaps, aes(x = num_peaks, y = peaks_overlapping_promoters, 
                         col = gene_type)) +
         geom_point() +
         geom_abline(slope = 1, linetype="dashed") +
  geom_smooth(method = "lm", se = FALSE, formula = "y ~ x") +
  scale_color_manual(values = c("darkorchid4", "cadetblue2"))+
  xlab("Peaks per DBP") +
  ylab("Peaks Overlapping Promoters") +
  ggtitle("DBP Peaks to Promoter Overlaps")
```

![](Final_files/figure-gfm/lnc%20&mrna%20overlaps-1.png)<!-- -->

``` r
# saving
#ggsave("figures/peaks_overlaps_by_gene_type.pdf")
```

\#results \#\# mrna starts to plateu after reaching a peak DBP coverage
before 50,000. While promoters overlapping lncrna per DBP hals plateues,
it’s not as dramatic as the mrna promoters.

\#Lets look at superbinders

``` r
# Make a list of genes that are "super binders" (peak feats counts peaks/feature counts and occurence [peaks that have more than 200 bound to promoter(could be more up to u)])
#super_binder_threshold <- 200
super_binder <- filter(promoter_peak_occurrence_df, number_of_dbp>= 200)

#write_csv(super_binder, "results/super_binders.csv")

#lncrna superbinders
lncrna_super_binders <- filter(super_binder, gene_type == "lncRNA" )
length(lncrna_super_binders)
```

    ## [1] 7

``` r
#mrna superbinders
mrna_super_binders <- filter(super_binder, gene_type == "protein_coding" )
length(mrna_super_binders)
```

    ## [1] 7

``` r
#hist of superbinders
ggplot(super_binder, aes(number_of_dbp))+
  ggtitle("all dbps")+
  geom_histogram(bins = 70, color = "black")
```

![](Final_files/figure-gfm/superbinder%20plots-1.png)<!-- -->

``` r
#ggsave("hist_all_superbinders.pdf")

#lncrna
ggplot(lncrna_super_binders, aes(number_of_dbp))+
  ggtitle("lncrna dbps")+
  geom_histogram(bins = 70, fill = "darkorchid4", color = "black")
```

![](Final_files/figure-gfm/superbinder%20plots-2.png)<!-- -->

``` r
#ggsave("hist_lncrna_superbinders.pdf")

#mrna
ggplot(mrna_super_binders, aes(number_of_dbp))+
  ggtitle("mrna dbps")+
  geom_histogram(bins = 70, fill = "cyan4", color = "black")
```

![](Final_files/figure-gfm/superbinder%20plots-3.png)<!-- -->

``` r
#ggsave("hist_mrna_superbinders.pdf")
```

\#results: (explain your findings of each section) \#\#the super binder
threshold is 200 DBPs \#\#Do lncRNAs also have super-binding promoters?
\#\#\#\#yes \#\#lncrna DBPs are more distributed. However mrna DBPs are
more localized to the DBPs they bind to and more are binding on the
DBPs.

# Clustering

``` r
# dbps on promoters object
DBPs_on_promoter <- lncrna_mrna_promoters %>%
                    as.data.frame() %>%
  dplyr::select(gene_id, gene_name)

#df for promoter dbps
promoter_dbps <- promoter_peak_occurrence_matrix %>%
  as.data.frame() %>%
  rownames_to_column("dbp") %>%
pivot_longer(2:ncol(.), names_to = "gene_id", values_to = "occurrence") %>%
  filter(occurrence == 1) %>%
  dplyr::select(-occurrence) %>%
  left_join(DBPs_on_promoter)
```

    ## Joining with `by = join_by(gene_id)`

``` r
write.table(promoter_dbps, "results/promoter_dbps.tsv")

#save as df too
promoter_dbps_df <- promoter_dbps %>% as.data.frame()

write.csv(promoter_dbps, "results/promoter_dbps.csv")

#ARID3A promoter dbps
arid_promoter <- promoter_dbps %>%
  filter(gene_name == "ARID3A")
```

``` r
# creating distance matrix
peak_occurrence_dist <- dist(promoter_peak_occurrence_matrix, method = "binary")

#cluster distance matrix
bin_hier <- hclust(peak_occurrence_dist, method = "complete")

promoter_dbps <- promoter_peak_occurrence_matrix %>%
  as.data.frame() %>%
  rownames_to_column("dbp") %>%
pivot_longer(2:ncol(.), names_to = "gene_id", values_to = "occurrence") %>%
  filter(occurrence == 1) %>%
  dplyr::select(-occurrence) %>%
  left_join(DBPs_on_promoter)
```

    ## Joining with `by = join_by(gene_id)`

``` r
#dendrogram chart of binding profies via promoters
ggdendro::ggdendrogram(bin_hier, rotate = FALSE,  size = 3,
                       theme_dendro = TRUE) +
   coord_flip() +
   scale_y_continuous() +
   scale_x_continuous(position = "top") +
   scale_x_continuous(breaks = seq_along(bin_hier$labels[bin_hier$order]),
             labels = bin_hier$labels[bin_hier$order], position = "top",
             expand = c(0,0)) +
   theme(axis.text.x = element_text(angle = 90, hjust  = 1)) +
   theme(axis.text.y = element_text(angle = 0,hjust = 1)) +
   scale_y_reverse(expand = c(0.01, 0)) +
   theme(
     plot.background = element_blank(),
     panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
     panel.border = element_blank()
   )
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

![](Final_files/figure-gfm/Distance%20Matrix-1.png)<!-- -->

``` r
#ggsave("figures/promoter_overlap_dendrogram.pdf", limitsize = FALSE, height = 55, width = 12)



# if we cluster by lncRNA and mRNA !!!seperately!!! what are some similarities and differences?
```

``` r
# make lnc and mrna promoters objects
lncrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type %in% "lncRNA"] 

mrna_promoters <- lncrna_mrna_promoters[lncrna_mrna_promoters$gene_type %in% "protein_coding"] 

# separate peak_occurrence_matrix into the 2 rna types
lncrna_peak_occurrence <- promoter_peak_occurrence_matrix[,lncrna_promoters$gene_id]

mrna_peak_occurrence <- promoter_peak_occurrence_matrix[,mrna_promoters$gene_id]

# make individual distance matrix
bin_hier_lncrna <- hclust(dist(lncrna_peak_occurrence, method = "binary"))

bin_hier_mrna <- hclust(dist(mrna_peak_occurrence, method = "binary"))

# plot with ggdendro (lncrna)
ggdendro::ggdendrogram(bin_hier_lncrna, rotate = T,  size = 3)
```

![](Final_files/figure-gfm/RNA%20Clusters-1.png)<!-- -->

``` r
#ggsave("figures/promoter_lncrna_overlap_dendrogram.pdf", limitsize = FALSE, height = 55, width = 12)

# plot with ggdendro (mrna)
ggdendro::ggdendrogram(bin_hier, rotate = TRUE,  size = 3)
```

![](Final_files/figure-gfm/RNA%20Clusters-2.png)<!-- -->

``` r
#ggsave("figures/promoter_mrna_overlap_dendrogram.pdf", limitsize = FALSE, height = 55, width = 12)
```

\#Results \#\#General Clusters: \#\#\# I was suprised to see that ARIDs
didnt cluster togeteher at all considering they have overlap in the
cancers/tumors they help to regulate. \#\#\# ZNFs are scattered through
out as is expected with their varying functions and many unknown ones.
\#\#RNA clusters: \#\#\# In both lnc and mRNA ARIDs stay in the same
general cluster group (i.e. ARID3A clusters with JUN and JUND), their
closeness does change slightly however.

\#metaplots \`\`\`{\#r DBPs on promoters} \#Very common plots to make

# Let’s look at the metaplot for all DBPs on lncRNA and mRNA promoters seperately (hint facet wrap).

metaplot\_df &lt;- data.frame(x = integer(), dens = numeric(), dbp =
character())

# for loop to populate DF

\#suppressWarnings(for(i in 1:length(filtered\_consensus\_list)) { \#
tmp\_df &lt;- profile\_tss(filtered\_consensus\_list\[\[i\]\],
lncrna\_mrna\_promoters) \#tmp\_df\$dbp &lt;-
names(filtered\_consensus\_list)\[\[i\]\] \#metaplot\_df &lt;-
bind\_rows(metaplot\_df, tmp\_df) \#})

# saving

\#write\_rds(metaplot\_df, “results/metaplot\_df.rds”)

\#write\_csv(metaplot\_df, “results/metaplot\_df.csv”)

metaplot\_df &lt;- read\_rds(“results/metaplot\_df.rds”)

metaplot\_filtered\_matrix &lt;- metaplot\_df %&gt;%
pivot\_wider(names\_from = x, values\_from = dens) %&gt;%
column\_to\_rownames(“dbp”) %&gt;% as.matrix() mm\_scaled &lt;-
t(scale(t(metaplot\_filtered\_matrix)))

\#metaplot matrix filtered metaplot\_filtered\_matrix &lt;- metaplot\_df
%&gt;% pivot\_wider(names\_from = x, values\_from = dens) %&gt;%
column\_to\_rownames(“dbp”) %&gt;% as.matrix() metaplot\_scaled &lt;-
t(scale(t(metaplot\_filtered\_matrix)))

write\_rds(metaplot\_scaled, “results/metaplot\_scaled.rds”)
metaplot\_scaled &lt;- read\_rds(“results/metaplot\_scaled.rds”)

metaplot\_hclust &lt;- hclust(dist(metaplot\_filtered\_matrix, method =
“complete”))

\#plot plot(hclust(dist(metaplot\_scaled), method = “complete”)) par(cex
= 0.4)

\#ggsave(“figures/dbp\_cluster.pdf”)

    # Results
    ##ARID3A has different cluster mates from the genes it clusters with in lnc and mRNA clusters. I wonder if it could be due to the difference in tpm, since these take that into account. 

    # RNAseq expression

    ```r
    #need samplesheet adn counts table
    samplesheet <- read_rds("../05_R_analyses/05_RNAseq/01_differential_expression/results/final_samplesheet.rds")

    # load rnaseq salmon
    salmon_tpm <- read.csv("/scratch/Shares/rinnclass/CLASS_2023/atg/Analysis_of_genome_2023/CLASSES/05_R_analyses/05_RNAseq/00_RNAseq_download_NF_core_pipeline/my_seq/results/salmon/salmon_merged_gene_tpm.csv")

    #tpm table with same order as sample sheet
    tpm <- salmon_tpm %>% 
      pivot_longer(cols = 2:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
      merge(samplesheet) %>%
      group_by(gene_id, condition) %>%
      summarize(tpm = mean(tpm, na.rm = T)) %>%
      pivot_wider(names_from = condition, values_from = tpm, names_prefix = "tpm_")

    ## `summarise()` has grouped output by 'gene_id'. You can override using the
    ## `.groups` argument.

``` r
#as matrix
tpm_matrix <- tpm %>% 
  column_to_rownames("gene_id") %>%
  as.matrix()
tpm_scaled <- t(scale(t(tpm_matrix)))
tpm_scaled <- tpm_scaled[complete.cases(tpm_scaled),]
```

``` r
# promoter features = promoter peak occurrence df
promoter_features_df <- promoter_peak_occurrence_df

#merge tpms
promoter_features_df <- promoter_features_df %>%
  left_join(tpm)
```

    ## Joining with `by = join_by(gene_id)`

``` r
#Total expression heat map
#pheatmap::pheatmap(tpm_scaled, show_rownames = FALSE)
pheatmap(tpm_scaled, cluster_columns = FALSE,  border = TRUE, border_color = "black", 
        show_column_names = FALSE,
        show_rownames = FALSE,
        use_raster = TRUE,
        column_gap = unit(0, "mm"))
```

![](Final_files/figure-gfm/tpm%20plots-1.png)<!-- -->

``` r
graphics.off()

#ggsave("figures/total_expression_heatmap.pdf", height = 40, width = 30)

#Let's make a heatmap of nuclear versus cytoplasmic expression

# Filter tpm_scaled to nuc v cyto
tpm_scaled_nuc_cyto <- tpm_scaled[,colnames(tpm_scaled) == "tpm_homo_sapiens_cytosolic_fraction"
                                  | colnames(tpm_scaled) == "tpm_homo_sapiens_nuclear_fraction"]

# plot
pheatmap::pheatmap(tpm_scaled_nuc_cyto, show_rownames = FALSE)

#ggsave("figures/nuc_cyto_heatmap.pdf", height =49, width = 12)

#lncRNA and mRNA expression in a density plot
ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_hepg2 + 0.01), color = gene_type, fill = gene_type))+
  geom_density(alpha = 0.2) +
    scale_color_manual(values = c("darkorchid4", "cyan4"), name = "Gene type")
```

    ## Warning: Removed 102 rows containing non-finite values (`stat_density()`).

``` r
#ggsave("figures/lnc_mrna_expression_density.pdf")

# nuclear RNAs density plot
ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_nuclear_fraction + 0.01), color = gene_type, fill = gene_type))+
  geom_density(alpha = 0.2)+
  scale_color_manual(values = c("darkorchid4", "cyan4"), name = "Gene type")
```

    ## Warning: Removed 102 rows containing non-finite values (`stat_density()`).

``` r
#ggsave("figures/lnc_mrna_nuc_expression_density.pdf")
```

\#Results \#\#Most of the RNA is nuclear, which makse since, becuase
both mRNA and lncRNA can be found in the nucleus, but lncRNA is not
abundantly found in the cytoplasm. \#\# The density plots shows
bimodality of lnc and mrna. Nuclear density has more expression in both
lnc and mrna, and lncRNAs peaks in density right shifted.

``` r
# scatter of dbps binding vs total RNA expression
ggplot(promoter_features_df, 
            aes(y = log2(tpm_homo_sapiens_hepg2 + 0.001), x = number_of_dbp, color = gene_type)) + 
geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_hepg2 > 0.001),
             shape = 17, alpha = 0.7) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("darkorchid4", "cyan4"), name = "Gene type")+
ggtitle("Expression vs. promoter binding events") + 
  xlab(expression('Number of TFs')) +
  ylab(expression(log[2](TPM))) 
```

    ## Warning: Removed 102 rows containing non-finite values (`stat_smooth()`).

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 102 rows containing non-finite values (`stat_smooth()`).

![](Final_files/figure-gfm/DBPs%20vs%20RNA%20expression-1.png)<!-- -->

``` r
# save
#ggsave("figures/dbp_v_rna_expresion.pdf")
```

\#results \#\# There is a positive trend when looking at number of
transcripts and DBPs bound the a site for both lnc and mRNA, which is
expected. Also as seen in previous results there is a sort of plateaue
where the DBPs bound decreases and levelks out.

``` r
#plot nuclear RNA abundance v number of DBPs bound to promoter
ggplot(promoter_features_df, 
            aes(y = log2(tpm_homo_sapiens_nuclear_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
  geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_nuclear_fraction > 0.001),
             shape = 17, alpha = 0.7) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("darkorchid4", "cyan4"), name = "Gene type") + 
  ggtitle("Nuclear Expression vs. RNA bound to promoters") + 
  xlab(expression('Number of DBPs')) +
  ylab(expression(log[2](TPM))) 
```

    ## Warning: Removed 102 rows containing non-finite values (`stat_smooth()`).

![](Final_files/figure-gfm/bound%20v%20nuc%20and%20cytosolic%20rna%20expression-1.png)<!-- -->

``` r
# saving figure
#ggsave("figures/nuc_expression_v_promoter_bound.pdf")

# plot cyto RNA abundance v # DBP on promoters
ggplot(promoter_features_df, 
            aes(y = log2(tpm_homo_sapiens_cytosolic_fraction + 0.001), x = number_of_dbp, color = gene_type)) + 
  geom_point(data = promoter_features_df %>% filter(tpm_homo_sapiens_cytosolic_fraction > 0.001),
             shape = 17, alpha = 0.7) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("darkorchid4", "cyan4"), name = "Gene type") + 
  ggtitle("Cytoplasmic Expression vs. promoter binding events") + 
  xlab(expression('Number of DBPs')) +
  ylab(expression(log[2](TPM))) 
```

    ## Warning: Removed 102 rows containing non-finite values (`stat_smooth()`).

![](Final_files/figure-gfm/bound%20v%20nuc%20and%20cytosolic%20rna%20expression-2.png)<!-- -->

``` r
# saving figure
#ggsave("figures/cyto_expression_v_promoter_bound.pdf")
```

\#Results \#\#Cytoplasmic expression is in less abundance than nuclear
expression, and the regression line has a slight drop that is seen on
the nuclear expression. These results confirm previous results made
about cytoplasmic and nuclear expression.

``` r
# If we zoom in on high binding promoters (> 200 DBPs) are there any that don't have any expression?

promoter_features_df$hepg2_reservoir <- 
  (promoter_features_df$number_of_dbp > 200 & 
               promoter_features_df$tpm_homo_sapiens_hepg2 < 0.001)

# filtering super binders
promoter_features_df <- promoter_features_df %>%
  mutate(superbinder = promoter_features_df$number_of_dbp > 200)

# setting column of superbinders
promoter_features_df <- promoter_features_df %>% 
  mutate(superbinder = ifelse(promoter_features_df$hepg2_reservoir ==T, "super_binder", "not_super_binder"))

# make density plot superbinder vs non superbinder expression in total RNA
ggplot(promoter_features_df, aes(x = log2(tpm_homo_sapiens_hepg2 + 0.01), color = superbinder, fill = superbinder))+
  geom_density(alpha = 0.2) +
  scale_color_manual(values = c("darkorchid4", "cyan4"))
```

    ## Warning: Removed 102 rows containing non-finite values (`stat_density()`).

![](Final_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# saving figure
#ggsave("figures/reservoir_v_not_rna_density.pdf")

# if we zoom in on high binding promoters (> 200 DBPs) are there any that don't have any expression?
no_expression_tpm <- filter(promoter_features_df, tpm_homo_sapiens_hepg2 < 0.1)

no_exp_high_binder_tpm <- filter(no_expression_tpm, superbinder == TRUE)

# save df as 
write.csv(no_exp_high_binder_tpm, "results/no_exp_high_binder_tpm.csv", row.names=TRUE)

lengthnoexp_highbinder <- nrow(no_exp_high_binder_tpm)

#print(paste0("There are ", lengthnoexp_highbinder, " super binders with no expression (TPM value of < 0.1)."))
```

\#Results \#\# All of the high binding promoters have expressions.
\#\#The RNA density of the Superbinders is very left shifted. However,
eventhough the nonsuperbinders are also leftward shifted, density is
trimodal.

```
# Let's make a heatmap of genes that are variable across samples 
pheatmap(metaplot_scaled, cluster_columns = FALSE,  border = TRUE, border_color = "black", 
        show_column_names = TRUE,
        use_raster = TRUE,
        column_gap = unit(0, "mm"))

ggsave("figures/heatmap_variable_genes.pdf",  height = 49, width = 40)
```
