---
title: "Wheat Rhizobacteria paper - Bioinformatic processing"
author: "Tessa Reid"
date: "28/10/2022"
---
# Code for Bioinformatic Processing of Amplicon Sequence Data


# 1. Processing of culture-independent (CI) 16S rRNA gene dataset from raw data
### In Qiime2
### Import fastq files in Qiime2
```{r}
source activate qiime2-2019.7

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 'Data to import' \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path CI-demux-paired-end.qza
```

## Visualize and verify sequence quality (total number of sequences = 15 264 143). Open qzv file on Qiime2View online
```{r}
qiime demux summarize \
  --i-data CI-demux-paired-end.qza \
  --o-visualization CI-demux-paired-end.qzv
```


## Join forward and reverse reads with the Vsearch command
```{r}
qiime vsearch join-pairs \
  --i-demultiplexed-seqs CI-demux-paired-end.qza \
  --p-truncqual 0 \
  --p-minlen 150 \
  --p-allowmergestagger \
  --p-maxdiffs 15 \
  --p-minmergelen 200 \
  --p-maxee 10 \
  --o-joined-sequences CI-demux-joined.qza \
  --verbose
```

## Visualize and verify resulting reads (after merging, number of sequences = 6 722 536). Open qzv file on Qiime2View online
```{r}
qiime demux summarize \
  --i-data CI-demux-joined.qza \
  --o-visualization CI-demux-joined.qzv
```
  
## Export merged sequences 
```{r}
qiime tools export \
  --input-path CI-demux-joined.qza \
  --output-path CI-demux-joined
```

## Import merged sequences as single-end in Qiime2
```{r}
qiime tools import \
  --type   'SampleData[SequencesWithQuality]' \
  --input-path CI-demux-joined \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path CI-demux-single-end.qza
```


## Sequence Denoising with DADA2 (after DADA2, number of ASVs = 21 690 and total number of sequences = 5 354 119)
```{r}
qiime dada2 denoise-single \
  --i-demultiplexed-seqs CI-demux-single-end.qza \
  --p-trunc-len 230 \
  --p-trim-left 5 \
  --o-representative-sequences CI-rep-seqs.qza \
  --o-table CI-table.qza \
  --o-denoising-stats CI-stats.qza
```

## Visualize the outputs of DADA2 in Qiime2View online
```{r}
qiime metadata tabulate \
  --m-input-file CI-denoising-stats.qza \
  --o-visualization CI-denoising-stats.qzv

qiime feature-table summarize \
  --i-table CI-table.qza \
  --o-visualization CI-table.qzv \
  --m-sample-metadata-file CI-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data CI-rep-seqs.qza \
  --o-visualization CI-rep-seqs.qzv
```

## Remove samples with low read frequencies (<100 reads) and ASVs with low abundance (less than 10 total observation count) from ASV table and rep-seq file


## Assign taxonomy with Qiime2 using SILVA 138 database
```{r}
qiime feature-classifier classify-consensus-vsearch \
    --i-query CI-rep-seqs-filtered.qza \
    --i-reference-reads silva_138_ASVs.qza \
    --i-reference-taxonomy silva_138_taxonomy.qza \
    --o-classification CI-taxonomy.qza \
    --p-threads 30
```

## Visualise taxa results to check for non-bacterial ASVs. Open qzv file on Qiime2View online
```{r}
qiime taxa barplot \
  --i-table CI-table-filtered.qza \
  --i-taxonomy CI-taxonomy.qza \
  --m-metadata-file CI-metadata.tsv \
  --o-visualization CI-taxa-bar-plots.qzv
```

## Remove Unassigned, Archaeal ASVs from ASV table and rep-seqs file
## Re-assign taxonomy with filtered rep-seqs file to produce a taxonomy table with only bacterial ASVs


## Make phylogenic tree in Qiime2 from filtered rep-seqs file (non-bacterial, low abundance sequences removed)
```{r}
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences CI-rep-seqs-filtered.qza \
  --o-alignment CI-aligned-rep-seqs.qza \
  --o-masked-alignment CI-masked-aligned-rep-seqs.qza \
  --o-tree CI-unrooted-tree.qza \
  --o-rooted-tree CI-rooted-tree.qza
```


## Export ASV table, taxonomy table, and tree for analysis in phyloseq

### Taxonomy Table
```{r}
qiime tools export \
  --input-path CI-taxonomy-16S.qza \
  --output-path CI-phyloseq
```

### Phylogenetic Tree
```{r}
qiime tools export \
  --input-path CI-rooted-tree.qza \
  --output-path CI-phyloseq
```

### ASV Table
```{r}
qiime tools export \
  --input-path CI-filtered-table-16S.qza \
  --output-path CI-phyloseq
```

### Convert Biom files to tsv
```{r}
biom convert \
  -i CI-phyloseq/feature-table.biom \
  -o CI-phyloseq/otu_table.txt \
  --to-tsv
```





######################################################################################

# 2. Processing of culture-dependent (CD) 16S rRNA gene dataset from raw data
### In Qiime2
### Import fastq files in Qiime2
```{r}
source activate qiime2-2019.7

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 'Data to import' \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path CD-demux-paired-end.qza
```

## Visualize and verify sequence quality (total number of sequences = 21 606 517). Open qzv file on Qiime2View online
```{r}
qiime demux summarize \
  --i-data CD-demux-paired-end.qza \
  --o-visualization CD-demux-paired-end.qzv
```


## Sequence Merging and Denoising with DADA2 (after DADA2, number of ASVs = 29 516 and total number of sequences = 1 847 197)
```{r}
qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux-paired-end.qza \
--p-trim-left-f 13 \
--p-trim-left-r 13 \
--p-trunc-len-f 250 \
--p-trunc-len-r 250 \
--o-table CD-table.qza \
--o-representative-sequences CD-rep-seqs.qza \
--o-denoising-stats CD-denoising-stats.qza \
--verbose \
--p-n-threads 32
```

## Visualize the outputs of DADA2 in Qiime2View online
```{r}
qiime metadata tabulate \
  --m-input-file CD-denoising-stats.qza \
  --o-visualization CD-denoising-stats.qzv

qiime feature-table summarize \
  --i-table CD-table.qza \
  --o-visualization CD-table.qzv \
  --m-sample-metadata-file CD-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data CD-rep-seqs.qza \
  --o-visualization CD-rep-seqs.qzv
```

## Remove samples with low read frequencies (<100 reads) and ASVs with low abundance (less than 10 total observation count) from ASV table and rep-seq file


## Assign taxonomy with Qiime2 using SILVA 138 database
```{r}
qiime feature-classifier classify-consensus-vsearch \
    --i-query CD-rep-seqs-filtered.qza \
    --i-reference-reads silva_138_ASVs.qza \
    --i-reference-taxonomy silva_138_taxonomy.qza \
    --o-classification CD-taxonomy.qza \
    --p-threads 30
```

## Visualise taxa results to check for non-bacterial ASVs. Open qzv file on Qiime2View online
```{r}
qiime taxa barplot \
  --i-table CD-table-filtered.qza \
  --i-taxonomy CD-taxonomy.qza \
  --m-metadata-file CD-metadata.tsv \
  --o-visualization CD-taxa-bar-plots.qzv
```

## Remove Unassigned, Archaeal ASVs from ASV table and rep-seqs file
## Re-assign taxonomy with filtered rep-seqs file to produce a taxonomy table with only bacterial ASVs


## Make phylogenic tree in Qiime2 from filtered rep-seqs file (non-bacterial, low abundance sequences removed)
```{r}
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences CD-rep-seqs-filtered.qza \
  --o-alignment CD-aligned-rep-seqs.qza \
  --o-masked-alignment CD-masked-aligned-rep-seqs.qza \
  --o-tree CD-unrooted-tree.qza \
  --o-rooted-tree CD-rooted-tree.qza
```


## Export ASV table, taxonomy table, and tree for analysis in phyloseq

### Taxonomy Table
```{r}
qiime tools export \
  --input-path CD-taxonomy-16S.qza \
  --output-path CD-phyloseq
```

### Phylogenetic Tree
```{r}
qiime tools export \
  --input-path CD-rooted-tree.qza \
  --output-path CD-phyloseq
```

### ASV Table
```{r}
qiime tools export \
  --input-path CD-filtered-table-16S.qza \
  --output-path CD-phyloseq
```

### Convert Biom files to tsv
```{r}
biom convert \
  -i CD-phyloseq/feature-table.biom \
  -o CD-phyloseq/otu_table.txt \
  --to-tsv
```
