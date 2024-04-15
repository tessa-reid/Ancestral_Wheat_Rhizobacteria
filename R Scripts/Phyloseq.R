############################################################
#R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################
#R 4.2.2
#setwd()

library(phyloseq) #v1.42.0
library(rncl) #v0.8.7

####### 1. Phyloseq object for culture-independent amplicon sequences #######


#Input files: Files produced by qiime2 (otu_table_CI.txt, taxonomy_CI.txt, tree_CI.nwk) and sample metadata (metadata_CI.txt).



#1.1. Importing the files into R:

otu <- read.table("asv_abundances_CI.txt", header=TRUE, row.names=1)
taxonomy <- read.table(file="asv_taxonomy_CI.txt", sep="\t", header=TRUE, fill=TRUE)
sample <- read.table(file="sample_data_CI.txt", sep="\t",header=TRUE, row.names = 1)
tree  <- read_newick_phylo(file="tree_CI.nwk",simplify=FALSE)



#1.2. Create the phyloseq object:

#otu table
otu <- as.matrix(otu)
otu.table <- otu_table(otu,taxa_are_rows=TRUE,errorIfNULL = TRUE)

#taxonomy table
row.names(taxonomy) <- taxonomy$ASVs
taxonomy <- as.matrix(taxonomy)
taxonomy.table <- tax_table(taxonomy)

#sample table
sample=sample_data(sample,errorIfNULL = TRUE)
#define factors and reorder levels
sample$SampleType = factor(sample$SampleType)
sample$Blanks = factor(sample$Blanks)
sample$Fertilization = factor(sample$Fertilization, levels = c("Non-fertilized", "Fertilized"))
sample$Soil = factor(sample$Soil, levels = c("Unplanted", "Rhizosphere", "Rhizoplane"))
sample$Ancestral_class = factor(sample$Ancestral_class, levels = c("Unplanted","Wild ancestor", "Traditional cultivar", "Commercial cross"))
sample$Ploidy = factor(sample$Ploidy, levels = c("Unplanted", "Diploid", "Tetraploid", "Hexaploid"))
sample$Genome = factor(sample$Genome, levels = c("Unplanted", "AA", "AABB", "AABBDD", "BB", "DD"))
sample$Plant_species = factor(sample$Plant_species)
sample$Plant = factor(sample$Plant, levels = c("Unplanted", "ASpeltoides", "ATauschii", "Turartu", "Tmonococcum", "Tdicoccoides", "Tcarthilicum", "Tpolonicum", "Tturanicum", "Tmacha", "Tspelta", "Chidham", "Red Llammas", "Cadenza", "Victor", "Avalon", "Hereward", "Malacca", "Gallant", "Crusoe"))
sample$PlantFertilizationNiche = factor(sample$PlantFertilizationNiche)
sample$Block = factor(sample$Block)


#final object
physeq <- phyloseq(otu.table, taxonomy.table, sample, tree)
physeq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10298 taxa and 285 samples ]
#sample_data() Sample Data:       [ 285 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 10298 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10298 tips and 10022 internal nodes ]


#Rarefaction curve to check data
library(vegan)
otu.rarecurve = rarecurve(as.data.frame(t(otu_table(physeq))), step = 10000, label = F)
#Here, there are samples that only have one observed species, which can be attributed to the sequencing run.
#These samples will be removed.


#1.3. Remove samples with only one observed species:
min.lib<-min(sample_sums(physeq)) #81
physeq.rar<-rarefy_even_depth(physeq, min.lib, replace=TRUE, rngseed = 9242)
#`set.seed(9242)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(9242); .Random.seed` for the full vector
#...
#7110OTUs were removed because they are no longer 
#present in any sample after random subsampling
#
#...
observed.species <- microbiome::alpha(physeq.rar, index = c("observed"))
obs <- sample_data(observed.species)
physeq1 <- merge_phyloseq(physeq, obs)
physeq1 <- subset_samples(physeq1, observed>1)
physeq1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10298 taxa and 248 samples ]
#sample_data() Sample Data:       [ 248 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 10298 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10298 tips and 10022 internal nodes ]



#1.4. Filter blank samples from phyloseq object
sample_data(physeq1)$is.neg <- sample_data(physeq1)$Blanks == "Sample"
library(decontam) #v1.16.0
contamdf.prev <- isContaminant(physeq1, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
#FALSE  TRUE 
#10205    93
physeq2 <- prune_taxa(!contamdf.prev$contaminant, physeq1)
physeq2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10205 taxa and 248 samples ]
#sample_data() Sample Data:       [ 248 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 10205 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10205 tips and 9936 internal nodes ]


# Remove all non-samples (mocks and blanks) from main phyloseq object
physeq3 <- subset_samples(physeq2, SampleType == "Sample")
physeq3 <- filter_taxa(physeq3, function(x) sum(x) > 0, TRUE)
physeq3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10193 taxa and 244 samples ]
#sample_data() Sample Data:       [ 244 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 10193 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10193 tips and 9924 internal nodes ]




#1.5. Remove reads not presence in at least 3 replicates

source("filter_by_presence.R", local = TRUE)

physeq4 <- filter_by_presence(physeq3, "PlantFertilizationNiche", 75)
physeq4
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1326 taxa and 244 samples ]
#sample_data() Sample Data:       [ 244 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1326 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1326 tips and 1313 internal nodes ]



#1.6. Check and filter ASVs unclassified at phylum level:

#get phylum levels
get_taxa_unique(physeq4, "Phylum")
#[1] "Proteobacteria"    "Actinobacteriota"  "Cyanobacteria"     "Bdellovibrionota"  "Acidobacteriota"  
#[6] "Bacteroidota"      "Firmicutes"        "Elusimicrobiota"   "Armatimonadota"    "Chloroflexi"      
#[11] "Spirochaetota"     ""                  "Latescibacterota"  "FCPU426"           "WPS_2"            
#[16] "Methylomirabilota" "Planctomycetota"   "Verrucomicrobiota" "Gemmatimonadota"   "Nitrospirota"     
#[21] "Patescibacteria"   "Dependentiae"      "Deinococcota"      "Myxococcota"       "RCP2_54"          
#[26] "Desulfobacterota"  "Fibrobacterota"  

#Remove unclassified reads at Phylum level:
physeq5 <- subset_taxa(physeq4, !is.na(Phylum) & !Phylum %in% c(""))
physeq5
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1320 taxa and 244 samples ]
#sample_data() Sample Data:       [ 244 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1320 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1320 tips and 1307 internal nodes ]
saveRDS(physeq5, "physeq_CI.rds")



#1.7. Normalize feature table for beta diversity calculations using DESeq2, 
library(DESeq2) #v1.36.0

#Convert phyloseq object to deseq object
dds = phyloseq_to_deseq2(physeq5, ~  Fertilization+Soil)

#Normalize with regularized logarithmic transformation
dds <- dds[ rowSums(counts(dds)) > 3, ] #filter low count rows
#estimate size factors
cts <- counts(dds) 
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)

#Variance stabilizing transformation
NormVST<- t(assay(varianceStabilizingTransformation(dds, blind=TRUE)))
#Set negative values to zero
NormVST[NormVST < 0] <- 0.0


#Normalised ASV table
physeq.norm = physeq5
otu_table(physeq.norm) <- otu_table(NormVST, FALSE)
otu.norm = as(otu_table(physeq.norm), "matrix")
otu.norm = t(otu.norm) # transpose otu table
otu.norm.table <- otu_table(otu.norm,taxa_are_rows=TRUE,errorIfNULL = TRUE)


#Merge normalised ASV table to phyloseq
physeq.norm <- phyloseq(otu.norm.table, taxonomy.table, sample, tree)
physeq.norm
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1320 taxa and 244 samples ]
#sample_data() Sample Data:       [ 244 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1320 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1320 tips and 1307 internal nodes ]
saveRDS(physeq.norm, "physeq_norm_CI.rds")


# Create table, number of features for each phylum
table(tax_table(physeq.norm)[, "Phylum"], exclude = NULL)
#Acidobacteriota  Actinobacteriota    Armatimonadota      Bacteroidota  Bdellovibrionota 
#            158               110                 3               233                12 
#Chloroflexi     Cyanobacteria      Deinococcota      Dependentiae  Desulfobacterota 
#        110                29                 1                 2                12 
#Elusimicrobiota           FCPU426    Fibrobacterota        Firmicutes   Gemmatimonadota 
#              7                 2                 6                39                20 
#Latescibacterota Methylomirabilota       Myxococcota      Nitrospirota   Patescibacteria 
#               7                 6                19                 9                46 
#Planctomycetota    Proteobacteria           RCP2_54     Spirochaetota Verrucomicrobiota 
#              3               415                 2                 3                54 
#WPS_2 
#   12 






####### 2. Phyloseq object for culture-dependent amplicon sequences #######


#Input files: Files produced by qiime2 (otu_table_CD.txt, taxonomy_CD.txt, tree_CD.nwk) and sample metadata (metadata_CD.txt).



#2.1. Importing the files into R:
otu <- read.table("asv_abundances_CD.txt", header=TRUE, row.names=1)
taxonomy <- read.table(file="asv_taxonomy_CD.txt", sep="\t", header=TRUE, fill=TRUE)
sample <- read.table(file="sample_data_CD.txt", sep="\t",header=TRUE, row.names = 1)
tree  <- read_newick_phylo(file="tree_CD.nwk",simplify=FALSE)



#2.2. Creating the phyloseq object:

#asv table
otu <- as.matrix(otu)
otu.table <- otu_table(otu,taxa_are_rows=TRUE,errorIfNULL = TRUE)

#taxonomy table
row.names(taxonomy) <- taxonomy$ASVs
taxonomy <- as.matrix(taxonomy)
taxonomy.table <- tax_table(taxonomy)

#sample table
sample=sample_data(sample,errorIfNULL = TRUE)
#define factors and reorder levels
sample$Fertilization = factor(sample$Fertilization, levels = c("Non-fertilized", "Fertilized"))
sample$Soil = factor(sample$Soil, levels = c("Unplanted", "Rhizosphere", "Rhizoplane"))
sample$Isolate_function = factor(sample$Isolate_function, levels = c("non-PGPR", "PGPR"))
sample$Ancestral_class = factor(sample$Ancestral_class, levels = c("Unplanted","Wild ancestor", "Traditional cultivar", "Commercial cross"))
sample$Ploidy = factor(sample$Ploidy, levels = c("Unplanted", "Diploid", "Tetraploid", "Hexaploid"))
sample$Genome = factor(sample$Genome, levels = c("Unplanted", "AA","BB", "DD" ,"AABB", "AABBDD"))
sample$Plant_species = factor(sample$Plant_species)
sample$Plant = factor(sample$Plant)
sample$PlantFertilizationPGPR = factor(sample$PlantFertilizationPGPR)


#final object
physeq <- phyloseq(otu.table, taxonomy.table, sample, tree)
physeq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6440 taxa and 196 samples ]
#sample_data() Sample Data:       [ 196 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 6440 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 6440 tips and 5678 internal nodes ]


#Rarefaction curve to check data
library(vegan)
otu.rarecurve = rarecurve(as.data.frame(t(otu_table(physeq))), step = 10000, label = F)


#2.3. Filtering blank samples for phyloseq object
sample_data(physeq)$is.neg <- sample_data(physeq)$Blanks == "Sample"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
#FALSE  TRUE 
#6426    14 
physeq1 <- prune_taxa(!contamdf.prev$contaminant, physeq)
physeq1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6426 taxa and 196 samples ]
#sample_data() Sample Data:       [ 196 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 6426 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 6426 tips and 5673 internal nodes ]


# Remove all non samples (mocks and blanks) from main phyloseq object
physeq2 <- subset_samples(physeq1, SampleType == "Sample")
physeq2 <- filter_taxa(physeq2, function(x) sum(x) > 0, TRUE)
physeq2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6406 taxa and 193 samples ]
#sample_data() Sample Data:       [ 193 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 6406 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 6406 tips and 5653 internal nodes ]



#PGPR vs. non-PGPR ASVs
library(ggplot2)
library(MiscMetabar)
library(grid)
venn <- venn_phyloseq(physeq2,
                      "Isolate_function",
                      min_nb_seq = 0,
                      print_values = T) +
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("#edae49", "#66a182"))
#non-PGPR = 3037
#nonPGPR&PGPR = 1639
#PGPR = 1730


#2.4. Remove reads not presence in at least 2 out of 3 replicates

#filter by presence 
physeq3 <- filter_by_presence(physeq2, "PlantFertilizationPGPR", 66)
physeq3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1339 taxa and 193 samples ]
#sample_data() Sample Data:       [ 193 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1339 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1339 tips and 853 internal nodes ]
saveRDS(physeq3, "physeq_CD.rds")


#get phylum levels
get_taxa_unique(physeq3, "Phylum")
#[1] "Spirochaetota"    "Nitrospirota"     "Chloroflexi"      "Actinobacteriota"
#[5] "Acidobacteriota"  "Desulfobacterota" "Firmicutes"       "Proteobacteria"  
#[9] "Bacteroidota" 



# Create table, number of features for each phylum
table(tax_table(physeq3)[, "Phylum"], exclude = NULL)
#Acidobacteriota Actinobacteriota     Bacteroidota      Chloroflexi Desulfobacterota 
#              1              228              119                1                2 
#Firmicutes     Nitrospirota   Proteobacteria    Spirochaetota 
#       161                2              824                1 


#Class abundances
physeq.CD <- readRDS("physeq_CD.rds")
library("metagMisc") #v.0.0.4
taxa.abd <- phyloseq_to_df(physeq.CD, addtax = T, addtot = F, addmaxrank = F,sorting = "abundance")
writexl::write_xlsx(taxa.abd, "taxa.abd.CD.raw.xlsx")



################### END ###################