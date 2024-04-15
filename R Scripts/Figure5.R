############################################################
#R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################
#R 4.2.2
#setwd()

library(DESeq2); packageVersion("DESeq2") #v1.36.0
library(phyloseq); packageVersion("phyloseq") #v1.42.0
library(ggplot2); packageVersion("ggplot2") #v3.4.2


#Input files: Files produced by phyloseq (phyloseq object)

#Directories:
dir.create("Figure5")
dir.create("Figure5/Data")
dir.create("Figure5/Plots")

#Custom functions
#tern_e function from: https://github.com/BulgarelliD-Lab/Microbiota_mapping/blob/main/QRMC-3HS_Fig1_SFig3/tern_e.R
#tern_e_modifed: modified for aesthetics including point shape and axis labels.

#Abbreviations:
## CD: culture-dependent
## PGPR: plant growth-promoting rhizobacteria
## Non-PGPR: isolates with no plant growth-promoting phenotype
## FRP: fertilized rhizoplane
## NFRP: non-fertilized rhizoplane
## ASV: amplicon sequence variant


# Load data
CD_ASV <- readRDS("physeq_CD.rds")
CD_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1339 taxa and 193 samples ]
#sample_data() Sample Data:       [ 193 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1339 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1339 tips and 853 internal nodes ]

#Four subsets:
## PGPR isolated from the rhizoplane of fertilized wheats
## PGPR isolated from the rhizoplane of non-fertilized wheats
## Non-PGPR isolated from the rhizoplane of fertilized wheats
## Non-PGPR isolated from the rhizoplane of non-fertilized wheats


# Subset data
CD_PGPR_FRP_ASV <- subset_samples(CD_ASV, Isolate_function =="PGPR" & Fertilization =="Fertilized" & Soil == "Rhizoplane")
CD_PGPR_FRP_ASV <- filter_taxa(CD_PGPR_FRP_ASV, function(x) sum(x) > 0, TRUE)
CD_PGPR_FRP_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1121 taxa and 32 samples ]
#sample_data() Sample Data:       [ 32 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1121 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1121 tips and 702 internal nodes ]

CD_PGPR_NFRP_ASV <- subset_samples(CD_ASV, Isolate_function =="PGPR" & Fertilization =="Non-fertilized" & Soil == "Rhizoplane")
CD_PGPR_NFRP_ASV <- filter_taxa(CD_PGPR_NFRP_ASV, function(x) sum(x) > 0, TRUE)
CD_PGPR_NFRP_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1105 taxa and 52 samples ]
#sample_data() Sample Data:       [ 52 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1105 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1105 tips and 670 internal nodes ]


CD_NONPGPR_FRP_ASV <- subset_samples(CD_ASV, Isolate_function =="non-PGPR" & Fertilization =="Fertilized" & Soil == "Rhizoplane")
CD_NONPGPR_FRP_ASV <- filter_taxa(CD_NONPGPR_FRP_ASV, function(x) sum(x) > 0, TRUE)
CD_NONPGPR_FRP_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1142 taxa and 52 samples ]
#sample_data() Sample Data:       [ 52 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1142 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1142 tips and 690 internal nodes ]


CD_NONPGPR_NFRP_ASV <- subset_samples(CD_ASV, Isolate_function =="non-PGPR" & Fertilization =="Non-fertilized" & Soil == "Rhizoplane")
CD_NONPGPR_NFRP_ASV <- filter_taxa(CD_NONPGPR_NFRP_ASV, function(x) sum(x) > 0, TRUE)
CD_NONPGPR_NFRP_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1152 taxa and 51 samples ]
#sample_data() Sample Data:       [ 51 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1152 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1152 tips and 689 internal nodes ]



#Transform data to account for zeros
CD_PGPR_FRP_ASV<-transform_sample_counts(CD_PGPR_FRP_ASV, function(x) (x+1))
CD_PGPR_NFRP_ASV<-transform_sample_counts(CD_PGPR_NFRP_ASV, function(x) (x+1))
CD_NONPGPR_FRP_ASV <- transform_sample_counts(CD_NONPGPR_FRP_ASV, function(x) (x+1))
CD_NONPGPR_NFRP_ASV <- transform_sample_counts(CD_NONPGPR_NFRP_ASV, function(x) (x+1))




#DESeq2
ddsTAXA_CD_PGPR_FRP <- phyloseq_to_deseq2(CD_PGPR_FRP_ASV, ~ Ploidy)
#converting counts to integer mode
ddsTAXA_CD_PGPR_NFRP <- phyloseq_to_deseq2(CD_PGPR_NFRP_ASV, ~ Ploidy)
#converting counts to integer mode
ddsTAXA_CD_NONPGPR_FRP <- phyloseq_to_deseq2(CD_NONPGPR_FRP_ASV, ~ Ploidy)
#converting counts to integer mode
ddsTAXA_CD_NONPGPR_NFRP <- phyloseq_to_deseq2(CD_NONPGPR_NFRP_ASV, ~ Ploidy)
#converting counts to integer mode


#Run deseq - local regression fit automatically substituted
ddsTAXA_CD_PGPR_FRP2 = DESeq(ddsTAXA_CD_PGPR_FRP, test="Wald", fitType="local")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 366 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing
ddsTAXA_CD_PGPR_NFRP2 = DESeq(ddsTAXA_CD_PGPR_NFRP, test="Wald", fitType="local")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 378 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing
ddsTAXA_CD_NONPGPR_FRP2 = DESeq(ddsTAXA_CD_NONPGPR_FRP, test="Wald", fitType="local")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 385 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing
ddsTAXA_CD_NONPGPR_NFRP2 = DESeq(ddsTAXA_CD_NONPGPR_NFRP, test="Wald", fitType="local")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 345 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing



#Check results
summary(ddsTAXA_CD_PGPR_FRP2)
#[1] "DESeqDataSet object of length 1121 with 27 metadata columns"
summary(ddsTAXA_CD_PGPR_NFRP2)
#[1] "DESeqDataSet object of length 1105 with 27 metadata columns"
summary(ddsTAXA_CD_NONPGPR_FRP2)
#[1] "DESeqDataSet object of length 1142 with 27 metadata columns"
summary(ddsTAXA_CD_NONPGPR_NFRP2)
#[1] "DESeqDataSet object of length 1152 with 27 metadata columns"



#Intercepts
resultsNames(ddsTAXA_CD_PGPR_FRP2)
#[1] "Intercept"                    "Ploidy_Tetraploid_vs_Diploid" "Ploidy_Hexaploid_vs_Diploid" 
resultsNames(ddsTAXA_CD_PGPR_NFRP2)
#[1] "Intercept"                    "Ploidy_Tetraploid_vs_Diploid" "Ploidy_Hexaploid_vs_Diploid" 
resultsNames(ddsTAXA_CD_NONPGPR_FRP2)
#[1] "Intercept"                    "Ploidy_Tetraploid_vs_Diploid" "Ploidy_Hexaploid_vs_Diploid" 
resultsNames(ddsTAXA_CD_NONPGPR_NFRP2)
#[1] "Intercept"                    "Ploidy_Tetraploid_vs_Diploid" "Ploidy_Hexaploid_vs_Diploid" 
 
###### Tetraploid vs Hexaploid still missing #####



#Export results
ddsTAXA2_CD_PGPR_FRP_res <- results(ddsTAXA_CD_PGPR_FRP2)
ddsTAXA2_CD_PGPR_NFRP_res <- results(ddsTAXA_CD_PGPR_NFRP2)
ddsTAXA2_CD_NONPGPR_FRP_res <- results(ddsTAXA_CD_NONPGPR_FRP2)
ddsTAXA2_CD_NONPGPR_NFRP_res <- results(ddsTAXA_CD_NONPGPR_NFRP2)


#check
summary(ddsTAXA2_CD_PGPR_FRP_res)
#
#out of 1121 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 40, 3.6%
#LFC < 0 (down)     : 151, 13%
#outliers [1]       : 93, 8.3%
#low counts [2]     : 152, 14%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#
summary(ddsTAXA2_CD_PGPR_NFRP_res)
#
#out of 1105 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 63, 5.7%
#LFC < 0 (down)     : 71, 6.4%
#outliers [1]       : 0, 0%
#low counts [2]     : 150, 14%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#
summary(ddsTAXA2_CD_NONPGPR_FRP_res)
#
#out of 1142 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 133, 12%
#LFC < 0 (down)     : 64, 5.6%
#outliers [1]       : 0, 0%
#low counts [2]     : 45, 3.9%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#
summary(ddsTAXA2_CD_NONPGPR_NFRP_res)
#
#out of 1152 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 52, 4.5%
#LFC < 0 (down)     : 54, 4.7%
#outliers [1]       : 0, 0%
#low counts [2]     : 90, 7.8%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#




# Perform contrasts between groups

#PGPGR F RP
ddsTAXA2_CD_PGPR_FRP_TvD <- results(ddsTAXA_CD_PGPR_FRP2, contrast = c("Ploidy", "Tetraploid", "Diploid"), alpha=0.05)
ddsTAXA2_CD_PGPR_FRP_TvB <- results(ddsTAXA_CD_PGPR_FRP2, contrast = c("Ploidy", "Hexaploid","Tetraploid"), alpha=0.05)
ddsTAXA2_CD_PGPR_FRP_DvB <- results(ddsTAXA_CD_PGPR_FRP2, contrast = c("Ploidy", "Hexaploid", "Diploid"), alpha=0.05)


#PGPR NF RP
ddsTAXA2_CD_PGPR_NFRP_TvD <- results(ddsTAXA_CD_PGPR_NFRP2, contrast = c("Ploidy", "Tetraploid", "Diploid"), alpha=0.05)
ddsTAXA2_CD_PGPR_NFRP_TvB <- results(ddsTAXA_CD_PGPR_NFRP2, contrast = c("Ploidy", "Hexaploid","Tetraploid"), alpha=0.05)
ddsTAXA2_CD_PGPR_NFRP_DvB <- results(ddsTAXA_CD_PGPR_NFRP2, contrast = c("Ploidy", "Hexaploid", "Diploid"), alpha=0.05)


#Non PGPGR F RP
ddsTAXA2_CD_NONPGPR_FRP_TvD <- results(ddsTAXA_CD_NONPGPR_FRP2, contrast = c("Ploidy", "Tetraploid", "Diploid"), alpha=0.05)
ddsTAXA2_CD_NONPGPR_FRP_TvB <- results(ddsTAXA_CD_NONPGPR_FRP2, contrast = c("Ploidy","Hexaploid","Tetraploid"), alpha=0.05)
ddsTAXA2_CD_NONPGPR_FRP_DvB <- results(ddsTAXA_CD_NONPGPR_FRP2, contrast = c("Ploidy", "Hexaploid", "Diploid"), alpha=0.05)


#Non PGPGR NF RP
ddsTAXA2_CD_NONPGPR_NFRP_TvD <- results(ddsTAXA_CD_NONPGPR_NFRP2, contrast = c("Ploidy", "Tetraploid", "Diploid"), alpha=0.05)
ddsTAXA2_CD_NONPGPR_NFRP_TvB <- results(ddsTAXA_CD_NONPGPR_NFRP2, contrast = c("Ploidy","Hexaploid","Tetraploid"), alpha=0.05)
ddsTAXA2_CD_NONPGPR_NFRP_DvB <- results(ddsTAXA_CD_NONPGPR_NFRP2, contrast = c("Ploidy", "Hexaploid", "Diploid"), alpha=0.05)



#Matrix of results

alpha <- 0.05

#CD PGPR FRP
allTAXA_CD_PGPR_FRP <- results(ddsTAXA_CD_PGPR_FRP2)
allTAXA_CD_PGPR_FRP = cbind(as(allTAXA_CD_PGPR_FRP, "data.frame"), as(tax_table(CD_PGPR_FRP_ASV)[rownames(allTAXA_CD_PGPR_FRP), ], "matrix"))

allTAXA_CD_PGPR_FRP_TvD = cbind(as(ddsTAXA2_CD_PGPR_FRP_TvD, "data.frame"), as(tax_table(CD_PGPR_FRP_ASV)[rownames(ddsTAXA2_CD_PGPR_FRP_TvD), ], "matrix"))
allTAXA_CD_PGPR_FRP_TvB = cbind(as(ddsTAXA2_CD_PGPR_FRP_TvB, "data.frame"), as(tax_table(CD_PGPR_FRP_ASV)[rownames(ddsTAXA2_CD_PGPR_FRP_TvB), ], "matrix"))
allTAXA_CD_PGPR_FRP_DvB = cbind(as(ddsTAXA2_CD_PGPR_FRP_DvB, "data.frame"), as(tax_table(CD_PGPR_FRP_ASV)[rownames(ddsTAXA2_CD_PGPR_FRP_DvB), ], "matrix"))

sigTAXA_CD_PGPR_FRP_TvD <- subset(allTAXA_CD_PGPR_FRP_TvD, padj < alpha)
sigTAXA_CD_PGPR_FRP_TvB <- subset(allTAXA_CD_PGPR_FRP_TvB, padj < alpha)
sigTAXA_CD_PGPR_FRP_DvB <- subset(allTAXA_CD_PGPR_FRP_DvB, padj < alpha)

writexl::write_xlsx(sigTAXA_CD_PGPR_FRP_TvD, "Figure5/Data/deseq-ternary-CD-PGPR-FRP-TvD.xlsx")
writexl::write_xlsx(sigTAXA_CD_PGPR_FRP_TvB, "Figure5/Data/deseq-ternary-CD-PGPR-FRP-TvB.xlsx")
writexl::write_xlsx(sigTAXA_CD_PGPR_FRP_DvB, "Figure5/Data/deseq-ternary-CD-PGPR-FRP-DvB.xlsx")


#CD PGPR NFRP
allTAXA_CD_PGPR_NFRP <- results(ddsTAXA_CD_PGPR_NFRP2)
allTAXA_CD_PGPR_NFRP = cbind(as(allTAXA_CD_PGPR_NFRP, "data.frame"), as(tax_table(CD_PGPR_NFRP_ASV)[rownames(allTAXA_CD_PGPR_NFRP), ], "matrix"))

allTAXA_CD_PGPR_NFRP_TvD = cbind(as(ddsTAXA2_CD_PGPR_NFRP_TvD, "data.frame"), as(tax_table(CD_PGPR_NFRP_ASV)[rownames(ddsTAXA2_CD_PGPR_NFRP_TvD), ], "matrix"))
allTAXA_CD_PGPR_NFRP_TvB = cbind(as(ddsTAXA2_CD_PGPR_NFRP_TvB, "data.frame"), as(tax_table(CD_PGPR_NFRP_ASV)[rownames(ddsTAXA2_CD_PGPR_NFRP_TvB), ], "matrix"))
allTAXA_CD_PGPR_NFRP_DvB = cbind(as(ddsTAXA2_CD_PGPR_NFRP_DvB, "data.frame"), as(tax_table(CD_PGPR_NFRP_ASV)[rownames(ddsTAXA2_CD_PGPR_NFRP_DvB), ], "matrix"))

sigTAXA_CD_PGPR_NFRP_TvD <- subset(allTAXA_CD_PGPR_NFRP_TvD, padj < alpha)
sigTAXA_CD_PGPR_NFRP_TvB <- subset(allTAXA_CD_PGPR_NFRP_TvB, padj < alpha)
sigTAXA_CD_PGPR_NFRP_DvB <- subset(allTAXA_CD_PGPR_NFRP_DvB, padj < alpha)

writexl::write_xlsx(sigTAXA_CD_PGPR_NFRP_TvD, "Figure5/Data/deseq-ternary-CD-PGPR-NFRP-TvD.xlsx")
writexl::write_xlsx(sigTAXA_CD_PGPR_NFRP_TvB, "Figure5/Data/deseq-ternary-CD-PGPR-NFRP-TvB.xlsx")
writexl::write_xlsx(sigTAXA_CD_PGPR_NFRP_DvB, "Figure5/Data/deseq-ternary-CD-PGPR-NFRP-DvB.xlsx")


#CD NONPGPR FRP
allTAXA_CD_NONPGPR_FRP <- results(ddsTAXA_CD_NONPGPR_FRP2)
allTAXA_CD_NONPGPR_FRP = cbind(as(allTAXA_CD_NONPGPR_FRP, "data.frame"), as(tax_table(CD_NONPGPR_FRP_ASV)[rownames(allTAXA_CD_NONPGPR_FRP), ], "matrix"))

allTAXA_CD_NONPGPR_FRP_TvD = cbind(as(ddsTAXA2_CD_NONPGPR_FRP_TvD, "data.frame"), as(tax_table(CD_NONPGPR_FRP_ASV)[rownames(ddsTAXA2_CD_NONPGPR_FRP_TvD), ], "matrix"))
allTAXA_CD_NONPGPR_FRP_TvB = cbind(as(ddsTAXA2_CD_NONPGPR_FRP_TvB, "data.frame"), as(tax_table(CD_NONPGPR_FRP_ASV)[rownames(ddsTAXA2_CD_NONPGPR_FRP_TvB), ], "matrix"))
allTAXA_CD_NONPGPR_FRP_DvB = cbind(as(ddsTAXA2_CD_NONPGPR_FRP_DvB, "data.frame"), as(tax_table(CD_NONPGPR_FRP_ASV)[rownames(ddsTAXA2_CD_NONPGPR_FRP_DvB), ], "matrix"))

sigTAXA_CD_NONPGPR_FRP_TvD <- subset(allTAXA_CD_NONPGPR_FRP_TvD, padj < alpha)
sigTAXA_CD_NONPGPR_FRP_TvB <- subset(allTAXA_CD_NONPGPR_FRP_TvB, padj < alpha)
sigTAXA_CD_NONPGPR_FRP_DvB <- subset(allTAXA_CD_NONPGPR_FRP_DvB, padj < alpha)

writexl::write_xlsx(sigTAXA_CD_NONPGPR_FRP_TvD, "Figure5/Data/deseq-ternary-CD-NONPGPR-FRP-TvD.xlsx")
writexl::write_xlsx(sigTAXA_CD_NONPGPR_FRP_TvB, "Figure5/Data/deseq-ternary-CD-NONPGPR-FRP-TvB.xlsx")
writexl::write_xlsx(sigTAXA_CD_NONPGPR_FRP_DvB, "Figure5/Data/deseq-ternary-CD-NONPGPR-FRP-DvB.xlsx")


#CD NONPGPR NFRP
allTAXA_CD_NONPGPR_NFRP <- results(ddsTAXA_CD_NONPGPR_NFRP2)
allTAXA_CD_NONPGPR_NFRP = cbind(as(allTAXA_CD_NONPGPR_NFRP, "data.frame"), as(tax_table(CD_NONPGPR_NFRP_ASV)[rownames(allTAXA_CD_NONPGPR_NFRP), ], "matrix"))

allTAXA_CD_NONPGPR_NFRP_TvD = cbind(as(ddsTAXA2_CD_NONPGPR_NFRP_TvD, "data.frame"), as(tax_table(CD_NONPGPR_NFRP_ASV)[rownames(ddsTAXA2_CD_NONPGPR_NFRP_TvD), ], "matrix"))
allTAXA_CD_NONPGPR_NFRP_TvB = cbind(as(ddsTAXA2_CD_NONPGPR_NFRP_TvB, "data.frame"), as(tax_table(CD_NONPGPR_NFRP_ASV)[rownames(ddsTAXA2_CD_NONPGPR_NFRP_TvB), ], "matrix"))
allTAXA_CD_NONPGPR_NFRP_DvB = cbind(as(ddsTAXA2_CD_NONPGPR_NFRP_DvB, "data.frame"), as(tax_table(CD_NONPGPR_NFRP_ASV)[rownames(ddsTAXA2_CD_NONPGPR_NFRP_DvB), ], "matrix"))

sigTAXA_CD_NONPGPR_NFRP_TvD <- subset(allTAXA_CD_NONPGPR_NFRP_TvD, padj < alpha)
sigTAXA_CD_NONPGPR_NFRP_TvB <- subset(allTAXA_CD_NONPGPR_NFRP_TvB, padj < alpha)
sigTAXA_CD_NONPGPR_NFRP_DvB <- subset(allTAXA_CD_NONPGPR_NFRP_DvB, padj < alpha)

writexl::write_xlsx(sigTAXA_CD_NONPGPR_NFRP_TvD, "Figure5/Data/deseq-ternary-CD-NONPGPR-NFRP-TvD.xlsx")
writexl::write_xlsx(sigTAXA_CD_NONPGPR_NFRP_TvB, "Figure5/Data/deseq-ternary-CD-NONPGPR-NFRP-TvB.xlsx")
writexl::write_xlsx(sigTAXA_CD_NONPGPR_NFRP_DvB, "Figure5/Data/deseq-ternary-CD-NONPGPR-NFRP-DvB.xlsx")



# Ternary ####

## CD PFPR FRP

# need to get the means of each compartment for the ternary plots 
TAXA2baseMeanPerLvl <- sapply( levels(ddsTAXA_CD_PGPR_FRP2$Ploidy), function(lvl) rowMeans( counts(ddsTAXA_CD_PGPR_FRP2,normalized=TRUE)[,ddsTAXA_CD_PGPR_FRP2$Ploidy == lvl, drop=F] ) )
head(TAXA2baseMeanPerLvl)
#                                  Diploid Tetraploid Hexaploid
#6377df2d1c652fa1e9e826dcac9fe2f8 1.151607   1.179785  1.867197
#a300787dbb86a1d0c7f1e1f6a189c6e1 1.151607   1.179785  2.817548
#51704f1b33d5221bd143b8ed286416f2 1.151607   1.179785  1.740483
#950c70d1350c2f44e0bcbe9c2b17bec9 1.151607   1.179785  2.057267
#ae42e69ef0c6f4763aaa358d5f544076 1.151607   1.179785  2.310694
#a4a832c8752233c8b3e0efb8efb41d32 1.151607   1.179785  3.458651
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
#   Diploid Tetraploid  Hexaploid 
#  10636.87   10224.64   10081.37 


#Define significant ASVs for ternary plot
TetraploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_PGPR_FRP_TvB, padj<0.05 & log2FoldChange<0)), 
  rownames(subset(ddsTAXA2_CD_PGPR_FRP_TvD, padj<0.05 & log2FoldChange>0))))

DiploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_PGPR_FRP_TvD, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CD_PGPR_FRP_DvB, padj<0.05 & log2FoldChange<0))))

HexaploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_PGPR_FRP_TvB, padj<0.05 & log2FoldChange>0)),
  rownames(subset(ddsTAXA2_CD_PGPR_FRP_DvB, padj<0.05 & log2FoldChange>0))))


#get phylum levels
get_taxa_unique(CD_PGPR_FRP_ASV, "Phylum")
#[1] "Nitrospirota"     "Chloroflexi"      "Acidobacteriota"  "Desulfobacterota" "Firmicutes"       "Actinobacteriota"
#[7] "Proteobacteria"   "Bacteroidota"



#Colour coding based on phyla
Acidobacteriota <- unique(rownames(subset(allTAXA_CD_PGPR_FRP, Phylum=="Acidobacteriota")))
Actinobacteria <- unique(rownames(subset(allTAXA_CD_PGPR_FRP, Phylum=="Actinobacteriota")))
Bacteroidetes <- unique(rownames(subset(allTAXA_CD_PGPR_FRP, Phylum=="Bacteroidota")))
Chloroflexi <- unique(rownames(subset(allTAXA_CD_PGPR_FRP, Phylum=="Chloroflexi")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CD_PGPR_FRP, Phylum=="Desulfobacterota")))
Firmicutes <- unique(rownames(subset(allTAXA_CD_PGPR_FRP, Phylum=="Firmicutes")))
Nitrospirota <- unique(rownames(subset(allTAXA_CD_PGPR_FRP, Phylum=="Nitrospirota")))
Proteobacteria <- unique(rownames(subset(allTAXA_CD_PGPR_FRP, Phylum=="Proteobacteria")))


#Colours
red<-"#E64B35FF"
blue<-"#4DBBD5B2"
turquoise<-"#6F99ADFF"
navy<-"#3C5488FF"
pale_red<-"#F39B7FFF"
light_blue<-"#8491B4FF"
brown<-"#7E6148FF"
dark_yellow<-"#edae49"



#Colour codes for plotting
TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% Acidobacteriota, pale_red,"#004d40")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_colors[Actinobacteria] <- red
TAXA2ternary_colors[Bacteroidetes] <- blue
TAXA2ternary_colors[Chloroflexi] <- light_blue
TAXA2ternary_colors[Firmicutes] <- turquoise
TAXA2ternary_colors[Nitrospirota] <- brown
TAXA2ternary_colors[Proteobacteria] <- navy
TAXA2ternary_colors[Desulfobacterota] <- dark_yellow



#Defining shape of points based on significant ASVs

#Based on ploidy
##TAXA2ternary_shapes <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% DiploidTaxa, 17, 19)
##names(TAXA2ternary_shapes) <- rownames(TAXA2baseMeanPerLvl)
##TAXA2ternary_shapes[TetraploidTaxa] <- 15
##TAXA2ternary_shapes[HexaploidTaxa] <- 18

#Based on significance
TAXA2ternary_shapes <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% DiploidTaxa, 17, 19)
names(TAXA2ternary_shapes) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_shapes[TetraploidTaxa] <- 17
TAXA2ternary_shapes[HexaploidTaxa] <- 17


### Plotting ternary with colours and shapes
library(grid); packageVersion("grid") #v4.2.2
source("tern_e_modified.R", local = TRUE)
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)

# Create plot and save as pdf
pdf(file = "Figure5/Plots/deseq-ternary-CD-PGPR-FRP.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)
dev.off()
#RStudioGD 
#2 

#####################

## CD PFPR NFRP

# need to get the means of each compartment for the ternary plots 
TAXA2baseMeanPerLvl <- sapply( levels(ddsTAXA_CD_PGPR_NFRP2$Ploidy), function(lvl) rowMeans( counts(ddsTAXA_CD_PGPR_NFRP2,normalized=TRUE)[,ddsTAXA_CD_PGPR_NFRP2$Ploidy == lvl, drop=F] ) )
head(TAXA2baseMeanPerLvl)
#                                  Diploid Tetraploid Hexaploid
#deee6403c3bbe581ab814bc36cc322f3 1.141647   1.142957  2.094486
#51704f1b33d5221bd143b8ed286416f2 1.141647   1.142957  2.350965
#bb9c50e936acb16f51f766fc8bcb168c 1.141647   1.142957  1.219402
#df55797fd7f0dd31b3c5db9d206387a8 1.141647   1.142957  2.372922
#1d8b86f875f974b6972a91dafef16d7b 1.141647   1.142957  1.537615
#bf2774d5c5f52c8c79f6ddf8b37db192 1.141647   1.142957  1.696721
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
#Hexaploid    Diploid Tetraploid 
# 9864.734   9275.189   9178.365


#Define significant ASVs for ternary plot
TetraploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_PGPR_NFRP_TvB, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CD_PGPR_NFRP_TvD, padj<0.05 & log2FoldChange>0))))

DiploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_PGPR_NFRP_TvD, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CD_PGPR_NFRP_DvB, padj<0.05 & log2FoldChange<0))))

HexaploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_PGPR_NFRP_TvB, padj<0.05 & log2FoldChange>0)),
  rownames(subset(ddsTAXA2_CD_PGPR_NFRP_DvB, padj<0.05 & log2FoldChange>0))))


#get phylum levels
get_taxa_unique(CD_PGPR_NFRP_ASV, "Phylum")
#[1] "Spirochaetota"    "Chloroflexi"      "Actinobacteriota" "Firmicutes"       "Proteobacteria"   "Bacteroidota" 


#Colour codes based on phyla
Actinobacteria <- unique(rownames(subset(allTAXA_CD_PGPR_NFRP, Phylum=="Actinobacteriota")))
Bacteroidetes <- unique(rownames(subset(allTAXA_CD_PGPR_NFRP, Phylum=="Bacteroidota")))
Chloroflexi <- unique(rownames(subset(allTAXA_CD_PGPR_NFRP, Phylum=="Chloroflexi")))
Firmicutes <- unique(rownames(subset(allTAXA_CD_PGPR_NFRP, Phylum=="Firmicutes")))
Proteobacteria <- unique(rownames(subset(allTAXA_CD_PGPR_NFRP, Phylum=="Proteobacteria")))
Spirochaetota <- unique(rownames(subset(allTAXA_CD_PGPR_NFRP, Phylum=="Spirochaetota")))

grey <- "grey37"

#Colour codes for plotting
TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% Actinobacteria, red, "#004d40")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_colors[Bacteroidetes] <- blue
TAXA2ternary_colors[Chloroflexi] <- light_blue
TAXA2ternary_colors[Firmicutes] <- turquoise
TAXA2ternary_colors[Proteobacteria] <- navy
TAXA2ternary_colors[Spirochaetota] <- grey


#Defining shape of points based on significant ASVs

#Based on ploidy
##TAXA2ternary_shapes <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% DiploidTaxa, 17, 19)
##names(TAXA2ternary_shapes) <- rownames(TAXA2baseMeanPerLvl)
##TAXA2ternary_shapes[TetraploidTaxa] <- 15
##TAXA2ternary_shapes[HexaploidTaxa] <- 18

#Based on significance
TAXA2ternary_shapes <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% DiploidTaxa, 17, 19)
names(TAXA2ternary_shapes) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_shapes[TetraploidTaxa] <- 17
TAXA2ternary_shapes[HexaploidTaxa] <- 17



### Plotting ternary with colours and shapes
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)

# Create plot and save as pdf
pdf(file = "Figure5/Plots/deseq-ternary-CD-PGPR-NFRP.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)
dev.off()
#RStudioGD 
#2 



#################

## CD NONPGPR FRP

# need to get the means of each compartment for the ternary plots 
TAXA2baseMeanPerLvl <- sapply( levels(ddsTAXA_CD_NONPGPR_FRP2$Ploidy), function(lvl) rowMeans( counts(ddsTAXA_CD_NONPGPR_FRP2,normalized=TRUE)[,ddsTAXA_CD_NONPGPR_FRP2$Ploidy == lvl, drop=F] ) )
head(TAXA2baseMeanPerLvl)
#                                   Diploid Tetraploid Hexaploid
#d79241eb6b2c1edf33e269087dae32b6  1.259879   1.243870  2.779890
#5901227df973943b4ad12c9149587329 20.808576   2.766226 10.842487
#5338557c9f92ef0d3fcc5e2b41bda35b 13.130115   4.357786  6.221074
#f185f16d68a0f3a3881d4fb2a0df23d6  6.228804   5.395515  7.158739
#ebfd46a060a2fffe3030ac70f15d58ac  6.081690   7.998323  5.571806
#68fef29365b75ba4d6d691356657c960  3.648426   1.243870  2.294646
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
# Hexaploid    Diploid Tetraploid 
#  11065.73   10685.35   10632.28


#Define significant ASVs for ternary plot
TetraploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_NONPGPR_FRP_TvB, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CD_NONPGPR_FRP_TvD, padj<0.05 & log2FoldChange>0))))

DiploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_NONPGPR_FRP_TvD, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CD_NONPGPR_FRP_DvB, padj<0.05 & log2FoldChange<0))))

HexaploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_NONPGPR_FRP_TvB, padj<0.05 & log2FoldChange>0)),
  rownames(subset(ddsTAXA2_CD_NONPGPR_FRP_DvB, padj<0.05 & log2FoldChange>0))))


#get phylum levels
get_taxa_unique(CD_NONPGPR_FRP_ASV, "Phylum")
#[1] "Firmicutes"       "Actinobacteriota" "Proteobacteria"   "Bacteroidota"  


#Colour codes based on phyla
Actinobacteria <- unique(rownames(subset(allTAXA_CD_NONPGPR_FRP, Phylum=="Actinobacteriota")))
Bacteroidetes <- unique(rownames(subset(allTAXA_CD_NONPGPR_FRP, Phylum=="Bacteroidota")))
Firmicutes <- unique(rownames(subset(allTAXA_CD_NONPGPR_FRP, Phylum=="Firmicutes")))
Proteobacteria <- unique(rownames(subset(allTAXA_CD_NONPGPR_FRP, Phylum=="Proteobacteria")))


#Colour codes for plotting
TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% Actinobacteria, red,"#004d40")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_colors[Bacteroidetes] <- blue
TAXA2ternary_colors[Firmicutes] <- turquoise
TAXA2ternary_colors[Proteobacteria] <- navy


#Defining shape of points based on significant ASVs

#Based on ploidy
##TAXA2ternary_shapes <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% DiploidTaxa, 17, 19)
##names(TAXA2ternary_shapes) <- rownames(TAXA2baseMeanPerLvl)
##TAXA2ternary_shapes[TetraploidTaxa] <- 15
##TAXA2ternary_shapes[HexaploidTaxa] <- 18

#Based on significance
TAXA2ternary_shapes <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% DiploidTaxa, 17, 19)
names(TAXA2ternary_shapes) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_shapes[TetraploidTaxa] <- 17
TAXA2ternary_shapes[HexaploidTaxa] <- 17



### Plotting ternary with colours and shapes
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)


# Create plot and save as pdf
pdf(file = "Figure5/Plots/deseq-ternary-CD-NONPGPR-FRP.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)
dev.off()
#RStudioGD 
#2 


#####################

## CD NONPFPR NFRP

# need to get the means of each compartment for the ternary plots 
TAXA2baseMeanPerLvl <- sapply( levels(ddsTAXA_CD_NONPGPR_NFRP2$Ploidy), function(lvl) rowMeans( counts(ddsTAXA_CD_NONPGPR_NFRP2,normalized=TRUE)[,ddsTAXA_CD_NONPGPR_NFRP2$Ploidy == lvl, drop=F] ) )
head(TAXA2baseMeanPerLvl)
#                                   Diploid Tetraploid Hexaploid
#bb9c50e936acb16f51f766fc8bcb168c  1.183644   1.176832  1.847078
#d79241eb6b2c1edf33e269087dae32b6  7.616763   1.176832  2.908258
#5901227df973943b4ad12c9149587329 12.039281  18.896524 11.146630
#5338557c9f92ef0d3fcc5e2b41bda35b 15.640476  11.237203 13.205257
#f185f16d68a0f3a3881d4fb2a0df23d6 19.792584  10.210519 10.764621
#ebfd46a060a2fffe3030ac70f15d58ac 14.640865   7.788295 11.512880
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
#  Diploid  Hexaploid Tetraploid 
#11495.546  10393.707   9917.901


#Define significant ASVs for ternary plot
TetraploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_NONPGPR_NFRP_TvB, padj<0.05 & log2FoldChange<0)), 
  rownames(subset(ddsTAXA2_CD_NONPGPR_NFRP_TvD, padj<0.05 & log2FoldChange>0))))

DiploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_NONPGPR_NFRP_TvD, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CD_NONPGPR_NFRP_DvB, padj<0.05 & log2FoldChange<0))))

HexaploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CD_NONPGPR_NFRP_TvB, padj<0.05 & log2FoldChange>0)), 
  rownames(subset(ddsTAXA2_CD_NONPGPR_NFRP_DvB, padj<0.05 & log2FoldChange>0)))) 


#get phylum levels
get_taxa_unique(CD_NONPGPR_NFRP_ASV, "Phylum")
#[1] "Actinobacteriota" "Firmicutes"       "Proteobacteria"   "Bacteroidota"  


#Colour codes based on phyla
Actinobacteria <- unique(rownames(subset(allTAXA_CD_NONPGPR_NFRP, Phylum=="Actinobacteriota")))
Bacteroidetes <- unique(rownames(subset(allTAXA_CD_NONPGPR_NFRP, Phylum=="Bacteroidota")))
Firmicutes <- unique(rownames(subset(allTAXA_CD_NONPGPR_NFRP, Phylum=="Firmicutes")))
Proteobacteria <- unique(rownames(subset(allTAXA_CD_NONPGPR_NFRP, Phylum=="Proteobacteria")))



#Colour codes for plotting
TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% Actinobacteria, red,"#004d40")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_colors[Bacteroidetes] <- blue
TAXA2ternary_colors[Firmicutes] <- turquoise
TAXA2ternary_colors[Proteobacteria] <- navy


#Defining shape of points based on significant ASVs

#Based on ploidy
##TAXA2ternary_shapes <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% DiploidTaxa, 17, 19)
##names(TAXA2ternary_shapes) <- rownames(TAXA2baseMeanPerLvl)
##TAXA2ternary_shapes[TetraploidTaxa] <- 15
##TAXA2ternary_shapes[HexaploidTaxa] <- 18

#Based on significance
TAXA2ternary_shapes <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% DiploidTaxa, 17, 19)
names(TAXA2ternary_shapes) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_shapes[TetraploidTaxa] <- 17
TAXA2ternary_shapes[HexaploidTaxa] <- 17



### Plotting ternary with colours and shapes
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)

# Create plot and save as pdf
pdf(file = "Figure5/Plots/deseq-ternary-CD-NONPGPR-NFRP.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)
dev.off()
#RStudioGD 
#2 


# make a legend for the ternary plots
library(stringr); packageVersion("stringr") #v1.5.0
library(tidyverse); packageVersion("tidyverse") #v2.0.0


### Colours
phyla <- c("Acidobacteriota", 
           "Actinobacteriota", 
           "Bacteroidota",
           "Chloroflexi", 
           "Desulfobacterota", 
           "Firmicutes", 
           "Nitrospirota",
           "Proteobacteria", 
           "Spirochaetota")

clrs <- c(pale_red,
          red,
          blue,
          light_blue,
          dark_yellow,
          turquoise,
          brown,
          navy,
          grey)


#Ploidy
##names <- c("Not significant", "Diploid", "Tetraploid", "Hexaploid")
##shape <- c(19, 17, 15, 18)


#Significance
names <- c("Not significant", 
           "Significant")

shape <- c(19,
           17)

# Check legends
plot.new()
legend("bottom",  legend = phyla, pch=16, cex=0.5, 
                bty='n', col = clrs, title = "Phyla", title.adj = 0)

legend("top",  legend = names, pch=shape, cex=0.5, 
       bty='n', title = "ASVs", title.adj = 0)


# Create legend and save as pdf
pdf(file = "Figure5/Plots/deseq-ternary-CD-legend.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot.new()
legend("bottom",  legend = phyla, pch=16, cex=0.5, 
       bty='n', col = clrs, title = "Phyla", title.adj = 0)
legend("top",  legend = names, pch=shape, cex=0.5, 
       bty='n', title = "ASVs", title.adj = 0)
dev.off()
#RStudioGD 
#2 


sessionInfo()
#R version 4.2.2 (2022-10-31)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Ventura 13.6.3
#
#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#
#locale:
#[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
#
#attached base packages:
#[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#[1] lubridate_1.9.2             forcats_1.0.0               dplyr_1.1.2                 purrr_1.0.1                
#[5] readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1                tidyverse_2.0.0            
#[9] stringr_1.5.0               ggplot2_3.4.2               phyloseq_1.42.0             DESeq2_1.36.0              
#[13] SummarizedExperiment_1.26.1 Biobase_2.58.0              MatrixGenerics_1.8.1        matrixStats_0.63.0         
#[17] GenomicRanges_1.48.0        GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.1           
#[21] BiocGenerics_0.44.0        
#
#loaded via a namespace (and not attached):
#[1] nlme_3.1-162           bitops_1.0-7           bit64_4.0.5            RColorBrewer_1.1-3    
#[5] httr_1.4.6             tools_4.2.2            utf8_1.2.3             R6_2.5.1              
#[9] vegan_2.6-4            DBI_1.1.3              mgcv_1.8-42            colorspace_2.1-0      
#[13] permute_0.9-7          rhdf5filters_1.10.0    ade4_1.7-22            withr_2.5.0           
#[17] tidyselect_1.2.0       bit_4.0.5              compiler_4.2.2         cli_3.6.1             
#[21] DelayedArray_0.22.0    scales_1.2.1           genefilter_1.78.0      digest_0.6.31         
#[25] XVector_0.38.0         pkgconfig_2.0.3        fastmap_1.1.1          rlang_1.1.1           
#[29] rstudioapi_0.14        RSQLite_2.3.1          generics_0.1.3         jsonlite_1.8.4        
#[33] BiocParallel_1.30.4    RCurl_1.98-1.12        magrittr_2.0.3         GenomeInfoDbData_1.2.9
#[37] biomformat_1.26.0      Matrix_1.5-4           Rcpp_1.0.10            munsell_0.5.0         
#[41] Rhdf5lib_1.20.0        fansi_1.0.4            ape_5.7-1              lifecycle_1.0.3       
#[45] stringi_1.7.12         MASS_7.3-60            zlibbioc_1.44.0        rhdf5_2.42.0          
#[49] plyr_1.8.8             blob_1.2.4             parallel_4.2.2         crayon_1.5.2          
#[53] lattice_0.21-8         Biostrings_2.66.0      splines_4.2.2          multtest_2.54.0       
#[57] annotate_1.74.0        hms_1.1.3              KEGGREST_1.36.3        locfit_1.5-9.7        
#[61] pillar_1.9.0           igraph_1.4.2           geneplotter_1.74.0     reshape2_1.4.4        
#[65] codetools_0.2-19       XML_3.99-0.14          glue_1.6.2             data.table_1.14.8     
#[69] png_0.1-8              vctrs_0.6.2            tzdb_0.3.0             foreach_1.5.2         
#[73] gtable_0.3.3           cachem_1.0.8           xtable_1.8-4           survival_3.5-5        
#[77] iterators_1.0.14       AnnotationDbi_1.58.0   memoise_2.0.1          writexl_1.4.2         
#[81] cluster_2.1.4          timechange_0.2.0 

#### Continued in Prism - Figure5.pzfx
##### END #####