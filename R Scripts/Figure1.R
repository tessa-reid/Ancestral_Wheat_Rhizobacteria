############################################################
#R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################
#R 4.2.2
#setwd()



#Pre-defined parameters for producing  plots:
red<-"#E64B35FF"
navy<-"#3C5488FF"

#Colour palettes
fertilizer.palette<-c(red, navy)


#Manually change colour scales on ggplot
library(ggplot2); packageVersion("ggplot2") #v3.4.2
colour.fill.fer<-scale_fill_manual(values = fertilizer.palette)
colour.point.fer<-scale_color_manual(values = fertilizer.palette)



#Set theme
theme_set(theme_bw())


#Input files: Files produced by phyloseq (phyloseq object, phyloseq object normalised by deseq)

#Directories:
dir.create("Figure1")
dir.create("Figure1/Data")
dir.create("Figure1/Plots")

suppressMessages(library(phyloseq)); packageVersion("phyloseq") #v1.42.0


#####Importing the files into R:
physeq.CI <- readRDS("physeq_CI.rds")
physeq.CI
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1320 taxa and 244 samples ]
#sample_data() Sample Data:       [ 244 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1320 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1320 tips and 1307 internal nodes ]
physeq.norm.CI <- readRDS("physeq_norm_CI.rds")
physeq.norm.CI
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1320 taxa and 244 samples ]
#sample_data() Sample Data:       [ 244 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1320 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1320 tips and 1307 internal nodes ]



#A. Alpha diversity analysis

#### Rarefication - determining cut off value
suppressMessages(library(MicrobiotaProcess)); packageVersion("MicrobiotaProcess") #v1.8.2
#rarefaction cure
rareres.CI <- get_rarecurve(obj=physeq.CI, chunks=400)
#There were 50 or more warnings (use warnings() to see the first 50)
prare1 <- ggrarecurve(obj=rareres.CI,
                          factorNames="Fertilization",
                          shadow=FALSE,
                          indexNames=c("Observe")) + 
  colour.point.fer +
  labs(x="Number of reads", y="Observed species") +
  theme(axis.text=element_text(size=8), 
        panel.grid=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
prare1


#Rarefy data
min.lib<-min(sample_sums(physeq.CI))
min.lib
#[1] 634
##too small to represent dataset

ps0.rar<-rarefy_even_depth(physeq.CI, 2000, replace=TRUE, rngseed = 9242)
#`set.seed(9242)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(9242); .Random.seed` for the full vector
#...
#2 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
#  
#Chidham.UnT.RP.2.22T.dicoccoides.T.RP.3.8
#...

rareres.CI.2000 <- get_rarecurve(obj=ps0.rar, chunks=400)

prare2 <- ggrarecurve(obj=rareres.CI.2000,
                        factorNames="Fertilization",
                        shadow=FALSE,
                        indexNames=c("Observe")) + 
  colour.point.fer +
  labs(x="Number of reads", y="Observed species") +
  theme(axis.text=element_text(size=8), 
        panel.grid=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
prare2
prare2 <- prare2 + ylim(0, 800)
library(ggpubr); packageVersion("ggpubr") #v0.6.0
rrare.plot <- ggarrange(prare1, prare2,common.legend = TRUE, legend = "right", labels = c("A)", "B)"), font.label = list(size=12, face="plain"))
rrare.plot

# Save as pdf
pdf(file = "Plots/rareplots1.pdf",   #file name
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches
rrare.plot
dev.off()
#RStudioGD 
#2 

rm(rareres.CI, rareres.CI.2000, prare1, prare2, rrare.plot) #save memory

detach("package:MicrobiotaProcess", unload = TRUE)
devtools::reload(pkgload::inst("phyloseq"), quiet = TRUE)

#### Alpha diversity metrics
library(microbiome); packageVersion("microbiome") #v1.19.1
hmp.div <- alpha(ps0.rar, index = c("all"))
hmp.meta <- meta(ps0.rar)
hmp.meta$sam_name <- rownames(hmp.meta)
hmp.div$sam_name <- rownames(hmp.div)
div.df <- merge(hmp.div, hmp.meta, by = "sam_name")
#Phylogenetic diversity
ps0.rar.asvtab <- as.data.frame(ps0.rar@otu_table)
ps0.rar.tree <- ps0.rar@phy_tree
ps0.rar@phy_tree
#Phylogenetic tree with 1320 tips and 1307 internal nodes.
#
#Tip labels:
#  c114f024375243308364cc8f91faeabf, f7df07512b80d985b360fc621a0893ea, c66c38d5daf4bd0a030031898490cc61, e3bdbe6e720b5021e9ae6358b9eae180, 6d76d0da1b52f9ee5e64554dab012c3e, 19069dcdac5f474ae061bade1059fbbc, ...
#Node labels:
#  0.868, 0.466, 0.752, 0.829, 0.560, 0.849, ...
#
#Rooted; includes branch lengths.
library(picante); packageVersion("picante") #v1.8.2
df.pd <- pd(t(ps0.rar.asvtab), ps0.rar.tree,include.root=T)
hmp.meta$Phylogenetic_Diversity <- df.pd$PD
div.df$Phylogenetic_Diversity <- df.pd$PD


#Organise data for plotting
div.df.1 <- div.df[, c("Soil", "Fertilization", "Block", "observed.x", "diversity_shannon", "evenness_simpson", "Phylogenetic_Diversity")]
writexl::write_xlsx(div.df.1, "Figure1/Data/alpha_diversity.xlsx")

colnames(div.df.1) <- c("Soil","Fertilization", "Block", "Observed", "Shannon index", "Simpson index", "Faith's PD")

div.df.melt <- reshape2::melt(div.df.1)
#Using Soil, Fertilization, Block as id variables


#Plot
alpha <- ggboxplot(div.df.melt, x = "Soil", y = "value",
                   fill = "Fertilization",
                   legend = "none") + 
  labs(title="",x="Soil", y = "Alpha diversity value") +
  scale_x_discrete(labels=c("UP", "RS", "RP")) +
  colour.fill.fer  +
  facet_wrap(~variable, nrow=1, ncol=4, scales = "free") +
  theme (panel.background = element_rect(colour = "black"),
         legend.position="none",
         axis.text.x=element_text(size=10))
alpha




#B. Beta diversity analysis
ord.norm<-ordinate(physeq.norm.CI, "PCoA", "bray")

pcoa.theme <- theme (
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position="right",
  axis.text.x=element_text(size=10) 
)

pcoa.bray<-plot_ordination(physeq.norm.CI, ord.norm, shape = "Soil" ,color = "Fertilization") +
  geom_point(size=4) +
  pcoa.theme +
  colour.point.fer  +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=Fertilization), show.legend=FALSE)
pcoa.bray #PC1 22.7% PC2 11.4%
#Too few points to calculate an ellipse
pcoa.bray <- pcoa.bray + labs(title="",x="PC1 (22.7%)", y = "PC2 (11.4%)")
pcoa.bray



#### Figures 1A and 1B
figures.1AB <- ggarrange(alpha, pcoa.bray,  ncol=2, nrow=1)
#Too few points to calculate an ellipse
figures.1AB

# Save as pdf
pdf(file = "Figure1/Plots/figures1AB.pdf",   #file name
    width = 12, # The width of the plot in inches
    height = 4) # The height of the plot in inches
figures.1AB
dev.off()
#RStudioGD 
#2 


#### Permutational ANOVA
relative.data.bray <- phyloseq::distance(physeq.norm.CI, method="bray")
metadata <- data.frame(sample_data(physeq.norm.CI))
perm <- how(nperm = 9999)
setBlocks(perm) <- with(metadata, Block)
#non-nested multifactorial
permanova<-adonis2(relative.data.bray ~ Fertilization*Soil*Ploidy*Ancestral_class*Genome*Plant_species, data = metadata, permutations = perm)
permanova
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative_data_bray ~ Fertilization * Soil * Ploidy * Ancestral_class * Genome * Plant_species, data = metadata, permutations = perm)
#                                   Df SumOfSqs      R2       F Pr(>F)    
#Fertilization                        1   12.099 0.17079 64.3384 0.0001 ***
#Soil                                 2    8.299 0.11715 22.0659 0.0001 ***
#Ploidy                               2    1.367 0.01930  3.6344 0.0001 ***
#Ancestral_class                      2    1.246 0.01758  3.3118 0.0001 ***
#Genome                               2    0.654 0.00923  1.7381 0.0149 *  
#Plant_species                        4    1.561 0.02204  2.0755 0.0005 ***
#Fertilization:Soil                   2    1.393 0.01966  3.7032 0.0001 ***
#Fertilization:Ploidy                 2    0.877 0.01237  2.3309 0.0020 ** 
#Soil:Ploidy                          2    0.472 0.00666  1.2545 0.1542    
#Fertilization:Ancestral_class        2    0.673 0.00950  1.7892 0.0165 *  
#Soil:Ancestral_class                 2    0.409 0.00577  1.0877 0.3065    
#Fertilization:Genome                 2    0.483 0.00681  1.2831 0.1322    
#Soil:Genome                          2    0.299 0.00422  0.7942 0.7635    
#Fertilization:Plant_species          4    1.198 0.01691  1.5930 0.0072 ** 
#Soil:Plant_species                   4    0.816 0.01151  1.0843 0.2818    
#Fertilization:Soil:Ploidy            2    0.406 0.00573  1.0795 0.3065    
#Fertilization:Soil:Ancestral_class   2    0.308 0.00434  0.8183 0.7289    
#Fertilization:Soil:Genome            2    0.350 0.00494  0.9297 0.5184    
#Fertilization:Soil:Plant_species     4    0.701 0.00990  0.9319 0.5894    
#Residual                           198   37.234 0.52559                   
#Total                              243   70.841 1.00000                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#nested multifactorial
permanova.nested<-adonis2(relative.data.bray ~ Fertilization*Soil*(Ploidy/Ancestral_class/Genome/Plant_species), data = metadata, permutations = perm)
permanova.nested
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative_data_bray ~ Fertilization * Soil * (Ploidy/Ancestral_class/Genome/Plant_species), data = metadata, permutations = perm)
#                                                                Df SumOfSqs      R2       F Pr(>F)    
#Fertilization                                                    1   12.099 0.17079 64.3384 0.0001 ***
#Soil                                                             2    8.299 0.11715 22.0659 0.0001 ***
#Ploidy                                                           2    1.367 0.01930  3.6344 0.0001 ***
#Fertilization:Soil                                               2    1.389 0.01960  3.6924 0.0001 ***
#Ploidy:Ancestral_class                                           3    1.681 0.02373  2.9804 0.0001 ***
#Fertilization:Ploidy                                             2    0.878 0.01239  2.3337 0.0012 ** 
#Soil:Ploidy                                                      2    0.477 0.00673  1.2676 0.1422    
#Ploidy:Ancestral_class:Genome                                    2    0.605 0.00854  1.6080 0.0315 *  
#Fertilization:Ploidy:Ancestral_class                             3    0.960 0.01355  1.7014 0.0084 ** 
#Soil:Ploidy:Ancestral_class                                      3    0.661 0.00933  1.1720 0.1970    
#Fertilization:Soil:Ploidy                                        2    0.398 0.00562  1.0578 0.3253    
#Ploidy:Ancestral_class:Genome:Plant_species                      3    1.161 0.01639  2.0582 0.0004 ***
#Fertilization:Ploidy:Ancestral_class:Genome                      2    0.457 0.00644  1.2139 0.1806    
#Soil:Ploidy:Ancestral_class:Genome                               2    0.300 0.00424  0.7985 0.7536    
#Fertilization:Soil:Ploidy:Ancestral_class                        3    0.528 0.00745  0.9357 0.5591    
#Fertilization:Ploidy:Ancestral_class:Genome:Plant_species        3    0.951 0.01342  1.6853 0.0086 ** 
#Soil:Ploidy:Ancestral_class:Genome:Plant_species                 3    0.573 0.00809  1.0154 0.4127    
#Fertilization:Soil:Ploidy:Ancestral_class:Genome                 2    0.337 0.00476  0.8966 0.5779    
#Fertilization:Soil:Ploidy:Ancestral_class:Genome:Plant_species   3    0.488 0.00689  0.8655 0.7046    
#Residual                                                       198   37.234 0.52559                   
#Total                                                          243   70.841 1.00000                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Checking dispersion
dispersion.f<-betadisper(relative_data_bray, group=metadata$Fertilization)
dispersion.s<-betadisper(relative_data_bray, group=metadata$Soil)
dispersion.p<-betadisper(relative_data_bray, group=metadata$Ploidy)
permutest(dispersion.f)
permutest(dispersion.s)
permutest(dispersion.p)
plot(dispersion.f, hull=FALSE, ellipse=TRUE)
plot(dispersion.s, hull=FALSE, ellipse=TRUE)
plot(dispersion.p, hull=FALSE, ellipse=TRUE)



# C. Differential abundance analysis

#### Subsetting data by Soil
physeq.da.up<-subset_samples(physeq.CI, Soil == "Unplanted")
physeq.da.up<-filter_taxa(physeq.da.up, function(x) sum(x) > 0, TRUE)

physeq.da.rs<-subset_samples(physeq.CI, Soil==c("Rhizosphere"))
physeq.da.rs<-filter_taxa(physeq.da.rs, function(x) sum(x) > 0, TRUE)

physeq.da.rp<-subset_samples(physeq.CI, Soil=="Rhizoplane")
physeq.da.rp<-filter_taxa(physeq.da.rp, function(x) sum(x) > 0, TRUE)


#transform data to account for zeros
physeq.da.up<-transform_sample_counts(physeq.da.up, function(x) (x+1))
physeq.da.rs<-transform_sample_counts(physeq.da.rs, function(x) (x+1))
physeq.da.rp<-transform_sample_counts(physeq.da.rp, function(x) (x+1))



#### Deseq2
library(DESeq2); packageVersion("DESeq2") #v1.36.0

#unplanted
diagdds.bs = phyloseq_to_deseq2(physeq.da.up, ~ Fertilization)
diagdds.bs = DESeq(diagdds.bs, test="Wald", fitType="local")

res.bs = results(diagdds.bs, cooksCutoff = FALSE)
alpha = 1.00
sigtab.bs = res.bs[which(res.bs$padj < alpha), ]
sigtab.bs = cbind(as(sigtab.bs, "data.frame"), as(tax_table(physeq.da.up)[rownames(sigtab.bs), ], "matrix"))


#rhizosphere
diagdds.rs = phyloseq_to_deseq2(physeq.da.rs, ~ Fertilization)
diagdds.rs = DESeq(diagdds.rs, test="Wald", fitType="local")


res.rs = results(diagdds.rs, cooksCutoff = FALSE)
sigtab.rs = res.rs[which(res.rs$padj < alpha), ]
sigtab.rs = cbind(as(sigtab.rs, "data.frame"), as(tax_table(physeq.da.rs)[rownames(sigtab.rs), ], "matrix"))


#rhizoplane
diagdds.rp = phyloseq_to_deseq2(physeq.da.rp, ~ Fertilization)
diagdds.rp = DESeq(diagdds.rp, test="Wald", fitType="local")


res.rp = results(diagdds.rp, cooksCutoff = FALSE)
sigtab.rp = res.rp[which(res.rp$padj < alpha), ]
sigtab.rp = cbind(as(sigtab.rp, "data.frame"), as(tax_table(physeq.da.rp)[rownames(sigtab.rp), ], "matrix"))



#Organise data for plotting
sigtab.bs$ASV <- "Not significant"
sigtab.bs$ASV[sigtab.bs$log2FoldChange > 0 & sigtab.bs$padj < 0.0001] <- "Fertilized"
sigtab.bs$ASV[sigtab.bs$log2FoldChange < 0 & sigtab.bs$padj < 0.0001] <- "Non-fertilized"
sigtab.bs$ASV = factor(sigtab.bs$ASV, levels = c("Non-fertilized", "Not significant" ,"Fertilized"))

sigtab.rs$ASV <- "Not significant"
sigtab.rs$ASV[sigtab.rs$log2FoldChange > 0 & sigtab.rs$padj < 0.0001] <- "Fertilized"
sigtab.rs$ASV[sigtab.rs$log2FoldChange < 0 & sigtab.rs$pvalue < 0.0001] <- "Non-fertilized"
sigtab.rs$ASV = factor(sigtab.rs$ASV, levels = c("Non-fertilized", "Not significant" ,"Fertilized"))

sigtab.rp$ASV <- "Not significant"
sigtab.rp$ASV[sigtab.rp$log2FoldChange > 0 & sigtab.rp$pvalue < 0.0001] <- "Fertilized"
sigtab.rp$ASV[sigtab.rp$log2FoldChange < 0 & sigtab.rp$pvalue < 0.0001] <- "Non-fertilized"
sigtab.rp$ASV = factor(sigtab.rp$ASV, levels = c("Non-fertilized", "Not significant" ,"Fertilized"))


#Relative abundance of ASVs
sigtab.bs$Abundance <- with(sigtab.bs,ave(baseMean, ASVs, FUN = function(x) 
  (sum(x)/sum(baseMean))*100))
sigtab.rs$Abundance <- with(sigtab.rs,ave(baseMean, ASVs, FUN = function(x) 
  (sum(x)/sum(baseMean))*100))
sigtab.rp$Abundance <- with(sigtab.rp,ave(baseMean, ASVs, FUN = function(x) 
  (sum(x)/sum(baseMean))*100))



#Export data
writexl::write_xlsx(sigtab.bs,"Figure1/Data/sigtab_bs.xlsx")
writexl::write_xlsx(sigtab.rs,"Figure1/Data/sigtab_rs.xlsx")
writexl::write_xlsx(sigtab.rp,"Figure1/Data/sigtab_rp.xlsx")



#### Volcano plots
plot.theme <- theme (
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position="right",
  axis.text.x=element_text(size=10) 
) 

y.axis.title <- my_y_title <- expression(paste("-log10(", italic("P"), " value)"))


volcano.bs <- ggplot(data=sigtab.bs,
                     aes(x=log2FoldChange, y=-log10(pvalue), col=ASV, size = Abundance)) +
  geom_point() + 
  plot.theme +
  scale_color_manual(values=c("black", navy)) +
  labs(title="Unplanted", x="", y=y.axis.title) +
  ylim(0, 100) +
  xlim(-10, 10) +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none") +
  scale_size(range=c(0.5,6),breaks=c(0.5, 1.5 ,3.0, 6.0) ,labels=c("0.5","1.5","3.0","6.0"),guide="legend")
volcano.bs


volcano.rs <- ggplot(data=sigtab.rs,
                     aes(x=log2FoldChange, y=-log10(pvalue), col=ASV, size = Abundance)) +
  geom_point() + 
  plot.theme +
  scale_color_manual(values=c(red, "black", navy)) +
  labs(title="Rhizosphere", x="log2FoldChange", y="") +
  ylim(0, 100) +
  xlim(-10, 10) +
  theme(axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none") +
  scale_size(range=c(0.5,6),breaks=c(0.5, 1.5 ,3.0, 6.0) ,labels=c("0.5","1.5","3.0","6.0"),guide="legend")
volcano.rs


volcano.rp <- ggplot(data=sigtab.rp,
                     aes(x=log2FoldChange, y=-log10(pvalue), col=ASV, size = Abundance)) +
  geom_point() + 
  plot.theme +
  scale_color_manual(values=c(red, "black", navy)) +
  labs(title="Rhizoplane", x="", y="") +
  ylim(0, 100) +
  xlim(-10, 11) +
  theme(axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12)) +
  scale_size(range=c(0.5,6),breaks=c(0.5, 1.5 ,3.0, 6.0) ,labels=c("0.5","1.5","3.0","6.0"),guide="legend") +
  labs(colour="ASV", size="Abundance (%)")
volcano.rp


#Figure1
p1 <- ggarrange(volcano.bs, volcano.rs, nrow = 1)
p1

figure.1C <- ggarrange(p1, volcano.rp, nrow = 1, ncol = 2, widths = c(2, 1.5))
figure.1C


# Save as pdf
pdf(file = "Figure1/Plots/figures1C.pdf",   #file name
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches
figure.1C
dev.off()
#RStudioGD 
#2 


#### Continued in Prism - Figure1.prism

################### END ###################
