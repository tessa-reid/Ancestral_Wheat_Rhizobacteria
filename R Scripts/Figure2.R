############################################################
#R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################
#R 4.2.2
#setwd("Figure2")



#Pre-defined parameters for producing  plots:
red<-"#E64B35FF"
navy<-"#3C5488FF"

#Colour palettes
fertiliser.palette<-c(red, navy)


#Manually change colour scales on ggplot
library(ggplot2); packageVersion("ggplot2") #v3.4.2
colour.fill.fer<-scale_fill_manual(values = fertiliser.palette)
colour.point.fer<-scale_color_manual(values = fertiliser.palette)



#Set theme
theme_set(theme_bw())



#Input files: Files produced by phyloseq (phyloseq object, phyloseq object normalised by deseq)


#Directories:
dir.create("Figure2")
dir.create("Figure2/Data")
dir.create("Figure2/Plots")

suppressMessages(library(phyloseq)); packageVersion("phyloseq") #v1.42.0


#####Importing the files into R:
physeq.CI <- readRDS("physeq_CI.rds")
physeq.norm.CI <- readRDS("physeq_norm_CI.rds")



##### Subsetting the data by Soil

physeq.rs<-subset_samples(physeq.CI, Soil=="Rhizosphere")
physeq.rs<-filter_taxa(physeq.rs, function(x) sum(x) > 0, TRUE)
physeq.rs
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1264 taxa and 127 samples ]
#sample_data() Sample Data:       [ 127 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1264 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1264 tips and 1256 internal nodes ]


physeq.rp<-subset_samples(physeq.CI, Soil=="Rhizoplane")
physeq.rp<-filter_taxa(physeq.rp, function(x) sum(x) > 0, TRUE)
physeq.rp
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1318 taxa and 110 samples ]
#sample_data() Sample Data:       [ 110 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1318 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1318 tips and 1305 internal nodes ]




#A. Alpha diversity analysis

#### Rarefication - determining cut off value
suppressMessages(library(MicrobiotaProcess)); packageVersion("MicrobiotaProcess") #v1.8.2
#rarefaction curve
rareres.rs <- get_rarecurve(obj=physeq.rs, chunks=400)
#There were 50 or more warnings (use warnings() to see the first 50)
rareres.rp <- get_rarecurve(obj=physeq.rp, chunks=400)
#There were 50 or more warnings (use warnings() to see the first 50)

#plots
rare.theme <- theme(axis.text=element_text(size=10), 
                    panel.grid=element_blank(),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    plot.title = element_text(hjust = 0.5))


prare1 <- ggrarecurve(obj=rareres.rs,
                      factorNames="Fertilization",
                      shadow=FALSE,
                      indexNames=c("Observe")) + 
  colour.point.fer +
  labs(title= "Rhizosphere", x="", y="") +
  rare.theme
prare1


prare2 <- ggrarecurve(obj=rareres.rp,
                      factorNames="Fertilization",
                      shadow=FALSE,
                      indexNames=c("Observe")) + 
  colour.point.fer +
  labs(title="Rhizoplane", x="", y="") +
  rare.theme
prare2


min.lib.rs<-min(sample_sums(physeq.rs)) #2125
min.lib.rp<-min(sample_sums(physeq.rp)) #634

#Rarefication
ps0.rar.rs<-rarefy_even_depth(physeq.rs, min.lib.rs, replace=TRUE, rngseed = 9242)
#`set.seed(9242)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(9242); .Random.seed` for the full vector
#...
#4OTUs were removed because they are no longer 
#present in any sample after random subsampling
#
#...

ps0.rar.rp<-rarefy_even_depth(physeq.rp, 2000, replace=TRUE, rngseed = 9242)
#`set.seed(9242)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(9242); .Random.seed` for the full vector
#...
#2 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
#  
#Chidham.UnT.RP.2.22T.dicoccoides.T.RP.3.8
#...
#7OTUs were removed because they are no longer 
#present in any sample after random subsampling
#
#...

rareres.rs.2125 <- get_rarecurve(obj=ps0.rar.rs, chunks=400)
#There were 50 or more warnings (use warnings() to see the first 50)
rareres.rp.2000 <- get_rarecurve(obj=ps0.rar.rp, chunks=400)


#plots
prare3 <- ggrarecurve(obj=rareres.rs.2125,
                      factorNames="Fertilization",
                      shadow=FALSE,
                      indexNames=c("Observe")) + 
  colour.point.fer +
  labs(title="" ,x="", y="") +
  rare.theme
prare3


prare4 <- ggrarecurve(obj=rareres.rp.2000,
                      factorNames="Fertilization",
                      shadow=FALSE,
                      indexNames=c("Observe")) + 
  colour.point.fer +
  labs(title="", x="", y="") +
  rare.theme
prare4


library(ggpubr) #v0.6.0
rrare.plot <- ggarrange(prare1 + ylim(0, 800),
                        prare2, 
                        prare3 + ylim(0, 800),
                        prare4+ ylim(0, 800),
                        common.legend = TRUE,
                        legend = "right",
                        labels = c("A)", "B)", "C)", "D)"),
                        font.label = list(size=12, face="plain"))
rrare.plot

rrare.plot2 <- annotate_figure(rrare.plot,
                               bottom = "Number of reads",
                               left = "Observed species")
rrare.plot2

# Save as pdf
pdf(file = "Figure2/Plots/rareplots2.pdf",   #file name
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
rrare.plot2
dev.off()
#RStudioGD 
#2 


detach("package:MicrobiotaProcess", unload = TRUE)
devtools::reload(pkgload::inst("phyloseq"), quiet = TRUE)


#### Alpha diversity metrics
hmp.div.rs <- microbiome::alpha(ps0.rar.rs, index = c("all"))
hmp.div.rp <- microbiome::alpha(ps0.rar.rp, index = c("all"))
library(microbiome) #v1.19.1
hmp.meta.rs <- meta(ps0.rar.rs)
hmp.meta.rp <- meta(ps0.rar.rp)
hmp.meta.rs$sam_name <- rownames(hmp.meta.rs)
hmp.meta.rp$sam_name <- rownames(hmp.meta.rp)
hmp.div.rs$sam_name <- rownames(hmp.div.rs)
hmp.div.rp$sam_name <- rownames(hmp.div.rp)
div.df.rs <- merge(hmp.div.rs, hmp.meta.rs, by = "sam_name")
div.df.rp <- merge(hmp.div.rp, hmp.meta.rp, by = "sam_name")
#Phylogenetic diversity
ps0.rar.rs.asvtab <- as.data.frame(ps0.rar.rs@otu_table)
ps0.rar.rp.asvtab <- as.data.frame(ps0.rar.rp@otu_table)
ps0.rar.rs.tree <- ps0.rar.rs@phy_tree
ps0.rar.rp.tree <- ps0.rar.rp@phy_tree
ps0.rar.rs@phy_tree
#
#Phylogenetic tree with 1260 tips and 1252 internal nodes.
#
#Tip labels:
#  c114f024375243308364cc8f91faeabf, f7df07512b80d985b360fc621a0893ea, c66c38d5daf4bd0a030031898490cc61, e3bdbe6e720b5021e9ae6358b9eae180, 6d76d0da1b52f9ee5e64554dab012c3e, 19069dcdac5f474ae061bade1059fbbc, ...
#Node labels:
#  0.868, 0.466, 0.752, 0.829, 0.560, 0.849, ...
#
#Rooted; includes branch lengths
ps0.rar.rp@phy_tree
#
#Phylogenetic tree with 1311 tips and 1298 internal nodes.
#
#Tip labels:
#  c114f024375243308364cc8f91faeabf, f7df07512b80d985b360fc621a0893ea, c66c38d5daf4bd0a030031898490cc61, e3bdbe6e720b5021e9ae6358b9eae180, 6d76d0da1b52f9ee5e64554dab012c3e, 19069dcdac5f474ae061bade1059fbbc, ...
#Node labels:
#  0.868, 0.466, 0.752, 0.829, 0.560, 0.849, ...
#
#Rooted; includes branch lengths.
library(picante) #v1.8.2
df.pd.rs <- pd(t(ps0.rar.rs.asvtab), ps0.rar.rs.tree,include.root=T)
df.pd.rp <- pd(t(ps0.rar.rp.asvtab), ps0.rar.rp.tree,include.root=T)
hmp.meta.rs$Phylogenetic_Diversity <- df.pd.rs$PD
hmp.meta.rp$Phylogenetic_Diversity <- df.pd.rp$PD
div.df.rs$Phylogenetic_Diversity <- df.pd.rs$PD
div.df.rp$Phylogenetic_Diversity <- df.pd.rp$PD



#### Organise data for plotting
div.df.rs.1 <- div.df.rs[, c("Soil", "Fertilization", "Ploidy", "Block", "observed.x", "diversity_shannon", "evenness_simpson", "Phylogenetic_Diversity")]
writexl::write_xlsx(div.df.rs.1, "Figure2/Data/alpha_diversity_rhizosphere.xlsx")

div.df.rp.1 <- div.df.rp[, c("Soil", "Fertilization", "Ploidy", "Block", "observed.x", "diversity_shannon", "evenness_simpson", "Phylogenetic_Diversity")]
writexl::write_xlsx(div.df.rp.1, "Figure2/Data/alpha_diversity_rhizoplane.xlsx")


colnames(div.df.rs.1) <- c("Soil", "Fertilization", "Ploidy", "Block", "Observed", "Shannon index", "Simpson index", "Faith's PD")
colnames(div.df.rp.1) <- c("Soil", "Fertilization", "Ploidy","Block", "Observed", "Shannon index", "Simpson index", "Faith's PD")


div.df.melt.rs <- reshape2::melt(div.df.rs.1)
div.df.melt.rp <- reshape2::melt(div.df.rp.1)



#Plots
alpha.rs <- ggboxplot(div.df.melt.rs, x = "Ploidy", y = "value",
                   fill = "Fertilization",
                   legend = "none") + 
  labs(title="Rhizosphere",x="Ploidy", y = "Alpha diversity value")+
  colour.fill.fer  +
  facet_wrap(~variable, nrow=1, ncol=4, scales = "free") +
  theme(panel.background = element_rect(colour = "black"),
         legend.position="none",
         axis.text.x=element_text(size=10, angle = 90, vjust = 0.5, hjust=1))
alpha.rs


alpha.rp <- ggboxplot(div.df.melt.rp, x = "Ploidy", y = "value",
                      fill = "Fertilization",
                      legend = "none") + 
  labs(title="Rhizoplane",x="Ploidy", y = "")+
  colour.fill.fer  +
  facet_wrap(~variable, nrow=1, ncol=4, scales = "free") +
  theme (panel.background = element_rect(colour = "black"),
         legend.position="right",
         axis.text.x=element_text(size=10, angle = 90, vjust = 0.5, hjust=1))
alpha.rp

figure.2A <- ggarrange(alpha.rs, alpha.rp, nrow = 1, ncol = 2, widths = c(1, 1.25))
figure.2A

# Save as pdf
pdf(file = "Figure2/Plots/alpha_diversity_plots_soil.pdf",   #file name
    width = 15, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches
figure.2A
dev.off()
#RStudioGD 
#2 


#Factors influencing beta diversity in rhizosphere and rhizoplane samples

#### Subsetting the data, abundances normalised by deseq VST, by Soil
physeq.norm.rs<-subset_samples(physeq.norm.CI, Soil=="Rhizosphere")
physeq.norm.rs<-filter_taxa(physeq.norm.rs, function(x) sum(x) > 0, TRUE)
physeq.norm.rs
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1264 taxa and 127 samples ]
#sample_data() Sample Data:       [ 127 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1264 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1264 tips and 1256 internal nodes ]


physeq.norm.rp<-subset_samples(physeq.norm.CI, Soil=="Rhizoplane")
physeq.norm.rp<-filter_taxa(physeq.norm.rp, function(x) sum(x) > 0, TRUE)
physeq.norm.rp
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1318 taxa and 110 samples ]
#sample_data() Sample Data:       [ 110 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1318 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1318 tips and 1305 internal nodes ]


#### Unconstrained ordination - Principal coordinates of analysis (PCoA)

ord.norm.rs<-ordinate(physeq.norm.rs, "PCoA", "bray")
ord.norm.rp<-ordinate(physeq.norm.rp, "PCoA", "bray")

pcoa.theme <- theme (
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position="right",
  axis.text.x=element_text(size=10) 
)

pcoa.bray.rs<-plot_ordination(physeq.norm.rs, ord.norm.rs, shape = "Ploidy" ,color = "Fertilization") +
  geom_point(size=4) +
  pcoa.theme +
  colour.point.fer  +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=Fertilization), show.legend=FALSE)
pcoa.bray.rs #PC1 22.7% PC2 11.4%

pcoa.bray.rp<-plot_ordination(physeq.norm.rp, ord.norm.rp, shape = "Ploidy" ,color = "Fertilization") +
  geom_point(size=4) +
  pcoa.theme +
  colour.point.fer  +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=Fertilization), show.legend=FALSE)
pcoa.bray.rp #PC1 32% PC2 7.6%


pcoa.bray.rs <- pcoa.bray.rs + labs(title="Rhizosphere",x="PC1 (22.7%)", y = "PC2 (11.4%)")
pcoa.bray.rp <- pcoa.bray.rp + labs(title="Rhizoplane",x="PC1 (32%)", y = "PC2 (7.6%)")

pcoa.bray.soil <- ggarrange(pcoa.bray.rs, pcoa.bray.rp, common.legend = TRUE, legend = "right", labels = c("A", "B"))
pcoa.bray.soil

# Save as pdf
pdf(file = "Figure2/Plots/pcoa_bray_soil.pdf",   #file name
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches
pcoa.bray.soil
dev.off()
#RStudioGD 
#2 


#### Permutational ANOVA

##Rhizosphere
relative.data.bray.rs <- phyloseq::distance(physeq.norm.rs, method="bray")
metadata.rs <- data.frame(sample_data(physeq.norm.rs))
perm <- how(nperm = 9999)
setBlocks(perm) <- with(metadata.rs, Block)

#non-nested multifactorial
permanova.rs<-adonis2(relative.data.bray.rs ~ Fertilization*Ploidy*Ancestral_class*Genome*Plant_species, data = metadata.rs, permutations = perm)
permanova.rs
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray ~ Fertilization * Ploidy * Ancestral_class * Genome * Plant_species, data = metadata, permutations = perm)
#                               Df SumOfSqs      R2       F Pr(>F)    
#Fertilization                   1    6.209 0.19295 32.3574 0.0001 ***
#Ploidy                          2    0.953 0.02962  2.4840 0.0008 ***
#Ancestral_class                 2    0.910 0.02828  2.3715 0.0018 ** 
#Genome                          2    0.445 0.01382  1.1591 0.2419    
#Plant_species                   4    1.016 0.03156  1.3230 0.0682 .  
#Fertilization:Ploidy            2    0.698 0.02168  1.8183 0.0124 *  
#Fertilization:Ancestral_class   2    0.506 0.01573  1.3190 0.1138    
#Fertilization:Genome            2    0.379 0.01179  0.9886 0.4693    
#Fertilization:Plant_species     4    0.916 0.02845  1.1928 0.1525    
#Residual                      105   20.150 0.62611                   
#Total                         126   32.182 1.00000                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#nested multifactorial
permanova.nested.rs<-adonis2(relative.data.bray.rs ~ Fertilization*(Ploidy/Ancestral_class/Genome/Plant_species), data = metadata.rs, permutations = perm)
permanova.nested.rs
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray ~ Fertilization * (Ploidy/Ancestral_class/Genome/Plant_species), data = metadata, permutations = perm)
#                                                           Df SumOfSqs      R2       F Pr(>F)    
#Fertilization                                               1    6.209 0.19295 32.3574 0.0001 ***
#Ploidy                                                      2    0.953 0.02962  2.4840 0.0012 ** 
#Ploidy:Ancestral_class                                      3    1.215 0.03774  2.1099 0.0005 ***
#Fertilization:Ploidy                                        2    0.695 0.02158  1.8098 0.0131 *  
#Ploidy:Ancestral_class:Genome                               2    0.420 0.01306  1.0952 0.3079    
#Fertilization:Ploidy:Ancestral_class                        3    0.734 0.02279  1.2742 0.1064    
#Ploidy:Ancestral_class:Genome:Plant_species                 3    0.737 0.02290  1.2800 0.1082    
#Fertilization:Ploidy:Ancestral_class:Genome                 2    0.367 0.01140  0.9561 0.5406    
#Fertilization:Ploidy:Ancestral_class:Genome:Plant_species   3    0.703 0.02184  1.2207 0.1486    
#Residual                                                  105   20.150 0.62611                   
#Total                                                     126   32.182 1.00000                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Checking dispersion
dispersion.rs<-betadisper(relative.data.bray.rs, group=metadata.rp$Fertilization)
permutest(dispersion.rs)
plot(dispersion.rs, hull=FALSE, ellipse=TRUE)



##Rhizoplane
relative.data.bray.rp <- phyloseq::distance(physeq.norm.rp, method="bray")
metadata.rp <- data.frame(sample_data(physeq.norm.rp))
perm <- how(nperm = 9999)
setBlocks(perm) <- with(metadata.rp, Block)

#non-nested multifactorial
permanova.rp<-adonis2(relative.data.bray.rp ~ Fertilization*Ploidy*Ancestral_class*Genome*Plant_species, data = metadata.rp, permutations = perm)
#permanova.rp
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rp, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rp ~ Fertilization * Ploidy * Ancestral_class * Genome * Plant_species, data = metadata.rp, permutations = perm)
#                               Df SumOfSqs      R2       F Pr(>F)    
#Fertilization                   1   7.1903 0.24106 38.0322 0.0001 ***
#Ploidy                          2   0.8887 0.02980  2.3505 0.0053 ** 
#Ancestral_class                 2   0.7543 0.02529  1.9948 0.0168 *  
#Genome                          2   0.5052 0.01694  1.3360 0.1120    
#Plant_species                   4   1.3827 0.04636  1.8284 0.0054 ** 
#Fertilization:Ploidy            2   0.5711 0.01915  1.5104 0.0714 .  
#Fertilization:Ancestral_class   2   0.4724 0.01584  1.2492 0.1738    
#Fertilization:Genome            2   0.4462 0.01496  1.1799 0.1985    
#Fertilization:Plant_species     4   0.9796 0.03284  1.2954 0.0976 .  
#Residual                       88  16.6372 0.55778                   
#Total                         109  29.8277 1.00000                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#nested multifactorial
permanova.nested.rp<-adonis2(relative.data.bray.rp ~ Fertilization*(Ploidy/Ancestral_class/Genome/Plant_species), data = metadata.rp, permutations = perm)
permanova.nested.rp
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray ~ Fertilization * (Ploidy/Ancestral_class/Genome/Plant_species), data = metadata, permutations = perm)
#                                                           Df SumOfSqs      R2       F Pr(>F)    
#Fertilization                                               1    6.209 0.19295 32.3574 0.0001 ***
#Ploidy                                                      2    0.953 0.02962  2.4840 0.0012 ** 
#Ploidy:Ancestral_class                                      3    1.215 0.03774  2.1099 0.0005 ***
#Fertilization:Ploidy                                        2    0.695 0.02158  1.8098 0.0131 *  
#Ploidy:Ancestral_class:Genome                               2    0.420 0.01306  1.0952 0.3079    
#Fertilization:Ploidy:Ancestral_class                        3    0.734 0.02279  1.2742 0.1064    
#Ploidy:Ancestral_class:Genome:Plant_species                 3    0.737 0.02290  1.2800 0.1082    
#Fertilization:Ploidy:Ancestral_class:Genome                 2    0.367 0.01140  0.9561 0.5406    
#Fertilization:Ploidy:Ancestral_class:Genome:Plant_species   3    0.703 0.02184  1.2207 0.1486    
#Residual                                                  105   20.150 0.62611                   
#Total                                                     126   32.182 1.00000                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Checking dispersion
dispersion.rp<-betadisper(relative.data.bray.rp, group=metadata.rp$Fertilization)
permutest(dispersion.rp)
plot(dispersion, hull=FALSE, ellipse=TRUE)


#B. Created in Prism from pseudo-F values - Figure2.pzfx


#C. Constrained ordination: Canonical Analysis of Principal coordinates

##### Convert factors to numerical values

#Rhizosphere
sample_data(physeq.norm.rs)$Fer <- as.numeric( sample_data(physeq.norm.rs)$Fertilization)
sample_data(physeq.norm.rs)$Ploi <- as.numeric( sample_data(physeq.norm.rs)$Ploidy)
sample_data(physeq.norm.rs)$Anc <- as.numeric( sample_data(physeq.norm.rs)$Ancestral_class)
sample_data(physeq.norm.rs)$Spe <- as.numeric( sample_data(physeq.norm.rs)$Plant_species)

#Rhizoplane
sample_data(physeq.norm.rp)$Fer <- as.numeric( sample_data(physeq.norm.rp)$Fertilization)
sample_data(physeq.norm.rp)$Ploi <- as.numeric( sample_data(physeq.norm.rp)$Ploidy)
sample_data(physeq.norm.rp)$Anc <- as.numeric( sample_data(physeq.norm.rp)$Ancestral_class)
sample_data(physeq.norm.rp)$Spe <- as.numeric( sample_data(physeq.norm.rp)$Plant_species)



#### Create ordination
cap.ord.rs = ordinate(physeq.norm.rs, method = "CAP", distance = "bray", formula = physeq.norm.rs ~ Fer * Ploi + Anc + Spe)
anova(cap.ord.rs)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = OTU ~ Fer * Ploi + Anc + Spe, data = data, distance = distance)
#          Df SumOfSqs      F Pr(>F)    
#Model      5   8.0821 7.9843  0.001 ***
#Residual 121  24.4964                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp = ordinate(physeq.norm.rp, method = "CAP", distance = "bray", formula = physeq.norm.rp ~ Fer * Ploi + Anc + Spe)
anova(cap.ord.rp)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = OTU ~ Fer * Ploi + Anc + Spe, data = data, distance = distance)
#          Df SumOfSqs      F Pr(>F)    
#Model      5    8.922 8.7916  0.001 ***
#Residual 104   21.108                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Add the environmental variables as arrows
arrowmat.rs <- vegan::scores(cap.ord.rs, display = "bp")
arrowmat.rp <- vegan::scores(cap.ord.rp, display = "bp")


#Add labels, make a data.frame
arrowdf.rs <- data.frame(labels = rownames(arrowmat.rs), arrowmat.rs)
arrowdf.rp <- data.frame(labels = rownames(arrowmat.rp), arrowmat.rp)



#Define the arrow aesthetic mapping
arrow.map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)


label.map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))



#### Plots
cap.plot.rs <- plot_ordination(
  physeq = physeq.norm.rs, 
  ordination = cap.ord.rs, 
  color = "Fertilization", 
  shape = "Ploidy",
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = fertiliser.palette) +
  theme (
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="right") +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="Rhizosphere bacterial community", x="CAP1 (19.1%)", y="CAP2 (3.1%)")
cap.plot.rs


cap.plot.rp <- plot_ordination(
  physeq = physeq.norm.rp, 
  ordination = cap.ord.rp, 
  color = "Fertilization", 
  shape = "Ploidy",
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = fertiliser.palette) +
  theme (
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="right") +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="Rhizoplane bacterial community", x="CAP1 (24.3%)", y="CAP2 (2.4%)")
cap.plot.rp



#Plots with arrow
cap.plot.rs <- cap.plot.rs + 
  geom_segment(
    mapping = arrow.map, 
    size = .5, 
    data = arrowdf.rp, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label.map, 
    size = 4,  
    data = arrowdf.rp, 
    show.legend = FALSE
  )
cap.plot.rs


cap.plot.rp<-cap.plot.rp + 
  geom_segment(
    mapping = arrow.map, 
    size = .5, 
    data = arrowdf.rp, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label.map, 
    size = 4,  
    data = arrowdf.rp, 
    show.legend = FALSE
  )
cap.plot.rp


#Final figure
figure.2C <- ggarrange(cap.plot.rs, cap.plot.rp, common.legend = TRUE, legend = "right")
figure.2C

# Save as pdf
pdf(file = "Figure2/Plots/cap_plots_soil.pdf",   #file name
    width = 8, # The width of the plot in inches
    height = 3) # The height of the plot in inches
figure.2C
dev.off()
#RStudioGD 
#2 


#### Continued in Prism - Figure2.pzfx

################### END ###################
