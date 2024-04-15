############################################################
#R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################
#R 4.2.2
#setwd()


#Input files: Files produced by phyloseq (phyloseq object normalised by deseq)

#Directories:
dir.create("Figure3")
dir.create("Data")
dir.create("Plots")

library(phyloseq); packageVersion("phyloseq") #v1.42.0



##Beta diversity analysis


#####Importing the files into R:
physeq.norm.CI <- readRDS("physeq_norm_CI.rds")
physeq.norm.CI
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1320 taxa and 244 samples ]
#sample_data() Sample Data:       [ 244 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1320 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1320 tips and 1307 internal nodes ]


#Four subsets:
## amplicon sequence variants from non-fertilized (nf) rhizosphere (rs)
## amplicon sequence variants from fertilized (f) rhizosphere (rs)
## amplicon sequence variants from non-fertilized (nf) rhizoplane (rp)
## amplicon sequence variants from fertilized (f) rhizoplane (rp)


##### Subsetting the data by Soil and Fertilization
physeq.norm.rs.nf<-subset_samples(physeq.norm.CI, Soil=="Rhizosphere" & Fertilization=="Non-fertilized")
physeq.norm.rs.nf<-filter_taxa(physeq.norm.rs.nf, function(x) sum(x) > 0, TRUE)
physeq.norm.rs.nf
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1114 taxa and 64 samples ]
#sample_data() Sample Data:       [ 64 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1114 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1114 tips and 1109 internal nodes ]


physeq.norm.rs.f<-subset_samples(physeq.norm.CI, Soil=="Rhizosphere" & Fertilization=="Fertilized")
physeq.norm.rs.f<-filter_taxa(physeq.norm.rs.f, function(x) sum(x) > 0, TRUE)
physeq.norm.rs.f
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1114 taxa and 63 samples ]
#sample_data() Sample Data:       [ 63 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1114 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1114 tips and 1107 internal nodes ]


physeq.norm.rp.nf<-subset_samples(physeq.norm.CI, Soil=="Rhizoplane" & Fertilization=="Non-fertilized")
physeq.norm.rp.nf<-filter_taxa(physeq.norm.rp.nf, function(x) sum(x) > 0, TRUE)
physeq.norm.rp.nf
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1253 taxa and 55 samples ]
#sample_data() Sample Data:       [ 55 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1253 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1253 tips and 1242 internal nodes ]


physeq.norm.rp.f<-subset_samples(physeq.norm.CI, Soil=="Rhizoplane" & Fertilization=="Fertilized")
physeq.norm.rp.f<-filter_taxa(physeq.norm.rp.f, function(x) sum(x) > 0, TRUE)
physeq.norm.rp.f
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1218 taxa and 55 samples ]
#sample_data() Sample Data:       [ 55 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1218 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1218 tips and 1205 internal nodes ]




####### Unconstrained ordination: Principal Coordinate Analysis (PCoA) Bray-Curtis ##########

#Dataframes
#Non-fertilized rhizosphere
ord.norm.rs.nf<-ordinate(physeq.norm.rs.nf, "PCoA", "bray")
pcoa.bray.rs.nf<-plot_ordination(physeq.norm.rs.nf, ord.norm.rs.nf)
pcoa.bray.rs.nf #PC1 10.3%, PC2 6.9%
pcoa.bray.rs.nf.df <- pcoa.bray.rs.nf$data
library(plyr); packageVersion("plyr") #v1.8.8
pcoa.bray.rs.nf.means <- ddply(pcoa.bray.rs.nf.df,
                               ~Ancestral_class+Ploidy+Plant_initials,
                               summarise,
                               PC1=mean(Axis.1),
                               x.se=sd(Axis.1)/sqrt((length(Axis.1))),
                               PC2=mean(Axis.2), 
                               y.se=sd(Axis.2)/sqrt((length(Axis.2)))
                               )


#Fertilized rhizosphere
ord.norm.rs.f<-ordinate(physeq.norm.rs.f, "PCoA", "bray")
pcoa.bray.rs.f<-plot_ordination(physeq.norm.rs.f, ord.norm.rs.f)
pcoa.bray.rs.f#PC1 21/4%, PC2 9.1%
pcoa.bray.rs.f.df <- pcoa.bray.rs.f$data
pcoa.bray.rs.f.means <- ddply(pcoa.bray.rs.f.df,
                              ~Ancestral_class+Ploidy+Plant_initials,
                              summarise,
                              PC1=mean(Axis.1),
                              x.se=sd(Axis.1)/sqrt((length(Axis.1))),
                              PC2=mean(Axis.2), 
                              y.se=sd(Axis.2)/sqrt((length(Axis.2)))
                              )


#Non-fertilized rhizoplane
ord.norm.rp.nf<-ordinate(physeq.norm.rp.nf, "PCoA", "bray")
pcoa.bray.rp.nf<-plot_ordination(physeq.norm.rp.nf, ord.norm.rp.nf)
pcoa.bray.rp.nf #PC1 16.5%, PC2 9.6%
pcoa.bray.rp.nf.df <- pcoa.bray.rp.nf$data
pcoa.bray.rp.nf.means <- ddply(pcoa.bray.rp.nf.df,
                               ~Ancestral_class+Ploidy+Plant_initials,
                               summarise,
                               PC1=mean(Axis.1),
                               x.se=sd(Axis.1)/sqrt((length(Axis.1))),
                               PC2=mean(Axis.2), 
                               y.se=sd(Axis.2)/sqrt((length(Axis.2)))
                               )


#Fertilized rhizoplane
ord.norm.rp.f<-ordinate(physeq.norm.rp.f, "PCoA", "bray")
pcoa.bray.rp.f<-plot_ordination(physeq.norm.rp.f, ord.norm.rp.f)
pcoa.bray.rp.f #PC1 16.6% PC2 12.4%
pcoa.bray.rp.f.df <- pcoa.bray.rp.f$data
pcoa.bray.rp.f.means <- ddply(pcoa.bray.rp.f.df,
                              ~Ancestral_class+Ploidy+Plant_initials,
                              summarise,
                              PC1=mean(Axis.1),
                              x.se=sd(Axis.1)/sqrt((length(Axis.1))),
                              PC2=mean(Axis.2), 
                              y.se=sd(Axis.2)/sqrt((length(Axis.2)))
                              )



#Pre-defined parameters for producing  plots:
pale.red<-"#F39B7FFF"
pale.blue<-"#8491B4FF"
pale.green<-"#91D1C2FF"


#Colour palettes
ploidy.palette<-c(pale.red, pale.blue, pale.green)


#Manually change colour scales on ggplot
library(ggplot2); packageVersion("ggplot2") #v3.4.2
colour.fill.ploi<-scale_fill_manual(values = ploidy.palette)
colour.point.ploi<-scale_color_manual(values = ploidy.palette)


#Set theme
theme.new <- theme (
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position="right",
  plot.title = element_text(hjust = 0.5)
) 

theme_set(theme_bw())



#Plots
pcoa.rs.nf.plot <- ggplot(pcoa.bray.rs.nf.means, 
                 aes(x=PC1, y=PC2, color=Ploidy, shape=Ancestral_class)) +
  geom_point(size=4) +
  theme.new +
  colour.point.ploi +
  scale_shape_manual(values=c(17,16,15)) +
  labs(title = "Rhizosphere", x="PC1 (10.7%)", y = "PC2 (6.9%)") +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  geom_errorbar(aes(ymin = PC2 - y.se, ymax = PC2 + y.se)) +
  geom_errorbarh(aes(xmax = PC1 + x.se, xmin = PC1 - x.se)) +
  geom_text(aes(label=Plant_initials), size=2, color="black", check_overlap = TRUE)
pcoa.rs.nf.plot


pcoa.rs.f.plot <- ggplot(pcoa.bray.rs.f.means, 
                 aes(x=PC1, y=PC2, color=Ploidy, shape=Ancestral_class)) +
  geom_point(size=4) +
  theme.new +
  colour.point.ploi +
  scale_shape_manual(values=c(17,16,15)) +
  labs(title = "Rhizosphere", x="PC1 (21.4%)", y = "PC2 (9.1%)") +
  geom_hline(yintercept = 0, size=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  geom_errorbar(aes(ymin = PC2 - y.se, ymax = PC2 + y.se)) +
  geom_errorbarh(aes(xmax = PC1 + x.se, xmin = PC1 - x.se)) +
  geom_text(aes(label=Plant_initials), size=2, color="black", check_overlap = TRUE)
pcoa.rs.f.plot


pcoa.rp.nf.plot <- ggplot(pcoa.bray.rp.nf.means, 
                 aes(x=PC1, y=PC2, color=Ploidy, shape=Ancestral_class)) +
  geom_point(size=4) +
  theme.new +
  colour.point.ploi +
  scale_shape_manual(values=c(17,16,15)) +
  labs(title = "Rhizoplane", x="PC1 (16.5%)", y = "PC2 (9.6%)") +
  geom_hline(yintercept = 0, size=0.1) + geom_vline(xintercept = 0, size=0.1) +
  geom_errorbar(aes(ymin = PC2 - y.se, ymax = PC2 + y.se)) +
  geom_errorbarh(aes(xmax = PC1 + x.se, xmin = PC1 - x.se)) +
  geom_text(aes(label=Plant_initials), size=2, color="black", check_overlap = TRUE)
pcoa.rp.nf.plot


pcoa.rp.f.plot <- ggplot(pcoa.bray.rp.f.means, 
                 aes(x=PC1, y=PC2, color=Ploidy, shape=Ancestral_class)) +
          geom_point(size=4) +
          theme.new +
          colour.point.ploi +
          scale_shape_manual(values=c(17,16,15)) +
          labs(title = "Rhizoplane", x="PC1 (16.6%)", y = "PC2 (12.4%)") +
          geom_hline(yintercept = 0, size=0.1) + 
          geom_vline(xintercept = 0, size=0.1) +
          geom_errorbar(aes(ymin = PC2 - y.se, ymax = PC2 + y.se)) +
          geom_errorbarh(aes(xmax = PC1 + x.se, xmin = PC1 - x.se)) +
          geom_text(aes(label=Plant_initials), size=2, color="black", check_overlap = TRUE)
pcoa.rp.f.plot


#Final figure
library(ggpubr); packageVersion("ggpubr") #v0.6.0
pcoa.plots <- ggarrange(pcoa.rs.nf.plot, pcoa.rp.nf.plot, pcoa.rs.f.plot, pcoa.rp.f.plot,
                        legend="right",
                        common.legend = TRUE,
                        nrow = 2,
                        ncol = 2,
                        align = "v",
                        labels = c("A", "", "B", ""))
pcoa.plots

# Save as pdf
pdf(file = "Figure3/Plots/pcoa-ploidy.pdf",   #file name
    width = 8, # The width of the plot in inches
    height = 6) # The height of the plot in inches
pcoa.plots
dev.off()
#RStudioGD 
#2 


#PERMANOVA test for the analysis of variation
#Non-fertilized rhizosphere
relative.data.bray.rs.nf <- phyloseq::distance(physeq.norm.rs.nf, method="bray")
metadata.rs.nf <- data.frame(sample_data(physeq.norm.rs.nf))
library(permute); packageVersion("permute") #v0.9.7
perm <- how(nperm = 9999)
setBlocks(perm) <- with(metadata.rs.nf, Block)


#non-nested multifactorial
library(vegan); packageVersion("vegan") #v2.6.4
permanova.rs.nf<-adonis2(relative.data.bray.rs.nf ~ Ploidy*Ancestral_class*Genome*Plant_species, data = metadata.rs.nf, permutations = perm)
permanova.rs.nf
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rs.nf, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rs.nf ~ Ploidy * Ancestral_class * Genome * Plant_species, data = metadata.rs.nf, permutations = perm)
#                Df SumOfSqs      R2      F Pr(>F)    
#Ploidy           2   0.6217 0.05152 1.7171 0.0002 ***
#Ancestral_class  2   0.5710 0.04731 1.5770 0.0015 ** 
#Genome           2   0.4035 0.03344 1.1144 0.2320    
#Plant_species    4   0.8767 0.07265 1.2106 0.0460 *  
#Residual        53   9.5956 0.79509                  
#Total           63  12.0686 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#nested multifactorial
permanova.rs.nf.nested<-adonis2(relative.data.bray.rs.nf ~ Ploidy/Ancestral_class/Genome/Plant_species, data = metadata.rs.nf, permutations = perm)
permanova.rs.nf.nested
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rs.nf, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rs.nf ~ Ploidy/Ancestral_class/Genome/Plant_species, data = metadata.rs.nf, permutations = perm)
#                                            Df SumOfSqs      R2      F Pr(>F)    
#Ploidy                                       2   0.6217 0.05152 1.7171 0.0006 ***
#Ploidy:Ancestral_class                       3   0.8720 0.07226 1.6055 0.0003 ***
#Ploidy:Ancestral_class:Genome                2   0.3690 0.03057 1.0190 0.4636    
#Ploidy:Ancestral_class:Genome:Plant_species  3   0.6103 0.05057 1.1236 0.1857    
#Residual                                    53   9.5956 0.79509                  
#Total                                       63  12.0686 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dispersion.rs.nf<-betadisper(relative.data.bray.rs.nf, group=metadata.rs.nf$Ploidy)
permutest(dispersion.rs.nf)
#
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.00976 0.0048793 0.8811    999  0.419
#Residuals 61 0.33780 0.0055376 
plot(dispersion.rs.nf, hull=FALSE, ellipse=TRUE)



#Fertilized Rhizosphere
relative.data.bray.rs.f <- phyloseq::distance(physeq.norm.rs.f, method="bray")
metadata.rs.f <- data.frame(sample_data(physeq.norm.rs.f))
perm <- how(nperm = 9999)
setBlocks(perm) <- with(metadata.rs.f, Block)


#non-nested multifactorial
permanova.rs.f<-adonis2(relative.data.bray.rs.f ~ Ploidy*Ancestral_class*Genome*Plant_species, data = metadata.rs.f, permutations = perm)
permanova.rs.f
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rs.f, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rs.f ~ Ploidy * Ancestral_class * Genome * Plant_species, data = metadata.rs.f, permutations = perm)
#                Df SumOfSqs      R2      F Pr(>F)    
#Ploidy           2   1.0272 0.07387 2.5304 0.0004 ***
#Ancestral_class  2   0.8630 0.06207 2.1261 0.0013 ** 
#Genome           2   0.4191 0.03015 1.0326 0.3904    
#Plant_species    4   1.0407 0.07485 1.2819 0.0656 .  
#Residual        52  10.5540 0.75906                  
#Total           62  13.9041 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#nested multifactorial
permanova.rs.f.nested<-adonis2(relative.data.bray.rs.f ~ Ploidy/Ancestral_class/Genome/Plant_species, data = metadata.rs.f, permutations = perm)
permanova.rs.f.nested
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rs.f, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rs.f ~ Ploidy/Ancestral_class/Genome/Plant_species, data = metadata.rs.f, permutations = perm)
#                                            Df SumOfSqs      R2      F Pr(>F)    
#Ploidy                                       2   1.0272 0.07387 2.5304 0.0004 ***
#Ploidy:Ancestral_class                       3   1.0785 0.07756 1.7712 0.0035 ** 
#Ploidy:Ancestral_class:Genome                2   0.4151 0.02985 1.0225 0.4171    
#Ploidy:Ancestral_class:Genome:Plant_species  3   0.8294 0.05965 1.3621 0.0475 *  
#Residual                                    52  10.5540 0.75906                  
#Total                                       62  13.9041 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


dispersion.rs.f<-betadisper(relative.data.bray.rs.f, group=metadata.rs.f$Ploidy)
permutest(dispersion.rs.f)
#
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)   
#Groups     2 0.030832 0.0154161 5.5206    999  0.009 **
#Residuals 60 0.167549 0.0027925                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
plot(dispersion.rs.f, hull=FALSE, ellipse=TRUE)



#Non-fertilized rhizoplane
relative.data.bray.rp.nf <- phyloseq::distance(physeq.norm.rp.nf, method="bray")
metadata.rp.nf <- data.frame(sample_data(physeq.norm.rp.nf))
perm <- how(nperm = 9999)
setBlocks(perm) <- with(metadata.rp.nf, Block)


#non-nested multifactorial
permanova.rp.nf<-adonis2(relative.data.bray.rp.nf ~ Ploidy*Ancestral_class*Genome*Plant_species, data = metadata.rp.nf, permutations = perm)
permanova.rp.nf
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rp.nf, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rp.nf ~ Ploidy * Ancestral_class * Genome * Plant_species, data = metadata.rp.nf, permutations = perm)
#                Df SumOfSqs      R2      F Pr(>F)    
#Ploidy           2   0.6964 0.06725 1.9766 0.0006 ***
#Ancestral_class  2   0.4884 0.04717 1.3864 0.0207 *  
#Genome           2   0.4100 0.03960 1.1639 0.1953    
#Plant_species    4   1.0100 0.09753 1.4334 0.0022 ** 
#Residual        44   7.7506 0.74846                  
#Total           54  10.3554 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#nested multifactorial
permanova.rp.nf.nested<-adonis2(relative.data.bray.rp.nf ~ Ploidy/Ancestral_class/Genome/Plant_species, data = metadata.rp.nf, permutations = perm)
permanova.rp.nf.nested
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rp.nf, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rp.nf ~ Ploidy/Ancestral_class/Genome/Plant_species, data = metadata.rp.nf, permutations = perm)
#                                            Df SumOfSqs      R2      F Pr(>F)    
#Ploidy                                       2   0.6964 0.06725 1.9766 0.0002 ***
#Ploidy:Ancestral_class                       3   0.7301 0.07050 1.3816 0.0203 *  
#Ploidy:Ancestral_class:Genome                2   0.4018 0.03880 1.1404 0.2140    
#Ploidy:Ancestral_class:Genome:Plant_species  3   0.7766 0.07499 1.4696 0.0048 ** 
#Residual                                    44   7.7506 0.74846                  
#Total                                       54  10.3554 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dispersion.rp.nf<-betadisper(relative.data.bray.rp.nf, group=metadata.rp.nf$Ploidy)
permutest(dispersion.rp.nf)
#
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.01470 0.0073486 0.4097    999  0.666
#Residuals 52 0.93263 0.0179353 
plot(dispersion.rp.nf, hull=FALSE, ellipse=TRUE)



#Fertilized Rhizoplane
relative.data.bray.rp.f <- phyloseq::distance(physeq.norm.rp.f, method="bray")
metadata.rp.f <- data.frame(sample_data(physeq.norm.rp.f))
perm <- how(nperm = 9999)
setBlocks(perm) <- with(metadata.rp.f, Block)


#non-nested multifactorial
permanova.rp.f<-adonis2(relative.data.bray.rp.f ~ Ploidy*Ancestral_class*Genome*Plant_species, data = metadata.rp.f, permutations = perm)
permanova.rp.f
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rp.f, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rp.f ~ Ploidy * Ancestral_class * Genome * Plant_species, data = metadata.rp.f, permutations = perm)
#                Df SumOfSqs      R2      F Pr(>F)   
#Ploidy           2   0.7709 0.06277 1.9085 0.0016 **
#Ancestral_class  2   0.7482 0.06092 1.8523 0.0025 **
#Genome           2   0.5431 0.04422 1.3445 0.0614 . 
#Plant_species    4   1.3332 0.10855 1.6502 0.0011 **
#Residual        44   8.8866 0.72355                 
#Total           54  12.2819 1.00000                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#nested multifactorial
permanova.rp.f.nested<-adonis2(relative.data.bray.rp.f ~ Ploidy/Ancestral_class/Genome/Plant_species, data = metadata.rp.f, permutations = perm)
permanova.rp.f.nested
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Blocks:  with(metadata.rp.f, Block) 
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = relative.data.bray.rp.f ~ Ploidy/Ancestral_class/Genome/Plant_species, data = metadata.rp.f, permutations = perm)
#                                            Df SumOfSqs      R2      F Pr(>F)    
#Ploidy                                       2   0.7709 0.06277 1.9085 0.0021 ** 
#Ploidy:Ancestral_class                       3   1.1718 0.09541 1.9339 0.0003 ***
#Ploidy:Ancestral_class:Genome                2   0.5051 0.04113 1.2506 0.1107    
#Ploidy:Ancestral_class:Genome:Plant_species  3   0.9475 0.07715 1.5638 0.0053 ** 
#Residual                                    44   8.8866 0.72355                  
#Total                                       54  12.2819 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dispersion.rp.f<-betadisper(relative.data.bray.rp.f, group=metadata.rp.f$Ploidy)
permutest(dispersion.rp.f)
#
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.00217 0.0010870 0.1591    999  0.877
#Residuals 52 0.35525 0.0068317
plot(dispersion.rp.f, hull=FALSE, ellipse=TRUE)




####### Constrained ordination: Canonical Analysis of Principal coordinates (CAP) Bray-Curtis ##########


#Bray-Curtis distances
bray.rs.nf <- phyloseq::distance(physeq = physeq.norm.rs.nf, method = "bray")
bray.rs.f <- phyloseq::distance(physeq = physeq.norm.rs.f, method = "bray")
bray.rp.nf <- phyloseq::distance(physeq = physeq.norm.rp.nf, method = "bray")
bray.rp.f <- phyloseq::distance(physeq = physeq.norm.rp.f, method = "bray")



# Ploidy
cap.ord.rs.nf <- ordinate(
  physeq = physeq.norm.rs.nf, 
  method = "CAP",
  distance = bray.rs.nf,
  formula = ~ Ploidy
)
anova(cap.ord.rs.nf)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ploidy, data = data)
#         Df SumOfSqs      F Pr(>F)    
#Model     2   0.6217 1.6566  0.001 ***
#Residual 61  11.4468                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rs.f <- ordinate(
  physeq = physeq.norm.rs.f, 
  method = "CAP",
  distance = bray.rs.f,
  formula = ~ Ploidy
)
anova(cap.ord.rs.f)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ploidy, data = data)
#         Df SumOfSqs      F Pr(>F)    
#Model     2   1.0272 2.3923  0.001 ***
#Residual 60  12.8811                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp.nf <- ordinate(
  physeq = physeq.norm.rp.nf, 
  method = "CAP",
  distance = bray.rp.nf,
  formula = ~ Ploidy
)
anova(cap.ord.rp.nf)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ploidy, data = data)
#         Df SumOfSqs      F Pr(>F)   
#Model     2   0.6964 1.8744  0.004 **
#Residual 52   9.6591                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp.f <- ordinate(
  physeq = physeq.norm.rp.f, 
  method = "CAP",
  distance = bray.rp.f,
  formula = ~ Ploidy
)
anova(cap.ord.rp.f)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ploidy, data = data)
#         Df SumOfSqs      F Pr(>F)   
#Model     2   0.7709 1.7413  0.006 **
#Residual 52  11.5110                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Plots
cap.plot.rs.nf <- plot_ordination(
  physeq = physeq.norm.rs.nf, 
  ordination = cap.ord.rs.nf, 
  color = "Ploidy", 
  axes = c(1,2)) +
  geom_point(size=4) +
  colour.point.ploi +
  theme.new +
  theme(plot.title=element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title = "Rhizosphere")
cap.plot.rs.nf


cap.plot.rs.f<-plot_ordination(
  physeq = physeq.norm.rs.f, 
  ordination = cap.ord.rs.f, 
  color = "Ploidy", 
  axes = c(1,2)) +
  geom_point(size=4) +
  colour.point.ploi +
  theme.new +
  theme(plot.title=element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title = "Rhizosphere")
cap.plot.rs.f


cap.plot.rp.nf<-plot_ordination(
  physeq = physeq.norm.rp.nf, 
  ordination = cap.ord.rp.nf, 
  color = "Ploidy", 
  axes = c(1,2)) +
  geom_point(size=4) +
  colour.point.ploi +
  theme.new +
  theme(plot.title=element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title = "Rhizoplane")
cap.plot.rp.nf


cap.plot.rp.f<-plot_ordination(
  physeq = physeq.norm.rp.f, 
  ordination = cap.ord.rp.f, 
  color = "Ploidy", 
  axes = c(1,2)) + 
  geom_point(size=4) +
  colour.point.ploi +
  theme.new +
  theme(plot.title=element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title = "Rhizoplane")
cap.plot.rp.f


#Figure
cap.ploidy<-ggarrange(cap.plot.rs.nf, cap.plot.rp.nf, cap.plot.rs.f, cap.plot.rp.f,
                      nrow=2,
                      ncol=2,
                      common.legend = TRUE,
                      legend = "right",
                      align = "v",
                      labels = c("A", "", "B", "")
                      )
cap.ploidy

# Save as pdf
pdf(file = "Figure3/Plots/cap-ploidy.pdf",   #file name
           width = 8, # The width of the plot in inches
           height = 6) # The height of the plot in inches
cap.ploidy
dev.off()
#RStudioGD 
#2 


# Ancestral class
cap.ord.rs.nf <- ordinate(
  physeq = physeq.norm.rs.nf, 
  method = "CAP",
  distance = bray.rs.nf,
  formula = ~ Ancestral_class
)
anova(cap.ord.rs.nf)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ancestral_class, data = data)
#         Df SumOfSqs      F Pr(>F)   
#Model     2   0.6244 1.6641  0.001 **
#Residual 61  11.4442                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rs.f <- ordinate(
  physeq = physeq.norm.rs.f, 
  method = "CAP",
  distance = bray.rs.f,
  formula = ~ Ancestral_class
)
anova(cap.ord.rs.f)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ancestral_class, data = data)
#         Df SumOfSqs      F Pr(>F)    
#Model     2   1.1585 2.7258  0.001 ***
#Residual 60  12.7498                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp.nf <- ordinate(
  physeq = physeq.norm.rp.nf, 
  method = "CAP",
  distance = bray.rp.nf,
  formula = ~ Ancestral_class
)
anova(cap.ord.rp.nf)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ancestral_class, data = data)
#         Df SumOfSqs      F Pr(>F)   
#Model     2   0.7050 1.8993  0.003 **
#Residual 52   9.6505                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp.f <- ordinate(
  physeq = physeq.norm.rp.f, 
  method = "CAP",
  distance = bray.rp.f,
  formula = ~ Ancestral_class
)
anova(cap.ord.rp.f)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ancestral_class, data = data)
#         Df SumOfSqs      F Pr(>F)   
#Model     2   0.8823 2.0124  0.001 **
#Residual 52  11.3996                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Plots
cap.plot.rs.nf <- plot_ordination(
  physeq = physeq.norm.rs.nf, 
  ordination = cap.ord.rs.nf, 
  color = "Ancestral_class", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = ploidy.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rs.nf

cap.plot.rs.f<-plot_ordination(
  physeq = physeq.norm.rs.f, 
  ordination = cap.ord.rs.f, 
  color = "Ancestral_class", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = ploidy.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rs.f

cap.plot.rp.nf<-plot_ordination(
  physeq = physeq.norm.rp.nf, 
  ordination = cap.ord.rp.nf, 
  color = "Ancestral_class", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = ploidy.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rp.nf

cap.plot.rp.f<-plot_ordination(
  physeq = physeq.norm.rp.f, 
  ordination = cap.ord.rp.f, 
  color = "Ancestral_class", 
  axes = c(1,2)) + 
  geom_point(size=4) +
  scale_color_manual(values = ploidy.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rp.f


#Figures Ancestry
cap.ancestry<-ggarrange(cap.plot.rs.nf, cap.plot.rp.nf, cap.plot.rs.f, cap.plot.rp.f,
                        nrow=1,
                        common.legend = TRUE,
                        legend = "right",
                        labels = c("A","B","C","D")
                        )
cap.ancestry



# Genome
cap.ord.rs.nf <- ordinate(
  physeq = physeq.norm.rs.nf, 
  method = "CAP",
  distance = bray.rs.nf,
  formula = ~ Genome
)
anova(cap.ord.rs.nf)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Ploidy, data = data)
#         Df SumOfSqs      F Pr(>F)    
#Model     2   0.6217 1.6566  0.004 ***
#Residual 61  11.4468                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rs.f <- ordinate(
  physeq = physeq.norm.rs.f, 
  method = "CAP",
  distance = bray.rs.f,
  formula = ~Genome
)
anova(cap.ord.rs.f)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Genome, data = data)
#        Df SumOfSqs      F Pr(>F)   
#Model     4   1.4461 1.6826  0.003 **
#Residual 58  12.4622                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp.nf <- ordinate(
  physeq = physeq.norm.rp.nf, 
  method = "CAP",
  distance = bray.rp.nf,
  formula = ~ Genome
)
anova(cap.ord.rp.nf)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Genome, data = data)
#         Df SumOfSqs      F Pr(>F)  
#Model     4   1.1493 1.5605  0.015 *
#Residual 50   9.2061                
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp.f <- ordinate(
  physeq = physeq.norm.rp.f, 
  method = "CAP",
  distance = bray.rp.f,
  formula = ~ Genome
)
anova(cap.ord.rp.f)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Genome, data = data)
#         Df SumOfSqs      F Pr(>F)   
#Model     4   1.2976 1.4767  0.014 *
#Residual 50  10.9843                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Plots

#colours
red<-"#E64B35FF"
navy<-"#3C5488FF"

genome.palette<-c(pale.red, pale.blue, pale.green, red, navy)

cap.plot.rs.nf <- plot_ordination(
  physeq = physeq.norm.rs.nf, 
  ordination = cap.ord.rs.nf, 
  color = "Genome", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = genome.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rs.nf


cap.plot.rs.f <- plot_ordination(
  physeq = physeq.norm.rs.f, 
  ordination = cap.ord.rs.f, 
  color = "Genome", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = genome.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rs.f


cap.plot.rp.nf<-plot_ordination(
  physeq = physeq.norm.rp.nf, 
  ordination = cap.ord.rp.nf, 
  color = "Genome", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = genome.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rp.nf


cap.plot.rp.f <- plot_ordination(
  physeq = physeq.norm.rp.f, 
  ordination = cap.ord.rp.f, 
  color = "Genome", 
  axes = c(1,2)) + 
  geom_point(size=4) +
  scale_color_manual(values = genome.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rp.f


#Figures genome
cap.genome<-ggarrange(cap.plot.rs.nf, cap.plot.rp.nf, cap.plot.rs.f, cap.plot.rp.f,
                      nrow=1,
                      common.legend = TRUE,
                      legend = "right",
                      labels = c("E","F","G","H")
                      )
cap.genome



# Plant species
cap.ord.rs.nf <- ordinate(
  physeq = physeq.norm.rs.nf, 
  method = "CAP",
  distance = bray.rs.nf,
  formula = ~ Plant_species
)
anova(cap.ord.rs.nf)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Plant_species, data = data)
#         Df SumOfSqs      F Pr(>F)    
#Model    10   2.4730 1.3659  0.001 ***
#Residual 53   9.5956                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rs.f <- ordinate(
  physeq = physeq.norm.rs.f, 
  method = "CAP",
  distance = bray.rs.f,
  formula = ~ Plant_species
)
anova(cap.ord.rs.f)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Plant_species, data = data)
#         Df SumOfSqs      F Pr(>F)    
#Model    10   3.3504 1.6501  0.001 ***
#Residual 52  10.5579                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp.nf <- ordinate(
  physeq = physeq.norm.rp.nf, 
  method = "CAP",
  distance = bray.rp.nf,
  formula = ~ Plant_species
)
anova(cap.ord.rp.nf)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Plant_species, data = data)
#         Df SumOfSqs      F Pr(>F)   
#Model    10   2.6048 1.4787  0.003 **
#Residual 44   7.7506                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cap.ord.rp.f <- ordinate(
  physeq = physeq.norm.rp.f, 
  method = "CAP",
  distance = bray.rp.f,
  formula = ~ Plant_species
)
anova(cap.ord.rp.f)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: capscale(formula = distance ~ Plant_species, data = data)
#         Df SumOfSqs      F Pr(>F)    
#Model    10   3.3953 1.6811  0.001 ***
#Residual 44   8.8866                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Plots

#colours
turquoise<-"#6F99ADFF"
bright.red<-"#DC0000FF"
brown<-"#7E6148FF"
light.brown<-"#B09C85FF"
purple<-"#7876B1FF"
yellow<-"#FFDC91FF"
blue<-"#4DBBD5B2"

plant.species.palette<-c(pale.red, pale.blue, pale.green, red, navy, turquoise, bright.red, brown, light.brown, purple, yellow, blue)

cap.plot.rs.nf <- plot_ordination(
  physeq = physeq.norm.rs.nf, 
  ordination = cap.ord.rs.nf, 
  color = "Plant_species", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = plant.species.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rs.nf


cap.plot.rs.f <- plot_ordination(
  physeq = physeq.norm.rs.f, 
  ordination = cap.ord.rs.f, 
  color = "Plant_species", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = plant.species.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rs.f


cap.plot.rp.nf <- plot_ordination(
  physeq = physeq.norm.rp.nf, 
  ordination = cap.ord.rp.nf, 
  color = "Plant_species", 
  axes = c(1,2)) +
  geom_point(size=4) +
  scale_color_manual(values = plant.species.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rp.nf


cap.plot.rp.f <- plot_ordination(
  physeq = physeq.norm.rp.f, 
  ordination = cap.ord.rp.f, 
  color = "Plant_species", 
  axes = c(1,2)) + 
  geom_point(size=4) +
  scale_color_manual(values = plant.species.palette) +
  theme.new +
  geom_hline(yintercept = 0, linewidth=0.1) + geom_vline(xintercept = 0, linewidth=0.1) +
  labs(title="")
cap.plot.rp.f


#Figures species
cap.plant<-ggarrange(cap.plot.rs.nf, cap.plot.rp.nf, cap.plot.rs.f, cap.plot.rp.f,
                     nrow=1,
                     common.legend = TRUE,
                     legend = "none",
                     labels = c("I","J","K","L")
                     )
cap.plant


#Final figure
cap.figure<-ggarrange(cap.ancestry,
                      cap.genome,
                      cap.plant,
                      nrow=3,
                      align = "v")
cap.figure


# Save as pdf
pdf(file = "Figure3/Plots/cap-figure.pdf",   #file name
    width = 7, # The width of the plot in inches
    height = 6) # The height of the plot in inches
cap.figure
dev.off()
#RStudioGD 
#2





################## Differential Abundance Analysis using DESeq2 ##################

library(DESeq2); packageVersion("DESeq2") #v1.36.0
library(phyloseq); packageVersion("phyloseq") #v1.42.0
library(ggplot2); packageVersion("ggplot2") #v3.4.2


#Input files: Files produced by phyloseq (phyloseq object)


# Load data
CI_ASV <- readRDS("physeq_CI.rds")
CI_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1320 taxa and 244 samples ]
#sample_data() Sample Data:       [ 244 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1320 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1320 tips and 1307 internal nodes ]



# Subset data
CI_FRS_ASV <- subset_samples(CI_ASV,  Fertilization =="Fertilized" & Soil == "Rhizosphere")
CI_FRS_ASV <- filter_taxa(CI_FRS_ASV, function(x) sum(x) > 0, TRUE)
CI_FRS_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1114 taxa and 63 samples ]
#sample_data() Sample Data:       [ 63 samples by 22 sample variables ]
#tax_table()   Taxonomy Table:    [ 1114 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1114 tips and 1107 internal nodes ]

CI_NFRS_ASV <- subset_samples(CI_ASV,  Fertilization =="Non-fertilized" & Soil == "Rhizosphere")
CI_NFRS_ASV <- filter_taxa(CI_NFRS_ASV, function(x) sum(x) > 0, TRUE)
CI_NFRS_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1114 taxa and 64 samples ]
#sample_data() Sample Data:       [ 64 samples by 22 sample variables ]
#tax_table()   Taxonomy Table:    [ 1114 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1114 tips and 1109 internal nodes ]
CI_NFRS_ASV = subset_samples(CI_NFRS_ASV,sample_names(CI_NFRS_ASV) != "T.spelta.UnT.RS.3.38")
CI_NFRS_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1114 taxa and 63 samples ]
#sample_data() Sample Data:       [ 63 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1114 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1114 tips and 1109 internal nodes ]


CI_FRP_ASV <- subset_samples(CI_ASV,  Fertilization =="Fertilized" & Soil == "Rhizoplane")
CI_FRP_ASV <- filter_taxa(CI_FRP_ASV, function(x) sum(x) > 0, TRUE)
CI_FRP_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1218 taxa and 55 samples ]
#sample_data() Sample Data:       [ 55 samples by 22 sample variables ]
#tax_table()   Taxonomy Table:    [ 1218 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1218 tips and 1205 internal nodes ]

CI_NFRP_ASV <- subset_samples(CI_ASV,  Fertilization =="Non-fertilized" & Soil == "Rhizoplane")
CI_NFRP_ASV <- filter_taxa(CI_NFRP_ASV, function(x) sum(x) > 0, TRUE)
CI_NFRP_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1253 taxa and 55 samples ]
#sample_data() Sample Data:       [ 55 samples by 22 sample variables ]
#tax_table()   Taxonomy Table:    [ 1253 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1253 tips and 1242 internal nodes ]

#Removing anomalous samples from ASVs
##This was determined previously whereby one ASV was very abundant but in only one sample, therefore skewing the results
OTU <-otu_table(CI_NFRP_ASV)
SAM <-sample_data(CI_NFRP_ASV)
TAX <-tax_table(CI_NFRP_ASV)
PHY <- phy_tree(CI_NFRP_ASV)

OTU["8d75ef69afd1ef924fce01320df67008","Crusoe.UnT.RP.1.20"] <-0
OTU["53e57c62a122683c12fd7e145485d79b","T.polonicum.UnT.RP.4.12"] <-0
CI_NFRP_ASV <-phyloseq(OTU, PHY, TAX, SAM)
CI_NFRP_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1253 taxa and 55 samples ]
#sample_data() Sample Data:       [ 55 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1253 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1253 tips and 1242 internal nodes ]


CI_NFRP_ASV = subset_samples(CI_NFRP_ASV,sample_names(CI_NFRP_ASV) != "Crusoe.UnT.RP.1.20" & sample_names(CI_NFRP_ASV) != "T.polonicum.UnT.RP.4.12")
CI_NFRP_ASV
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1253 taxa and 53 samples ]
#sample_data() Sample Data:       [ 53 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 1253 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1253 tips and 1242 internal nodes ]


#Transform data to account for zeros
CI_FRS_ASV<-transform_sample_counts(CI_FRS_ASV, function(x) (x+1))
CI_NFRS_ASV<-transform_sample_counts(CI_NFRS_ASV, function(x) (x+1))
CI_FRP_ASV<-transform_sample_counts(CI_FRP_ASV, function(x) (x+1))
CI_NFRP_ASV<-transform_sample_counts(CI_NFRP_ASV, function(x) (x+1))



#DESeq2
ddsTAXA_CI_FRS <- phyloseq_to_deseq2(CI_FRS_ASV, ~ Ploidy)
#converting counts to integer mode
ddsTAXA_CI_NFRS <- phyloseq_to_deseq2(CI_NFRS_ASV, ~ Ploidy)
#converting counts to integer mode
ddsTAXA_CI_FRP <- phyloseq_to_deseq2(CI_FRP_ASV, ~ Ploidy)
#converting counts to integer mode
ddsTAXA_CI_NFRP <- phyloseq_to_deseq2(CI_NFRP_ASV, ~ Ploidy)
#converting counts to integer mode



#Run deseq - local regression fit automatically substituted
ddsTAXA_CI_FRS2 = DESeq(ddsTAXA_CI_FRS, test="Wald", fitType="local")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 133 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing
ddsTAXA_CI_NFRS2 = DESeq(ddsTAXA_CI_NFRS, test="Wald", fitType="local")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 120 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing
ddsTAXA_CI_FRP2 = DESeq(ddsTAXA_CI_FRP, test="Wald", fitType="local")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 311 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing
ddsTAXA_CI_NFRP2 = DESeq(ddsTAXA_CI_NFRP, test="Wald", fitType="local")
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#-- replacing outliers and refitting for 171 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing


#Check results
summary(ddsTAXA_CI_FRS2)
#[1] "DESeqDataSet object of length 1114 with 27 metadata columns"
summary(ddsTAXA_CI_NFRS2)
#[1] "DESeqDataSet object of length 1114 with 27 metadata columns"
summary(ddsTAXA_CI_FRP2)
#[1] "DESeqDataSet object of length 1218 with 27 metadata columns"
summary(ddsTAXA_CI_NFRP2)
#[1] "DESeqDataSet object of length 1253 with 27 metadata columns"



#Intercepts
resultsNames(ddsTAXA_CI_FRS2)
#"Intercept"                      "Ploidy_Tetraploid_vs_Diploid" "Ploidy_Hexaploid_vs_Diploid" 
resultsNames(ddsTAXA_CI_NFRS2)
#[1] "Intercept"                    "Ploidy_Tetraploid_vs_Diploid"
#[3] "Ploidy_Hexaploid_vs_Diploid"
resultsNames(ddsTAXA_CI_FRP2)
#"Intercept"                      "Ploidy_Tetraploid_vs_Diploid" "Ploidy_Hexaploid_vs_Diploid" 
resultsNames(ddsTAXA_CI_NFRP2)
#[1] "Intercept"                    "Ploidy_Tetraploid_vs_Diploid"
#[3] "Ploidy_Hexaploid_vs_Diploid"

####Tetraploid vs Hexaploid still missing



#Export results
ddsTAXA2_CI_NFRS_res <- results(ddsTAXA_CI_NFRS2)
ddsTAXA2_CI_FRS_res <- results(ddsTAXA_CI_FRS2)
ddsTAXA2_CI_NFRP_res <- results(ddsTAXA_CI_NFRP2)
ddsTAXA2_CI_FRP_res <- results(ddsTAXA_CI_FRP2)


#View results summary
summary(ddsTAXA2_CI_NFRP_res)
#
#out of 1253 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 104, 8.3%
#LFC < 0 (down)     : 43, 3.4%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#
summary(ddsTAXA2_CI_FRS_res)
#
#out of 1114 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 107, 9.6%
#LFC < 0 (down)     : 46, 4.1%
#outliers [1]       : 0, 0%
#low counts [2]     : 324, 29%
#(mean count < 2)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#
summary(ddsTAXA2_CI_NFRP_res)
#
#out of 1253 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 104, 8.3%
#LFC < 0 (down)     : 43, 3.4%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#
summary(ddsTAXA2_CI_FRP_res)
#
#out of 1218 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 86, 7.1%
#LFC < 0 (down)     : 46, 3.8%
#outliers [1]       : 0, 0%
#low counts [2]     : 354, 29%
#(mean count < 2)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#



# Perform contrasts between groups

#CI F RS
ddsTAXA2_CI_FRS_TvD <- results(ddsTAXA_CI_FRS2, contrast = c("Ploidy", "Tetraploid", "Diploid"), alpha=0.05)
ddsTAXA2_CI_FRS_TvB <- results(ddsTAXA_CI_FRS2, contrast = c("Ploidy", "Hexaploid","Tetraploid"), alpha=0.05)
ddsTAXA2_CI_FRS_DvB <- results(ddsTAXA_CI_FRS2, contrast = c("Ploidy", "Hexaploid", "Diploid"), alpha=0.05)

#CI NF RS
ddsTAXA2_CI_NFRS_TvD <- results(ddsTAXA_CI_NFRS2, contrast = c("Ploidy", "Tetraploid", "Diploid"), alpha=0.05)
ddsTAXA2_CI_NFRS_TvB <- results(ddsTAXA_CI_NFRS2, contrast = c("Ploidy","Hexaploid","Tetraploid"), alpha=0.05)
ddsTAXA2_CI_NFRS_DvB <- results(ddsTAXA_CI_NFRS2, contrast = c("Ploidy", "Hexaploid", "Diploid"), alpha=0.05)

#CI F RP
ddsTAXA2_CI_FRP_TvD <- results(ddsTAXA_CI_FRP2, contrast = c("Ploidy", "Tetraploid", "Diploid"), alpha=0.05)
ddsTAXA2_CI_FRP_TvB <- results(ddsTAXA_CI_FRP2, contrast = c("Ploidy", "Hexaploid","Tetraploid"), alpha=0.05)
ddsTAXA2_CI_FRP_DvB <- results(ddsTAXA_CI_FRP2, contrast = c("Ploidy", "Hexaploid", "Diploid"), alpha=0.05)

#CI NF RP
ddsTAXA2_CI_NFRP_TvD <- results(ddsTAXA_CI_NFRP2, contrast = c("Ploidy", "Tetraploid", "Diploid"), alpha=0.05)
ddsTAXA2_CI_NFRP_TvB <- results(ddsTAXA_CI_NFRP2, contrast = c("Ploidy","Hexaploid","Tetraploid"), alpha=0.05)
ddsTAXA2_CI_NFRP_DvB <- results(ddsTAXA_CI_NFRP2, contrast = c("Ploidy", "Hexaploid", "Diploid"), alpha=0.05)



#Matrix of results - all taxa
allTAXA_CI_FRS <- results(ddsTAXA_CI_FRS2)
allTAXA_CI_FRS = cbind(as(allTAXA_CI_FRS, "data.frame"), as(tax_table(CI_FRS_ASV)[rownames(allTAXA_CI_FRS), ], "matrix"))

allTAXA_CI_NFRS <- results(ddsTAXA_CI_NFRS2)
allTAXA_CI_NFRS = cbind(as(allTAXA_CI_NFRS, "data.frame"), as(tax_table(CI_NFRS_ASV)[rownames(allTAXA_CI_NFRS), ], "matrix"))

allTAXA_CI_FRP <- results(ddsTAXA_CI_FRP2)
allTAXA_CI_FRP = cbind(as(allTAXA_CI_FRP, "data.frame"), as(tax_table(CI_FRP_ASV)[rownames(allTAXA_CI_FRP), ], "matrix"))

allTAXA_CI_NFRP <- results(ddsTAXA_CI_NFRP2)
allTAXA_CI_NFRP = cbind(as(allTAXA_CI_NFRP, "data.frame"), as(tax_table(CI_NFRP_ASV)[rownames(allTAXA_CI_NFRP), ], "matrix"))




#Matrix of results - all taxa contrasts

#NF RP
allTAXA_CI_NFRP_DvB = cbind(as(ddsTAXA2_CI_NFRP_DvB, "data.frame"), as(tax_table(CI_NFRP_ASV)[rownames(ddsTAXA2_CI_NFRP_DvB), ], "matrix"))
allTAXA_CI_NFRP_TvB = cbind(as(ddsTAXA2_CI_NFRP_TvB, "data.frame"), as(tax_table(CI_NFRP_ASV)[rownames(ddsTAXA2_CI_NFRP_TvB), ], "matrix"))
allTAXA_CI_NFRP_TvD = cbind(as(ddsTAXA2_CI_NFRP_TvD, "data.frame"), as(tax_table(CI_NFRP_ASV)[rownames(ddsTAXA2_CI_NFRP_TvD), ], "matrix"))


#F RP
allTAXA_CI_FRP_TvD = cbind(as(ddsTAXA2_CI_FRP_TvD, "data.frame"), as(tax_table(CI_FRP_ASV)[rownames(ddsTAXA2_CI_FRP_TvD), ], "matrix"))
allTAXA_CI_FRP_TvB = cbind(as(ddsTAXA2_CI_FRP_TvB, "data.frame"), as(tax_table(CI_FRP_ASV)[rownames(ddsTAXA2_CI_FRP_TvB), ], "matrix"))
allTAXA_CI_FRP_DvB = cbind(as(ddsTAXA2_CI_FRP_DvB, "data.frame"), as(tax_table(CI_FRP_ASV)[rownames(ddsTAXA2_CI_FRP_DvB), ], "matrix"))


#NF RS
allTAXA_CI_NFRS_DvB = cbind(as(ddsTAXA2_CI_NFRS_DvB, "data.frame"), as(tax_table(CI_NFRS_ASV)[rownames(ddsTAXA2_CI_NFRS_DvB), ], "matrix"))
allTAXA_CI_NFRS_TvB = cbind(as(ddsTAXA2_CI_NFRS_TvB, "data.frame"), as(tax_table(CI_NFRS_ASV)[rownames(ddsTAXA2_CI_NFRS_TvB), ], "matrix"))
allTAXA_CI_NFRS_TvD = cbind(as(ddsTAXA2_CI_NFRS_TvD, "data.frame"), as(tax_table(CI_NFRS_ASV)[rownames(ddsTAXA2_CI_NFRS_TvD), ], "matrix"))


#F RS
allTAXA_CI_FRS_TvD = cbind(as(ddsTAXA2_CI_FRS_TvD, "data.frame"), as(tax_table(CI_FRS_ASV)[rownames(ddsTAXA2_CI_FRS_TvD), ], "matrix"))
allTAXA_CI_FRS_TvB = cbind(as(ddsTAXA2_CI_FRS_TvB, "data.frame"), as(tax_table(CI_FRS_ASV)[rownames(ddsTAXA2_CI_FRS_TvB), ], "matrix"))
allTAXA_CI_FRS_DvB = cbind(as(ddsTAXA2_CI_FRS_DvB, "data.frame"), as(tax_table(CI_FRS_ASV)[rownames(ddsTAXA2_CI_FRS_DvB), ], "matrix"))



## Significant Taxa

alpha <- 0.05

#F RP
sigTAXA_CI_FRP_TvD <- subset(allTAXA_CI_FRP_TvD, padj < alpha)
sigTAXA_CI_FRP_TvB <- subset(allTAXA_CI_FRP_TvB, padj < alpha)
sigTAXA_CI_FRP_DvB <- subset(allTAXA_CI_FRP_DvB, padj < alpha)

writexl::write_xlsx(sigTAXA_CI_FRP_TvD, "Figure3/Data/deseq-ternary-CI-FRP-TvD.xlsx")
writexl::write_xlsx(sigTAXA_CI_FRP_TvB, "Figure3/Data/deseq-ternary-CI-FRP-TvB.xlsx")
writexl::write_xlsx(sigTAXA_CI_FRP_DvB, "Figure3/Data/deseq-ternary-CI-FRP-DvB.xlsx")
writexl::write_xlsx(allTAXA_CI_FRP, "Figure3/Data/deseq-ternary-CI-FRP.xlsx")


#NF RP
sigTAXA_CI_NFRP_TvD <- subset(allTAXA_CI_NFRP_TvD, padj < alpha)
sigTAXA_CI_NFRP_TvB <- subset(allTAXA_CI_NFRP_TvB, padj < alpha)
sigTAXA_CI_NFRP_DvB <- subset(allTAXA_CI_NFRP_DvB, padj < alpha)

writexl::write_xlsx(sigTAXA_CI_NFRP_TvD, "Figure3/Data/deseq-ternary-CI-NFRP-TvD.xlsx")
writexl::write_xlsx(sigTAXA_CI_NFRP_TvB, "Figure3/Data/deseq-ternary-CI-NFRP-TvB.xlsx")
writexl::write_xlsx(sigTAXA_CI_NFRP_DvB, "Figure3/Data/deseq-ternary-CI-NFRP-DvB.xlsx")
writexl::write_xlsx(allTAXA_CI_NFRP, "Figure3/Data/deseq-ternary-CI-NFRP.xlsx")


#F RS
sigTAXA_CI_FRS_TvD <- subset(allTAXA_CI_FRS_TvD, padj < alpha)
sigTAXA_CI_FRS_TvB <- subset(allTAXA_CI_FRS_TvB, padj < alpha)
sigTAXA_CI_FRS_DvB <- subset(allTAXA_CI_FRS_DvB, padj < alpha)

writexl::write_xlsx(sigTAXA_CI_FRS_TvD, "Figure3/Data/deseq-ternary-CI-FRS-TvD.xlsx")
writexl::write_xlsx(sigTAXA_CI_FRS_TvB, "Figure3/Data/deseq-ternary-CI-FRS-TvB.xlsx")
writexl::write_xlsx(sigTAXA_CI_FRS_DvB, "Figure3/Data/deseq-ternary-CI-FRS-DvB.xlsx")
writexl::write_xlsx(allTAXA_CI_FRS, "Figure3/Data/deseq-ternary-CI-FRS.xlsx")


#NF RS
sigTAXA_CI_NFRS_TvD <- subset(allTAXA_CI_NFRS_TvD, padj < alpha)
sigTAXA_CI_NFRS_TvB <- subset(allTAXA_CI_NFRS_TvB, padj < alpha)
sigTAXA_CI_NFRS_DvB <- subset(allTAXA_CI_NFRS_DvB, padj < alpha)

writexl::write_xlsx(sigTAXA_CI_NFRS_TvD, "Figure3/Data/deseq-ternary-CI-NFRS-TvD.xlsx")
writexl::write_xlsx(sigTAXA_CI_NFRS_TvB, "Figure3/Data/deseq-ternary-CI-NFRS-TvB.xlsx")
writexl::write_xlsx(sigTAXA_CI_NFRS_DvB, "Figure3/Data/deseq-ternary-CI-NFRS-DvB.xlsx")
writexl::write_xlsx(allTAXA_CI_NFRS, "Figure3/Data/deseq-ternary-CI-NFRS.xlsx")



# Ternary ####


### Fertilized RP ###

# need to get the means of each compartment for the ternary plots 
TAXA2baseMeanPerLvl <- sapply(levels(ddsTAXA_CI_FRP2$Ploidy),function(lvl) rowMeans( counts(ddsTAXA_CI_FRP2, normalized=TRUE)[,ddsTAXA_CI_FRP2$Ploidy == lvl, drop=F] ) )

head(TAXA2baseMeanPerLvl)
#                                  Diploid Tetraploid Hexaploid
#c114f024375243308364cc8f91faeabf 3.383694   2.673463  2.779368
#f7df07512b80d985b360fc621a0893ea 1.253085   2.603363  1.251886
#c66c38d5daf4bd0a030031898490cc61 1.253085   3.710876  1.251886
#e3bdbe6e720b5021e9ae6358b9eae180 1.253085   1.426074  1.251886
#6d76d0da1b52f9ee5e64554dab012c3e 2.577888   3.381112  1.837682
#19069dcdac5f474ae061bade1059fbbc 4.300133   6.391436  3.322629
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
#   Diploid  Hexaploid Tetraploid 
#  28205.53   25772.20   23048.33


#Defining significant ASVs
TetraploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_FRP_TvB, padj<0.05 & log2FoldChange<0)), 
  rownames(subset(ddsTAXA2_CI_FRP_TvD, padj<0.05 & log2FoldChange>0))))

DiploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_FRP_TvD, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CI_FRP_DvB, padj<0.05 & log2FoldChange<0))))

HexaploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_FRP_TvB, padj<0.05 & log2FoldChange>0)),
  rownames(subset(ddsTAXA2_CI_FRP_DvB, padj<0.05 & log2FoldChange>0))))


#get phylum levels
get_taxa_unique(CI_FRP_ASV, "Phylum")
#[1] "Proteobacteria"    "Actinobacteriota"  "Cyanobacteria"     "Bdellovibrionota" 
#[5] "Acidobacteriota"   "Bacteroidota"      "Elusimicrobiota"   "Armatimonadota"   
#[9] "Chloroflexi"       "Firmicutes"        "Spirochaetota"     "Latescibacterota" 
#[13] "FCPU426"           "WPS_2"             "Methylomirabilota" "Planctomycetota"  
#[17] "Verrucomicrobiota" "Gemmatimonadota"   "Nitrospirota"      "Patescibacteria"  
#[21] "Dependentiae"      "Deinococcota"      "Myxococcota"       "RCP2_54"          
#[25] "Desulfobacterota"  "Fibrobacterota" 

#Colour codes
Acidobacteriota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Acidobacteriota")))
Actinobacteria <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Actinobacteriota")))
Bacteroidetes <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Bacteroidota")))
Chloroflexi <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Chloroflexi")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Desulfobacterota")))
Firmicutes <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Firmicutes")))
Nitrospirota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Nitrospirota")))
Proteobacteria <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Proteobacteria")))
Armatimonadota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Armatimonadota")))
Bdellovibrionota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Bdellovibrionota")))
Cyanobacteria <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Cyanobacteria")))
Deinococcota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Deinococcota")))
Dependentiae <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Dependentiae")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Desulfobacterota")))
Elusimicrobiota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Elusimicrobiota")))
FCPU426 <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="FCPU426")))
Fibrobacterota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Fibrobacterota")))
Gemmatimonadota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Gemmatimonadota")))
Latescibacterota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Latescibacterota")))
Methylomirabilota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Methylomirabilota")))
Myxococcota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Myxococcota")))
Patescibacteria <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Patescibacteria")))
Planctomycetota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Planctomycetota")))
RCP2_54 <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="RCP2_54")))
Spirochaetota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Spirochaetota")))
Verrucomicrobiota <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="Verrucomicrobiota")))
WPS_2 <- unique(rownames(subset(allTAXA_CI_FRP, Phylum=="WPS_2")))


#Colours
red<-"#E64B35FF"
blue<-"#4DBBD5B2"
turquoise<-"#6F99ADFF"
navy<-"#3C5488FF"
pale_red<-"#F39B7FFF"
light_blue<-"#8491B4FF"
light_green<-"#91D1C2FF"
bright_red<-"#DC0000FF"
brown<-"#7E6148FF"
light_brown<-"#B09C85FF"
purple<-"#7876B1FF"
yellow<-"#FFDC91FF"
dark_yellow<-"#edae49"
green<-"#66a182"



#Colour codes for plotting
TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% Acidobacteriota, pale_red,"#004d40")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_colors[Actinobacteria] <- red
TAXA2ternary_colors[Bacteroidetes] <- blue
TAXA2ternary_colors[Chloroflexi] <- light_blue
TAXA2ternary_colors[Firmicutes] <- turquoise
TAXA2ternary_colors[Nitrospirota] <- brown
TAXA2ternary_colors[Proteobacteria] <- navy
TAXA2ternary_colors[Armatimonadota] <- bright_red
TAXA2ternary_colors[Bdellovibrionota] <- purple
TAXA2ternary_colors[Cyanobacteria] <- "skyblue4"
TAXA2ternary_colors[Deinococcota] <- "khaki2"
TAXA2ternary_colors[Dependentiae] <- "yellow2"
TAXA2ternary_colors[Desulfobacterota] <- dark_yellow
TAXA2ternary_colors[Elusimicrobiota] <- yellow
TAXA2ternary_colors[FCPU426] <- "yellowgreen"
TAXA2ternary_colors[Fibrobacterota] <- "wheat4"
TAXA2ternary_colors[Gemmatimonadota] <- "darkorchid"
TAXA2ternary_colors[Latescibacterota] <- "darkorange3"
TAXA2ternary_colors[Methylomirabilota] <- "salmon4"
TAXA2ternary_colors[Myxococcota] <- "saddlebrown"
TAXA2ternary_colors[Patescibacteria] <- light_brown
TAXA2ternary_colors[Planctomycetota] <- "steelblue4"
TAXA2ternary_colors[RCP2_54] <- "royalblue"
TAXA2ternary_colors[Spirochaetota] <- "grey37"
TAXA2ternary_colors[Verrucomicrobiota] <- "gray0"
TAXA2ternary_colors[WPS_2] <- "black"


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
library(grid)
source("tern_e_modified.R", local = TRUE)
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)


# Create plot and save as pdf
pdf(file = "Figure3/Plots/deseq-ternary-CI-FRP.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)
dev.off()
#RStudioGD 
#2 



##### Non-fertilized Rhizoplane

# need to get the means of each compartment for the ternary plots
TAXA2baseMeanPerLvl <- sapply(levels(ddsTAXA_CI_NFRP2$Ploidy),function(lvl) rowMeans( counts(ddsTAXA_CI_NFRP2, normalized=TRUE)[,ddsTAXA_CI_NFRP2$Ploidy == lvl, drop=F] ) )

head(TAXA2baseMeanPerLvl)
#                                   Diploid Tetraploid Hexaploid
#c114f024375243308364cc8f91faeabf  1.822457   3.363135  2.555440
#f7df07512b80d985b360fc621a0893ea  2.140071   3.674439  4.949886
#c66c38d5daf4bd0a030031898490cc61  1.579820   6.432817  7.253384
#e3bdbe6e720b5021e9ae6358b9eae180 12.373762   6.977500  8.639618
#6d76d0da1b52f9ee5e64554dab012c3e 14.837557   5.638880  9.550740
#19069dcdac5f474ae061bade1059fbbc 18.195468  33.455851 27.640100
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
#Tetraploid  Hexaploid    Diploid 
#  21776.51   21327.46   19952.10


#Define ASVs
TetraploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_NFRP_TvB, padj<0.05 & log2FoldChange<0)), 
  rownames(subset(ddsTAXA2_CI_NFRP_TvD, padj<0.05 & log2FoldChange>0))))

DiploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_NFRP_TvD, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CI_NFRP_DvB, padj<0.05 & log2FoldChange<0))))

HexaploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_NFRP_TvB, padj<0.05 & log2FoldChange>0)),
  rownames(subset(ddsTAXA2_CI_NFRP_DvB, padj<0.05 & log2FoldChange>0))))


#get phylum levels
get_taxa_unique(CI_NFRP_ASV, "Phylum")
#[1] "Proteobacteria"    "Actinobacteriota"  "Cyanobacteria"     "Bdellovibrionota" 
#[5] "Acidobacteriota"   "Bacteroidota"      "Firmicutes"        "Elusimicrobiota"  
#[9] "Armatimonadota"    "Chloroflexi"       "Spirochaetota"     "Latescibacterota" 
#[13] "FCPU426"           "WPS_2"             "Methylomirabilota" "Planctomycetota"  
#[17] "Verrucomicrobiota" "Gemmatimonadota"   "Nitrospirota"      "Patescibacteria"  
#[21] "Dependentiae"      "Deinococcota"      "Myxococcota"       "RCP2_54"          
#[25] "Desulfobacterota"  "Fibrobacterota" 


#Set colours
Acidobacteriota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Acidobacteriota")))
Actinobacteria <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Actinobacteriota")))
Bacteroidetes <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Bacteroidota")))
Chloroflexi <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Chloroflexi")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Desulfobacterota")))
Firmicutes <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Firmicutes")))
Nitrospirota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Nitrospirota")))
Proteobacteria <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Proteobacteria")))
Armatimonadota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Armatimonadota")))
Bdellovibrionota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Bdellovibrionota")))
Cyanobacteria <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Cyanobacteria")))
Deinococcota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Deinococcota")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Desulfobacterota")))
Elusimicrobiota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Elusimicrobiota")))
FCPU426 <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="FCPU426")))
Fibrobacterota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Fibrobacterota")))
Gemmatimonadota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Gemmatimonadota")))
Latescibacterota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Latescibacterota")))
Methylomirabilota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Methylomirabilota")))
Myxococcota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Myxococcota")))
Patescibacteria <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Patescibacteria")))
Planctomycetota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Planctomycetota")))
RCP2_54 <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="RCP2_54")))
Spirochaetota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Spirochaetota")))
Verrucomicrobiota <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="Verrucomicrobiota")))
WPS_2 <- unique(rownames(subset(allTAXA_CI_NFRP, Phylum=="WPS_2")))



#Colour codes for plotting
TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% Acidobacteriota, pale_red,"#004d40")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_colors[Actinobacteria] <- red
TAXA2ternary_colors[Bacteroidetes] <- blue
TAXA2ternary_colors[Chloroflexi] <- light_blue
TAXA2ternary_colors[Firmicutes] <- turquoise
TAXA2ternary_colors[Nitrospirota] <- brown
TAXA2ternary_colors[Proteobacteria] <- navy
TAXA2ternary_colors[Armatimonadota] <- bright_red
TAXA2ternary_colors[Bdellovibrionota] <- purple
TAXA2ternary_colors[Cyanobacteria] <- "skyblue4"
TAXA2ternary_colors[Deinococcota] <- "khaki2"
TAXA2ternary_colors[Desulfobacterota] <- dark_yellow
TAXA2ternary_colors[Elusimicrobiota] <- yellow
TAXA2ternary_colors[FCPU426] <- "yellowgreen"
TAXA2ternary_colors[Fibrobacterota] <- "wheat4"
TAXA2ternary_colors[Gemmatimonadota] <- "darkorchid"
TAXA2ternary_colors[Latescibacterota] <- "darkorange3"
TAXA2ternary_colors[Methylomirabilota] <- "salmon4"
TAXA2ternary_colors[Myxococcota] <- "saddlebrown"
TAXA2ternary_colors[Patescibacteria] <- light_brown
TAXA2ternary_colors[Planctomycetota] <- "steelblue4"
TAXA2ternary_colors[RCP2_54] <- "royalblue"
TAXA2ternary_colors[Spirochaetota] <- "grey37"
TAXA2ternary_colors[Verrucomicrobiota] <- "gray0"
TAXA2ternary_colors[WPS_2] <- "black"



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
pdf(file = "Figure3/Plots/deseq-ternary-CI-NFRP.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)
dev.off()
#RStudioGD 
#2 



### Fertilized RS ###

# need to get the means of each compartment for the ternary plots 
TAXA2baseMeanPerLvl <- sapply(levels(ddsTAXA_CI_FRS2$Ploidy),function(lvl) rowMeans( counts(ddsTAXA_CI_FRS2, normalized=TRUE)[,ddsTAXA_CI_FRS2$Ploidy == lvl, drop=F] ) )

head(TAXA2baseMeanPerLvl)
#                                  Diploid Tetraploid Hexaploid
#c114f024375243308364cc8f91faeabf 2.469682   2.917229  4.272441
#f7df07512b80d985b360fc621a0893ea 1.203691   2.642424  2.482573
#c66c38d5daf4bd0a030031898490cc61 8.283112   7.307997  3.813164
#6d76d0da1b52f9ee5e64554dab012c3e 1.203691   2.072560  2.937472
#19069dcdac5f474ae061bade1059fbbc 1.203691   1.190277  1.575141
#083a72b30daab20c9fc08a5cf2534bbf 1.203691   3.395985  5.499712
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
#  Hexaploid    Diploid Tetraploid 
#  10459.962  10217.541   8879.186 



#Define ASVs
TetraploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_FRS_TvB, padj<0.05 & log2FoldChange<0)), 
  rownames(subset(ddsTAXA2_CI_FRS_TvD, padj<0.05 & log2FoldChange>0))))

DiploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_FRS_TvD, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CI_FRS_DvB, padj<0.05 & log2FoldChange<0))))

HexaploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_FRS_TvB, padj<0.05 & log2FoldChange>0)),
  rownames(subset(ddsTAXA2_CI_FRS_DvB, padj<0.05 & log2FoldChange>0))))


#get phylum levels
get_taxa_unique(CI_FRS_ASV, "Phylum")
#[1] "Proteobacteria"    "Actinobacteriota"  "Cyanobacteria"     "Bdellovibrionota" 
#[5] "Acidobacteriota"   "Bacteroidota"      "Firmicutes"        "Elusimicrobiota"  
#[9] "Armatimonadota"    "Chloroflexi"       "Spirochaetota"     "Latescibacterota" 
#[13] "FCPU426"           "WPS_2"             "Methylomirabilota" "Planctomycetota"  
#[17] "Verrucomicrobiota" "Gemmatimonadota"   "Nitrospirota"      "Patescibacteria"  
#[21] "Dependentiae"      "Deinococcota"      "Myxococcota"       "RCP2_54"          
#[25] "Desulfobacterota"  "Fibrobacterota" 


#Colour codes based on phyla
Acidobacteriota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Acidobacteriota")))
Actinobacteria <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Actinobacteriota")))
Bacteroidetes <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Bacteroidota")))
Chloroflexi <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Chloroflexi")))
Deinococcota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Deinococcota")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Desulfobacterota")))
Firmicutes <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Firmicutes")))
Nitrospirota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Nitrospirota")))
Proteobacteria <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Proteobacteria")))
Armatimonadota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Armatimonadota")))
Bdellovibrionota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Bdellovibrionota")))
Cyanobacteria <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Cyanobacteria")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Desulfobacterota")))
Elusimicrobiota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Elusimicrobiota")))
FCPU426 <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="FCPU426")))
Fibrobacterota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Fibrobacterota")))
Gemmatimonadota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Gemmatimonadota")))
Latescibacterota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Latescibacterota")))
Methylomirabilota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Methylomirabilota")))
Myxococcota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Myxococcota")))
Patescibacteria <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Patescibacteria")))
Planctomycetota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Planctomycetota")))
RCP2_54 <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="RCP2_54")))
Spirochaetota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Spirochaetota")))
Verrucomicrobiota <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="Verrucomicrobiota")))
WPS_2 <- unique(rownames(subset(allTAXA_CI_FRS, Phylum=="WPS_2")))


#Colour codes for plotting
TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% Acidobacteriota, pale_red,"#004d40")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_colors[Actinobacteria] <- red
TAXA2ternary_colors[Bacteroidetes] <- blue
TAXA2ternary_colors[Chloroflexi] <- light_blue
TAXA2ternary_colors[Firmicutes] <- turquoise
TAXA2ternary_colors[Nitrospirota] <- brown
TAXA2ternary_colors[Proteobacteria] <- navy
TAXA2ternary_colors[Armatimonadota] <- bright_red
TAXA2ternary_colors[Bdellovibrionota] <- purple
TAXA2ternary_colors[Cyanobacteria] <- "skyblue4"
TAXA2ternary_colors[Deinococcota] <- "khaki2"
TAXA2ternary_colors[Desulfobacterota] <- dark_yellow
TAXA2ternary_colors[Elusimicrobiota] <- yellow
TAXA2ternary_colors[FCPU426] <- "yellowgreen"
TAXA2ternary_colors[Fibrobacterota] <- "wheat4"
TAXA2ternary_colors[Gemmatimonadota] <- "darkorchid"
TAXA2ternary_colors[Latescibacterota] <- "darkorange3"
TAXA2ternary_colors[Methylomirabilota] <- "salmon4"
TAXA2ternary_colors[Myxococcota] <- "saddlebrown"
TAXA2ternary_colors[Patescibacteria] <- light_brown
TAXA2ternary_colors[Planctomycetota] <- "steelblue4"
TAXA2ternary_colors[RCP2_54] <- "royalblue"
TAXA2ternary_colors[Spirochaetota] <- "grey37"
TAXA2ternary_colors[Verrucomicrobiota] <- "gray0"
TAXA2ternary_colors[WPS_2] <- "black"


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
pdf(file = "Figure3/Plots/deseq-ternary-CI-FRS.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tern_e_modified(TAXA2baseMeanPerLvl, col=TAXA2ternary_colors, shape=TAXA2ternary_shapes, grid = TRUE, labels = "outside", grid_color = "white", prop_size = 3)
dev.off()
#RStudioGD 
#2 


### NonFertilized RS ###

# need to get the means of each compartment for the ternary plots 
TAXA2baseMeanPerLvl <- sapply(levels(ddsTAXA_CI_NFRS2$Ploidy),function(lvl) rowMeans( counts(ddsTAXA_CI_NFRS2, normalized=TRUE)[,ddsTAXA_CI_NFRS2$Ploidy == lvl, drop=F] ) )

head(TAXA2baseMeanPerLvl)
#                                  Diploid Tetraploid Hexaploid
#c114f024375243308364cc8f91faeabf 1.215968   1.549133  4.395591
#f7df07512b80d985b360fc621a0893ea 7.417786   3.757900  7.378476
#c66c38d5daf4bd0a030031898490cc61 9.221895  10.282567  8.483635
#e3bdbe6e720b5021e9ae6358b9eae180 3.594048   1.245088  2.800300
#6d76d0da1b52f9ee5e64554dab012c3e 7.333094   7.894852  8.827008
#19069dcdac5f474ae061bade1059fbbc 5.649119  12.911017  9.827918
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
#  Hexaploid    Diploid Tetraploid 
#   9306.136   9171.448   8317.141


#Define ASVs
TetraploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_NFRS_TvB, padj<0.05 & log2FoldChange<0)), 
  rownames(subset(ddsTAXA2_CI_NFRS_TvD, padj<0.05 & log2FoldChange>0))))

DiploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_NFRS_TvD, padj<0.05 & log2FoldChange<0)),
  rownames(subset(ddsTAXA2_CI_NFRS_DvB, padj<0.05 & log2FoldChange<0))))

HexaploidTaxa <- unique(c(
  rownames(subset(ddsTAXA2_CI_NFRS_TvB, padj<0.05 & log2FoldChange>0)),
  rownames(subset(ddsTAXA2_CI_NFRS_DvB, padj<0.05 & log2FoldChange>0))))


#get phylum levels
get_taxa_unique(CI_NFRS_ASV, "Phylum")
#[1] "Proteobacteria"    "Actinobacteriota"  "Cyanobacteria"     "Bdellovibrionota" 
#[5] "Acidobacteriota"   "Bacteroidota"      "Firmicutes"        "Elusimicrobiota"  
#[9] "Armatimonadota"    "Chloroflexi"       "Spirochaetota"     "Latescibacterota" 
#[13] "FCPU426"           "WPS_2"             "Methylomirabilota" "Planctomycetota"  
#[17] "Verrucomicrobiota" "Gemmatimonadota"   "Nitrospirota"      "Patescibacteria"  
#[21] "Dependentiae"      "Deinococcota"      "Myxococcota"       "RCP2_54"          
#[25] "Desulfobacterota"  "Fibrobacterota" 


#Colour codes absed on phyla
Acidobacteriota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Acidobacteriota")))
Actinobacteria <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Actinobacteriota")))
Bacteroidetes <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Bacteroidota")))
Chloroflexi <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Chloroflexi")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Desulfobacterota")))
Firmicutes <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Firmicutes")))
Nitrospirota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Nitrospirota")))
Proteobacteria <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Proteobacteria")))
Armatimonadota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Armatimonadota")))
Bdellovibrionota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Bdellovibrionota")))
Cyanobacteria <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Cyanobacteria")))
Deinococcota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Deinococcota")))
Desulfobacterota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Desulfobacterota")))
Elusimicrobiota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Elusimicrobiota")))
FCPU426 <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="FCPU426")))
Fibrobacterota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Fibrobacterota")))
Gemmatimonadota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Gemmatimonadota")))
Latescibacterota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Latescibacterota")))
Methylomirabilota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Methylomirabilota")))
Myxococcota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Myxococcota")))
Patescibacteria <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Patescibacteria")))
Planctomycetota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Planctomycetota")))
RCP2_54 <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="RCP2_54")))
Spirochaetota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Spirochaetota")))
Verrucomicrobiota <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="Verrucomicrobiota")))
WPS_2 <- unique(rownames(subset(allTAXA_CI_NFRS, Phylum=="WPS_2")))


#Colour codes for plotting
TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvl) %in% Acidobacteriota, pale_red,"#004d40")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvl)
TAXA2ternary_colors[Actinobacteria] <- red
TAXA2ternary_colors[Bacteroidetes] <- blue
TAXA2ternary_colors[Chloroflexi] <- light_blue
TAXA2ternary_colors[Firmicutes] <- turquoise
TAXA2ternary_colors[Nitrospirota] <- brown
TAXA2ternary_colors[Proteobacteria] <- navy
TAXA2ternary_colors[Armatimonadota] <- bright_red
TAXA2ternary_colors[Bdellovibrionota] <- purple
TAXA2ternary_colors[Cyanobacteria] <- "skyblue4"
TAXA2ternary_colors[Deinococcota] <- "khaki2"
TAXA2ternary_colors[Desulfobacterota] <- dark_yellow
TAXA2ternary_colors[Elusimicrobiota] <- yellow
TAXA2ternary_colors[FCPU426] <- "yellowgreen"
TAXA2ternary_colors[Fibrobacterota] <- "wheat4"
TAXA2ternary_colors[Gemmatimonadota] <- "darkorchid"
TAXA2ternary_colors[Latescibacterota] <- "darkorange3"
TAXA2ternary_colors[Methylomirabilota] <- "salmon4"
TAXA2ternary_colors[Myxococcota] <- "saddlebrown"
TAXA2ternary_colors[Patescibacteria] <- light_brown
TAXA2ternary_colors[Planctomycetota] <- "steelblue4"
TAXA2ternary_colors[RCP2_54] <- "royalblue"
TAXA2ternary_colors[Spirochaetota] <- "grey37"
TAXA2ternary_colors[Verrucomicrobiota] <- "gray0"
TAXA2ternary_colors[WPS_2] <- "black"


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
pdf(file = "Figure3/Plots/deseq-ternary-CI-NFRS.pdf",   #file name
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


#Phyla
phyla <- c("Acidobacteriota", 
           "Actinobacteria",
           "Armatimonadota",
           "Bacteroidetes",
           "Bdellovibrionota",
           "Chloroflexi",
           "Cyanobacteria",
           "Deinococcota",
           "Dependentiae",
           "Desulfobacterota", 
           "Elusimicrobiota", 
           "FCPU426", 
           "Fibrobacterota",
           "Firmicutes", 
           "Gemmatimonadota", 
           "Latescibacterota", 
           "Methylomirabilota", 
           "Myxococcota", 
           "Nitrospirota", 
           "Patescibacteria", 
           "Planctomycetota", 
           "Proteobacteria", 
           "RCP2_54", 
           "Spirochaetota", 
           "Verrucomicrobiota", 
           "WPS_2")

clrs <- c(pale_red,
          red,
          bright_red,
          blue, purple,
          light_blue,
          "skyblue4",
          "khaki2",
          "yellow2",
          dark_yellow,
          yellow,
          "yellowgreen",
          "wheat4",
          turquoise,
          "darkorchid",
          "darkorange3",
          "salmon4",
          "saddlebrown",
          brown,
          light_brown,
          "steelblue4",
          navy,
          "royalblue",
          "grey37",
          "gray",
          "black")


#Ploidy
##names <- c("Not significant", "Diploid", "Tetraploid", "Hexaploid")
##shape <- c(19, 17, 15, 18)


#Abundance
names <- c("1.0", "0.75", "0.5", "0.25")

sizes <- c(3.0,
           3.0*0.75,
           3.0*0.5,
           3.0*0.25)


# Check legends

#Phyla
plot.new()
legend("left",  legend = phyla, pch=16, cex=0.5, 
       bty='n', col = clrs, title = "Phyla", title.adj = 0)

#ASVs
plot.new()
legend("center",  legend = names, pch=19, cex=0.5, y.intersp = 3, x.intersp = 3,
       bty='n', pt.cex = sizes, title = "Not significant", title.adj = 0)

legend("center",  legend = names, pch=17, cex=0.5, y.intersp = 3, x.intersp = 3,
       bty='n', pt.cex = sizes, title = "Significant", title.adj = 0)



# Create legend and save as pdf

#Phyla
pdf(file = "Figure3/Plots/deseq-ternary-CI-legend.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 8) # The height of the plot in inches
plot.new()
legend("bottom",  legend = phyla, pch=16, cex=0.5, 
       bty='n', col = clrs, title = "Phyla", title.adj = 0)
dev.off()
#RStudioGD 
#2 

#Shapes
pdf(file = "Figure3/Plots/deseq-ternary-CI-legend-shapes.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot.new()
legend("right",  legend = names, pch=19, cex=0.5, y.intersp = 3, x.intersp = 3,
       bty='n', pt.cex = sizes, title = "Not significant", title.adj = 0)
legend("center",  legend = names, pch=17, cex=0.5, y.intersp = 3, x.intersp = 3,
       bty='n', pt.cex = sizes, title = "Significant", title.adj = 0)
dev.off()
#RStudioGD 
#2 

#Abundances
pdf(file = "Figure3/Plots/deseq-ternary-CI-legend-shapes.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot.new()
legend("left",  legend = names, pch=19, 
       bty='n', pt.cex = sizes, title = "Not significant", title.adj = 0)
legend("center",  legend = names, pch=17, 
       bty='n', pt.cex = sizes, title = "Significant", title.adj = 0)
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
#[1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#[1] vegan_2.6-4     lattice_0.21-8  permute_0.9-7   plyr_1.8.8      phyloseq_1.42.0 ggplot2_3.4.2  
#
#loaded via a namespace (and not attached):
#[1] nlme_3.1-162           bitops_1.0-7           fs_1.6.2               usethis_2.1.6         
#[5] devtools_2.4.5         RColorBrewer_1.1-3     GenomeInfoDb_1.34.9    tools_4.2.2           
#[9] profvis_0.3.8          utf8_1.2.3             R6_2.5.1               DBI_1.1.3             
#[13] BiocGenerics_0.44.0    mgcv_1.8-42            colorspace_2.1-0       rhdf5filters_1.10.0   
#[17] ade4_1.7-22            withr_2.5.0            urlchecker_1.0.1       tidyselect_1.2.0      
#[21] prettyunits_1.1.1      processx_3.8.1         compiler_4.2.2         microbiome_1.19.1     
#[25] cli_3.6.1              Biobase_2.58.0         labeling_0.4.2         scales_1.2.1          
#[29] callr_3.7.3            stringr_1.5.0          digest_0.6.31          XVector_0.38.0        
#[33] pkgconfig_2.0.3        htmltools_0.5.5        sessioninfo_1.2.2      fastmap_1.1.1         
#[37] htmlwidgets_1.6.2      rlang_1.1.1            rstudioapi_0.14        shiny_1.7.4           
#[41] farver_2.1.1           generics_0.1.3         jsonlite_1.8.4         dplyr_1.1.2           
#[45] RCurl_1.98-1.12        magrittr_2.0.3         GenomeInfoDbData_1.2.9 biomformat_1.26.0     
#[49] Matrix_1.5-4           Rcpp_1.0.10            munsell_0.5.0          S4Vectors_0.36.1      
#[53] Rhdf5lib_1.20.0        fansi_1.0.4            ape_5.7-1              lifecycle_1.0.3       
#[57] stringi_1.7.12         MASS_7.3-60            zlibbioc_1.44.0        Rtsne_0.16            
#[61] rhdf5_2.42.0           pkgbuild_1.4.0         grid_4.2.2             parallel_4.2.2        
#[65] promises_1.2.0.1       crayon_1.5.2           miniUI_0.1.1.1         Biostrings_2.66.0     
#[69] splines_4.2.2          multtest_2.54.0        ps_1.7.5               pillar_1.9.0          
#[73] igraph_1.4.2           reshape2_1.4.4         codetools_0.2-19       stats4_4.2.2          
#[77] pkgload_1.3.2          glue_1.6.2             data.table_1.14.8      remotes_2.4.2         
#[81] vctrs_0.6.2            httpuv_1.6.10          foreach_1.5.2          tidyr_1.3.0           
#[85] gtable_0.3.3           purrr_1.0.1            cachem_1.0.8           mime_0.12             
#[89] xtable_1.8-4           later_1.3.1            survival_3.5-5         tibble_3.2.1          
#[93] iterators_1.0.14       memoise_2.0.1          IRanges_2.32.0         writexl_1.4.2         
#[97] cluster_2.1.4          ellipsis_0.3.2 


#### Continued in Prism - Figure3.prism
######### END #########