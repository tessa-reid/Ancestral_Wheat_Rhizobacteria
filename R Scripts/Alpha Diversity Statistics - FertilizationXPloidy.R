############################################################
#R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################

#https://www.cfholbert.com/blog/nonparametric_two_way_anova/
#https://bookdown.org/dereksonderegger/571/6-two-way-anova.html

#R 4.2.2
#setwd()


################### RHIZOSPHERE ###################

#### 1. Load data
div.df.rs.1 <- readxl::read_xlsx('Figure2/Data/alpha_diversity_rhizosphere.xlsx')
head(div.df.rs.1)
## A tibble: 6 × 8
#  Soil        Fertilization  Ploidy  Block observed.x diversity_shannon evenness_simpson Phylogenetic_Diversity
#  <chr>       <chr>          <chr>   <chr>      <dbl>             <dbl>            <dbl>                  <dbl>
#1 Rhizosphere Fertilized     Diploid 2            197              4.94            0.482                   23.4
#2 Rhizosphere Fertilized     Diploid 3            255              4.71            0.189                   25.2
#3 Rhizosphere Non-fertilized Diploid 1            251              5.28            0.639                   31.6
#4 Rhizosphere Non-fertilized Diploid 2            264              5.27            0.595                   29.9
#5 Rhizosphere Non-fertilized Diploid 3            233              5.22            0.665                   28.2
#6 Rhizosphere Fertilized     Diploid 1            184              4.45            0.154                   21.8


#### 2. Define factors

# Make fertilization variable ordered factor
div.df.rs.1$Fertilization <- factor(
  div.df.rs.1$Fertilization, ordered = TRUE,
  levels = c('Non-fertilized', 'Fertilized')
)

# Make ploidy variable ordered factor
div.df.rs.1$Ploidy <- factor(
  div.df.rs.1$Ploidy, ordered = TRUE,
  levels = c('Diploid', 'Tetraploid', 'Hexaploid')
)


str(div.df.rs.1)
#tibble [127 × 8] (S3: tbl_df/tbl/data.frame)
#$ Soil                  : chr [1:127] "Rhizosphere" "Rhizosphere" "Rhizosphere" "Rhizosphere" ...
#$ Fertilization         : Ord.factor w/ 2 levels "Non-fertilized"<..: 2 2 1 1 1 2 2 2 1 1 ...
#$ Ploidy                : Ord.factor w/ 3 levels "Diploid"<"Tetraploid"<..: 1 1 1 1 1 1 1 1 1 1 ...
#$ Block                 : chr [1:127] "2" "3" "1" "2" ...
#$ observed.x            : num [1:127] 197 255 251 264 233 184 176 179 356 236 ...
#$ diversity_shannon     : num [1:127] 4.94 4.71 5.28 5.27 5.22 ...
#$ evenness_simpson      : num [1:127] 0.482 0.189 0.639 0.595 0.665 ...
#$ Phylogenetic_Diversity: num [1:127] 23.4 25.2 31.6 29.9 28.2 ...




#### 3. Type 1 or Type 3 least sum of squares in two-factor anova model? Balanced design (equal n in groups) vs unbalanced design (unequal n in groups)


# Create frequency table
with(div.df.rs.1, table(Ploidy, Fertilization))
#Fertilization
#Ploidy       Non-fertilized Fertilized
#Diploid                14         11
#Tetraploid             13         13
#Hexaploid              37         39

### Sample numbers are uneven so a type iii two-way ANOVA will be employed


############################## OBSERVED SPECIES ################

#### 4. Visualise data

# Compute summary statistics
library(dplyr)
div.df.rs.1 %>%
  group_by(Fertilization, Ploidy) %>%
  summarise(
    count = n(),
    mean = round(mean(observed.x, na.rm = TRUE), 2),
    median = round(median(observed.x, na.rm = TRUE), 2),
    sd = round(sd(observed.x, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Fertilization'. You can override using the `.groups` argument.
## A tibble: 6 × 7
#  Fertilization  Ploidy     count  mean median    sd    cv
#  <ord>          <ord>      <int> <dbl>  <dbl> <dbl> <dbl>
#1 Non-fertilized Diploid       14  269.    263  52.5  0.2 
#2 Non-fertilized Tetraploid    13  245.    254  35.8  0.15
#3 Non-fertilized Hexaploid     37  267.    274  79.2  0.3 
#4 Fertilized     Diploid       11  204.    184  50.9  0.25
#5 Fertilized     Tetraploid    13  219.    220  40.1  0.18
#6 Fertilized     Hexaploid     39  205.    209  37.3  0.18



# Box plots and interaction plots

library(ggpubr)
library(cowplot)


p1 <- ggboxplot(
  div.df.rs.1, x = 'Ploidy', y = 'observed.x', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p2 <- ggline(
  div.df.rs.1, x = 'Ploidy', y = 'observed.x', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
)

plot_grid(p1, p2, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  observed.x ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
library(car)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: observed.x
#                      Sum Sq  Df   F value    Pr(>F)    
#(Intercept)          5380558   1 1735.5860 < 2.2e-16 ***
#Fertilization          63663   1   20.5356 1.384e-05 ***
#Ploidy                   267   2    0.0430    0.9579    
#Fertilization:Ploidy    6907   2    1.1139    0.3316    
#Residuals             375117 121                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(observed.x) ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(observed.x)
#                      Sum Sq  Df    F value    Pr(>F)    
#(Intercept)          2878.01   1 50838.0293 < 2.2e-16 ***
#Fertilization           1.09   1    19.2572 2.457e-05 ***
#Ploidy                  0.00   2     0.0402    0.9606    
#Fertilization:Ploidy    0.10   2     0.8473    0.4311    
#Residuals               6.85 121                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(observed.x) ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(observed.x)
#                     Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          399430   1 382.1032 < 2.2e-16 ***
#Fertilization         30997   1  29.6528 2.755e-07 ***
#Ploidy                   77   2   0.0369    0.9638    
#Fertilization:Ploidy   2234   2   1.0685    0.3468    
#Residuals            126487 121                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#         
#235.4724 
#
# Fertilization 
#    Non-fertilized Fertilized
#             262.9      207.6
#rep           64.0       63.0
#
# Ploidy 
#    Diploid Tetraploid Hexaploid
#        237      232.6       236
#rep      25       26.0        76
#
# Fertilization:Ploidy 
#                Ploidy
#Fertilization    Diploid Tetraploid Hexaploid
#  Non-fertilized 268.79  245.46     266.76   
#  rep             14.00   13.00      37.00   
#  Fertilized     203.64  219.23     204.90   
#  rep             11.00   13.00      39.00   


### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
library(ggfortify)
autoplot(aov.org, which=c(1,2)) + geom_point( aes(color=Fertilization:Ploidy))
autoplot(aov.log, which=c(1,2)) + geom_point( aes(color=Fertilization:Ploidy))
autoplot(aov.rnk, which=c(1,2)) + geom_point( aes(color=Fertilization:Ploidy))


#Normality test
res.org = aov.org$resid #extract residuals
shapiro.test(x = res.org)

res.log = aov.log$resid #extract residuals
shapiro.test(x = res.log)

res.rnk = aov.rnk$resid #extract residuals
shapiro.test(x = res.rnk)


#Equal variance
leveneTest(observed.x ~ Fertilization * Ploidy, data = div.df.rs.1)
leveneTest(log(observed.x) ~ Fertilization * Ploidy, data = div.df.rs.1)
leveneTest(rank(observed.x) ~ Fertilization * Ploidy, data = div.df.rs.1)





# Compute estimated marginal means for factor combinations - pairwise comparisons
library(emmeans)
emmeans(aov.rnk, pairwise ~ Fertilization | Ploidy)



#### Shannon ####
#### 4. Visualise data

# Compute summary statistics
div.df.rs.1 %>%
  group_by(Fertilization, Ploidy) %>%
  summarise(
    count = n(),
    mean = round(mean(diversity_shannon, na.rm = TRUE), 2),
    median = round(median(diversity_shannon, na.rm = TRUE), 2),
    sd = round(sd(diversity_shannon, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Fertilization'. You can override using the `.groups`
#argument.
## A tibble: 6 × 7
# Fertilization  Ploidy     count  mean median    sd    cv
# <ord>          <ord>      <int> <dbl>  <dbl> <dbl> <dbl>
#1 Non-fertilized Diploid       14  5.3    5.28  0.2   0.04
#2 Non-fertilized Tetraploid    13  5.23   5.29  0.18  0.03
#3 Non-fertilized Hexaploid     37  5.25   5.34  0.38  0.07
#4 Fertilized     Diploid       11  4.68   4.71  0.34  0.07
#5 Fertilized     Tetraploid    13  4.9    4.95  0.25  0.05
#6 Fertilized     Hexaploid     39  4.62   4.73  0.46  0.1 



# Box plots and interaction plots

p3 <- ggboxplot(
  div.df.rs.1, x = 'Ploidy', y = 'diversity_shannon', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p4 <- ggline(
  div.df.rs.1, x = 'Ploidy', y = 'diversity_shannon', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF'))

plot_grid(p3, p4, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  diversity_shannon ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: diversity_shannon
#                      Sum Sq  Df   F value    Pr(>F)    
#(Intercept)          2437.11   1 18614.399 < 2.2e-16 ***
#Fertilization           6.81   1    51.977 5.285e-11 ***
#Ploidy                  0.34   2     1.315    0.2723    
#Fertilization:Ploidy    0.44   2     1.694    0.1881    
#Residuals              15.84 121                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(diversity_shannon) ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(diversity_shannon)
#                      Sum Sq  Df    F value    Pr(>F)    
#(Intercept)          251.461   1 39242.7294 < 2.2e-16 ***
#Fertilization          0.283   1    44.2115 8.991e-10 ***
#Ploidy                 0.019   2     1.5050    0.2261    
#Fertilization:Ploidy   0.020   2     1.5948    0.2072    
#Residuals              0.775 121                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(diversity_shannon) ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(diversity_shannon)
#                     Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          411172   1 568.5778 < 2.2e-16 ***
#Fertilization         56401   1  77.9926 9.696e-15 ***
#Ploidy                 1035   2   0.7153    0.4911    
#Fertilization:Ploidy   2308   2   1.5959    0.2070    
#Residuals             87502 121                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#         
#4.974315 
#
# Fertilization 
#    Non-fertilized Fertilized
#             5.257      4.687
#rep         64.000     63.000
#
# Ploidy 
#    Diploid Tetraploid Hexaploid
#      4.999      5.066     4.935
#rep  25.000     26.000    76.000
#
# Fertilization:Ploidy 
#                Ploidy
#Fertilization    Diploid Tetraploid Hexaploid
#  Non-fertilized  5.30    5.23       5.25    
#  rep            14.00   13.00      37.00    
#  Fertilized      4.68    4.90       4.62    
#  rep            11.00   13.00      39.00    


### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
autoplot(aov.org, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.log, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.rnk, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))


#Normality test
res.org = aov.org$resid #extract residuals
shapiro.test(x = res.org)

res.log = aov.log$resid #extract residuals
shapiro.test(x = res.log)

res.rnk = aov.rnk$resid #extract residuals
shapiro.test(x = res.rnk)


#Equal variance
leveneTest(diversity_shannon ~ Fertilization * Ploidy, data = div.df.rs.1)
leveneTest(log(diversity_shannon) ~ Fertilization * Ploidy, data = div.df.rs.1)
leveneTest(rank(diversity_shannon) ~ Fertilization * Ploidy, data = div.df.rs.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Ploidy)





#### Simpson ####
#### 4. Visualise data

# Compute summary statistics
div.df.rs.1 %>%
  group_by(Fertilization, Ploidy) %>%
  summarise(
    count = n(),
    mean = round(mean(evenness_simpson, na.rm = TRUE), 2),
    median = round(median(evenness_simpson, na.rm = TRUE), 2),
    sd = round(sd(evenness_simpson, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Fertilization'. You can override using the `.groups`
#argument.
## A tibble: 6 × 7
# Fertilization  Ploidy     count  mean median    sd    cv
# <ord>          <ord>      <int> <dbl>  <dbl> <dbl> <dbl>
#1 Non-fertilized Diploid       14  0.61   0.6   0.05  0.08
#2 Non-fertilized Tetraploid    13  0.62   0.62  0.05  0.08
#3 Non-fertilized Hexaploid     37  0.59   0.62  0.09  0.15
#4 Fertilized     Diploid       11  0.29   0.23  0.17  0.59
#5 Fertilized     Tetraploid    13  0.39   0.42  0.2   0.51
#6 Fertilized     Hexaploid     39  0.28   0.3   0.15  0.54



# Box plots and interaction plots

p5 <- ggboxplot(
  div.df.rs.1, x = 'Ploidy', y = 'evenness_simpson', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p6 <- ggline(
  div.df.rs.1, x = 'Ploidy', y = 'evenness_simpson', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
)

plot_grid(p5, p6, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  evenness_simpson ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: evenness_simpson
#                      Sum Sq  Df   F value  Pr(>F)    
#(Intercept)          21.0037   1 1300.0078 < 2e-16 ***
#Fertilization         2.0527   1  127.0530 < 2e-16 ***
#Ploidy                0.0994   2    3.0768 0.04973 *  
#Fertilization:Ploidy  0.0338   2    1.0468 0.35423    
#Residuals             1.9549 121                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(evenness_simpson)
#                     Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          82.323   1 329.5439 < 2.2e-16 ***
#Fertilization        16.554   1  66.2655 4.016e-13 ***
#Ploidy                0.930   2   1.8605    0.1600    
#Fertilization:Ploidy  0.429   2   0.8592    0.4261    
#Residuals            30.227 121                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(evenness_simpson)
#                     Sum Sq  Df  F value  Pr(>F)    
#(Intercept)          426949   1 795.8357 < 2e-16 ***
#Fertilization         72764   1 135.6332 < 2e-16 ***
#Ploidy                 4666   2   4.3485 0.01501 *  
#Fertilization:Ploidy    505   2   0.4707 0.62570    
#Residuals             64914 121                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#          
#0.4540534 
#
# Fertilization 
#    Non-fertilized Fertilized
#            0.6033     0.3025
#rep        64.0000    63.0000
#
# Ploidy 
#    Diploid Tetraploid Hexaploid
#     0.4515     0.5081    0.4364
#rep 25.0000    26.0000   76.0000
#
# Fertilization:Ploidy 
#                Ploidy
#Fertilization    Diploid Tetraploid Hexaploid
#  Non-fertilized  0.61    0.62       0.59    
#  rep            14.00   13.00      37.00    
#  Fertilized      0.29    0.39       0.28    
#  rep            11.00   13.00      39.00    



### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
autoplot(aov.org, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.log, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.rnk, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))


#Normality test
res.org = aov.org$resid #extract residuals
shapiro.test(x = res.org)

res.log = aov.log$resid #extract residuals
shapiro.test(x = res.log)

res.rnk = aov.rnk$resid #extract residuals
shapiro.test(x = res.rnk)


#Equal variance
leveneTest(evenness_simpson ~ Fertilization * Ploidy, data = div.df.rs.1)
leveneTest(log(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rs.1)
leveneTest(rank(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rs.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Ploidy)
emmeans(aov.rnk, pairwise ~ Ploidy | Fertilization)





#### Faith's PD ####
#### 4. Visualise data

# Compute summary statistics
div.df.rs.1 %>%
  group_by(Fertilization, Ploidy) %>%
  summarise(
    count = n(),
    mean = round(mean(Phylogenetic_Diversity, na.rm = TRUE), 2),
    median = round(median(Phylogenetic_Diversity, na.rm = TRUE), 2),
    sd = round(sd(Phylogenetic_Diversity, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Fertilization'. You can override using the `.groups`
#argument.
## A tibble: 6 × 7
#  Fertilization  Ploidy     count  mean median    sd    cv
#  <ord>          <ord>      <int> <dbl>  <dbl> <dbl> <dbl>
#1 Non-fertilized Diploid       14  29.5   29.8  3.48  0.12
#2 Non-fertilized Tetraploid    13  27.7   28.0  2.04  0.07
#3 Non-fertilized Hexaploid     37  30.0   30.6  4.72  0.16
#4 Fertilized     Diploid       11  22.6   22.0  2.55  0.11
#5 Fertilized     Tetraploid    13  24.6   24.5  2.66  0.11
#6 Fertilized     Hexaploid     39  23.0   23.1  2.93  0.13


# Box plots and interaction plots

p7 <- ggboxplot(
  div.df.rs.1, x = 'Ploidy', y = 'Phylogenetic_Diversity', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p8 <- ggline(
  div.df.rs.1, x = 'Ploidy', y = 'Phylogenetic_Diversity', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
)

plot_grid(p7, p8, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  Phylogenetic_Diversity ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: Phylogenetic_Diversity
#                     Sum Sq  Df   F value    Pr(>F)    
#(Intercept)           67229   1 5483.4320 < 2.2e-16 ***
#Fertilization           785   1   64.0587 8.315e-13 ***
#Ploidy                    7   2    0.2654   0.76736    
#Fertilization:Ploidy     80   2    3.2740   0.04123 *  
#Residuals              1484 121                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(Phylogenetic_Diversity) ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(Phylogenetic_Diversity)
#                      Sum Sq  Df    F value    Pr(>F)    
#(Intercept)          1033.50   1 57028.4557 < 2.2e-16 ***
#Fertilization           1.14   1    62.7304 1.294e-12 ***
#Ploidy                  0.00   2     0.1212   0.88597    
#Fertilization:Ploidy    0.11   2     2.9059   0.05853 .  
#Residuals               2.19 121                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(Phylogenetic_Diversity) ~ Fertilization * Ploidy, data = div.df.rs.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(Phylogenetic_Diversity)
#                     Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          388175   1 547.0812 < 2.2e-16 ***
#Fertilization         54899   1  77.3734 1.174e-14 ***
#Ploidy                  374   2   0.2638    0.7686    
#Fertilization:Ploidy   4422   2   3.1163    0.0479 *  
#Residuals             85854 121                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#         
#26.39451 
#
# Fertilization 
#    Non-fertilized Fertilized
#             29.45      23.29
#rep          64.00      63.00
#
# Ploidy 
#    Diploid Tetraploid Hexaploid
#      26.11      26.17     26.56
#rep   25.00      26.00     76.00
#
# Fertilization:Ploidy 
#                Ploidy
#Fertilization    Diploid Tetraploid Hexaploid
#  Non-fertilized 29.52   27.67      30.05    
#  rep            14.00   13.00      37.00    
#  Fertilized     22.56   24.62      23.05    
#  rep            11.00   13.00      39.00    



### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
autoplot(aov.org, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.log, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.rnk, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))


#Normality test
res.org = aov.org$resid #extract residuals
shapiro.test(x = res.org)

res.log = aov.log$resid #extract residuals
shapiro.test(x = res.log)

res.rnk = aov.rnk$resid #extract residuals
shapiro.test(x = res.rnk)


#Equal variance
leveneTest(evenness_simpson ~ Fertilization * Ploidy, data = div.df.rs.1)
leveneTest(log(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rs.1)
leveneTest(rank(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rs.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Ploidy)
emmeans(aov.rnk, pairwise ~ Ploidy | Fertilization)




################### RHIZOPLANE ###################


#### 1. Load data
div.df.rp.1 <- readxl::read_xlsx('Figure2/Data/alpha_diversity_rhizoplane.xlsx')
head(div.df.rp.1)
## A tibble: 6 × 8
#  Soil       Fertilization  Ploidy  Block observed.x diversity_shannon evenness_simpson Phylogenetic_Diversity
#  <chr>      <chr>          <chr>   <chr>      <dbl>             <dbl>            <dbl>                  <dbl>
#1 Rhizoplane Fertilized     Diploid 2            221              4.29            0.121                   17.9
#2 Rhizoplane Fertilized     Diploid 3            131              3.49            0.110                   11.1
#3 Rhizoplane Non-fertilized Diploid 1            308              5.23            0.292                   30.0
#4 Rhizoplane Non-fertilized Diploid 2            129              4.54            0.541                   16.2
#5 Rhizoplane Non-fertilized Diploid 3            422              5.30            0.186                   32.4
#6 Rhizoplane Fertilized     Diploid 1            308              5.47            0.622                   31.1



#### 2. Define factors

# Make fertilization variable ordered factor
div.df.rp.1$Fertilization <- factor(
  div.df.rp.1$Fertilization, ordered = TRUE,
  levels = c('Non-fertilized', 'Fertilized')
)

# Make ploidy variable ordered factor
div.df.rp.1$Ploidy <- factor(
  div.df.rp.1$Ploidy, ordered = TRUE,
  levels = c('Diploid', 'Tetraploid', 'Hexaploid')
)


str(div.df.rp.1)
#tibble [108 × 8] (S3: tbl_df/tbl/data.frame)
#$ Soil                  : chr [1:108] "Rhizoplane" "Rhizoplane" "Rhizoplane" "Rhizoplane" ...
#$ Fertilization         : Ord.factor w/ 2 levels "Non-fertilized"<..: 2 2 1 1 1 2 2 1 1 2 ...
#$ Ploidy                : Ord.factor w/ 3 levels "Diploid"<"Tetraploid"<..: 1 1 1 1 1 1 1 1 1 3 ...
#$ Block                 : chr [1:108] "2" "3" "1" "2" ...
#$ observed.x            : num [1:108] 221 131 308 129 422 308 113 255 254 108 ...
#$ diversity_shannon     : num [1:108] 4.29 3.49 5.23 4.54 5.3 ...
#$ evenness_simpson      : num [1:108] 0.121 0.11 0.292 0.541 0.186 ...
#$ Phylogenetic_Diversity: num [1:108] 17.9 11.1 30 16.2 32.4 ...



#### 3. Type 1 or Type 3 least sum of squares in two-factor anova model? Balanced design (equal n in groups) vs unbalanced design (unequal n in groups)

# Create frequency table
with(div.df.rp.1, table(Ploidy, Fertilization))
#            Fertilization
#Ploidy       Non-fertilized Fertilized
#Diploid                10          9
#Tetraploid             13         13
#Hexaploid              31         32


### Sample numbers are uneven so a type iii two-way ANOVA will be employed


############################## OBSERVED SPECIES ################

#### 4. Visualise data

# Compute summary statistics
div.df.rp.1 %>%
  group_by(Fertilization, Ploidy) %>%
  summarise(
    count = n(),
    mean = round(mean(observed.x, na.rm = TRUE), 2),
    median = round(median(observed.x, na.rm = TRUE), 2),
    sd = round(sd(observed.x, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Fertilization'. You can override using the `.groups` argument.
## A tibble: 6 × 7
#  Fertilization  Ploidy     count  mean median    sd    cv
#  <ord>          <ord>      <int> <dbl>  <dbl> <dbl> <dbl>
#1 Non-fertilized Diploid       10  267.   254.  96.6  0.36
#2 Non-fertilized Tetraploid    13  356    352  140.   0.39
#3 Non-fertilized Hexaploid     31  359.   372  107.   0.3 
#4 Fertilized     Diploid        9  168.   148   68.6  0.41
#5 Fertilized     Tetraploid    13  231.   172  120.   0.52
#6 Fertilized     Hexaploid     32  165.   152.  65.5  0.4 


# Box plots and interaction plots

p9 <- ggboxplot(
  div.df.rp.1, x = 'Ploidy', y = 'observed.x', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p10 <- ggline(
  div.df.rp.1, x = 'Ploidy', y = 'observed.x', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
)

plot_grid(p9, p10, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  observed.x ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: observed.x
#                      Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          5577711   1 565.3199 < 2.2e-16 ***
#Fertilization         406973   1  41.2481 4.317e-09 ***
#Ploidy                 62548   2   3.1697   0.04618 *  
#Fertilization:Ploidy   42979   2   2.1780   0.11850    
#Residuals            1006380 102                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(observed.x) ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(observed.x)
#                      Sum Sq  Df    F value    Pr(>F)    
#(Intercept)          2475.45   1 16113.2205 < 2.2e-16 ***
#Fertilization           6.97   1    45.3827 9.801e-10 ***
#Ploidy                  0.79   2     2.5751   0.08109 .  
#Fertilization:Ploidy    0.71   2     2.3000   0.10544    
#Residuals              15.67 102                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(observed.x) ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(observed.x)
#                     Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          242969   1 431.0082 < 2.2e-16 ***
#Fertilization         25025   1  44.3929 1.392e-09 ***
#Ploidy                 2989   2   2.6511   0.07543 .  
#Fertilization:Ploidy   2585   2   2.2927   0.10618    
#Residuals             57500 102                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#       
#261.25 
#
# Fertilization 
#    Non-fertilized Fertilized
#               341      181.5
#rep             54       54.0
#
# Ploidy 
#    Diploid Tetraploid Hexaploid
#      216.2      293.3     261.6
#rep    19.0       26.0      63.0
#
# Fertilization:Ploidy 
#                Ploidy
#Fertilization    Diploid Tetraploid Hexaploid
#  Non-fertilized 267.2   356.0      358.5    
#  rep             10.0    13.0       31.0    
#  Fertilized     168.3   230.6      165.2    
#  rep              9.0    13.0       32.0 


### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
library(ggfortify)
autoplot(aov.org, which=c(1,2)) + geom_point( aes(color=Fertilization:Ploidy))
autoplot(aov.log, which=c(1,2)) + geom_point( aes(color=Fertilization:Ploidy))
autoplot(aov.rnk, which=c(1,2)) + geom_point( aes(color=Fertilization:Ploidy))


#Normality test
res.org = aov.org$resid #extract residuals
shapiro.test(x = res.org)

res.log = aov.log$resid #extract residuals
shapiro.test(x = res.log)

res.rnk = aov.rnk$resid #extract residuals
shapiro.test(x = res.rnk)


#Equal variance
leveneTest(observed.x ~ Fertilization * Ploidy, data = div.df.rp.1)
leveneTest(log(observed.x) ~ Fertilization * Ploidy, data = div.df.rp.1)
leveneTest(rank(observed.x) ~ Fertilization * Ploidy, data = div.df.rp.1)



# Compute estimated marginal means for factor combinations - pairwise comparisons
library(emmeans)
emmeans(aov.rnk, pairwise ~ Fertilization | Ploidy)



#### Shannon ####
#### 4. Visualise data

# Compute summary statistics
div.df.rp.1 %>%
  group_by(Fertilization, Ploidy) %>%
  summarise(
    count = n(),
    mean = round(mean(diversity_shannon, na.rm = TRUE), 2),
    median = round(median(diversity_shannon, na.rm = TRUE), 2),
    sd = round(sd(diversity_shannon, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Fertilization'. You can override using the `.groups`
#argument.
## A tibble: 6 × 7
# Fertilization  Ploidy     count  mean median    sd    cv
#  <ord>          <ord>      <int> <dbl>  <dbl> <dbl> <dbl>
#1 Non-fertilized Diploid       10  5.02   5.14  0.38  0.08
#2 Non-fertilized Tetraploid    13  5.18   5.25  0.67  0.13
#3 Non-fertilized Hexaploid     31  5.22   5.46  0.63  0.12
#4 Fertilized     Diploid        9  3.65   3.81  1.09  0.3 
#5 Fertilized     Tetraploid    13  3.96   3.94  0.76  0.19
#6 Fertilized     Hexaploid     32  3.57   3.53  0.63  0.18



# Box plots and interaction plots

p11 <- ggboxplot(
  div.df.rp.1, x = 'Ploidy', y = 'diversity_shannon', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p12 <- ggline(
  div.df.rp.1, x = 'Ploidy', y = 'diversity_shannon', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF'))

plot_grid(p11, p12, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  diversity_shannon ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: diversity_shannon
#                      Sum Sq  Df   F value    Pr(>F)    
#(Intercept)          1650.87   1 3563.0823 < 2.2e-16 ***
#Fertilization          41.99   1   90.6314 9.402e-16 ***
#Ploidy                  0.73   2    0.7891    0.4570    
#Fertilization:Ploidy    0.97   2    1.0442    0.3557    
#Residuals              47.26 102                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(diversity_shannon) ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(diversity_shannon)
#                      Sum Sq  Df   F value    Pr(>F)    
#(Intercept)          178.996   1 5692.0926 < 2.2e-16 ***
#Fertilization          2.441   1   77.6310 3.456e-14 ***
#Ploidy                 0.057   2    0.9065    0.4072    
#Fertilization:Ploidy   0.054   2    0.8658    0.4238    
#Residuals              3.208 102                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(diversity_shannon) ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(diversity_shannon)
#                     Sum Sq  Df  F value Pr(>F)    
#(Intercept)          250703   1 598.8182 <2e-16 ***
#Fertilization         40438   1  96.5879 <2e-16 ***
#Ploidy                  673   2   0.8040 0.4503    
#Fertilization:Ploidy   1155   2   1.3796 0.2563    
#Residuals             42704 102                    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#         
#4.425459 
#
# Fertilization 
#    Non-fertilized Fertilized
#             5.175      3.676
#rep         54.000     54.000
#
# Ploidy 
#    Diploid Tetraploid Hexaploid
#       4.33      4.567     4.396
#rep   19.00     26.000    63.000
#
# Fertilization:Ploidy 
#                Ploidy
#Fertilization    Diploid Tetraploid Hexaploid
#  Non-fertilized  5.02    5.18       5.22    
#  rep            10.00   13.00      31.00    
#  Fertilized      3.65    3.96       3.57    
#  rep             9.00   13.00      32.00    


### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
autoplot(aov.org, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.log, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.rnk, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))


#Normality test
res.org = aov.org$resid #extract residuals
shapiro.test(x = res.org)

res.log = aov.log$resid #extract residuals
shapiro.test(x = res.log)

res.rnk = aov.rnk$resid #extract residuals
shapiro.test(x = res.rnk)


#Equal variance
leveneTest(diversity_shannon ~ Fertilization * Ploidy, data = div.df.rp.1)
leveneTest(log(diversity_shannon) ~ Fertilization * Ploidy, data = div.df.rp.1)
leveneTest(rank(diversity_shannon) ~ Fertilization * Ploidy, data = div.df.rp.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Ploidy)





#### Simpson ####
#### 4. Visualise data

# Compute summary statistics
div.df.rp.1 %>%
  group_by(Fertilization, Ploidy) %>%
  summarise(
    count = n(),
    mean = round(mean(evenness_simpson, na.rm = TRUE), 2),
    median = round(median(evenness_simpson, na.rm = TRUE), 2),
    sd = round(sd(evenness_simpson, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Fertilization'. You can override using the `.groups`
#argument.
## A tibble: 6 × 7
# Fertilization  Ploidy     count  mean median    sd    cv
#  <ord>          <ord>      <int> <dbl>  <dbl> <dbl> <dbl>
#1 Non-fertilized Diploid       10  0.37   0.41  0.16  0.43
#2 Non-fertilized Tetraploid    13  0.33   0.33  0.16  0.48
#3 Non-fertilized Hexaploid     31  0.32   0.32  0.13  0.41
#4 Fertilized     Diploid        9  0.17   0.11  0.19  1.12
#5 Fertilized     Tetraploid    13  0.09   0.06  0.09  1   
#6 Fertilized     Hexaploid     32  0.1    0.08  0.08  0.8 



# Box plots and interaction plots

p13 <- ggboxplot(
  div.df.rp.1, x = 'Ploidy', y = 'evenness_simpson', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p14 <- ggline(
  div.df.rp.1, x = 'Ploidy', y = 'evenness_simpson', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
)

plot_grid(p13, p14, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  evenness_simpson ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: evenness_simpson
#                     Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          4.3723   1 274.2902 < 2.2e-16 ***
#Fertilization        0.9886   1  62.0209 3.801e-12 ***
#Ploidy               0.0558   2   1.7491    0.1791    
#Fertilization:Ploidy 0.0040   2   0.1253    0.8823    
#Residuals            1.6259 102                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(evenness_simpson)
#                      Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          292.932   1 552.3186 < 2.2e-16 ***
#Fertilization         31.395   1  59.1941 9.354e-12 ***
#Ploidy                 1.130   2   1.0652    0.3485    
#Fertilization:Ploidy   0.045   2   0.0425    0.9584    
#Residuals             54.097 102                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(evenness_simpson)
#                     Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          263762   1 458.3784 < 2.2e-16 ***
#Fertilization         34306   1  59.6180 8.164e-12 ***
#Ploidy                 1664   2   1.4458    0.2403    
#Fertilization:Ploidy     54   2   0.0468    0.9543    
#Residuals             58693 102                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#          
#0.2184628 
#
# Fertilization 
#    Non-fertilized Fertilized
#            0.3283     0.1087
#rep        54.0000    54.0000
#
# Ploidy 
#    Diploid Tetraploid Hexaploid
#     0.2672     0.2089    0.2077
#rep 19.0000    26.0000   63.0000
#
# Fertilization:Ploidy 
#                Ploidy
#Fertilization    Diploid Tetraploid Hexaploid
#  Non-fertilized  0.37    0.33       0.32    
#  rep            10.00   13.00      31.00    
#  Fertilized      0.17    0.09       0.10    
#  rep             9.00   13.00      32.00    


### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
autoplot(aov.org, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.log, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.rnk, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))


#Normality test
res.org = aov.org$resid #extract residuals
shapiro.test(x = res.org)

res.log = aov.log$resid #extract residuals
shapiro.test(x = res.log)

res.rnk = aov.rnk$resid #extract residuals
shapiro.test(x = res.rnk)


#Equal variance
leveneTest(evenness_simpson ~ Fertilization * Ploidy, data = div.df.rp.1)
leveneTest(log(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rp.1)
leveneTest(rank(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rp.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Ploidy)





#### Faith's PD ####
#### 4. Visualise data

# Compute summary statistics
div.df.rp.1 %>%
  group_by(Fertilization, Ploidy) %>%
  summarise(
    count = n(),
    mean = round(mean(Phylogenetic_Diversity, na.rm = TRUE), 2),
    median = round(median(Phylogenetic_Diversity, na.rm = TRUE), 2),
    sd = round(sd(Phylogenetic_Diversity, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Fertilization'. You can override using the `.groups`
#argument.
## A tibble: 6 × 7
# Fertilization  Ploidy     count  mean median    sd    cv
#  <ord>          <ord>      <int> <dbl>  <dbl> <dbl> <dbl>
#1 Non-fertilized Diploid       10  24.7   25.0  5.33  0.22
#2 Non-fertilized Tetraploid    13  27.7   26.4  7.11  0.26
#3 Non-fertilized Hexaploid     31  30.1   30.1  5.19  0.17
#4 Fertilized     Diploid        9  16.0   13.8  6.17  0.39
#5 Fertilized     Tetraploid    13  20.5   18.0  6.29  0.31
#6 Fertilized     Hexaploid     32  16.3   15.6  4.1   0.25


# Box plots and interaction plots

p15 <- ggboxplot(
  div.df.rp.1, x = 'Ploidy', y = 'Phylogenetic_Diversity', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p16 <- ggline(
  div.df.rp.1, x = 'Ploidy', y = 'Phylogenetic_Diversity', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
)

plot_grid(p15, p16, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  Phylogenetic_Diversity ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: Phylogenetic_Diversity
#                     Sum Sq  Df   F value    Pr(>F)    
#(Intercept)           42713   1 1466.6165 < 2.2e-16 ***
#Fertilization          2050   1   70.3980 2.888e-13 ***
#Ploidy                  168   2    2.8794   0.06075 .  
#Fertilization:Ploidy    243   2    4.1644   0.01826 *  
#Residuals              2971 102                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(Phylogenetic_Diversity) ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(Phylogenetic_Diversity)
#                     Sum Sq  Df    F value    Pr(>F)    
#(Intercept)          783.60   1 13086.6738 < 2.2e-16 ***
#Fertilization          4.54   1    75.8680 5.751e-14 ***
#Ploidy                 0.37   2     3.0849   0.05003 .  
#Fertilization:Ploidy   0.48   2     3.9755   0.02175 *  
#Residuals              6.11 102                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(Phylogenetic_Diversity) ~ Fertilization * Ploidy, data = div.df.rp.1,
  contrasts = list(
    Fertilization = 'contr.sum',
    Ploidy = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(Phylogenetic_Diversity)
#                     Sum Sq  Df  F value    Pr(>F)    
#(Intercept)          238269   1 523.6810 < 2.2e-16 ***
#Fertilization         31822   1  69.9397 3.315e-13 ***
#Ploidy                 2806   2   3.0839   0.05007 .  
#Fertilization:Ploidy   3575   2   3.9287   0.02272 *  
#Residuals             46409 102                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#         
#22.88727 
#
# Fertilization 
#    Non-fertilized Fertilized
#             28.51      17.26
#rep          54.00      54.00
#
# Ploidy 
#    Diploid Tetraploid Hexaploid
#      20.27      24.11     23.17
#rep   19.00      26.00     63.00
#
# Fertilization:Ploidy 
#                Ploidy
#Fertilization    Diploid Tetraploid Hexaploid
#  Non-fertilized 24.66   27.70      30.10    
#  rep            10.00   13.00      31.00    
#  Fertilized     16.02   20.51      16.29    
#  rep             9.00   13.00      32.00 


### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
autoplot(aov.org, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.log, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))
autoplot(aov.rnk, which=c(1,2)) + geom_point( aes(color=Soil:Fertilization))


#Normality test
res.org = aov.org$resid #extract residuals
shapiro.test(x = res.org)

res.log = aov.log$resid #extract residuals
shapiro.test(x = res.log)

res.rnk = aov.rnk$resid #extract residuals
shapiro.test(x = res.rnk)


#Equal variance
leveneTest(evenness_simpson ~ Fertilization * Ploidy, data = div.df.rp.1)
leveneTest(log(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rp.1)
leveneTest(rank(evenness_simpson) ~ Fertilization * Ploidy, data = div.df.rp.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Ploidy)
emmeans(aov.rnk, pairwise ~ Ploidy | Fertilization)


######################### END #########################