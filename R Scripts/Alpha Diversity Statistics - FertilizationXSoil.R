############################################################
R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################

https://www.cfholbert.com/blog/nonparametric_two_way_anova/
https://bookdown.org/dereksonderegger/571/6-two-way-anova.html

#R 4.2.2
#setwd()



#### 1. Load data
div.df.1 <- readxl::read_xlsx('Figure1/Data/alpha_diversity.xlsx')
head(div.df.1)
## A tibble: 6 × 7
#  Soil        Fertilization  Block observed.x diversity_shannon evenness_simpson Phylogenetic_Diversity
#  <chr>       <chr>          <chr>      <dbl>             <dbl>            <dbl>                  <dbl>
#1 Rhizoplane  Fertilized     2            218              4.28            0.122                   18.3
#2 Rhizoplane  Fertilized     3            131              3.48            0.109                   10.8
#3 Rhizosphere Fertilized     2            197              4.91            0.449                   22.5
#4 Rhizosphere Fertilized     3            251              4.71            0.192                   23.6
#5 Rhizoplane  Non-fertilized 1            309              5.24            0.297                   30.4
#6 Rhizoplane  Non-fertilized 2            127              4.50            0.535                   15.6



#### 2. Define factors

# Make fertilization variable ordered factor
div.df.1$Fertilization <- factor(
  div.df.1$Fertilization, ordered = TRUE,
  levels = c('Non-fertilized', 'Fertilized')
)

# Make soil variable ordered factor
div.df.1$Soil <- factor(
  div.df.1$Soil, ordered = TRUE,
  levels = c('Unplanted', 'Rhizosphere', 'Rhizoplane')
)


str(div.df.1)
#tibble [242 × 7] (S3: tbl_df/tbl/data.frame)
#$ Soil                  : Ord.factor w/ 3 levels "Unplanted"<"Rhizosphere"<..: 3 3 2 2 3 3 3 2 2 2 ...
#$ Fertilization         : Ord.factor w/ 2 levels "Non-fertilized"<..: 2 2 2 2 1 1 1 1 1 1 ...
#$ Block                 : chr [1:242] "2" "3" "2" "3" ...
#$ observed.x            : num [1:242] 218 131 197 251 309 127 423 252 260 233 ...
#$ diversity_shannon     : num [1:242] 4.28 3.48 4.91 4.71 5.24 ...
#$ evenness_simpson      : num [1:242] 0.122 0.109 0.449 0.192 0.297 ...
#$ Phylogenetic_Diversity: num [1:242] 18.3 10.8 22.5 23.6 30.4 ...




#### 3. Type 1 or Type 3 least sum of squares in two-factor anova model? Balanced design (equal n in groups) vs unbalanced design (unequal n in groups)

# Create frequency table
with(div.df.1, table(Soil, Fertilization))
#Fertilization
#Soil          Non-fertilized Fertilized
#Unplanted                3          4
#Rhizosphere             64         63
#Rhizoplane              54         54

### Sample numbers are uneven so a type iii two-way ANOVA will be employed


############################## OBSERVED SPECIES ################

#### 4. Visualise data

# Compute summary statistics
library(dplyr)
div.df.1 %>%
  group_by(Soil, Fertilization) %>%
  summarise(
    count = n(),
    mean = round(mean(observed.x, na.rm = TRUE), 2),
    median = round(median(observed.x, na.rm = TRUE), 2),
    sd = round(sd(observed.x, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Soil'. You can override using the `.groups` argument.
## A tibble: 6 × 7
#  Soil        Fertilization  count  mean median    sd    cv
#  <ord>       <ord>          <int> <dbl>  <dbl> <dbl> <dbl>
#1 Unplanted   Non-fertilized     3  289.   279   39.5  0.14
#2 Unplanted   Fertilized         4  233.   240.  40.0  0.17
#3 Rhizosphere Non-fertilized    64  262.   263   67.1  0.26
#4 Rhizosphere Fertilized        63  206.   207   39.5  0.19
#5 Rhizoplane  Non-fertilized    54  343.   332  119.   0.35
#6 Rhizoplane  Fertilized        54  184.   155   88.0  0.48




# Box plots and interaction plots

library(ggpubr)
library(cowplot)


p1 <- ggboxplot(
  div.df.1, x = 'Soil', y = 'observed.x', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p2 <- ggline(
  div.df.1, x = 'Soil', y = 'observed.x', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
  )

plot_grid(p1, p2, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  observed.x ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
library(car)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: observed.x
#                    Sum Sq  Df  F value    Pr(>F)    
#(Intercept)        3530881   1 536.0848 < 2.2e-16 ***
#Soil                 50344   2   3.8218   0.02326 *  
#Fertilization       113461   1  17.2265 4.635e-05 ***
#Soil:Fertilization  160569   2  12.1894 9.165e-06 ***
#Residuals          1554396 236                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(observed.x) ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(observed.x)
#                    Sum Sq  Df    F value    Pr(>F)    
#(Intercept)        1653.96   1 16012.7906 < 2.2e-16 ***
#Soil                  0.12   2     0.5702    0.5662    
#Fertilization         1.85   1    17.8804 3.366e-05 ***
#Soil:Fertilization    2.74   2    13.2562 3.499e-06 ***
#Residuals            24.38 236                         
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(observed.x) ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(observed.x)
#                   Sum Sq  Df  F value    Pr(>F)    
#(Intercept)        937119   1 289.2519 < 2.2e-16 ***
#Soil                 6287   2   0.9703  0.380491    
#Fertilization       69538   1  21.4636 5.964e-06 ***
#Soil:Fertilization  47920   2   7.3955  0.000767 ***
#Residuals          764593 236                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#         
#247.9752 
#
# Soil 
#    Unplanted Rhizosphere Rhizoplane
#          257       234.5      263.3
#rep         7       127.0      108.0
#
# Fertilization 
#    Non-fertilized Fertilized
#               299      196.9
#rep            121      121.0
#
# Soil:Fertilization 
#             Fertilization
#Soil          Non-fertilized Fertilized
#  Unplanted   289.3          232.8     
#  rep           3.0            4.0     
#  Rhizosphere 262.2          206.3     
#  rep          64.0           63.0     
#  Rhizoplane  343.0          183.5     
#  rep          54.0           54.0  


### Checking normality assumption and homogeneity of variances on original, log-transformed, and rank-transformed data
library(ggfortify)
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
leveneTest(observed.x ~ Soil * Fertilization, data = div.df.1)
leveneTest(log(observed.x) ~ Soil * Fertilization, data = div.df.1)
leveneTest(rank(observed.x) ~ Soil * Fertilization, data = div.df.1)





# Compute estimated marginal means for factor combinations - pairwise comparisons
library(emmeans)
emmeans(aov.rnk, pairwise ~ Fertilization | Soil)
emmeans(aov.rnk, pairwise ~ Soil | Fertilization)



#### Shannon ####
#### 4. Visualise data

# Compute summary statistics
div.df.1 %>%
  group_by(Soil, Fertilization) %>%
  summarise(
    count = n(),
    mean = round(mean(diversity_shannon, na.rm = TRUE), 2),
    median = round(median(diversity_shannon, na.rm = TRUE), 2),
    sd = round(sd(diversity_shannon, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Soil'. You can override using the `.groups` argument.
## A tibble: 6 × 7
#Soil        Fertilization  count  mean median    sd    cv
#<ord>       <ord>          <int> <dbl>  <dbl> <dbl> <dbl>
#1 Unplanted   Non-fertilized     3  5.34   5.27  0.16  0.03
#2 Unplanted   Fertilized         4  5.13   5.16  0.17  0.03
#3 Rhizosphere Non-fertilized    64  5.25   5.31  0.31  0.06
#4 Rhizosphere Fertilized        63  4.69   4.74  0.41  0.09
#5 Rhizoplane  Non-fertilized    54  5.18   5.27  0.59  0.11
#6 Rhizoplane  Fertilized        54  3.68   3.65  0.74  0.2




# Box plots and interaction plots

p3 <- ggboxplot(
  div.df.1, x = 'Soil', y = 'diversity_shannon', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p4 <- ggline(
  div.df.1, x = 'Soil', y = 'diversity_shannon', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF'))

plot_grid(p3, p4, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  diversity_shannon ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: diversity_shannon
#Sum Sq  Df  F value    Pr(>F)    
#(Intercept)        1314.87   1 4871.744 < 2.2e-16 ***
#Soil                 18.82   2   34.862 5.429e-14 ***
#Fertilization         7.96   1   29.498 1.390e-07 ***
#Soil:Fertilization   13.63   2   25.251 1.155e-10 ***
#Residuals            63.70 236                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(diversity_shannon) ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(diversity_shannon)
#                    Sum Sq  Df  F value    Pr(>F)    
#(Intercept)        136.373   1 8188.035 < 2.2e-16 ***
#Soil                 1.222   2   36.698 1.327e-14 ***
#Fertilization        0.401   1   24.087 1.718e-06 ***
#Soil:Fertilization   0.878   2   26.351 4.682e-11 ***
#Residuals            3.931 236                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(diversity_shannon) ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(diversity_shannon)
#                   Sum Sq  Df  F value    Pr(>F)    
#(Intercept)        995010   1 454.4806 < 2.2e-16 ***
#Soil                72007   2  16.4449 2.060e-07 ***
#Fertilization       84239   1  38.4771 2.461e-09 ***
#Soil:Fertilization  41902   2   9.5696 0.0001009 ***
#Residuals          516683 236                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#         
#4.738247 
#
# Soil 
#    Unplanted Rhizosphere Rhizoplane
#        5.222       4.973       4.43
#rep     7.000     127.000     108.00
#
# Fertilization 
#    Non-fertilized Fertilized
#             5.224      4.252
#rep        121.000    121.000
#
# Soil:Fertilization 
#             Fertilization
#Soil          Non-fertilized Fertilized
#  Unplanted    5.34           5.13     
#  rep          3.00           4.00     
#  Rhizosphere  5.25           4.69     
#  rep         64.00          63.00     
#  Rhizoplane   5.18           3.68     
#  rep         54.00          54.00 


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
leveneTest(diversity_shannon ~ Soil * Fertilization, data = div.df.1)
leveneTest(log(diversity_shannon) ~ Soil * Fertilization, data = div.df.1)
leveneTest(rank(diversity_shannon) ~ Soil * Fertilization, data = div.df.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Soil)
emmeans(aov.rnk, pairwise ~ Soil | Fertilization)





#### Simpson ####
#### 4. Visualise data

# Compute summary statistics
div.df.1 %>%
  group_by(Soil, Fertilization) %>%
  summarise(
    count = n(),
    mean = round(mean(evenness_simpson, na.rm = TRUE), 2),
    median = round(median(evenness_simpson, na.rm = TRUE), 2),
    sd = round(sd(evenness_simpson, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Soil'. You can override using the `.groups` argument.
## A tibble: 6 × 7
#Soil        Fertilization  count  mean median    sd    cv
#<ord>       <ord>          <int> <dbl>  <dbl> <dbl> <dbl>
#1 Unplanted   Non-fertilized     3  0.55   0.56  0.05  0.09
#2 Unplanted   Fertilized         4  0.58   0.58  0.02  0.03
#3 Rhizosphere Non-fertilized    64  0.6    0.62  0.07  0.12
#4 Rhizosphere Fertilized        63  0.3    0.29  0.16  0.53
#5 Rhizoplane  Non-fertilized    54  0.33   0.34  0.14  0.42
#6 Rhizoplane  Fertilized        54  0.11   0.08  0.11  1   



# Box plots and interaction plots

p5 <- ggboxplot(
  div.df.1, x = 'Soil', y = 'evenness_simpson', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p6 <- ggline(
  div.df.1, x = 'Soil', y = 'evenness_simpson', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
)

plot_grid(p5, p6, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  evenness_simpson ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: evenness_simpson
#                   Sum Sq  Df F value    Pr(>F)    
#(Intercept)        9.3205   1 601.151 < 2.2e-16 ***
#Soil               3.5277   2 113.764 < 2.2e-16 ***
#Fertilization      0.3687   1  23.778 1.987e-06 ***
#Soil:Fertilization 0.2477   2   7.989 0.0004394 ***
#Residuals          3.6590 236                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(evenness_simpson) ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(evenness_simpson)
#                   Sum Sq  Df  F value    Pr(>F)    
#(Intercept)        72.359   1 201.0207 < 2.2e-16 ***
#Soil               57.954   2  80.5013 < 2.2e-16 ***
#Fertilization       6.413   1  17.8166 3.472e-05 ***
#Soil:Fertilization  4.118   2   5.7205  0.003749 ** 
#Residuals          84.950 236                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(evenness_simpson) ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(evenness_simpson)
#                    Sum Sq  Df  F value    Pr(>F)    
#(Intercept)        1055197   1 669.3092 < 2.2e-16 ***
#Soil                375250   2 119.0101 < 2.2e-16 ***
#Fertilization        37339   1  23.6839 2.077e-06 ***
#Soil:Fertilization   27408   2   8.6925 0.0002278 ***
#Residuals           372065 236                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#          
#0.3514655 
#
# Soil 
#    Unplanted Rhizosphere Rhizoplane
#       0.5647      0.4535     0.2177
#rep    7.0000    127.0000   108.0000
#
# Fertilization 
#    Non-fertilized Fertilized
#            0.4792     0.2237
#rep       121.0000   121.0000
#
# Soil:Fertilization 
#             Fertilization
#Soil          Non-fertilized Fertilized
#  Unplanted    0.55           0.58     
#  rep          3.00           4.00     
#  Rhizosphere  0.60           0.30     
#  rep         64.00          63.00     
#  Rhizoplane   0.33           0.11     
#  rep         54.00          54.00 
  


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
leveneTest(evenness_simpson ~ Soil * Fertilization, data = div.df.1)
leveneTest(log(evenness_simpson) ~ Soil * Fertilization, data = div.df.1)
leveneTest(rank(evenness_simpson) ~ Soil * Fertilization, data = div.df.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Soil)
emmeans(aov.rnk, pairwise ~ Soil | Fertilization)





#### Faith's PD ####
#### 4. Visualise data

# Compute summary statistics
div.df.1 %>%
  group_by(Soil, Fertilization) %>%
  summarise(
    count = n(),
    mean = round(mean(Phylogenetic_Diversity, na.rm = TRUE), 2),
    median = round(median(Phylogenetic_Diversity, na.rm = TRUE), 2),
    sd = round(sd(Phylogenetic_Diversity, na.rm = TRUE), 2),
    cv = round(sd/mean, 2),
  ) %>%
  ungroup()
#`summarise()` has grouped output by 'Soil'. You can override using the `.groups` argument.
## A tibble: 6 × 7
#  Soil        Fertilization  count  mean median    sd    cv
#  <ord>       <ord>          <int> <dbl>  <dbl> <dbl> <dbl>
#1 Unplanted   Non-fertilized     3  28.5   27.1  3.15  0.11
#2 Unplanted   Fertilized         4  25.2   25.7  2.26  0.09
#3 Rhizosphere Non-fertilized    64  29.4   29.7  4.09  0.14
#4 Rhizosphere Fertilized        63  23.1   22.5  2.83  0.12
#5 Rhizoplane  Non-fertilized    54  28.7   28.9  6.04  0.21
#6 Rhizoplane  Fertilized        54  17.4   16    5.55  0.32



# Box plots and interaction plots

p7 <- ggboxplot(
  div.df.1, x = 'Soil', y = 'Phylogenetic_Diversity', fill = 'Fertilization',
  palette = c('#E64B35FF', '#3C5488FF')
)

p8 <- ggline(
  div.df.1, x = 'Soil', y = 'Phylogenetic_Diversity', color = 'Fertilization',
  add = c('mean_se', 'dotplot'), size = 1,
  palette = c('#E64B35FF', '#3C5488FF')
)

plot_grid(p7, p8, ncol = 1, align = 'v')



#### 5. Two-factor ANOVA type 111

#ANOVA procedure on original data
aov.org <- aov(
  Phylogenetic_Diversity ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.org, type = 'III')
#Anova Table (Type III tests)
#
#Response: Phylogenetic_Diversity
#                   Sum Sq  Df   F value    Pr(>F)    
#(Intercept)         35549   1 1628.4331 < 2.2e-16 ***
#Soil                  624   2   14.2823 1.396e-06 ***
#Fertilization         668   1   30.5919 8.428e-08 ***
#Soil:Fertilization    409   2    9.3579 0.0001227 ***
#Residuals            5152 236                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on log-transformed data
aov.log <- aov(
  log(Phylogenetic_Diversity) ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.log, type = 'III')
#Anova Table (Type III tests)
#
#Response: log(Phylogenetic_Diversity)
#                   Sum Sq  Df   F value    Pr(>F)    
#(Intercept)        566.44   1 13976.385 < 2.2e-16 ***
#Soil                 1.95   2    24.016 3.208e-10 ***
#Fertilization        1.18   1    29.230 1.572e-07 ***
#Soil:Fertilization   1.25   2    15.395 5.197e-07 ***
#Residuals            9.56 236                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA procedure on rank-transformed data
aov.rnk <- aov(
  rank(Phylogenetic_Diversity) ~ Soil * Fertilization, data = div.df.1,
  contrasts = list(
    Soil = 'contr.sum',
    Fertilization = 'contr.sum'
  )
)
Anova(aov.rnk, type = 'III')
#Anova Table (Type III tests)
#
#Response: rank(Phylogenetic_Diversity)
#                   Sum Sq  Df  F value    Pr(>F)    
#(Intercept)        901229   1 359.1097 < 2.2e-16 ***
#Soil                56169   2  11.1908 2.274e-05 ***
#Fertilization       80922   1  32.2448 3.977e-08 ***
#Soil:Fertilization  23870   2   4.7557  0.009445 ** 
#Residuals          592270 236                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model summary
model.tables(aov.org, type = "means", se = TRUE)
#Design is unbalanced - use se.contrast() for se's
#Tables of means
#Grand mean
#         
#24.81776 
#
# Soil 
#    Unplanted Rhizosphere Rhizoplane
#        26.63       26.24      23.03
#rep      7.00      127.00     108.00
#
# Fertilization 
#    Non-fertilized Fertilized
#             29.03      20.61
#rep         121.00     121.00
#
# Soil:Fertilization 
#             Fertilization
#Soil          Non-fertilized Fertilized
#  Unplanted   28.53          25.22     
#  rep          3.00           4.00     
#  Rhizosphere 29.36          23.07     
#  rep         64.00          63.00     
#  Rhizoplane  28.66          17.39     
#  rep         54.00          54.00 



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
leveneTest(evenness_simpson ~ Soil * Fertilization, data = div.df.1)
leveneTest(log(evenness_simpson) ~ Soil * Fertilization, data = div.df.1)
leveneTest(rank(evenness_simpson) ~ Soil * Fertilization, data = div.df.1)


# Compute estimated marginal means for factor combinations - pairwise comparisons
emmeans(aov.rnk, pairwise ~ Fertilization | Soil)
emmeans(aov.rnk, pairwise ~ Soil | Fertilization)




######################### END #########################