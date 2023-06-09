library(vegan)
library(ggplot2)
library(ggsci)
library(dplyr)
library(SoDA)


#VPA analysis uses the partial RDA principle, and the results of the two calculations are consistent

#Read file
gcf_abundance <- read.csv("gcf_abundance.csv",header = T,row.names = 1)
environment_factor <- read.csv("environment.csv",header = T,row.names = 1)
taxa <- read.csv("16s_phylum_abundance.csv",header = T,row.names = 1)
taxa1 <- as.data.frame(t(taxa))
taxa_dominant <-taxa1[,1:15]

#Sorted by sample name
environment_factor_s <- environment_factor[order(rownames(environment_factor)),]
gcf_abundance_s <- gcf_abundance[order(rownames(gcf_abundance)),]
gcf_abundance_s.h <- decostand(gcf_abundance_s, "hellinger")
taxa_s <- taxa_dominant[order(rownames(taxa_dominant)),]

#Data normalization of environmental variables，mean value 0，standard deviation 1
env.scaled<-data.frame(scale(environment_factor_s,center = F)) #数据标准化，均值为0，标准差为1

# Coordinates to pcnm
geo <- subset(environment_factor_s,select = c(longitude,latitude))
geo.xy <- geoXY(geo$longitude, geo$latitude, unit = 1)
pcnm1 <- pcnm(dist(geo.xy))
df_pcnm <- data.frame(pcnm1$vectors)
rownames(df_pcnm) <- rownames(environment_factor)


#In multi-factor analysis, such as RDA analysis and multiple linear regression, any addition of new variables will make the model seem to fit better, but the actual situation is obviously not the case. In order to solve the overfitting problem, it is a common index to adjust R square
# Create the full model
rda_full <- rda(gcf_abundance_s.h~., data = env.scaled)
rda_back <- ordistep(rda_full, direction = 'backward',trace = 0)

# Create the zero model
rda_null <- rda(gcf_abundance_s.h~1, data = env.scaled)
rda_frwd <- ordistep(rda_null, formula(rda_full), direction = 'forward',trace = 0)
rda_both <- ordistep(rda_null, formula(rda_full), direction = 'both',trace = 0)

rda_back
rda_frwd
rda_both

#Based on the above results, pH, MC,NO3,TN,DOC, AK, MAT, MAP were selected


##Variance partition analysis
df_clim <- env.scaled %>% select(MAT,MAP)
df_soil <- env.scaled %>% select(pH, MC,NO3,TN,DOC, AK)
df_geo <- cbind(df_clim,df_pcnm)
df_taxa <- taxa_s

vpt <- varpart(gcf_abundance_s.h, df_geo,df_taxa, df_soil,chisquare=FALSE)
vpt

#Figure plot
plot(
  vpt,
  bg = c("#6A6599FF","#DF8F44FF","#00A1D5FF"),
  alpha=150,
  id.size = 1.1,
  cex = 1.2,
  Xnames = c('Geographic+Climatic', 'Microbial', 'Edaphic')
)



#Significance test
formula_soil <- formula(clfs_surface_s.h ~ as.matrix(df_soil) +Condition(as.matrix(df_geo))+Condition(as.matrix(df_taxa)))
formula_pcnm_clim <- formula(clfs_surface_s.h ~ as.matrix(df_geo) +Condition(as.matrix(df_soil))+Condition(as.matrix(df_taxa)))
formula_taxa <- formula(clfs_surface_s.h ~ as.matrix(df_taxa) +Condition(as.matrix(df_geo))+Condition(as.matrix(df_soil)))

anova(rda(formula_taxa))
anova(rda(formula_soil))
anova(rda(formula_pcnm_clim))
