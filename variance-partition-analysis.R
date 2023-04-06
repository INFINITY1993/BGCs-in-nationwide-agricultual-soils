library(vegan)
library(ggplot2)
library(ggsci)
library(dplyr)
library(SoDA)


#VPA分析就是采用的partial RDA原理，两者计算的结果是一致的


gcf_abundance <- read.csv("gcf_abundance.csv",header = T,row.names = 1)
environment_factor <- read.csv("environment.csv",header = T,row.names = 1)
taxa <- read.csv("16s_phylum_abundance.csv",header = T,row.names = 1)
taxa1 <- as.data.frame(t(taxa))
taxa_dominant <-taxa1[,1:15]

#按样本名排序
environment_factor_s <- environment_factor[order(rownames(environment_factor)),]
gcf_abundance_s <- gcf_abundance[order(rownames(gcf_abundance)),]
gcf_abundance_s.h <- decostand(gcf_abundance_s, "hellinger")
taxa_s <- taxa_dominant[order(rownames(taxa_dominant)),]

#环境因子数据标准化，均值为0，标准差为1
env.scaled<-data.frame(scale(environment_factor_s,center = F)) #数据标准化，均值为0，标准差为1

#经纬度转换pcnm
geo <- subset(environment_factor_s,select = c(longitude,latitude))
geo.xy <- geoXY(geo$longitude, geo$latitude, unit = 1)
pcnm1 <- pcnm(dist(geo.xy))
df_pcnm <- data.frame(pcnm1$vectors)
rownames(df_pcnm) <- rownames(environment_factor)


#在多因素的分析中，例如 RDA 分析和多元线性回归，
#任何增加新的变量，都会使得模型看起来拟合得更好，然而实际情况显然并非如此，为了解决过拟合问题调整R方是一个很常用的指标

#全模型逐渐消减变量
# Create the full model
rda_full <- rda(gcf_abundance_s.h~., data = env.scaled)
rda_back <- ordistep(rda_full, direction = 'backward',trace = 0)

#从零模型逐渐增加变量以及双向选择变量
# Create the zero model
rda_null <- rda(gcf_abundance_s.h~1, data = env.scaled)
rda_frwd <- ordistep(rda_null, formula(rda_full), direction = 'forward',trace = 0)
rda_both <- ordistep(rda_null, formula(rda_full), direction = 'both',trace = 0)

rda_back
rda_frwd
rda_both

##根据以上三个运行结果，选择筛选后的解释变量pH, MC,NO3,TN,DOC, AK, MAT,MAP


##方差分解分析
df_clim <- env.scaled %>% select(MAT,MAP)
df_soil <- env.scaled %>% select(pH, MC,NO3,TN,DOC, AK)
df_geo <- cbind(df_clim,df_pcnm)
df_taxa <- taxa_s

vpt <- varpart(gcf_abundance_s.h, df_geo,df_taxa, df_soil,chisquare=FALSE)
vpt

#绘制图表
plot(
  vpt,
  bg = c("#6A6599FF","#DF8F44FF","#00A1D5FF"),
  alpha=150,
  id.size = 1.1,
  cex = 1.2,
  Xnames = c('Geographic+Climatic', 'Microbial', 'Edaphic')
)



#显著性检验

formula_soil <- formula(clfs_surface_s.h ~ as.matrix(df_soil) +Condition(as.matrix(df_geo))+Condition(as.matrix(df_taxa)))
formula_pcnm_clim <- formula(clfs_surface_s.h ~ as.matrix(df_geo) +Condition(as.matrix(df_soil))+Condition(as.matrix(df_taxa)))
formula_taxa <- formula(clfs_surface_s.h ~ as.matrix(df_taxa) +Condition(as.matrix(df_geo))+Condition(as.matrix(df_soil)))

anova(rda(formula_taxa))
anova(rda(formula_soil))
anova(rda(formula_pcnm_clim))
