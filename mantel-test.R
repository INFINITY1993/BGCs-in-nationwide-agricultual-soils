library(corrplot) 
library(vegan) 
library(ggcor)
library(ggplot2)
library(dplyr)


environment_factor <- read.csv("environment_factor.csv",header = T)
gcf_abundance <- read.csv("gcf_abundance.csv",header = T)
df <- merge(environment_factor,gcf_abundance,by="sample")

df_gcf <- df[ ,14:ncol(df)]
df_env <- df[,2:13]

mantel <- mantel_test(df_gcf, df_env,
                      spec.select = list(NRPS=1:2431,PKSI = 3838:4161,
                                         Hybrids = 4162:4249,
                                         PKSother = 4250:4881,RiPPs=4882:6698,Terpene=6703:8303,Others=2432:3837)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

quickcor(df_env, type = "upper",legend.position="right") +
  geom_circle2() +
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  scale_fill_gradient2(low="red",mid="white",high="navy")+
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
