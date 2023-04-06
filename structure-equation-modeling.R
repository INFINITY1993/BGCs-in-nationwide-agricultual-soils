library(lavaan)

sem_data <- read.csv("SEM_analysis_data.csv",header = T,row.names = 1)


sem_data.scaled<-data.frame(scale(sem_data,center = F)) #数据标准化，均值为0，标准差为1


model <- '
# regressions

PC1~pH+NO3+MC+MAT
BGC_richness~PC1
BGC_abundance~MAT+MAP+pH+PC1
'
Fit <- lavaan::sem(model, data=sem_data.scaled)
summary(Fit, rsquare=T, standardized=T,fit.measures=TRUE)
residuals(Fit, type="cor")
modificationIndices(Fit,standardized=F)

fitmeasures(Fit,fit.measures="all")


#SEM plot
library(semPlot)

semPaths(Fit)
semPaths(Fit,"std",nCharNodes=0,layout="spring",edge.label.cex = 1,sizeMan =6,sizeMan2=6,minimum=0.1,exoCov=FALSE,residuals=FALSE)

