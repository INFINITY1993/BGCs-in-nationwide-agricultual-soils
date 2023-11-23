library(lavaan)

#read file
sem_data <- read.csv("SEM_analysis_data.csv",header = T,row.names = 1)

#data normalization，mean value 0，standard deviation 1
sem_data.scaled<-data.frame(scale(sem_data,center = F)) 

#model construction
model <- '
# regressions

MiCo~pH+NO3+MAT+MC+AP
BGC.Co~pH+MAP+MiCo
BGC_D~MiCo
'
Fit <- lavaan::sem(model, data=sem_data.scaled)
summary(Fit, rsquare=T, standardized=T,fit.measures=TRUE)

modificationIndices(Fit,standardized=F)

fitMeasures(Fit,c("chisq","df","gfi","agfi","cfi","nfi","ifi","srmr","rmsea"))


#SEM plot
library(semPlot)

semPaths(Fit)
semPaths(Fit,"std",nCharNodes=0,layout="spring",edge.label.cex = 1,sizeMan =6,sizeMan2=6,minimum=0.1,exoCov=FALSE,residuals=FALSE)

