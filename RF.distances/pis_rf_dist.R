#This script is written by Lukas J Musher. Please cite:
#Musher, L.J., et al. ... TINAMOU PAPER UPDATE

library(ggplot2)
library(gplots)
library(phytools)
library(ape)
library(phangorn)
library(ips)
library(ggpubr)

#must have alignments and trees--if there are 20 genes there must be 20 genes--and must be in same directory 
#and named appropriately so you know you have the right tree with each gene (they must be in the exact same order)
#ideally you will have alignments and trees named exactly the same e.g. uce-10.nexus uce-10.tre 

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/CDS_complete_nexus/")

# files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
# trees <- list.files(path="./",pattern="*.tre",full.names=F,recursive=FALSE)
# 
# spec_tree <- root(read.tree(file = "../../Trees/CDS.complete.ASTRAL.PP.tre"), "Rhea_penn", resolve.root = T)
# plot(spec_tree)
# 
# pis<-function(x){
#   len<-length(x[[1]])
#   x<-read.nexus.data(x)
#   site<-c()
#   counts=0
#   for(j in 1:length(x[[1]])){
#     for(i in 1:length(x)){
#       site[i]<-(x[[i]][j])
#     }
#     if(!("?" %in% unique(site))&!("-" %in% unique(site))){
#       if (length(unique(site))>=2){
#         counts=counts+1
#       }
#     }
#     else if(!("?" %in% unique(site))){
#       if (length(unique(site))>=3){
#         counts=counts+1
#       }
#     }
#     else if(!("-" %in% unique(site))){
#       if (length(unique(site))>=3){
#         counts=counts+1
#       }
#     }
#     else if (length(unique(site))>4){
#       #print(unique(site))
#       counts=counts+1
#     }
#     else counts=counts
#   }
#   return(counts)
# }
# 
# #how many bp sites in the alignment
# site_count<-function(x){
#   x<-read.nexus.data(x)
#   len<-length(x[[1]])
#   return(len)
# }
# 
# aln<-c()
# pars.site<-c()
# tot.sites<-c()
# rf.distance<-c()
# pars.site.prop<-c()
# 
# for (i in 1:length(files)){
#   print(c(i,files[i]))
#   aln[i]<-files[i]
#   pars.site[i]<-as.numeric(pis(files[i]))
#   tot.sites[i]<-as.numeric(site_count(files[i]))
#   pars.site.prop[i]<-as.numeric(pars.site[i]/tot.sites[i])
#   skip_to_next <- FALSE
#   
#   # Note that print(b) fails since b doesn't exist
#   
#   tryCatch(rf.distance[i]<-as.numeric(RF.dist(read.tree(trees[i]),spec_tree), error = function(e) { skip_to_next <<- TRUE}))
#   
#   if(skip_to_next) { next }
#   
# }
# 
# df.CDS<-as.data.frame(cbind(aln,pars.site,tot.sites,pars.site.prop,rf.distance))
# 
# write.csv(df.CDS,"../CDS.csv")

df.CDS<-read.csv("../CDS.csv")

df.CDS$pars.site<-as.numeric(df.CDS$pars.site)
df.CDS$tot.sites<-as.numeric(df.CDS$tot.sites)
df.CDS$pars.site.prop<-as.numeric(df.CDS$pars.site.prop)
df.CDS$rf.distance<-as.numeric(df.CDS$rf.distance)

mean(df.CDS$pars.site) #922.5146
sd(df.CDS$pars.site) #690.084

mean(df.CDS$pars.site.prop) #0.3717782
sd(df.CDS$pars.site.prop) #0.1033898

mean(df.CDS$rf.distance) #58.52413
sd(df.CDS$rf.distance) #21.11341

model<-lm(rf.distance~pars.site.prop, data = df.CDS)

CDS.prop.plot<-ggplot(df.CDS,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("CDS", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot

summary(model)
# Residual standard error: 20.96 on 2505 degrees of freedom
# Multiple R-squared:  0.01487,	Adjusted R-squared:  0.01448 
# F-statistic: 37.81 on 1 and 2505 DF,  p-value: 9.026e-10

model<-lm(rf.distance~pars.site, data = df.CDS)

CDS.pis.plot<-ggplot(df.CDS,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("CDS", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")

CDS.pis.plot

# Residual standard error: 20.68 on 2505 degrees of freedom
# Multiple R-squared:  0.04124,	Adjusted R-squared:  0.04086 
# F-statistic: 107.8 on 1 and 2505 DF,  p-value: < 2.2e-16

model<-lm(rf.distance~tot.sites, data = df.CDS)

CDS.tot.plot<-ggplot(df.CDS,aes(x=tot.sites,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("CDS", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Locus length (bp)") + ylab("RF Distance")

CDS.tot.plot

# Residual standard error: 20.36 on 2505 degrees of freedom
# Multiple R-squared:  0.07008,	Adjusted R-squared:  0.06971 
# F-statistic: 188.8 on 1 and 2505 DF,  p-value: < 2.2e-16

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/alignments_uce-1000flank-complete-nexus/")

# files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
# trees <- list.files(path="./",pattern="*.tre",full.names=F,recursive=FALSE)
# 
# spec_tree <- root(read.tree(file = "../../Trees/UCEs.1000flank.complete.ASTRAL.PP.tre"), "Rhea_penn_GCA_003342835", resolve.root = T)
# plot(spec_tree)
# 
# aln<-c()
# pars.site<-c()
# tot.sites<-c()
# rf.distance<-c()
# pars.site.prop<-c()
# 
# for (i in 1:length(files)){
#   print(c(i,files[i]))
#   aln[i]<-files[i]
#   pars.site[i]<-pis(files[i])
#   tot.sites[i]<-site_count(files[i])
#   pars.site.prop[i]<-pars.site[i]/tot.sites[i]
#   skip_to_next <- FALSE
#   
#   # Note that print(b) fails since b doesn't exist
#   
#   tryCatch(rf.distance[i]<-RF.dist(read.tree(trees[i]),spec_tree), error = function(e) { skip_to_next <<- TRUE})
#   
#   if(skip_to_next) { next }
#   
# }
# 
# df.uce1000<-as.data.frame(cbind(aln,pars.site,tot.sites,pars.site.prop,rf.distance))
# 
# df.uce1000$pars.site<-as.numeric(df.uce1000$pars.site)
# df.uce1000$tot.sites<-as.numeric(df.uce1000$tot.sites)
# df.uce1000$pars.site.prop<-as.numeric(df.uce1000$pars.site.prop)
# df.uce1000$rf.distance<-as.numeric(df.uce1000$rf.distance)
# 
# write.csv(df.uce1000,"../uce1000.csv")

df.uce1000<-read.csv("../uce1000.csv")
model<-lm(rf.distance~pars.site.prop, data = df.uce1000)

mean(df.uce1000$pars.site) #1035.713
sd(df.uce1000$pars.site) #250.5732

mean(df.uce1000$pars.site.prop) #0.4178852
sd(df.uce1000$pars.site.prop) #0.08261385

mean(df.uce1000$rf.distance) #29.86748
sd(df.uce1000$rf.distance) #8.468606

uce1000.prop.plot<-ggplot(df.uce1000,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

uce1000.prop.plot

# Residual standard error: 8.311 on 2790 degrees of freedom
# Multiple R-squared:  0.03733,	Adjusted R-squared:  0.03698 
# F-statistic: 108.2 on 1 and 2790 DF,  p-value: < 2.2e-16

model<-lm(rf.distance~pars.site, data=df.uce1000)

uce1000.pis.plot<-ggplot(df.uce1000,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")

uce1000.pis.plot

# Residual standard error: 8.425 on 2790 degrees of freedom
# Multiple R-squared:  0.01067,	Adjusted R-squared:  0.01031 
# F-statistic: 30.08 on 1 and 2790 DF,  p-value: 4.51e-08

model<-lm(rf.distance~tot.sites, data=df.uce1000)

uce1000.tot.plot<-ggplot(df.uce1000,aes(x=tot.sites,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Locus length (bp)") + ylab("RF Distance")

uce1000.tot.plot

# Residual standard error: 8.429 on 2790 degrees of freedom
# Multiple R-squared:  0.009766,	Adjusted R-squared:  0.009411 
# F-statistic: 27.52 on 1 and 2790 DF,  p-value: 1.675e-07

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/alignments_uce-300flank-complete-nexus/")

# files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
# trees <- list.files(path="./",pattern="*.tre",full.names=F,recursive=FALSE)
# 
# aln<-c()
# pars.site<-c()
# tot.sites<-c()
# rf.distance<-c()
# pars.site.prop<-c()
# 
# for (i in 1:length(files)){
#   print(c(i,files[i]))
#   aln[i]<-files[i]
#   pars.site[i]<-pis(files[i])
#   tot.sites[i]<-site_count(files[i])
#   pars.site.prop[i]<-pars.site[i]/tot.sites[i]
#   skip_to_next <- FALSE
#   
#   # Note that print(b) fails since b doesn't exist
#   
#   tryCatch(rf.distance[i]<-RF.dist(read.tree(trees[i]),spec_tree), error = function(e) { skip_to_next <<- TRUE})
#   
#   if(skip_to_next) { next }
#   
# }
# 
# df.uce300<-as.data.frame(cbind(aln,pars.site,tot.sites,pars.site.prop,rf.distance))
# 
# df.uce300$pars.site<-as.numeric(df.uce300$pars.site)
# df.uce300$tot.sites<-as.numeric(df.uce300$tot.sites)
# df.uce300$pars.site.prop<-as.numeric(df.uce300$pars.site.prop)
# df.uce300$rf.distance<-as.numeric(df.uce300$rf.distance)
# 
# write.csv(df.uce300,"../uce300.csv")

df.uce300<-read.csv("../uce300.csv")

mean(df.uce300$pars.site) #206.8736
sd(df.uce300$pars.site) #88.36566

mean(df.uce300$pars.site.prop) #0.2648982
sd(df.uce300$pars.site.prop) #0.1048726

mean(df.uce300$rf.distance) #59.16151
sd(df.uce300$rf.distance) #18.40835

model<-lm(rf.distance~pars.site.prop, data=df.uce300)

uce300.prop.plot<-ggplot(df.uce300,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

uce300.prop.plot

# Residual standard error: 12.26 on 2877 degrees of freedom
# Multiple R-squared:  0.5569,	Adjusted R-squared:  0.5568 
# F-statistic:  3616 on 1 and 2877 DF,  p-value: < 2.2e-16

model<-lm(rf.distance~pars.site, data=df.uce300)

uce300.pis.plot<-ggplot(df.uce300,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")

uce300.pis.plot

# Residual standard error: 12.88 on 2877 degrees of freedom
# Multiple R-squared:  0.5106,	Adjusted R-squared:  0.5104 
# F-statistic:  3001 on 1 and 2877 DF,  p-value: < 2.2e-16

model<-lm(rf.distance~tot.sites, data=df.uce300)

uce300.tot.plot<-ggplot(df.uce300,aes(x=tot.sites,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Locus length (bp)") + ylab("RF Distance")

uce300.tot.plot

# Residual standard error: 18.05 on 2877 degrees of freedom
# Multiple R-squared:  0.03857,	Adjusted R-squared:  0.03823 
# F-statistic: 115.4 on 1 and 2877 DF,  p-value: < 2.2e-16

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/alignments_uce-100flank-complete-nexus/")

# files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
# trees <- list.files(path="./",pattern="*.tre",full.names=F,recursive=FALSE)
# 
# aln<-c()
# pars.site<-c()
# tot.sites<-c()
# rf.distance<-c()
# pars.site.prop<-c()
# 
# for (i in 1:length(files)){
#   print(c(i,files[i]))
#   aln[i]<-files[i]
#   pars.site[i]<-pis(files[i])
#   tot.sites[i]<-site_count(files[i])
#   pars.site.prop[i]<-pars.site[i]/tot.sites[i]
#   skip_to_next <- FALSE
#   
#   # Note that print(b) fails since b doesn't exist
#   
#   tryCatch(rf.distance[i]<-RF.dist(read.tree(trees[i]),spec_tree), error = function(e) { skip_to_next <<- TRUE})
#   
#   if(skip_to_next) { next }
#   
# }
# 
# df.uce100<-as.data.frame(cbind(aln,pars.site,tot.sites,pars.site.prop,rf.distance))
# df.uce100<-na.omit(df.uce100)
# df.uce100$pars.site<-as.numeric(df.uce100$pars.site)
# df.uce100$tot.sites<-as.numeric(df.uce100$tot.sites)
# df.uce100$pars.site.prop<-as.numeric(df.uce100$pars.site.prop)
# df.uce100$rf.distance<-as.numeric(df.uce100$rf.distance)
# 
# write.csv(df.uce100,"../uce100.csv")

df.uce100<-read.csv("../uce100.csv")

mean(df.uce100$pars.site) #45.20298
sd(df.uce100$pars.site) #28.84223

mean(df.uce100$pars.site.prop) #0.1347156
sd(df.uce100$pars.site.prop) #0.0828647

mean(na.omit(df.uce100$rf.distance)) #106.5005
sd(na.omit(df.uce100$rf.distance)) #20.10551

model<-lm(rf.distance~pars.site.prop, data=df.uce100)

uce100.prop.plot<-ggplot(df.uce100,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

uce100.prop.plot

# Residual standard error: 19.73 on 2879 degrees of freedom
# (6 observations deleted due to missingness)
# Multiple R-squared:  0.03752,	Adjusted R-squared:  0.03718 
# F-statistic: 112.2 on 1 and 2879 DF,  p-value: < 2.2e-16

model<-lm(rf.distance~pars.site, data=df.uce100)

uce100.pis.plot<-ggplot(df.uce100,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")

uce100.pis.plot

# Residual standard error: 19.68 on 2879 degrees of freedom
# (6 observations deleted due to missingness)
# Multiple R-squared:  0.04178,	Adjusted R-squared:  0.04145 
# F-statistic: 125.5 on 1 and 2879 DF,  p-value: < 2.2e-16

model<-lm(rf.distance~tot.sites, data=df.uce100)

uce100.tot.plot<-ggplot(df.uce100,aes(x=tot.sites,y=rf.distance)) +
  geom_point(alpha=0.6) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(model)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(model))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Locus length (bp)") + ylab("RF Distance")

uce100.tot.plot

# Residual standard error: 19.68 on 2879 degrees of freedom
# (6 observations deleted due to missingness)
# Multiple R-squared:  0.04178,	Adjusted R-squared:  0.04145 
# F-statistic: 125.5 on 1 and 2879 DF,  p-value: < 2.2e-16

pdf(file = "../regressions_genetrees.pdf", width = 14, height = 7)
ggarrange(CDS.pis.plot,uce100.pis.plot,uce300.pis.plot,uce1000.pis.plot,
          CDS.prop.plot,uce100.prop.plot,uce300.prop.plot,uce1000.prop.plot,
          CDS.tot.plot,uce100.tot.plot,uce300.tot.plot,uce1000.tot.plot,
          ncol = 4, nrow = 3)
dev.off()

pdf(file = "../regressions_genetrees2.pdf", width = 14, height = 7)
ggarrange(CDS.pis.plot,uce100.pis.plot,uce300.pis.plot,uce1000.pis.plot,
          CDS.prop.plot,uce100.prop.plot,uce300.prop.plot,uce1000.prop.plot,
          ncol = 4, nrow = 2)
dev.off()

###########################
###Kruskall Wallace Test###
###########################

dataset<-c(rep("BUSCOs",length(df.CDS[,1])),rep("UCE100Flank",length(df.uce100[,1])),rep("UCE300Flank",length(df.uce300[,1])),rep("UCE1000Flank",length(df.uce1000[,1])))
PIS<-as.numeric(c(df.CDS$pars.site,df.uce100$pars.site,df.uce300$pars.site,df.uce1000$pars.site))
PIS.prop<-as.numeric(c(df.CDS$pars.site.prop,df.uce100$pars.site.prop,df.uce300$pars.site.prop,df.uce1000$pars.site.prop))
RF<-as.numeric(c(df.CDS$rf.distance,df.uce100$rf.distance,df.uce300$rf.distance,df.uce1000$rf.distance))

df.PIS.RF<-data.frame(cbind(dataset,PIS,PIS.prop,RF))

kruskal.test(PIS~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  PIS by dataset
# Kruskal-Wallis chi-squared = 578, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(as.numeric(df.PIS.RF$PIS), df.PIS.RF$dataset, p.adjust.method = "hochberg")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  as.numeric(df.PIS.RF$PIS) and df.PIS.RF$dataset 
# 
# BUSCOs UCE1000Flank UCE100Flank
# UCE1000Flank <2e-16 -            -          
#   UCE100Flank  <2e-16 <2e-16       -          
#   UCE300Flank  <2e-16 <2e-16       <2e-16     

kruskal.test(PIS.prop~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  PIS.prop by dataset
# Kruskal-Wallis chi-squared = 6408.9, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(as.numeric(df.PIS.RF$PIS.prop), df.PIS.RF$dataset, p.adjust.method = "hochberg")


# BUSCOs UCE1000Flank UCE100Flank
# UCE1000Flank <2e-16 -            -          
#   UCE100Flank  <2e-16 <2e-16       -          
#   UCE300Flank  <2e-16 <2e-16       <2e-16     
# 
# P value adjustment method: bonferroni 

kruskal.test(RF~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  RF by dataset
# Kruskal-Wallis chi-squared = 2740.1, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(as.numeric(df.PIS.RF$RF), df.PIS.RF$dataset, p.adjust.method = "hochberg")

# BUSCOs UCE1000Flank UCE100Flank
# UCE1000Flank <2e-16 -            -          
#   UCE100Flank  <2e-16 <2e-16       -          
#   UCE300Flank  0.25      <2e-16       <2e-16     
# 
# P value adjustment method: bonferroni

###########################
#FITTING NON LINEAR MODELS#
###########################

df.uce100<-read.csv("../uce100.csv")
for (i in 1:length(df.uce100$pars.site.prop)){
  if (df.uce100$pars.site.prop[i]==0){
    df.uce100$pars.site.prop[i]<-NA
  }
}
df.uce100<-na.omit(df.uce100)
df.uce100$pars.site.prop<-as.numeric(df.uce100$pars.site.prop)

fit.lm <- lm(rf.distance~pars.site, data = df.uce100)
fit.log <- lm(rf.distance~log(pars.site), data = df.uce100)

summary(fit.lm) #AIC: 25321
summary(fit.log) #AIC: 25304 Lowest AIC

# Call:
#   lm(formula = rf.distance ~ log(pars.site), data = df.uce100)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -70.216 -12.948   2.969  15.143  40.419 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    128.0628     1.8614   68.80   <2e-16 ***
#   log(pars.site)  -6.0136     0.5082  -11.83   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 19.62 on 2876 degrees of freedom
# Multiple R-squared:  0.04642,	Adjusted R-squared:  0.04609 
# F-statistic:   140 on 1 and 2876 DF,  p-value: < 2.2e-16

uce100.pis.plot2<-ggplot(df.uce100,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")

uce100.pis.plot2

fit.lm <- lm(rf.distance~pars.site.prop, data = df.uce100)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.uce100)

summary(fit.lm) #AIC: 25334
summary(fit.log) #AIC: 25304 Lowest AIC

# Call:
#   lm(formula = rf.distance ~ log(pars.site.prop), data = df.uce100)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -70.023 -12.831   3.013  15.186  40.610 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          93.0237     1.1936   77.93   <2e-16 ***
#   log(pars.site.prop)  -6.0591     0.5122  -11.83   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 19.62 on 2876 degrees of freedom
# Multiple R-squared:  0.04641,	Adjusted R-squared:  0.04608 
# F-statistic:   140 on 1 and 2876 DF,  p-value: < 2.2e-16

uce100.prop.plot2<-ggplot(df.uce100,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

uce100.prop.plot2

df.uce300<-read.csv("../uce300.csv")
for (i in 1:length(df.uce300$pars.site.prop)){
  if (df.uce300$pars.site.prop[i]==0){
    df.uce300$pars.site.prop[i]<-NA
  }
}
df.uce300<-na.omit(df.uce300)
df.uce300$pars.site.prop<-as.numeric(df.uce300$pars.site.prop)

fit.lm <- lm(rf.distance~pars.site, data = df.uce300)
fit.log <- lm(rf.distance~log(pars.site), data = df.uce300)

summary(fit.lm) #AIC: 22890
summary(fit.log) #AIC: 21945 Lowest AIC

# Call:
#   lm(formula = rf.distance ~ log(pars.site), data = df.uce300)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -31.090  -7.627  -0.711   6.741  56.860 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    212.2052     2.1152  100.32   <2e-16 ***
#   log(pars.site) -29.2997     0.4031  -72.69   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 10.93 on 2877 degrees of freedom
# Multiple R-squared:  0.6475,	Adjusted R-squared:  0.6474 
# F-statistic:  5284 on 1 and 2877 DF,  p-value: < 2.2e-16

uce300.pis.plot2<-ggplot(df.uce300,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")

uce300.pis.plot2

fit.lm <- lm(rf.distance~pars.site.prop, data = df.uce300)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.uce300)

summary(fit.lm) #AIC: 22604
summary(fit.log) #AIC: 21797 Lowest AIC

# Call:
#   lm(formula = rf.distance ~ log(pars.site.prop), data = df.uce300)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -32.403  -7.296  -0.612   6.712  54.765 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          14.5544     0.6225   23.38   <2e-16 ***
#   log(pars.site.prop) -31.2986     0.4140  -75.60   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 10.65 on 2877 degrees of freedom
# Multiple R-squared:  0.6652,	Adjusted R-squared:  0.6651 
# F-statistic:  5716 on 1 and 2877 DF,  p-value: < 2.2e-16

uce300.prop.plot2<-ggplot(df.uce300,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

uce300.prop.plot2

df.uce1000<-read.csv("../uce1000.csv")
for (i in 1:length(df.uce1000$pars.site.prop)){
  if (df.uce1000$pars.site.prop[i]==0){
    df.uce1000$pars.site.prop[i]<-NA
  }
}

df.uce1000<-na.omit(df.uce1000)
df.uce1000$pars.site.prop<-as.numeric(df.uce1000$pars.site.prop)

fit.lm <- lm(rf.distance~pars.site.prop, data = df.uce1000)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.uce1000)

summary(fit.lm) #AIC: 19752
summary(fit.log) #AIC: 19710 Lowest AIC

# Call:
#   lm(formula = rf.distance ~ log(pars.site.prop), data = df.uce1000)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -20.926  -5.795  -0.813   4.814  60.568 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          22.2254     0.6407   34.69   <2e-16 ***
#   log(pars.site.prop)  -8.5341     0.6939  -12.30   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.249 on 2790 degrees of freedom
# Multiple R-squared:  0.05142,	Adjusted R-squared:  0.05108 
# F-statistic: 151.2 on 1 and 2790 DF,  p-value: < 2.2e-16

uce1000.prop.plot2<-ggplot(df.uce1000,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

uce1000.prop.plot2

fit.lm <- lm(rf.distance~pars.site, data = df.uce1000)
fit.log <- lm(rf.distance~log(pars.site), data = df.uce1000)

summary(fit.lm) #AIC: 19828
summary(fit.log) #AIC: 19785 Lowest AIC

# Call:
#   lm(formula = rf.distance ~ log(pars.site), data = df.uce1000)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -20.827  -5.842  -0.780   4.887  62.023 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     64.6431     4.0553  15.940   <2e-16 ***
#   log(pars.site)  -5.0330     0.5865  -8.582   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.36 on 2790 degrees of freedom
# Multiple R-squared:  0.02572,	Adjusted R-squared:  0.02537 
# F-statistic: 73.65 on 1 and 2790 DF,  p-value: < 2.2e-16

uce1000.pis.plot2<-ggplot(df.uce1000,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")

uce1000.pis.plot2

df.CDS<-read.csv("../CDS.csv")
for (i in 1:length(df.CDS$pars.site.prop)){
  if (df.CDS$pars.site.prop[i]==0){
    df.CDS$pars.site.prop[i]<-NA
  }
}
df.CDS<-na.omit(df.CDS)
df.CDS$pars.site.prop<-as.numeric(df.CDS$pars.site.prop)
df.CDS$pars.site<-as.numeric(df.CDS$pars.site)
df.CDS$rf.distance<-as.numeric(df.CDS$rf.distance)

fit.lm <- lm(rf.distance~pars.site.prop, data = df.CDS)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.CDS)

summary(fit.lm) #AIC: 22374 Lowest AIC
summary(fit.log) #AIC: 22385

# Call:
#   lm(formula = rf.distance ~ pars.site.prop, data = df.CDS)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -46.37 -15.98  -1.78  14.41  77.40 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      49.266      1.563  31.526  < 2e-16 ***
#   pars.site.prop   24.903      4.050   6.149 9.03e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 20.96 on 2505 degrees of freedom
# Multiple R-squared:  0.01487,	Adjusted R-squared:  0.01448 
# F-statistic: 37.81 on 1 and 2505 DF,  p-value: 9.026e-10

CDS.prop.plot2<-ggplot(df.CDS,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("CDS", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot2

fit.lm <- lm(rf.distance~pars.site, data = df.CDS)
fit.log <- lm(rf.distance~log(pars.site), data = df.CDS)

summary(fit.lm) #AIC: 22306
summary(fit.log) #AIC: 22286 Lowest AIC
# 
# Call:
#   lm(formula = rf.distance ~ log(pars.site), data = df.CDS)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -42.775 -16.040  -2.472  13.961  73.677 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    103.2857     3.9573   26.10   <2e-16 ***
#   log(pars.site)  -6.7871     0.5968  -11.37   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 20.59 on 2505 degrees of freedom
# Multiple R-squared:  0.0491,	Adjusted R-squared:  0.04872 
# F-statistic: 129.3 on 1 and 2505 DF,  p-value: < 2.2e-16

CDS.pis.plot2<-ggplot(df.CDS,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("CDS", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.pis.plot2


####Combine UCEs###
dataset<-c(rep("uce 100", length(df.uce100[,1])),rep("uce 300", length(df.uce300[,1])),rep("uce 1000", length(df.uce1000[,1])))

df.uce<-rbind(df.uce100,df.uce300,df.uce1000)
df.uce<-cbind(dataset,df.uce)

fit.lm <- lm(rf.distance~pars.site.prop, data = df.uce)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.uce)

summary(fit.lm) #AIC: 76369
summary(fit.log) #AIC: 77413

# Call:
#   lm(formula = rf.distance ~ pars.site.prop, data = df.uce)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -68.409 -14.706  -2.148  12.980  96.411 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     118.5185     0.4785   247.7   <2e-16 ***
#   pars.site.prop -195.5041     1.5521  -126.0   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 21.06 on 8547 degrees of freedom
# Multiple R-squared:  0.6499,	Adjusted R-squared:  0.6499 
# F-statistic: 1.587e+04 on 1 and 8547 DF,  p-value: < 2.2e-16

uce.prop.plot2<-ggplot(df.uce,aes(x=pars.site.prop,y=rf.distance, col=dataset)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values = c("#ca7dcc",
                                "#1b98e0",
                                "#353436"))+
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 3 Datasets", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+
  theme(legend.position = "none")

uce.prop.plot2

fit.lm <- lm(rf.distance~pars.site, data = df.uce)
fit.log <- lm(rf.distance~log(pars.site), data = df.uce)

summary(fit.lm) #AIC: 77987
summary(fit.log) #AIC: 72165 Lowest AIC

# Call:
#   lm(formula = rf.distance ~ log(pars.site), data = df.uce)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -93.005  -9.806  -1.036   8.944  65.219 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    179.0050     0.6651   269.2   <2e-16 ***
#   log(pars.site) -21.7212     0.1226  -177.1   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 16.47 on 8547 degrees of freedom
# Multiple R-squared:  0.7859,	Adjusted R-squared:  0.7858 
# F-statistic: 3.137e+04 on 1 and 8547 DF,  p-value: < 2.2e-16

uce.pis.plot2<-ggplot(df.uce,aes(x=pars.site,y=rf.distance, col=dataset)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  scale_color_manual(values = c("#ca7dcc",
                                "#1b98e0",
                                "#353436"))+
  ggtitle("UCE 3 Datasets", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+
  theme(legend.position = c(0.75,0.75),legend.background = element_rect(fill="white",
                                                                        size=0.5, linetype="solid", 
                                                                        colour ="black"))

uce.pis.plot2

#plot(Figure 6)

pdf(file = "../../figs/Figure6b.pdf", width = 14, height = 7)
ggarrange(CDS.pis.plot2,uce100.pis.plot2,uce300.pis.plot2,uce1000.pis.plot2,uce.pis.plot2,
          CDS.prop.plot2,uce100.prop.plot2,uce300.prop.plot2,uce1000.prop.plot2,uce.prop.plot2,
          ncol = 5, nrow = 2)
dev.off()


####################
#####RF Heatmap#####
####################

uce.100.rax<-read.tree("../../Trees/UCE.100flank.crypt.raxml.boots.tre")
uce.100.ast<-read.tree("../../Trees/UCE.100flank.crypt.ASTRAL.PP.tre")

uce.300.rax<-read.tree("../../Trees/UCE.300flank.crypt.raxml.boots.tre")
uce.300.ast<-read.tree("../../Trees/UCE.300flank.crypt.ASTRAL.PP.tre")

uce.1000.rax<-read.tree("../../Trees/UCE.1000flank.crypt.raxml.boots.tre")
uce.1000.ast<-read.tree("../../Trees/UCEs.1000flank.crypt.ASTRAL.PP.tre")

cds.rax<-read.tree("../../Trees/CDS.crypt.raxml.boots.tre")
cds.ast<-read.tree("../../Trees/CDS.crypt.ASTRAL.PP.tre")

trees<-c(uce.100.rax,uce.100.ast,uce.300.rax,uce.300.ast,uce.1000.rax,uce.1000.ast,cds.rax,cds.ast)
trees.list<-c("uce.100.rax","uce.100.ast","uce.300.rax","uce.300.ast","uce.1000.rax","uce.1000.ast","cds.rax","cds.ast")
par(mfrow=c(1,2))
plot(cds.rax)
plot(cds.ast)
trees1<-c()
trees2<-c()
RF<-c()
count=0

for (i in 1:8){
  for (j in 1:8){
    count=count+1
    trees1[count]<-trees.list[i]
    trees2[count]<-trees.list[j]
    print(trees1[count])
    print(trees2[count])
    RF[count]<-RF.dist(tree1=trees[i][[1]],tree2=trees[j][[1]])
  }
}


tab<-data.frame(cbind(trees1,trees2,as.numeric(RF)))

tab$trees1<-factor(trees1, levels = c("cds.rax","uce.100.rax","uce.300.rax","uce.1000.rax","cds.ast","uce.100.ast","uce.300.ast","uce.1000.ast"))
tab$trees2<-factor(trees2, levels = c("cds.rax","uce.100.rax","uce.300.rax","uce.1000.rax","cds.ast","uce.100.ast","uce.300.ast","uce.1000.ast"))

library(hrbrthemes)


##Plot Figure 3 heatmap

ggplot(tab, aes(trees1, trees2, fill= as.numeric(RF))) + 
  geom_tile() +
  scale_fill_gradient2(low = "steelblue3", high = "tomato3", mid = "ivory", 
                       midpoint = max(RF)/2, limit = c(min(RF),max(RF)), space = "Lab", 
                       name="RF Distance") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+xlab("") + ylab("")+ 
  geom_text(aes(trees1, trees2, label = RF), color = "black", size = 3) +
  
  ggtitle(label = "Species Tree Differences for Crypt. Clade A")

