#This script is written by Lukas J Musher. Please cite:
#Musher, L.J., et al. ... TINAMOU PAPER UPDATE

library(ggplot2)
library(gplots)
library(phytools)
library(ape)
library(phangorn)
library(ips)
library(ggpubr)

#must have alignments and trees--if there are 20 alignments there must be 20 trees and 20 trees with collapsed low 
#support nodes--and must be in same directory and named appropriately so you know you have the right tree with each 
#gene (they must be in the exact same order) ideally you will have alignments and trees named exactly the same e.g. 
#uce-10.nexus uce-10.treefile uce-10.collapsed.tre

#function to get number of parsimony informative sites
pis<-function(x){
  len<-length(x[[1]])
  x<-read.nexus.data(x)
  site<-c()
  counts=0
  for(j in 1:length(x[[1]])){
    for(i in 1:length(x)){
      site[i]<-(x[[i]][j])
    }
    if(!("?" %in% unique(site))&!("-" %in% unique(site))){
      if (length(unique(site))>=2){
        counts=counts+1
      }
    }
    else if(!("?" %in% unique(site))){
      if (length(unique(site))>=3){
        counts=counts+1
      }
    }
    else if(!("-" %in% unique(site))){
      if (length(unique(site))>=3){
        counts=counts+1
      }
    }
    else if (length(unique(site))>4){
      #print(unique(site))
      counts=counts+1
    }
    else counts=counts
  }
  return(counts)
}

#how many bp sites in the alignment
site_count<-function(x){
  x<-read.nexus.data(x)
  len<-length(x[[1]])
  return(len)
}

#function to get gc content proportion
get_gc_content<-function(x){
  x<-read.nexus.data(x)
  site<-c()
  counts=0
  for(j in 1:length(x[[1]])){
    for(i in 1:length(x)){
      counts=counts+1
      site[counts]<-(x[[i]][j])
    }
  }
  g<-sum(site=="g")
  c<-sum(site=="c")
  t<-sum(site=="t")
  a<-sum(site=="a")
  gc.cont<-(g+c)/(g+c+a+t)
  return(gc.cont)
  
}

#function to build table of variables for genes and gene trees including rf distances for both bifurcating and collapsed low-support nodes
get_rf_gc_pis_table<-function(files,trees.coll,trees,spec_tree,out.path){
  aln<-c()
  pars.site<-c()
  tot.sites<-c()
  rf.distance<-c()
  rf.distance.coll<-c()
  pars.site.prop<-c()
  gc.cont<-c()
  
  for (i in 1:length(files)){
    print(c(i,files[i]))
    aln[i]<-files[i]
    pars.site[i]<-as.numeric(pis(files[i]))
    tot.sites[i]<-as.numeric(site_count(files[i]))
    pars.site.prop[i]<-as.numeric(pars.site[i]/tot.sites[i])
    gc.cont[i]<-get_gc_content(files[i])
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch(rf.distance[i]<-as.numeric(RF.dist(read.tree(trees[i]),spec_tree), error = function(e) { skip_to_next <<- TRUE}))
    
    if(skip_to_next) { next }
    
    tryCatch(rf.distance.coll[i]<-as.numeric(RF.dist(read.tree(trees.coll[i]),spec_tree), error = function(e) { skip_to_next <<- TRUE}))
    
    if(skip_to_next) { next }
    
  }
  
  df<-as.data.frame(cbind(aln,pars.site,tot.sites,pars.site.prop,rf.distance,rf.distance.coll, gc.cont))
  write.csv(df,out.path)
  return(df)
  
}

#set working directory
setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/CDS_complete_nexus/")

#define files, trees, and collapsed trees
files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees <- list.files(path="./",pattern="*treefile",full.names=F,recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*collapsed.tre",full.names=F,recursive=FALSE)

#check to make sure same numbers of trees and alignments
length(files)
length(trees)
length(trees.coll)

#define species tree 
spec_tree <- root(read.tree(file = "../../Trees/CDS.complete.ASTRAL.PP.tre"), "Rhea_penn", resolve.root = T)
plot(spec_tree)

#get results for BUSCOs
#df.CDS<-get_rf_gc_pis_table(files,trees.coll,trees,spec_tree,out.path="../CDS.collapsed.csv")
df.CDS<-read.csv("../CDS.collapsed.csv")

df.CDS$pars.site<-as.numeric(df.CDS$pars.site)
df.CDS$tot.sites<-as.numeric(df.CDS$tot.sites)
df.CDS$pars.site.prop<-as.numeric(df.CDS$pars.site.prop)
df.CDS$rf.distance<-as.numeric(df.CDS$rf.distance)
df.CDS$rf.distance.coll<-as.numeric(df.CDS$rf.distance.coll)
df.CDS$gc.cont<-as.numeric(df.CDS$gc.cont)

mean(df.CDS$pars.site) #922.5146
sd(df.CDS$pars.site) #690.084

mean(df.CDS$pars.site.prop) #0.3717782
sd(df.CDS$pars.site.prop) #0.1033898

mean(df.CDS$rf.distance) #58.52413
sd(df.CDS$rf.distance) #21.11341

mean(df.CDS$rf.distance.coll) #48.90307
sd(df.CDS$rf.distance.coll) #16.69681

mean(df.CDS$gc.cont) #0.5088858
sd(df.CDS$gc.cont) #0.0618374

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/alignments_uce-1000flank-complete-nexus/")

files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*.collapsed.tre",full.names=F,recursive=FALSE)
trees <- list.files(path="./",pattern="*.treefile",full.names=F,recursive=FALSE)

length(trees.coll)
length(files)
length(trees)

spec_tree <- root(read.tree(file = "../../Trees/UCEs.1000flank.complete.ASTRAL.PP.tre"), "Rhea_penn_GCA_003342835", resolve.root = T)
plot(spec_tree)

#df.uce1000<-get_rf_gc_pis_table(files,trees.coll,trees,spec_tree,out.path="../uce1000.collapsed.csv")
df.uce1000<-read.csv("../uce1000.collapsed.csv")

mean(df.uce1000$pars.site) #1035.713
sd(df.uce1000$pars.site) #250.5732
mean(df.uce1000$pars.site.prop) #0.4178852
sd(df.uce1000$pars.site.prop) #0.08261385
mean(df.uce1000$rf.distance) #29.86748
sd(df.uce1000$rf.distance) #8.468606
mean(df.uce1000$rf.distance.coll) #26.66082
sd(df.uce1000$rf.distance.coll) #7.134672
mean(df.uce1000$gc.cont) #0.4017327
sd(df.uce1000$gc.cont) #0.05414255

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/alignments_uce-300flank-complete-nexus/")

files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*.collapsed.tre",full.names=F,recursive=FALSE)
trees <- list.files(path="./",pattern="*.treefile",full.names=F,recursive=FALSE)

length(trees.coll)
length(files)
length(trees)

#df.uce300<-get_rf_gc_pis_table(files,trees.coll,trees,spec_tree,out.path="../uce300.collapsed.csv")
df.uce300<-read.csv("../uce300.collapsed.csv")

mean(df.uce300$pars.site) #206.8736
sd(df.uce300$pars.site) #88.36566
mean(df.uce300$pars.site.prop) #0.2648982
sd(df.uce300$pars.site.prop) #0.1048726
mean(df.uce300$rf.distance) #59.16151
sd(df.uce300$rf.distance) #18.40835
mean(df.uce300$rf.distance.coll) #46.59257
sd(df.uce300$rf.distance.coll) #13.05453
mean(df.uce300$gc.cont) #0.3927467
sd(df.uce300$gc.cont) #0.05663905

plot(df.uce300$pars.site,df.uce300$rf.distance.coll)

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/alignments_uce-100flank-complete-nexus/")
files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*collapsed.tre",full.names=F,recursive=FALSE)
trees <- list.files(path="./",pattern="*.treefile",full.names=F,recursive=FALSE)

length(trees.coll)
length(files)
length(trees)

plot(read.tree("uce-1003.treefile"))

#df.uce100<-get_rf_gc_pis_table(files,trees.coll,trees,spec_tree,out.path="../uce100.collapsed.csv")

df.uce100<-read.csv("../uce100.collapsed.csv")

mean(df.uce100$pars.site) #45.20298
sd(df.uce100$pars.site) #28.84223
mean(df.uce100$pars.site.prop) #0.1347156
sd(df.uce100$pars.site.prop) #0.0828647
mean(na.omit(df.uce100$rf.distance)) #106.5005
sd(na.omit(df.uce100$rf.distance)) #20.10551
mean(na.omit(df.uce100$rf.distance.coll)) #107.9626
sd(na.omit(df.uce100$rf.distance.coll)) #20.1613
mean(na.omit(df.uce100$gc.cont)) #0.3899228
sd(na.omit(df.uce100$gc.cont)) #0.05682132

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/alignments_uce-100flank-complete-nexus/")
files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*collapsed.tre",full.names=F,recursive=FALSE)
trees <- list.files(path="./",pattern="*.treefile",full.names=F,recursive=FALSE)

length(trees.coll)
length(files)
length(trees)

plot(read.tree("uce-1003.treefile"))

df.uce100<-get_rf_gc_pis_table(files,trees.coll,trees,spec_tree,out.path="../uce100.collapsed.csv")

df.uce100<-read.csv("../uce100.collapsed.csv")

mean(df.uce100$pars.site) #45.20298
sd(df.uce100$pars.site) #28.84223
mean(df.uce100$pars.site.prop) #0.1347156
sd(df.uce100$pars.site.prop) #0.0828647
mean(na.omit(df.uce100$rf.distance)) #106.5005
sd(na.omit(df.uce100$rf.distance)) #20.10551
mean(na.omit(df.uce100$rf.distance.coll)) #106.5005
sd(na.omit(df.uce100$rf.distance.coll)) #20.10551
mean(na.omit(df.uce100$gc.cont)) #106.5005
sd(na.omit(df.uce100$gc.cont)) #20.10551

###########################
###Kruskall Wallace Test###
###########################

dataset<-c(rep("BUSCOs",length(df.CDS[,1])),rep("UCE100Flank",length(df.uce100[,1])),rep("UCE300Flank",length(df.uce300[,1])),rep("UCE1000Flank",length(df.uce1000[,1])))
PIS<-as.numeric(c(df.CDS$pars.site,df.uce100$pars.site,df.uce300$pars.site,df.uce1000$pars.site))
PIS.prop<-as.numeric(c(df.CDS$pars.site.prop,df.uce100$pars.site.prop,df.uce300$pars.site.prop,df.uce1000$pars.site.prop))
RF<-as.numeric(c(df.CDS$rf.distance,df.uce100$rf.distance,df.uce300$rf.distance,df.uce1000$rf.distance))
RF.coll<-as.numeric(c(df.CDS$rf.distance.coll,df.uce100$rf.distance.coll,df.uce300$rf.distance.coll,df.uce1000$rf.distance.coll))
GC<-as.numeric(c(df.CDS$gc.cont,df.uce100$gc.cont,df.uce300$gc.cont,df.uce1000$gc.cont))

df.PIS.RF<-data.frame(cbind(dataset,PIS,PIS.prop,RF,RF.coll,GC))

kruskal.test(PIS~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  PIS by dataset
# Kruskal-Wallis chi-squared = 580.06, df = 3, p-value < 2.2e-16

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
# Kruskal-Wallis chi-squared = 6406.1, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(as.numeric(df.PIS.RF$PIS.prop), df.PIS.RF$dataset, p.adjust.method = "hochberg")


# BUSCOs UCE1000Flank UCE100Flank
# UCE1000Flank <2e-16 -            -          
#   UCE100Flank  <2e-16 <2e-16       -          
#   UCE300Flank  <2e-16 <2e-16       <2e-16     
# 
# P value adjustment method: hochberg 

kruskal.test(RF.coll~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  RF by dataset
# Kruskal-Wallis chi-squared = 2898.6, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(as.numeric(df.PIS.RF$RF), df.PIS.RF$dataset, p.adjust.method = "hochberg")

# BUSCOs UCE1000Flank UCE100Flank
# UCE1000Flank <2e-16 -            -          
#   UCE100Flank  <2e-16 <2e-16       -          
#   UCE300Flank  0.0075      <2e-16       <2e-16     
# 
# P value adjustment method: hochberg

kruskal.test(RF.coll~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  RF by dataset
# Kruskal-Wallis chi-squared = 3051.8, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(as.numeric(df.PIS.RF$RF.coll), df.PIS.RF$dataset, p.adjust.method = "hochberg")

# BUSCOs UCE1000Flank UCE100Flank
# UCE1000Flank <2e-16 -            -          
#   UCE100Flank  <2e-16 <2e-16       -          
#   UCE300Flank  8.6e-07      <2e-16       <2e-16     
# 
# P value adjustment method: hochberg

kruskal.test(GC~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  RF by dataset
# Kruskal-Wallis chi-squared = 4079.4, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(as.numeric(df.PIS.RF$GC), df.PIS.RF$dataset, p.adjust.method = "hochberg")

#             BUSCOs  UCE1000Flank UCE100Flank
# UCE1000Flank < 2e-16 -            -          
# UCE100Flank  < 2e-16 3.8e-13      -          
# UCE300Flank  < 2e-16 8.0e-12      0.39   
# P value adjustment method: hochberg

###########################
#FITTING NON LINEAR MODELS#
###########################

z.uces<-read.table("../../z.chrom.list.txt")$V1

df.uce100<-read.csv("../uce100.collapsed.csv")
for (i in 1:length(df.uce100$pars.site.prop)){
  if (df.uce100$pars.site.prop[i]==0){
    df.uce100$pars.site.prop[i]<-NA
  }
}

sum(df.uce100$aln %in% z.uces)

chrom<-c()
for (i in 1:length(df.uce100$pars.site.prop)){
  if (df.uce100$aln[i] %in% z.uces){
    chrom[i]<-"Z-chromosome"
    print("found Z")
  }
  if (!(df.uce100$aln[i] %in% z.uces)){
    chrom[i]<-"Autosome"
  }
  
}

df.uce100<-cbind(df.uce100,chrom)
sum(df.uce100$chrom=="Z-chromosome")

df.uce100<-na.omit(df.uce100)
df.uce100$pars.site.prop<-as.numeric(df.uce100$pars.site.prop)

fit.glm <- glm(rf.distance.coll~pars.site, data = df.uce100)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.uce100)

summary(fit.glm) #AIC: 21684 Lowest AIC
summary(fit.glog) #AIC: 22041

fit.lm <- lm(rf.distance.coll~pars.site, data = df.uce100)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.uce100)
summary(fit.lm) 
summary(fit.log)

uce100.pis.plot2<-ggplot(df.uce100,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+
  geom_point(data=df.uce100[df.uce100$chrom=="Z-chromosome",], mapping =aes(x=pars.site, y=rf.distance.coll), color="yellow", size=1)

uce100.pis.plot2

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.uce100)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.uce100)

summary(fit.glm) #AIC: 21525 Lowest AIC
summary(fit.glog) #AIC: 22119

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.uce100)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.uce100)
summary(fit.lm)
summary(fit.log)

uce100.prop.plot2<-ggplot(df.uce100,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+
  geom_point(data=df.uce100[df.uce100$chrom=="Z-chromosome",], mapping =aes(x=pars.site.prop, y=rf.distance.coll), color="yellow", size=1)

uce100.prop.plot2

df.uce300<-read.csv("../uce300.collapsed.csv")
for (i in 1:length(df.uce300$pars.site.prop)){
  if (df.uce300$pars.site.prop[i]==0){
    df.uce300$pars.site.prop[i]<-NA
  }
}

sum(df.uce300$aln %in% z.uces)

chrom<-c()
for (i in 1:length(df.uce300$pars.site.prop)){
  if (df.uce300$aln[i] %in% z.uces){
    chrom[i]<-"Z-chromosome"
    print("found Z")
  }
  if (!(df.uce300$aln[i] %in% z.uces)){
    chrom[i]<-"Autosome"
  }
  
}

df.uce300<-cbind(df.uce300,chrom)
sum(df.uce300$chrom=="Z-chromosome")

df.uce300<-na.omit(df.uce300)
df.uce300$pars.site.prop<-as.numeric(df.uce300$pars.site.prop)

fit.glm <- glm(rf.distance.coll~pars.site, data = df.uce300)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.uce300)

summary(fit.glm) #AIC: 21345
summary(fit.glog) #AIC: 20216 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site, data = df.uce300)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.uce300)
summary(fit.lm) 
summary(fit.log)

uce300.pis.plot2<-ggplot(df.uce300,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+
  geom_point(data=df.uce300[df.uce300$chrom=="Z-chromosome",], mapping =aes(x=pars.site, y=rf.distance.coll), color="yellow", size=1)

uce300.pis.plot2

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.uce300)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.uce300)

summary(fit.glm) #AIC: 21092
summary(fit.glog) #AIC: 20047 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.uce300)
fit.log <- lm(rf.distance.coll~log(pars.site.prop), data = df.uce300)
summary(fit.lm) #AIC: 21092
summary(fit.log) #AIC: 20047 Lowest AIC

uce300.prop.plot2<-ggplot(df.uce300,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+
  geom_point(data=df.uce300[df.uce300$chrom=="Z-chromosome",], mapping =aes(x=pars.site.prop, y=rf.distance.coll), color="yellow", size=1)

uce300.prop.plot2

fit.glm <- glm(rf.distance.coll~gc.cont, data = df.uce300)
fit.glog <- glm(rf.distance.coll~log(gc.cont), data = df.uce300)

summary(fit.glm) #AIC: 22941 Lowest AIC
summary(fit.glog) #AIC: 22942

fit.lm <- lm(rf.distance.coll~gc.cont, data = df.uce300)
#fit.log <- lm(rf.distance.coll~log(gc.cont), data = df.uce300)

uce300.gc.plot2<-ggplot(df.uce300,aes(x=gc.cont,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("GC Content") + ylab("RF Distance")+
  geom_point(data=df.uce300[df.uce300$chrom=="Z-chromosome",], mapping =aes(x=gc.cont, y=rf.distance.coll), color="yellow", size=1)

uce300.gc.plot2

df.uce1000<-read.csv("../uce1000.collapsed.csv")
for (i in 1:length(df.uce1000$pars.site.prop)){
  if (df.uce1000$pars.site.prop[i]==0){
    df.uce1000$pars.site.prop[i]<-NA
  }
}

sum(df.uce1000$aln %in% z.uces)

chrom<-c()
for (i in 1:length(df.uce1000$pars.site.prop)){
  if (df.uce1000$aln[i] %in% z.uces){
    chrom[i]<-"Z-chromosome"
    print("found Z")
  }
  if (!(df.uce1000$aln[i] %in% z.uces)){
    chrom[i]<-"Autosome"
  }
  
}

df.uce1000<-cbind(df.uce1000,chrom)
sum(df.uce1000$chrom=="Z-chromosome")

plot(as.factor(df.uce1000$chrom),df.uce1000$rf.distance.coll)

df.uce1000<-na.omit(df.uce1000)
df.uce1000$pars.site.prop<-as.numeric(df.uce1000$pars.site.prop)

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.uce1000)
fit.glog <- glm(rf.distance~log(pars.site.prop), data = df.uce1000)

summary(fit.glm) #AIC: 18843 Lowest AIC
summary(fit.glog) #AIC: 19711

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.uce1000)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.uce1000)
summary(fit.lm) 
summary(fit.log)

uce1000.prop.plot2<-ggplot(df.uce1000,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+
  geom_point(data=df.uce1000[df.uce1000$chrom=="Z-chromosome",], mapping =aes(x=pars.site.prop, y=rf.distance.coll), color="yellow", size=1)

uce1000.prop.plot2

fit.glm <- glm(rf.distance.coll~pars.site, data = df.uce1000)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.uce1000)

summary(fit.glm) #AIC: 18893
summary(fit.glog) #AIC: 18870 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site, data = df.uce1000)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.uce1000)
summary(fit.lm)
summary(fit.log)

uce1000.pis.plot2<-ggplot(df.uce1000,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+
  geom_point(data=df.uce1000[df.uce1000$chrom=="Z-chromosome",], mapping =aes(x=pars.site, y=rf.distance.coll), color="yellow", size=1)

uce1000.pis.plot2

fit.glm <- glm(rf.distance.coll~gc.cont, data = df.uce1000)
fit.glog <- glm(rf.distance.coll~log(gc.cont), data = df.uce1000)

summary(fit.glm) #AIC: 18884 Lowest AIC
summary(fit.glog) #AIC: 18884

fit.lm <- lm(rf.distance.coll~gc.cont, data = df.uce1000)
#fit.log <- lm(rf.distance.coll~log(gc.cont), data = df.uce300)

uce1000.gc.plot2<-ggplot(df.uce1000,aes(x=gc.cont,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("GC Content") + ylab("RF Distance")+
  geom_point(data=df.uce1000[df.uce1000$chrom=="Z-chromosome",], mapping =aes(x=gc.cont, y=rf.distance.coll), color="yellow", size=1)

uce1000.gc.plot2

df.CDS<-read.csv("../CDS.collapsed.csv")
for (i in 1:length(df.CDS$pars.site.prop)){
  if (df.CDS$pars.site.prop[i]==0){
    df.CDS$pars.site.prop[i]<-NA
  }
}
df.CDS<-na.omit(df.CDS)
df.CDS$pars.site.prop<-as.numeric(df.CDS$pars.site.prop)
df.CDS$pars.site<-as.numeric(df.CDS$pars.site)
df.CDS$rf.distance<-as.numeric(df.CDS$rf.distance)

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.CDS)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.CDS)

summary(fit.glm) #AIC: 21177 Lowest AIC
summary(fit.glog) #AIC: 21189

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.CDS)
fit.log <- lm(rf.distance.coll~log(pars.site.prop), data = df.CDS)
summary(fit.lm) 
summary(fit.log)

CDS.prop.plot2<-ggplot(df.CDS,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("BUSCOs", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot2

fit.glm <- glm(rf.distance.coll~pars.site, data = df.CDS)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.CDS)

summary(fit.glm) #AIC: 21179 Lowest AIC
summary(fit.glog) #AIC: 21179

fit.lm <- lm(rf.distance.coll~pars.site, data = df.CDS)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.CDS)
summary(fit.lm) 
summary(fit.log)

CDS.pis.plot2<-ggplot(df.CDS,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("BUSCOs", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# of informative sites") + ylab("RF Distance")

CDS.pis.plot2

fit.glm <- glm(rf.distance.coll~gc.cont, data = df.CDS)
fit.glog <- glm(rf.distance.coll~log(gc.cont), data = df.CDS)

summary(fit.glm) #AIC: 21234 Lowest AIC
summary(fit.glog) #AIC: 21234

fit.lm <- lm(rf.distance.coll~gc.cont, data = df.CDS)
#fit.log <- lm(rf.distance.coll~log(gc.cont), data = df.CDS)

CDS.gc.plot2<-ggplot(df.CDS,aes(x=gc.cont,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("CDS", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("GC Content") + ylab("RF Distance")

CDS.gc.plot2

####Combine UCEs###
dataset<-c(rep("uce 100", length(df.uce100[,1])),rep("uce 300", length(df.uce300[,1])),rep("uce 1000", length(df.uce1000[,1])))

df.uce<-rbind(df.uce100,df.uce300,df.uce1000)
df.uce<-cbind(dataset,df.uce)

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.uce)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.uce)

summary(fit.glm) #AIC: 76049
summary(fit.glog) #AIC: 74480 Lowest AIC

fit.lm <-lm(rf.distance.coll~pars.site.prop, data = df.uce)
fit.log <-lm(rf.distance.coll~log(pars.site.prop), data = df.uce)
summary(fit.lm) 
summary(fit.log)

uce.prop.plot2<-ggplot(df.uce,aes(x=pars.site.prop,y=rf.distance.coll, col=dataset)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values = c("#ca7dcc",
                                "#1b98e0",
                                "#353436"))+
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 3 Datasets", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+
  theme(legend.position = "none")+
  geom_point(data=df.uce[df.uce$chrom=="Z-chromosome",], mapping =aes(x=pars.site.prop, y=rf.distance.coll), color="Yellow", size=1)

uce.prop.plot2

fit.glm <- glm(rf.distance.coll~pars.site, data = df.uce)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.uce)

summary(fit.glm) #AIC: 80245
summary(fit.glog) #AIC: 69350 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site, data = df.uce)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.uce)
summary(fit.lm) 
summary(fit.log)

uce.pis.plot2<-ggplot(df.uce,aes(x=pars.site,y=rf.distance.coll, col=dataset)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", se = T ,formula = y~log(x), color="red")+
  scale_color_manual(values = c("#ca7dcc",
                                "#1b98e0",
                                "#353436"))+
  ggtitle("UCE 3 Datasets", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+
  theme(legend.position = c(0.75,0.75),legend.background = element_rect(fill="white",
                                                                        size=0.5, linetype="solid", 
                                                                        colour ="black"))+
  geom_point(data=df.uce[df.uce$chrom=="Z-chromosome",], mapping =aes(x=pars.site, y=rf.distance.coll), color="yellow", size=1)

uce.pis.plot2

#plot(Figure 6)

pdf(file = "../../figs/Figure6b.pdf", width = 14, height = 7)
ggarrange(CDS.pis.plot2,uce100.pis.plot2,uce300.pis.plot2,uce1000.pis.plot2,uce.pis.plot2,
          CDS.prop.plot2,uce100.prop.plot2,uce300.prop.plot2,uce1000.prop.plot2,uce.prop.plot2,
          ncol = 5, nrow = 2)
dev.off()

#boxplots
chrom<-rep(NA,length(df.CDS[,1]))
df.CDS<-cbind(df.CDS,chrom)
dataset<-c(rep("uce 100", length(df.uce100[,1])),rep("uce 300", length(df.uce300[,1])),rep("uce 1000", length(df.uce1000[,1])),rep("BUSCO", length(df.CDS[,1])))

df.all<-rbind(df.uce100,df.uce300,df.uce1000,df.CDS)
df.all<-cbind(dataset,df.all)

df.all$dataset<-factor(df.all$dataset, levels = c("BUSCO", "uce 100", "uce 300", "uce 1000"))

a<-ggplot(df.all, aes(x=dataset, y=rf.distance.coll)) + 
  geom_violin()+ geom_boxplot(width=0.1)+ggtitle("RF Distance by dataset")+
  theme_bw() + xlab("") + ylab("RF Distance")

b<-ggplot(df.all, aes(x=dataset, y=pars.site)) + 
  geom_violin()+ geom_boxplot(width=0.1)+ggtitle("PIS by dataset")+
  theme_bw() + xlab("") + ylab("Number of PIS")

c<-ggplot(df.all, aes(x=dataset, y=pars.site.prop*100)) + 
  geom_violin()+ geom_boxplot(width=0.1)+ggtitle("Percent PIS by dataset")+
  theme_bw() + xlab("") + ylab("Percent PIS")

d<-ggplot(df.all, aes(x=dataset, y=gc.cont)) + 
  geom_violin()+ geom_boxplot(width=0.1)+ggtitle("GC content by dataset")+
  theme_bw() + xlab("") + ylab("GC content")

names(df.uce)

#does RF of Z-chromosome UCEs differ from RF of Autosomal UCEs?
df.uce$dataset<-factor(df.uce$dataset, levels = c("uce 100", "uce 300", "uce 1000"))
e<-ggplot(df.uce, aes(x=dataset, y=rf.distance.coll)) + 
  geom_violin(aes(color = chrom), trim = FALSE, position = position_dodge(0.9) ) +
  geom_boxplot(aes(color = chrom), width = 0.1, position = position_dodge(0.9)) +
  scale_color_manual(values = c("#353436","#1b98e0"))+
  theme_bw() + xlab("") + ylab("RF Distance")+ggtitle("RF Distance by chromosome")

e

mean(df.uce1000$rf.distance.coll[df.uce1000$chrom=="Autosome"]) #26.77393
sd(df.uce1000$rf.distance.coll[df.uce1000$chrom=="Autosome"]) #7.053628

mean(df.uce1000$rf.distance.coll[df.uce1000$chrom=="Z-chromosome"]) #22.68421
sd(df.uce1000$rf.distance.coll[df.uce1000$chrom=="Z-chromosome"]) #8.771485

kruskal.test(rf.distance.coll~chrom,data=df.uce1000)

# Kruskal-Wallis rank sum test
# 
# data:  rf.distance.coll by chrom
# Kruskal-Wallis chi-squared = 29.299, df = 1, p-value = 6.202e-08

kruskal.test(rf.distance.coll~chrom,data=df.uce300)

# Kruskal-Wallis rank sum test
# 
# data:  rf.distance.coll by chrom
# Kruskal-Wallis chi-squared = 0.26771, df = 1, p-value = 0.6049

kruskal.test(rf.distance.coll~chrom,data=df.uce100)

# Kruskal-Wallis rank sum test
# 
# data:  rf.distance.coll by chrom
# Kruskal-Wallis chi-squared = 0.35855, df = 1, p-value = 0.5493

###Concordance factors boxplots
setwd("~/Documents/ANSDU/Tinamous/concordance_factors/")

cfs<-na.omit(read.csv("CFS.all.datasets.csv"))
cfs$Dataset<-factor(cfs$Dataset, levels = c("BUSCOs","UCE100Flank", "UCE300Flank", "UCE1000Flank"))

kruskal.test(cfs$gCF[cfs$Method=="astral"]~cfs$Dataset[cfs$Method=="astral"])
kruskal.test(cfs$gCF[cfs$Dataset=="UCE1000Flank"]~cfs$Method[cfs$Dataset=="UCE1000Flank"])

mean(cfs$gCF[cfs$Dataset=="UCE1000Flank" & cfs$Method=="iqtree"])
mean(na.omit(cfs$gCF[cfs$Dataset=="UCE1000Flank" & cfs$Method=="astral"]))

pairwise.wilcox.test(cfs$gCF[cfs$Method=="astral"], cfs$Dataset[cfs$Method=="astral"], p.adjust.method = "hochberg")

f<-ggplot(cfs, aes(x=Dataset, y=gCF)) + 
  geom_violin(aes(color = Method), trim = FALSE, position = position_dodge(0.9) ) +
  geom_boxplot(aes(color = Method), width = 0.1, position = position_dodge(0.9)) +
  scale_color_manual(values = c("#353436","#ca7dcc"))+
  theme_bw() + xlab("") + ylab("RF Distance")+ggtitle("RF Distance by chromosome")

f

pdf(file = "../figs/Figure5b.pdf", width = 18, height = 7)
ggarrange(b,c,d,a,e,f,
          ncol = 3, nrow = 2)
dev.off()



#Look at CDS results for UCE-like subsamples
a<-mean(df.uce1000$pars.site.prop)+(2*sd(df.uce1000$pars.site.prop))
b<-mean(df.uce1000$pars.site.prop)-(2*sd(df.uce1000$pars.site.prop))
c<-mean(df.uce300$pars.site.prop)+(2*sd(df.uce300$pars.site.prop))
d<-mean(df.uce300$pars.site.prop)-(2*sd(df.uce300$pars.site.prop))
e<-mean(df.uce100$pars.site.prop)+(2*sd(df.uce100$pars.site.prop))
f<-mean(df.uce100$pars.site.prop)-(2*sd(df.uce100$pars.site.prop))

g<-mean(df.CDS$gc.cont)+(2*sd(df.CDS$gc.cont))
h<-mean(df.CDS$gc.cont)-(2*sd(df.CDS$gc.cont))

#filter out extremes of GC
df.CDS.filt.gc<-df.CDS[df.CDS$gc.cont<=g & df.CDS$gc.cont>=h,]

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.CDS.filt.gc)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.CDS.filt.gc)

summary(fit.glm) #AIC: 20449 Lowest AIC
summary(fit.glog) #AIC: 20461

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.CDS.filt.gc)

CDS.prop.plot3<-ggplot(df.CDS.filt.gc,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("BUSCOs filtered by GC content", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot3

fit.glm <- glm(rf.distance~pars.site, data = df.CDS.filt.gc)
fit.glog <- glm(rf.distance~log(pars.site), data = df.CDS.filt.gc)

summary(fit.glm) #AIC: 21603
summary(fit.glog) #AIC: 21578 Lowest AIC

fit.log <- lm(rf.distance~log(pars.site), data = df.CDS.filt.gc)

CDS.pis.plot4<-ggplot(df.CDS.filt.gc,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs filtered by GC content", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# of informative sites") + ylab("RF Distance")

CDS.pis.plot4

#CDS with comparable prop variable sites to all UCEs datasets
df.CDS.filt1<-df.CDS[df.CDS$pars.site.prop>=f & df.CDS$pars.site.prop<=a,]

fit.glm <- glm(rf.distance~pars.site.prop, data = df.CDS.filt1)
fit.glog <- glm(rf.distance~log(pars.site.prop), data = df.CDS.filt1)

summary(fit.glm) #AIC: 22017 Lowest AIC
summary(fit.glog) #AIC: 22024

fit.lm <- lm(rf.distance~pars.site.prop, data = df.CDS.filt1)

CDS.prop.plot6<-ggplot(df.CDS.filt1,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("BUSCOs: UCE-like filtered by PIS", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot6

fit.glm <- glm(rf.distance~pars.site, data = df.CDS.filt1)
fit.glog <- glm(rf.distance~log(pars.site), data = df.CDS.filt1)

summary(fit.glm) #AIC: 21905
summary(fit.glog) #AIC: 21883 Lowest AIC

fit.log <- lm(rf.distance~log(pars.site), data = df.CDS.filt1)

CDS.pis.plot7<-ggplot(df.CDS.filt1,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs: UCE-like filtered by PIS", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# of informative sites") + ylab("RF Distance")

CDS.pis.plot7

#filter by UCE and gc content
df.CDS.filt2<-df.CDS.filt1[df.CDS.filt1$gc.cont<=g & df.CDS.filt1$gc.cont>=h,]

fit.glm <- glm(rf.distance~pars.site.prop, data = df.CDS.filt2)
fit.glog <- glm(rf.distance~log(pars.site.prop), data = df.CDS.filt2)

summary(fit.glm) #AIC: 21264 Lowest AIC
summary(fit.glog) #AIC: 21270

fit.lm <- lm(rf.distance~pars.site.prop, data = df.CDS.filt2)

CDS.prop.plot8<-ggplot(df.CDS.filt2,aes(x=pars.site.prop,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("BUSCOs: UCE-like filtered by GC & PIS", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot8

fit.glm <- glm(rf.distance~pars.site, data = df.CDS.filt2)
fit.glog <- glm(rf.distance~log(pars.site), data = df.CDS.filt2)

summary(fit.glm) #AIC: 21153
summary(fit.glog) #AIC: 21130 Lowest AIC

fit.log <- lm(rf.distance~log(pars.site), data = df.CDS.filt2)

CDS.pis.plot9<-ggplot(df.CDS.filt2,aes(x=pars.site,y=rf.distance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs: UCE-like filtered by GC & PIS", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# of informative sites") + ylab("RF Distance")

CDS.pis.plot9


pdf(file = "../figs/Filtered_CDS_plots.pdf", width = 18, height = 7)
ggarrange(CDS.prop.plot3,CDS.prop.plot6,CDS.prop.plot8,CDS.pis.plot4,CDS.pis.plot7,CDS.pis.plot9,
          ncol = 3, nrow = 2)
dev.off()

####################
#####RF Heatmap#####
####################
setwd("~/Documents/ANSDU/Tinamous/Trees/")
uce.100.iqt<-read.tree("cryp.UCE100.iqtree.tre")
uce.100.ast<-read.tree("cryp.UCE100.astral.tre")
plot(uce.100.iqt)
uce.300.iqt<-read.tree("cryp.UCE300.iqtree.tre")
uce.300.ast<-read.tree("cryp.UCE300.astral.tre")

uce.1000.iqt<-read.tree("cryp.UCE1000.iqtree.tre")
uce.1000.ast<-read.tree("cryp.UCE1000.astral.tre")

uce.1000.ChrZ.iqt<-read.tree("ChrZ.UCE1000.crypt.iqtree.tre")
uce.1000.ChrZ.ast<-read.tree("ChrZ.UCE1000.crypt.astral.rooted.tre")
plot(uce.1000.ChrZ.ast)
cds.iqt<-read.tree("cryp.BUSCOs.iqtree.tre")
cds.ast<-read.tree("cryp.BUSCOs.astral.tre")

cds.iqt$tip.label[2]<-"CryCin_GCA_003342915"
cds.iqt$tip.label[4]<-"CryUnd_GCA_013389825"
cds.ast$tip.label[8]<-"CryCin_GCA_003342915"
cds.ast$tip.label[12]<-"CryUnd_GCA_013389825"

trees<-c(uce.100.iqt,uce.100.ast,uce.300.iqt,uce.300.ast,uce.1000.iqt,uce.1000.ast,uce.1000.ChrZ.iqt,uce.1000.ChrZ.ast,cds.iqt,cds.ast)
trees.list<-c("uce.100.iqt","uce.100.ast","uce.300.iqt","uce.300.ast","uce.1000.iqt","uce.1000.ast","uce.1000.ChrZ.iqt","uce.1000.ChrZ.ast","cds.iqt","cds.ast")

trees1<-c()
trees2<-c()
RF<-c()
count=0

for (i in 1:10){
  for (j in 1:10){
    count=count+1
    trees1[count]<-trees.list[i]
    trees2[count]<-trees.list[j]
    print(trees1[count])
    print(trees2[count])
    RF[count]<-RF.dist(tree1=trees[i][[1]],tree2=trees[j][[1]])
  }
}


tab<-data.frame(cbind(trees1,trees2,as.numeric(RF)))

tab$trees1<-factor(trees1, levels = c("cds.iqt","uce.100.iqt","uce.300.iqt","uce.1000.iqt","uce.1000.ChrZ.iqt","cds.ast","uce.100.ast","uce.300.ast","uce.1000.ast","uce.1000.ChrZ.ast"))
tab$trees2<-factor(trees2, levels = c("cds.iqt","uce.100.iqt","uce.300.iqt","uce.1000.iqt","uce.1000.ChrZ.iqt","cds.ast","uce.100.ast","uce.300.ast","uce.1000.ast","uce.1000.ChrZ.ast"))

library(hrbrthemes)


##Plot Figure 3 heatmap

pdf(file = "../figs/heatmap.raw.pdf", width = 7, height = 7)
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

dev.off()

#plot figure 3 trees
pdf(file = "../figs/Cryp.trees.pdf", width = 14, height = 7)
par(mfrow=c(2,4))
plot.phylo(ladderize(uce.300.iqt), use.edge.length = F)
plot.phylo(ladderize(uce.300.ast), use.edge.length = F)
plot.phylo(ladderize(uce.1000.iqt), use.edge.length = F)
plot.phylo(ladderize(uce.1000.ast), use.edge.length = F)
plot.phylo(ladderize(uce.1000.ChrZ.iqt), use.edge.length = F)
plot.phylo(ladderize(uce.1000.ChrZ.ast), use.edge.length = F)
plot.phylo(ladderize(cds.iqt), use.edge.length = F)
plot.phylo(ladderize(cds.ast), use.edge.length = F)
dev.off()
