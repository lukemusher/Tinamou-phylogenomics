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

#function to get variance in gc content proportion 

get_gc_content<-function(x){
  x<-read.nexus.data(x)
  site<-c()
  counts=0
  gc.cont<-c()
  for(i in 1:length(x)){
    g<-sum(x[[i]]=="g")
    c<-sum(x[[i]]=="c")
    t<-sum(x[[i]]=="t")
    a<-sum(x[[i]]=="a")
    gc.cont[i]<-(g+c)/(g+c+a+t)
  }
  gc.cont.var<-var(gc.cont)
  return(gc.cont.var)
  
}

#function to build table of variables for genes and gene trees including rf distances for both bifurcating and collapsed low-support nodes
get_rf_gc_pis_table<-function(files,trees.coll,trees,spec_tree,out.path){
  aln<-c()
  pars.site<-c()
  tot.sites<-c()
  #rf.distance<-c()
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
    
    #tryCatch(rf.distance[i]<-as.numeric(RF.dist(read.tree(trees[i]),spec_tree), error = function(e) { skip_to_next <<- TRUE}))
    
    #if(skip_to_next) { next }
    
    tryCatch(rf.distance.coll[i]<-as.numeric(RF.dist(read.tree(trees.coll[i]),spec_tree), error = function(e) { skip_to_next <<- TRUE}))
    
    if(skip_to_next) { next }
    
  }
  
  #df<-as.data.frame(cbind(aln,pars.site,tot.sites,pars.site.prop,rf.distance,rf.distance.coll, gc.cont))
  df<-as.data.frame(cbind(aln,pars.site,tot.sites,pars.site.prop,rf.distance.coll, gc.cont))
  write.csv(df,out.path)
  return(df)
  
}

#set working directory
setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/cds-100p-gts-collapsed/")

#define files, trees, and collapsed trees
files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*treefile",full.names=F,recursive=FALSE)

#check to make sure same numbers of trees and alignments
length(files)
length(trees.coll)

#define species tree 
spec_tree <- root(read.tree(file = "../../Trees/cds.100p.astral.tre"), "rhea_penn", resolve.root = T)
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

mean(df.CDS$tot.sites) #1364.601
sd(df.CDS$tot.sites) #1011.409

mean(df.CDS$pars.site) #529.5792
sd(df.CDS$pars.site) #393.2861

mean(df.CDS$pars.site.prop) #0.3982965
sd(df.CDS$pars.site.prop) #0.1211765

mean(df.CDS$rf.distance.coll) #59.23043
sd(df.CDS$rf.distance.coll) #20.66592

mean(df.CDS$gc.cont) #0.0001931018
sd(df.CDS$gc.cont) #0.0003469101

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/uce1000-100p-gts-collapsed/")

files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*.treefile",full.names=F,recursive=FALSE)

length(trees.coll)
length(files)

#####commented code below shows that sample of gene-to-gene-tree comparisons is similar to gene-to-species-tree comparisons
# rf.uce1000<-c()
# counts=0
# for (i in 1:1000){
#   for (j in 1001:1300){
#     if(i!=j){
#       counts=counts+1
#       tree1<-read.tree(trees.coll[i])
#       tree2<-read.tree(trees.coll[j])
#       rf.uce1000[counts]<-RF.dist(tree1,tree2)
#     }
#   }
# }
# 
# length(rf)
# mean(rf)
# sd(rf)

spec_tree <- root(read.tree(file = "../../Trees/uce1000.full.75p.astral.tre"), "rhea_penn", resolve.root = T)
plot(spec_tree)

#df.uce1000<-get_rf_gc_pis_table(files,trees.coll,trees,spec_tree,out.path="../uce1000.collapsed.csv")
df.uce1000<-read.csv("../uce1000.collapsed.csv")

mean(df.uce1000$pars.site) #682.5251
sd(df.uce1000$pars.site) #164.1699
mean(df.uce1000$pars.site.prop) #0.3892043
sd(df.uce1000$pars.site.prop) #0.08824309
mean(df.uce1000$rf.distance.coll) #38.45586
sd(df.uce1000$rf.distance.coll) #9.698341
mean(df.uce1000$gc.cont) #8.983829e-05
sd(df.uce1000$gc.cont) #0.0001145198

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/uce300-100p-gts-collapsed/")

files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*.treefile",full.names=F,recursive=FALSE)

length(trees.coll)
length(files)

spec_tree <- root(read.tree(file = "../../Trees/uce300.full.75p.astral.tre"), "rhea_penn", resolve.root = T)
plot(spec_tree)

#df.uce300<-get_rf_gc_pis_table(files,trees.coll,trees,spec_tree,out.path="../uce300.collapsed.csv")
df.uce300<-read.csv("../uce300.collapsed.csv")

mean(df.uce300$pars.site) #165
sd(df.uce300$pars.site) #62.24948
mean(df.uce300$pars.site.prop) #0.2419654
sd(df.uce300$pars.site.prop) #0.09586637
mean(df.uce300$rf.distance.coll) #72.6188
sd(df.uce300$rf.distance.coll) #19.8526
mean(df.uce300$gc.cont) #6.935002e-05
sd(df.uce300$gc.cont) #9.554192e-05

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/RF.distances/uce100-100p-gts-collapsed/")
files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees.coll <- list.files(path="./",pattern="*.treefile",full.names=F,recursive=FALSE)

length(trees.coll)
length(files)

spec_tree <- root(read.tree(file = "../../Trees/uce1000.full.100p.astral.tre"), "rhea_penn", resolve.root = T)
plot(spec_tree)

#df.uce100<-get_rf_gc_pis_table(files,trees.coll,trees,spec_tree,out.path="../uce100.collapsed.csv")
df.uce100<-read.csv("../uce100.collapsed.csv")

mean(df.uce100$pars.site) #41.07301
sd(df.uce100$pars.site) #23.55901
mean(df.uce100$pars.site.prop) #0.1285544
sd(df.uce100$pars.site.prop) #0.07501452
mean(na.omit(df.uce100$rf.distance.coll)) #123.174
sd(na.omit(df.uce100$rf.distance.coll)) #20.24933
mean(na.omit(df.uce100$gc.cont)) #4.887498e-05
sd(na.omit(df.uce100$gc.cont)) #8.040409e-05

###########################
###Kruskall Wallace Test###
###########################
library(FSA)
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
# Kruskal-Wallis chi-squared = 3510.1, df = 3, p-value < 2.2e-16

dunnTest(as.numeric(PIS)~dataset,data=df.PIS.RF)
# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Holm method.
# 
# Comparison         Z       P.unadj         P.adj
# 1      BUSCOs - UCE1000Flank -16.29205  1.123945e-59  1.123945e-59
# 2       BUSCOs - UCE100Flank  65.42887  0.000000e+00  0.000000e+00
# 3 UCE1000Flank - UCE100Flank  85.01284  0.000000e+00  0.000000e+00
# 4       BUSCOs - UCE300Flank  32.08076 8.180616e-226 1.636123e-225
# 5 UCE1000Flank - UCE300Flank  50.35382  0.000000e+00  0.000000e+00
# 6  UCE100Flank - UCE300Flank -34.85645 3.400439e-266 1.020132e-265
# Warning message:
#   dataset was coerced to a factor. 

kruskal.test(PIS.prop~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  PIS.prop by dataset
# Kruskal-Wallis chi-squared = 6222.1, df = 3, p-value < 2.2e-16

dunnTest(as.numeric(PIS.prop)~dataset,data=df.PIS.RF)
# Comparison           Z       P.unadj         P.adj
# 1      BUSCOs - UCE1000Flank  -0.8777843  3.800608e-01  3.800608e-01
# 2       BUSCOs - UCE100Flank  64.2621533  0.000000e+00  0.000000e+00
# 3 UCE1000Flank - UCE100Flank  67.6725898  0.000000e+00  0.000000e+00
# 4       BUSCOs - UCE300Flank  36.3192892 8.026995e-289 2.408098e-288
# 5 UCE1000Flank - UCE300Flank  38.6427948  0.000000e+00  0.000000e+00
# 6  UCE100Flank - UCE300Flank -29.1983902 2.032652e-187 4.065304e-187
# Warning message:
#   dataset was coerced to a factor.

kruskal.test(RF.coll~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  RF.coll by dataset
# Kruskal-Wallis chi-squared = 3830.5, df = 3, p-value < 2.2e-16

dunnTest(as.numeric(RF.coll)~dataset,data=df.PIS.RF)

# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Holm method.
# 
# Comparison          Z   P.unadj     P.adj
# 1      BUSCOs - UCE1000Flank  28.09590 1.099270e-173 2.198540e-173
# 2       BUSCOs - UCE100Flank -53.86013  0.000000e+00  0.000000e+00
# 3 UCE1000Flank - UCE100Flank -85.34613  0.000000e+00  0.000000e+00
# 4       BUSCOs - UCE300Flank -16.40912  1.645828e-60  1.645828e-60
# 5 UCE1000Flank - UCE300Flank -46.41495  0.000000e+00  0.000000e+00
# 6  UCE100Flank - UCE300Flank  39.16224  0.000000e+00  0.000000e+00
# Warning message:
#   dataset was coerced to a factor. 

kruskal.test(RF.coll~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  RF by dataset
# Kruskal-Wallis chi-squared = 3830.5, df = 3, p-value < 2.2e-16

dunnTest(as.numeric(RF.coll)~dataset,data=df.PIS.RF)

# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Holm method.
# 
# Comparison         Z       P.unadj         P.adj
# 1      BUSCOs - UCE1000Flank  28.09590 1.099270e-173 2.198540e-173
# 2       BUSCOs - UCE100Flank -53.86013  0.000000e+00  0.000000e+00
# 3 UCE1000Flank - UCE100Flank -85.34613  0.000000e+00  0.000000e+00
# 4       BUSCOs - UCE300Flank -16.40912  1.645828e-60  1.645828e-60
# 5 UCE1000Flank - UCE300Flank -46.41495  0.000000e+00  0.000000e+00
# 6  UCE100Flank - UCE300Flank  39.16224  0.000000e+00  0.000000e+00
# Warning message:
#   dataset was coerced to a factor. 

kruskal.test(GC~dataset,data=df.PIS.RF)

# Kruskal-Wallis rank sum test
# 
# data:  RF by dataset
# Kruskal-Wallis chi-squared = 229.01, df = 3, p-value < 2.2e-16

dunnTest(as.numeric(GC)~dataset,data=df.PIS.RF)
# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Holm method.
# 
# Comparison         Z       P.unadj         P.adj
# 1      BUSCOs - UCE1000Flank  16.93067  2.673314e-64  8.019942e-64
# 2       BUSCOs - UCE100Flank  40.50733  0.000000e+00  0.000000e+00
# 3 UCE1000Flank - UCE100Flank  24.36320 4.201024e-131 1.680409e-130
# 4       BUSCOs - UCE300Flank  27.72719 3.283363e-169 1.641681e-168
# 5 UCE1000Flank - UCE300Flank  11.10166  1.231360e-28  1.231360e-28
# 6  UCE100Flank - UCE300Flank -13.34457  1.274118e-40  2.548236e-40
# Warning message:
#   dataset was coerced to a factor.

###########################
#FITTING LINEAR MODELS#
###########################

z.uces<-read.table("~/Documents/ANSDU/Tinamous/z.chrom.list.txt")$V1

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

summary(fit.glm) #AIC: 19897 Lowest AIC
summary(fit.glog) #AIC: 20794

fit.lm <- lm(rf.distance.coll~pars.site, data = df.uce100)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.uce100)
summary(fit.lm) 
summary(fit.log)

uce100.pis.plot2<-ggplot(df.uce100,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+ylim(0,155)+
  geom_point(data=df.uce100[df.uce100$chrom=="Z-chromosome",], mapping =aes(x=pars.site, y=rf.distance.coll), color="yellow", size=1)

uce100.pis.plot2

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.uce100)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.uce100)

summary(fit.glm) #AIC: 19959 Lowest AIC
summary(fit.glog) #AIC: 20843

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.uce100)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.uce100)
summary(fit.lm)
summary(fit.log)

uce100.prop.plot2<-ggplot(df.uce100,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 100 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+ylim(0,155)+
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
summary(fit.glm) #AIC: 21092
summary(fit.glog) #AIC: 20331 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site, data = df.uce300)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.uce300)
summary(fit.lm) 
summary(fit.log)

uce300.pis.plot2<-ggplot(df.uce300,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+ylim(0,155)+
  geom_point(data=df.uce300[df.uce300$chrom=="Z-chromosome",], mapping =aes(x=pars.site, y=rf.distance.coll), color="yellow", size=1)

uce300.pis.plot2

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.uce300)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.uce300)

summary(fit.glm) #AIC: 21320
summary(fit.glog) #AIC: 20523 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.uce300)
fit.log <- lm(rf.distance.coll~log(pars.site.prop), data = df.uce300)
summary(fit.lm) #AIC: 21092
summary(fit.log) #AIC: 20047 Lowest AIC

uce300.prop.plot2<-ggplot(df.uce300,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+ylim(0,155)+
  geom_point(data=df.uce300[df.uce300$chrom=="Z-chromosome",], mapping =aes(x=pars.site.prop, y=rf.distance.coll), color="yellow", size=1)

uce300.prop.plot2

fit.glm <- glm(rf.distance.coll~gc.cont, data = df.uce300)
fit.glog <- glm(rf.distance.coll~log(gc.cont), data = df.uce300)

summary(fit.glm) #AIC: 23704
summary(fit.glog) #AIC: 23120 lowest

fit.lm <- lm(rf.distance.coll~gc.cont, data = df.uce300)
fit.log <- lm(rf.distance.coll~log(gc.cont), data = df.uce300)

uce300.gc.plot2<-ggplot(df.uce300,aes(x=gc.cont,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 300 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Variance in GC Content") + ylab("RF Distance")+ylim(0,155)+
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
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.uce1000)

summary(fit.glm) #AIC: 19035
summary(fit.glog) #AIC: 18927 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.uce1000)
fit.log <- lm(rf.distance~log(pars.site.prop), data = df.uce1000)
summary(fit.lm) 
summary(fit.log)

uce1000.prop.plot2<-ggplot(df.uce1000,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+ylim(0,155)+
  geom_point(data=df.uce1000[df.uce1000$chrom=="Z-chromosome",], mapping =aes(x=pars.site.prop, y=rf.distance.coll), color="yellow", size=1)

uce1000.prop.plot2

fit.glm <- glm(rf.distance.coll~pars.site, data = df.uce1000)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.uce1000)

summary(fit.glm) #AIC: 18683
summary(fit.glog) #AIC: 18506 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site, data = df.uce1000)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.uce1000)
summary(fit.lm)
summary(fit.log)

uce1000.pis.plot2<-ggplot(df.uce1000,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+ylim(0,155)+
  geom_point(data=df.uce1000[df.uce1000$chrom=="Z-chromosome",], mapping =aes(x=pars.site, y=rf.distance.coll), color="yellow", size=1)

uce1000.pis.plot2

fit.glm <- glm(rf.distance.coll~gc.cont, data = df.uce1000)
fit.glog <- glm(rf.distance.coll~log(gc.cont), data = df.uce1000)

summary(fit.glm) #AIC: 19358 Lowest AIC
summary(fit.glog) #AIC: 19391

fit.lm <- lm(rf.distance.coll~gc.cont, data = df.uce1000)
#fit.log <- lm(rf.distance.coll~log(gc.cont), data = df.uce300)

uce1000.gc.plot2<-ggplot(df.uce1000,aes(x=gc.cont,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("UCE 1000 Flank", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("GC Content") + ylab("RF Distance")+ylim(0,155)+
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
df.CDS$rf.distance.coll<-as.numeric(df.CDS$rf.distance.coll)

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.CDS)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.CDS)

summary(fit.glm) #AIC: 20228
summary(fit.glog) #AIC: 20208 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site.prop, data = df.CDS)
fit.log <- lm(rf.distance.coll~log(pars.site.prop), data = df.CDS)
summary(fit.lm) 
summary(fit.log)

CDS.prop.plot2<-ggplot(df.CDS,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+ylim(0,155)

CDS.prop.plot2

fit.glm <- glm(rf.distance.coll~pars.site, data = df.CDS)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.CDS)

summary(fit.glm) #AIC: 19291
summary(fit.glog) #AIC: 18669 Lowest AIC

fit.lm <- lm(rf.distance.coll~pars.site, data = df.CDS)
fit.log <- lm(rf.distance.coll~log(pars.site), data = df.CDS)
summary(fit.lm) 
summary(fit.log)

CDS.pis.plot2<-ggplot(df.CDS,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# of informative sites") + ylab("RF Distance")+ylim(0,155)

CDS.pis.plot2

fit.glm <- glm(rf.distance.coll~gc.cont, data = df.CDS)
fit.glog <- glm(rf.distance.coll~log(gc.cont), data = df.CDS)

summary(fit.glm) #AIC: 20197 Lowest AIC
summary(fit.glog) #AIC: 20210

fit.lm <- lm(rf.distance.coll~gc.cont, data = df.CDS)
fit.log <- lm(rf.distance.coll~log(gc.cont), data = df.CDS)

CDS.gc.plot2<-ggplot(df.CDS,aes(x=gc.cont,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~x, color="red")+
  ggtitle("CDS", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("GC Content") + ylab("RF Distance")+ylim(0,155)

CDS.gc.plot2

####Combine UCEs###
dataset<-c(rep("uce 100", length(df.uce100[,1])),rep("uce 300", length(df.uce300[,1])),rep("uce 1000", length(df.uce1000[,1])))

df.uce<-rbind(df.uce100,df.uce300,df.uce1000)
df.uce<-cbind(dataset,df.uce)
write.csv(df.uce,"../uces.combined.csv")

####for each duplicate value in column "aln" take random alignment (i.e., for each locus use only 100, 300, or 1000 rather than keeping all three)
library(plyr)
subsampled_data <- ddply(df.uce,.(aln),
                         function(x) {
                           x[sample(nrow(x),size=1),]
                         })

df.uce<-subsampled_data
df.uce[df.uce$aln=="uce-1177.nexus",] #confirm only represented once

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.uce)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.uce)

summary(fit.glm) #AIC: 23630
summary(fit.glog) #AIC: 23311 Lowest AIC

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
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")+ylim(0,155)+
  theme(legend.position = "none")+
  geom_point(data=df.uce[df.uce$chrom=="Z-chromosome",], mapping =aes(x=pars.site.prop, y=rf.distance.coll), color="Yellow", size=1)

uce.prop.plot2

fit.glm <- glm(rf.distance.coll~pars.site, data = df.uce)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.uce)

summary(fit.glm) #AIC: 24592
summary(fit.glog) #AIC: 20767 Lowest AIC

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
  theme_bw() + xlab("# informative sites") + ylab("RF Distance")+ylim(0,155)+
  theme(legend.position = c(0.75,0.75),legend.background = element_rect(fill="white",
                                                                        size=0.5, linetype="solid", 
                                                                        colour ="black"))+
  geom_point(data=df.uce[df.uce$chrom=="Z-chromosome",], mapping =aes(x=pars.site, y=rf.distance.coll), color="yellow", size=1)

uce.pis.plot2


#plot(Figure 6)

pdf(file = "../../figs_30Oct24/cors.pdf", width = 14, height = 7)
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

b<-ggplot(df.all, aes(x=dataset, y=log(pars.site))) + 
  geom_violin()+ geom_boxplot(width=0.1)+ggtitle("PIS by dataset")+
  theme_bw() + xlab("") + ylab("Number of PIS")

c<-ggplot(df.all, aes(x=dataset, y=pars.site.prop)) + 
  geom_violin()+ geom_boxplot(width=0.1)+ggtitle("Percent PIS by dataset")+
  theme_bw() + xlab("") + ylab("Proportion PIS")

d<-ggplot(df.all, aes(x=dataset, y=log(gc.cont))) + 
  geom_violin()+ geom_boxplot(width=0.1)+ggtitle("GC content variation by dataset")+
  theme_bw() + xlab("") + ylab("Variance in GC content")

names(df.uce)

#does RF of Z-chromosome UCEs differ from RF of Autosomal UCEs?
df.uce$dataset<-factor(df.uce$dataset, levels = c("uce 100", "uce 300", "uce 1000"))
e<-ggplot(df.all[df.all$dataset!="BUSCO",], aes(x=dataset, y=rf.distance.coll)) + 
  geom_violin(aes(color = chrom), trim = FALSE, position = position_dodge(0.9) ) +
  geom_boxplot(aes(color = chrom), width = 0.1, position = position_dodge(0.9)) +
  scale_color_manual(values = c("#353436","#1b98e0"))+
  theme_bw() + xlab("") + ylab("RF Distance")+ggtitle("RF Distance by chromosome")

e

mean(df.uce1000$rf.distance.coll[df.uce1000$chrom=="Autosome"]) #38.58408
sd(df.uce1000$rf.distance.coll[df.uce1000$chrom=="Autosome"]) #9.673809

mean(df.uce1000$rf.distance.coll[df.uce1000$chrom=="Z-chromosome"]) #32.22642
sd(df.uce1000$rf.distance.coll[df.uce1000$chrom=="Z-chromosome"]) #8.889582

kruskal.test(rf.distance.coll~chrom,data=df.uce1000)

# Kruskal-Wallis rank sum test
# 
# data:  rf.distance.coll by chrom
# Kruskal-Wallis chi-squared = 26.645, df = 1, p-value = 2.445e-07
fligner.test(rf.distance.coll~chrom,data=df.uce1000)

kruskal.test(rf.distance.coll~chrom,data=df.uce300)

# Kruskal-Wallis rank sum test
# 
# data:  rf.distance.coll by chrom
# Kruskal-Wallis chi-squared = 0.36594, df = 1, p-value = 0.5452

kruskal.test(rf.distance.coll~chrom,data=df.uce100)

# Kruskal-Wallis rank sum test
# 
# data:  rf.distance.coll by chrom
# Kruskal-Wallis chi-squared = 0.02013, df = 1, p-value = 0.8872

###Concordance factors boxplots
setwd("~/Documents/ANSDU/Tinamous/concordance_factors/")

cfs<-na.omit(read.csv("CFS.all.datasets.csv"))
unique(cfs$Dataset)
cfs$Dataset<-factor(cfs$Dataset, levels = c("BUSCOs","UCE100", "UCE300", "UCE1000-Autosomes", "UCE1000-ChrZ"))

kruskal.test(cfs$gCF[cfs$Method=="astral"]~cfs$Dataset[cfs$Method=="astral"])
kruskal.test(cfs$gCF[cfs$Dataset=="UCE1000-Autosomes"]~cfs$Method[cfs$Dataset=="UCE1000-Autosomes"])

mean(cfs$gCF[cfs$Dataset=="UCE1000-Autosomes" & cfs$Method=="iqtree"])
mean(na.omit(cfs$gCF[cfs$Dataset=="UCE1000-Autosomes" & cfs$Method=="astral"]))

dunnTest(cfs$gCF[cfs$Method=="astral"]~cfs$Dataset[cfs$Method=="astral"])

# f<-ggplot(cfs[cfs$Method=="astral",], aes(x=Dataset, y=gCF)) + 
#   geom_violin(aes(color = Method), trim = FALSE, position = position_dodge(0.9) ) +
#   geom_boxplot(aes(color = Method), width = 0.1, position = position_dodge(0.9)) +
#   scale_color_manual(values = c("#353436","#ca7dcc"))+
#   theme_bw() + xlab("") + ylab("gCF")+ggtitle("Condcordance factors")


f<-ggplot(cfs[cfs$Method=="astral",], aes(x=Dataset, y=gCF)) + 
  geom_violin() +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  theme_bw() + xlab("") + ylab("gCF")+ggtitle("Condcordance factors")

f
pdf(file = "../figs_30Oct24/violin_plots.pdf", width = 13, height = 8)
ggarrange(b,c,d,a,e,f,
          ncol = 2, nrow = 3)
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

summary(fit.glm) #AIC: 19712
summary(fit.glog) #AIC: 19688 Lowest AIC

fit.lm <- lm(rf.distance.coll~log(pars.site.prop), data = df.CDS.filt.gc)

CDS.prop.plot3<-ggplot(df.CDS.filt.gc,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs filtered by GC content", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot3

fit.glm <- glm(rf.distance.coll~pars.site, data = df.CDS.filt.gc)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.CDS.filt.gc)

summary(fit.glm) #AIC: 18778
summary(fit.glog) #AIC: 18140 Lowest AIC

fit.log <- lm(rf.distance.coll~log(pars.site), data = df.CDS.filt.gc)

CDS.pis.plot4<-ggplot(df.CDS.filt.gc,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs filtered by GC content", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# of informative sites") + ylab("RF Distance")

CDS.pis.plot4

#CDS with comparable prop variable sites to all UCEs datasets
df.CDS.filt1<-df.CDS[df.CDS$pars.site.prop>=f & df.CDS$pars.site.prop<=a,]

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.CDS.filt1)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.CDS.filt1)

summary(fit.glm) #AIC: 18274
summary(fit.glog) #AIC: 18247 Lowest AIC

fit.lm <- lm(rf.distance~log(pars.site.prop), data = df.CDS.filt1)

CDS.prop.plot6<-ggplot(df.CDS.filt1,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs: UCE-like filtered by PIS", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot6

fit.glm <- glm(rf.distance.coll~pars.site, data = df.CDS.filt1)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.CDS.filt1)

summary(fit.glm) #AIC: 17342
summary(fit.glog) #AIC: 16702 Lowest AIC

fit.log <- lm(rf.distance.coll~log(pars.site), data = df.CDS.filt1)

CDS.pis.plot7<-ggplot(df.CDS.filt1,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs: UCE-like filtered by PIS", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# of informative sites") + ylab("RF Distance")

CDS.pis.plot7

#filter by UCE and gc content
df.CDS.filt2<-df.CDS.filt1[df.CDS.filt1$gc.cont<=g & df.CDS.filt1$gc.cont>=h,]

fit.glm <- glm(rf.distance.coll~pars.site.prop, data = df.CDS.filt2)
fit.glog <- glm(rf.distance.coll~log(pars.site.prop), data = df.CDS.filt2)

summary(fit.glm) #AIC: 18007
summary(fit.glog) #AIC: 17980 Lowest AIC

fit.lm <- lm(rf.distance.coll~log(pars.site.prop), data = df.CDS.filt2)

CDS.prop.plot8<-ggplot(df.CDS.filt2,aes(x=pars.site.prop,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs: UCE-like filtered by GC & PIS", subtitle = paste("Adj. R-squared = ",round(summary(fit.lm)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.lm))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("Proportion of informative sites") + ylab("RF Distance")

CDS.prop.plot8

fit.glm <- glm(rf.distance.coll~pars.site, data = df.CDS.filt2)
fit.glog <- glm(rf.distance.coll~log(pars.site), data = df.CDS.filt2)

summary(fit.glm) #AIC: 17096
summary(fit.glog) #AIC: 16459 Lowest AIC

fit.log <- lm(rf.distance.coll~log(pars.site), data = df.CDS.filt2)

CDS.pis.plot9<-ggplot(df.CDS.filt2,aes(x=pars.site,y=rf.distance.coll)) +
  geom_point(alpha=0.6) +
  geom_smooth(method=lm, se = T ,formula = y~log(x), color="red")+
  ggtitle("BUSCOs: UCE-like filtered by GC & PIS", subtitle = paste("Adj. R-squared = ",round(summary(fit.log)$adj.r.squared,digits = 2), ", P = ",round(coef(summary(fit.log))[2,4],digits = 4),sep = ""))+
  theme_bw() + xlab("# of informative sites") + ylab("RF Distance")

CDS.pis.plot9


pdf(file = "../figs_30Oct24/Filtered_CDS_plots.pdf", width = 18, height = 7)
ggarrange(CDS.prop.plot3,CDS.prop.plot6,CDS.prop.plot8,CDS.pis.plot4,CDS.pis.plot7,CDS.pis.plot9,
          ncol = 3, nrow = 2)
dev.off()
par(mfrow=c(1,1))
####################
#####RF Heatmap#####
####################
setwd("~/Documents/ANSDU/Tinamous/Trees/")
spec<-c("cuyap42937","cuver35355","cusim651488")
uce.100.iqt<-drop.tip(read.tree("crypt.uce100.full.100p.iqtree.treefile"),spec)
uce.100.ast<-drop.tip(read.tree("crypt.uce100.full.100p.astral.tre"),spec)
uce.100.iqt.75<-drop.tip(read.tree("crypt.uce100.full.75p.iqtree.treefile"),spec)
uce.100.ast.75<-drop.tip(read.tree("crypt.uce100.full.75p.astral.tre"),spec)
plot(uce.100.iqt.75)
uce.300.iqt<-drop.tip(read.tree("crypt.uce300.full.100p.iqtree.treefile"),spec)
uce.300.ast<-drop.tip(read.tree("crypt.uce300.full.100p.astral.tre"),spec)
uce.300.iqt.75<-drop.tip(read.tree("crypt.uce300.full.75p.iqtree.treefile"),spec)
uce.300.ast.75<-drop.tip(read.tree("crypt.uce300.full.75p.astral.tre"),spec)
plot(uce.300.iqt)
uce.1000.iqt<-drop.tip(read.tree("crypt.uce1000.full.100p.iqtree.treefile"),spec)
uce.1000.ast<-drop.tip(read.tree("crypt.uce1000.full.100p.astral.tre"),spec)
uce.1000.iqt.75<-drop.tip(read.tree("crypt.uce1000.full.75p.iqtree.treefile"),spec)
uce.1000.ast.75<-drop.tip(read.tree("crypt.uce1000.full.75p.astral.tre"),spec)
plot(uce.1000.iqt.75)
uce.1000.ChrZ.iqt<-drop.tip(read.tree("crypt.uce1000.chrz.75p.iqtree.treefile"),spec)
uce.1000.ChrZ.ast<-drop.tip(read.tree("crypt.uce1000.ChrZ.75p.astral.tre"),spec)
plot(uce.1000.ChrZ.iqt)
uce.1000.autosomes.iqt<-drop.tip(read.tree("crypt.uce1000.autosomes.75p.iqtree.treefile"),spec)
uce.1000.autosomes.ast<-drop.tip(read.tree("crypt.uce1000.autosomes.75p.astral.tre"),spec)
plot(uce.1000.autosomes.iqt)

cds.iqt<-drop.tip(read.tree("crypt.cds.100p.iqtree.treefile"),spec)
cds.ast<-drop.tip(read.tree("crypt.cds.100p.astral.tre"),spec)
cds.iqt.75<-drop.tip(read.tree("crypt.cds.75p.iqtree.treefile"),spec)
cds.ast.75<-drop.tip(read.tree("crypt.cds.75p.astral.tre"),spec)
plot(cds.ast.75)
cds.ast.75$tip.label
RF.dist(cds.ast.75,uce.1000.iqt.75)
#full
# trees<-c(uce.100.iqt,uce.100.ast,uce.100.iqt.75,uce.100.ast.75,uce.300.iqt,uce.300.ast,uce.300.iqt.75,uce.300.ast.75,uce.1000.iqt,uce.1000.ast,uce.1000.iqt.75,uce.1000.ast.75,uce.1000.ChrZ.iqt,uce.1000.ChrZ.ast,cds.iqt,cds.ast,cds.iqt.75,cds.ast.75)
# trees.list<-c("uce.100.iqt","uce.100.ast","uce.100.iqt.75","uce.100.ast.75","uce.300.iqt","uce.300.ast","uce.300.iqt.75","uce.300.ast.75","uce.1000.iqt","uce.1000.ast","uce.1000.iqt.75","uce.1000.ast.75","uce.1000.ChrZ.iqt","uce.1000.ChrZ.ast","cds.iqt","cds.ast","cds.iqt.75","cds.ast.75")

#75% only
trees<-c(uce.100.iqt.75,uce.100.ast.75,uce.300.iqt.75,uce.300.ast.75,uce.1000.iqt.75,uce.1000.ast.75,uce.1000.autosomes.iqt,uce.1000.autosomes.ast,uce.1000.ChrZ.iqt,uce.1000.ChrZ.ast,cds.iqt.75,cds.ast.75)
trees.list<-c("uce.100.iqt.75","uce.100.ast.75","uce.300.iqt.75","uce.300.ast.75","uce.1000.iqt.75","uce.1000.ast.75","uce.1000.autosomes.iqt","uce.1000.autosomes.ast","uce.1000.ChrZ.iqt","uce.1000.ChrZ.ast","cds.iqt.75","cds.ast.75")

trees1<-c()
trees2<-c()
RF<-c()
count=0

for (i in 1:length(trees)){
  for (j in 1:length(trees)){
    count=count+1
    trees1[count]<-trees.list[i]
    trees2[count]<-trees.list[j]
    print(trees1[count])
    print(trees2[count])
    RF[count]<-RF.dist(tree1=trees[i][[1]],tree2=trees[j][[1]])
  }
}


tab<-data.frame(cbind(trees1,trees2,as.numeric(RF)))

tab$trees1<-factor(trees1, levels = c("cds.iqt.75","uce.100.iqt.75","uce.300.iqt.75","uce.1000.iqt.75","uce.1000.autosomes.iqt","uce.1000.ChrZ.iqt","cds.ast.75","uce.100.ast.75","uce.300.ast.75","uce.1000.ast.75","uce.1000.autosomes.ast","uce.1000.ChrZ.ast"))
tab$trees2<-factor(trees2, levels = c("cds.iqt.75","uce.100.iqt.75","uce.300.iqt.75","uce.1000.iqt.75","uce.1000.autosomes.iqt","uce.1000.ChrZ.iqt","cds.ast.75","uce.100.ast.75","uce.300.ast.75","uce.1000.ast.75","uce.1000.autosomes.ast","uce.1000.ChrZ.ast"))

library(hrbrthemes)


##Plot RF distances heatmap

pdf(file = "../figs_30Oct24/heatmap.raw.pdf", width = 7, height = 7)
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

support.tab<-na.omit(read.csv("../RF.distances/Topologies_Support.csv"))
unique(support.tab$Metric)

a<-(support.tab$Value[support.tab$Metric=="Bootstrap/PP"]-min(support.tab$Value[support.tab$Metric=="Bootstrap/PP"]))/(max(support.tab$Value[support.tab$Metric=="Bootstrap/PP"])-min(support.tab$Value[support.tab$Metric=="Bootstrap/PP"]))*100
b<-(support.tab$Value[support.tab$Metric=="gCF"]-min(support.tab$Value[support.tab$Metric=="gCF"]))/(max(support.tab$Value[support.tab$Metric=="gCF"])-min(support.tab$Value[support.tab$Metric=="gCF"]))*100
c<-100-(support.tab$Value[support.tab$Metric=="100-rf.t1"]-min(support.tab$Value[support.tab$Metric=="100-rf.t1"]))/(max(support.tab$Value[support.tab$Metric=="100-rf.t1"])-min(support.tab$Value[support.tab$Metric=="100-rf.t1"]))*100
d<-100-(support.tab$Value[support.tab$Metric=="100-rf.t2"]-min(support.tab$Value[support.tab$Metric=="100-rf.t2"]))/(max(support.tab$Value[support.tab$Metric=="100-rf.t2"])-min(support.tab$Value[support.tab$Metric=="100-rf.t2"]))*100
e<-(support.tab$Value[support.tab$Metric=="gCF.cladeA"]-min(support.tab$Value[support.tab$Metric=="gCF.cladeA"]))/(max(support.tab$Value[support.tab$Metric=="gCF.cladeA"])-min(support.tab$Value[support.tab$Metric=="gCF.cladeA"]))*100
f<-(support.tab$Value[support.tab$Metric=="Bootstrap/PP.cladeA"]-min(support.tab$Value[support.tab$Metric=="Bootstrap/PP.cladeA"]))/(max(support.tab$Value[support.tab$Metric=="Bootstrap/PP.cladeA"])-min(support.tab$Value[support.tab$Metric=="Bootstrap/PP.cladeA"]))*100
Value2<-round(c(a,b,c,d,e,f),digits=0)
support.tab$Value2<-Value2
support.tab$Value[support.tab$Metric=="100-rf.t1"]<-100-support.tab$Value[support.tab$Metric=="100-rf.t1"]
support.tab$Value[support.tab$Metric=="100-rf.t2"]<-100-support.tab$Value[support.tab$Metric=="100-rf.t2"]
support.tab$Value[support.tab$Metric=="gCF"]<-round(support.tab$Value[support.tab$Metric=="gCF"],digits = 0)
support.tab$Value[support.tab$Metric=="gCF.cladeA"]<-round(support.tab$Value[support.tab$Metric=="gCF.cladeA"],digits = 0)

support.tab$Tree<-factor(support.tab$Tree, levels = c("uce.100.iqt","uce.100.iqt.75","uce.300.iqt","uce.300.iqt.75","uce.1000.iqt","uce.1000.iqt.75","uce.1000.autosomes.iqt","uce.1000.ChrZ.iqt","cds.iqt","cds.iqt.75","uce.100.ast","uce.100.ast.75","uce.300.ast","uce.300.ast.75","uce.1000.ast","uce.1000.ast.75","uce.1000.autosomes.ast","uce.1000.ChrZ.ast","cds.ast","cds.ast.75"))
support.tab$Metric<-factor(support.tab$Metric, levels = c("Bootstrap/PP","Bootstrap/PP.cladeA","gCF","gCF.cladeA","100-rf.t1","100-rf.t2"))

pdf(file = "../figs_30Oct24/heatmap2.pdf", width = 16, height = 4)
ggplot(support.tab, aes(Tree, Metric, fill= as.numeric(Value2))) + 
  geom_tile() +
  scale_fill_gradient2(low = "steelblue3", high = "tomato3" , 
                       midpoint = max(support.tab$Value2)/2, limit = c(min(support.tab$Value2),max(support.tab$Value2)), space = "Lab", 
                       name="normalized score") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+xlab("") + ylab("")+ 
  geom_text(aes(Tree, Metric, label = Value), color = "black", size = 3) +
  
  ggtitle(label = "Species Tree Support for Crypt. Clade A")
dev.off()

###only plot trees/datasets for main text

filt.tree<-factor(c("uce.100.iqt","uce.100.ast","uce.300.iqt","uce.300.ast","uce.100.iqt.75","uce.100.ast.75","uce.300.iqt.75","uce.300.ast.75"))
support.tab2<-support.tab[!(support.tab$Tree%in%filt.tree),]

a<-(support.tab2$Value[support.tab2$Metric=="Bootstrap/PP"]-min(support.tab2$Value[support.tab2$Metric=="Bootstrap/PP"]))/(max(support.tab2$Value[support.tab2$Metric=="Bootstrap/PP"])-min(support.tab2$Value[support.tab2$Metric=="Bootstrap/PP"]))*100
b<-(support.tab2$Value[support.tab2$Metric=="gCF"]-min(support.tab2$Value[support.tab2$Metric=="gCF"]))/(max(support.tab2$Value[support.tab2$Metric=="gCF"])-min(support.tab2$Value[support.tab2$Metric=="gCF"]))*100
c<-100-(support.tab2$Value[support.tab2$Metric=="100-rf.t1"]-min(support.tab2$Value[support.tab2$Metric=="100-rf.t1"]))/(max(support.tab2$Value[support.tab2$Metric=="100-rf.t1"])-min(support.tab2$Value[support.tab2$Metric=="100-rf.t1"]))*100
d<-100-(support.tab2$Value[support.tab2$Metric=="100-rf.t2"]-min(support.tab2$Value[support.tab2$Metric=="100-rf.t2"]))/(max(support.tab2$Value[support.tab2$Metric=="100-rf.t2"])-min(support.tab2$Value[support.tab2$Metric=="100-rf.t2"]))*100
e<-(support.tab2$Value[support.tab2$Metric=="gCF.cladeA"]-min(support.tab2$Value[support.tab2$Metric=="gCF.cladeA"]))/(max(support.tab2$Value[support.tab2$Metric=="gCF.cladeA"])-min(support.tab2$Value[support.tab2$Metric=="gCF.cladeA"]))*100
f<-(support.tab2$Value[support.tab2$Metric=="Bootstrap/PP.cladeA"]-min(support.tab2$Value[support.tab2$Metric=="Bootstrap/PP.cladeA"]))/(max(support.tab2$Value[support.tab2$Metric=="Bootstrap/PP.cladeA"])-min(support.tab2$Value[support.tab2$Metric=="Bootstrap/PP.cladeA"]))*100
Value2<-round(c(a,b,c,d,e,f),digits=0)
support.tab2$Value2<-Value2


pdf(file = "../Figs_August2025/heatmap.pdf", width = 16, height = 4)
ggplot(support.tab2, aes(Tree, Metric, fill= as.numeric(Value2))) + 
  geom_tile() +
  scale_fill_gradient2(low = "steelblue3", high = "tomato3" , 
                       midpoint = max(support.tab$Value2)/2, limit = c(min(support.tab$Value2),max(support.tab$Value2)), space = "Lab", 
                       name="normalized score") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+xlab("") + ylab("")+ 
  geom_text(aes(Tree, Metric, label = Value), color = "black", size = 3) +
  
  ggtitle(label = "Species Tree Support for Crypt. Clade A")
dev.off()

###compare node heights from Z and Autosomes

setwd("~/Documents/ANSDU/Tinamous/ASTRAL/")
outs<-c("cvar27099","crker147057","cryund","crybou","crdui419412","ceery21146","ccgol2213","crnoc59478","cuund34614","cuyap42937","cbar42396","cbcos53413","cscap613347","cssou586295","csnig6048","crysou_gca_013389845","csalbgur081","ccin21636","cpta650391","crber182233","eeele16893","eepat56833","eudele_gca_003342815","effor12668","tpen7597","tiing186371","nbor79764","ndaga68765","nmann614501","nmnig16884","nomin344282","tanan133116ab","nccin2861","rhmac145784","rrruf13855","nobra61359","noorn61422","noros645439","nocur425359","notper_gca_003342845","nppen6601","ntac27158","notorn_gca_013398335","npnie62872","notpen_gca_013398315","rhea_penn_gca_003342835","anodid_gca_006937325","nbbon33138","nncad473927","not_nig_gca_013398345","njul32810","notjul_gca_013398735","tgut85944","tingut_gca_000705375","tiosg311159","tsol25977","ttkle11280","tmmaj20810","tmser456416","tmper40426","cocas11195","coobs25896","cooch113498","copun86899","cparmnt800","cttat35476","cbre625102","crcas434024")
trees<-read.tree("autosomes-cladeA-100p.trees")

plot(drop.tip(trees[[1]],outs))
t.ht<-c()
for (i in 1:length(trees)){
  trees[[i]]<-drop.tip(root(trees[[i]],"cvar27099",resolve.root = T),outs)
  tree<-trees[[i]]
  t.ht[i]<-max(nodeHeights(tree))
}

plot(trees[[7]])

trees2<-read.tree("chrz-cladeA-100p.trees")

t.ht.Z<-c()

for (i in 1:length(trees2)){
  tree<-trees2[[i]]
  tree<-drop.tip(root(tree,"cvar27099",resolve.root = T),outs)
  t.ht.Z[i]<-max(nodeHeights(tree))
}
plot(trees[[1]])

hts<-as.numeric(append(t.ht,t.ht.Z))
chrom<-c(rep("autosome",length(trees)),rep("Z-chrom",length(trees2)))
df.hts<-data.frame(cbind(chrom,hts))
df.hts$hts<-as.numeric(df.hts$hts)
df.hts$chrom<-factor(df.hts$chrom,levels = c("autosome","Z-chrom"))


pdf(file = "../figs_30Oct24/TreeHeights.pdf", width = 7, height = 7)
bp.hts<-ggplot(df.hts, aes(x=chrom, y=hts)) + 
  geom_boxplot(outliers = F)+ggtitle("Tree heights by chromosome type")+
  theme_bw() + xlab("") + ylab("Tree height (subst. per site)")

bp.hts
dev.off()

mean(t.ht.Z)/mean(t.ht)

kruskal.test(hts~chrom,data=df.hts)

outs<-c("cvar27099","crybou","cryund_gca_013389825","crycin_gca_003342915","ccgol2213","crnoc59478","cusim651488","cuver35355","cuund34614","cuyap42937","cbar42396","cbcos53413","cscap613347","cssou586295","csnig6048","crysou_gca_013389845","csalbgur081","ccin21636","cpta650391","crber182233","eeele16893","eepat56833","eudele_gca_003342815","effor12668","tpen7597","tiing186371","nbor79764","ndaga68765","nmann614501","nmnig16884","nomin344282","tanan133116ab","nccin2861","rhmac145784","rrruf13855","nobra61359","noorn61422","noros645439","nocur425359","notper_gca_003342845","nppen6601","ntac27158","notorn_gca_013398335","npnie62872","notpen_gca_013398315","rhea_penn_gca_003342835","anodid_gca_006937325","nbbon33138","nncad473927","not_nig_gca_013398345","njul32810","notjul_gca_013398735","tgut85944","tingut_gca_000705375","tiosg311159","tsol25977","ttkle11280","tmmaj20810","tmser456416","tmper40426","cocas11195","coobs25896","cooch113498","copun86899","cparmnt800","cttat35476","cbre625102","crcas434024")
ins<-c("crker147057","crdui419412","ceery21146","cstr9577")
trees<-root(read.tree("autosomes-cladeA-100p.trees"),"cvar27099")
plot(trees[[1]])

aut.mon.T2<-c()
aut.mon.T3<-c()
aut.mon.T1<-c()
counts=0
for(j in 1:3){
  for (i in 1:length(trees)){
    counts=counts+1
    rand<-sample(length(ins),1)
    tax1<-ins[rand]
    ins.new<-ins[-rand]
    tree<-drop.tip(trees[[i]],c(outs,ins.new))
    trees[[counts]]<-tree
    aut.mon.T2[counts]<-is.monophyletic(tree,c("ctra66583",tax1)) #T2/T4
    aut.mon.T3[counts]<-is.monophyletic(tree,c("ctra66583","cratr320360")) #T3
    aut.mon.T1[counts]<-is.monophyletic(tree,c(tax1,"cratr320360")) #T1
    
  }
}
sum(aut.mon.T2) #T2/T4
sum(aut.mon.T3) #T3
sum(aut.mon.T1) #T1
plot(trees[[1]])

hts.aut.mon.T2<-c()
for (i in 1:length(trees[aut.mon.T2])){
  tree<-drop.tip(trees[aut.mon.T2][[i]],"cratr320360")
  hts.aut.mon.T2[i]<-max(nodeHeights(tree))
}

mean(hts.aut.mon.T2)

hts.aut.mon.T3<-c()
for (i in 1:length(trees[aut.mon.T3])){
  tree<-drop.tip(trees[aut.mon.T3][[i]],ins)
  hts.aut.mon.T3[i]<-max(nodeHeights(tree))
}
mean(hts.aut.mon.T3)

hts.aut.mon.T1<-c()
for (i in 1:length(trees[aut.mon.T1])){
  tree<-drop.tip(trees[aut.mon.T1][[i]],"ctra66583")
  hts.aut.mon.T1[i]<-max(nodeHeights(tree))
}
mean(hts.aut.mon.T1)
tree<-c(rep("T2",length(hts.aut.mon.T2)),rep("T3",length(hts.aut.mon.T3)),rep("T1", length(hts.aut.mon.T1)))
hts<-c(hts.aut.mon.T2,hts.aut.mon.T3,hts.aut.mon.T1)
df.hts2<-data.frame(cbind(tree,hts))
df.hts2$tree<-factor(df.hts2$tree, levels=c("T2","T3","T1"))
df.hts2$hts<-as.numeric(df.hts2$hts)

kruskal.test(hts~tree, data=df.hts2)

# Kruskal-Wallis rank sum test
# 
# data:  hts by tree
# Kruskal-Wallis chi-squared = 41.982, df = 2, p-value = 7.653e-10

dunnTest(hts~tree, data=df.hts2)

# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Holm method.
# 
# Comparison         Z      P.unadj        P.adj
# 1    T1 - T2 -6.277236 3.446457e-10 1.033937e-09
# 2    T1 - T3 -1.583528 1.133011e-01 1.133011e-01
# 3    T2 - T3  4.489111 7.152100e-06 1.430420e-05

bp.hts<-ggplot(df.hts2, aes(x=tree, y=hts)) + 
  geom_boxplot(outliers = F)+ggtitle("Tree heights by chromosome type")+
  theme_bw() + xlab("") + ylab("Tree height (subst. per site)")

bp.hts

trees2<-read.tree("chrz-cladeA-100p.trees")
Z.mon.T2<-c()
Z.mon.T3<-c()
Z.mon.T1<-c()
for (i in 1:length(trees2)){
  rand<-sample(length(ins),1)
  tax1<-ins[rand]
  ins.new<-ins[-rand]
  tree<-drop.tip(trees2[[i]],c(outs,ins.new))
  trees2[[i]]<-tree
  Z.mon.T2[i]<-is.monophyletic(tree,c("ctra66583",tax1))
  Z.mon.T3[i]<-is.monophyletic(tree,c("ctra66583","cratr320360"))
  Z.mon.T1[i]<-is.monophyletic(tree,c(tax1,"cratr320360"))
  
}
sum(aut.mon.T1)
sum(Z.mon.T1)
sum(aut.mon.T2)
sum(Z.mon.T2)
sum(aut.mon.T3)
sum(Z.mon.T3)

hts.Z.mon.T2<-c()
for (i in 1:length(trees2[Z.mon.T2])){
  tree<-drop.tip(trees2[Z.mon.T2][[i]],"cratr320360")
  hts.Z.mon.T2[i]<-max(nodeHeights(tree))
}

mean(hts.Z.mon.T2)

hts.Z.mon.T3<-c()
for (i in 1:length(trees2[Z.mon.T3])){
  tree<-drop.tip(trees2[Z.mon.T3][[i]],"cstr9577")
  hts.Z.mon.T3[i]<-max(nodeHeights(tree))
}
mean(hts.Z.mon.T3)

hts.Z.mon.T1<-c()
for (i in 1:length(trees2[Z.mon.T1])){
  tree<-drop.tip(trees2[Z.mon.T1][[i]],"ctra66583")
  hts.Z.mon.T1[i]<-max(nodeHeights(tree))
}
mean(hts.Z.mon.T1)
plot(tree)
tree<-c(rep("T2",length(hts.Z.mon.T2)),rep("T3",length(hts.Z.mon.T3)),rep("T1", length(hts.Z.mon.T1)))
hts<-c(hts.Z.mon.T2,hts.Z.mon.T3,hts.Z.mon.T1)
df.hts3<-data.frame(cbind(tree,hts))
df.hts3$tree<-factor(df.hts3$tree, levels=c("T2","T3","T1"))
df.hts3$hts<-as.numeric(df.hts3$hts)

dataset<-c(rep("Autosomes",length(df.hts2$hts)),rep("Z-chromosome", length(df.hts3$hts)))

df.hts4<-rbind(df.hts2,df.hts3)
df.hts4<-data.frame(cbind(dataset,df.hts4))

df.hts4$tree<-factor(df.hts4$tree, levels = c("T1","T2","T3"))

#pdf(file = "../figs_30Oct24/TreeHeights_by_topology.pdf", width = 7, height = 5)
bp.hts<-ggplot(df.hts4, aes(x=tree, y=hts)) + 
  #geom_violin(aes(color = dataset), trim = FALSE, position = position_dodge(0.9) ) +
  geom_boxplot(aes(color = dataset), width = 0.75, position = position_dodge(0.9), outliers = F) +
  #geom_point(aes(color = dataset),position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),alpha=0.3)+
  scale_color_manual(values = c("#353436","#1b98e0"))+
  theme_bw() + xlab("") + ylab("Node Height (subst. per site)")+ggtitle("Tree Heights for alternative topologies")

bp.hts
#dev.off()

kruskal.test(hts~tree, data=df.hts2)

# Kruskal-Wallis rank sum test
# 
# data:  hts by tree
# Kruskal-Wallis chi-squared = 41.982, df = 2, p-value = 7.653e-10

dunnTest(hts~tree, data=df.hts2)

# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Holm method.
# 
# Comparison         Z      P.unadj        P.adj
# 1    T1 - T2 -6.277236 3.446457e-10 1.033937e-09
# 2    T1 - T3 -1.583528 1.133011e-01 1.133011e-01
# 3    T2 - T3  4.489111 7.152100e-06 1.430420e-05

kruskal.test(hts~tree, data=df.hts3)

# Kruskal-Wallis rank sum test
# 
# data:  hts by tree
# Kruskal-Wallis chi-squared = 26.92, df = 2, p-value = 1.427e-06

dunnTest(hts~tree, data=df.hts3)

# Comparison          Z      P.unadj        P.adj
# 1    T1 - T2 -3.6097121 3.065370e-04 6.130740e-04
# 2    T1 - T3 -4.4982125 6.852719e-06 2.055816e-05
# 3    T2 - T3 -0.9433557 3.454990e-01 3.454990e-01

kruskal.test(hts~dataset, data=df.hts4[df.hts4$tree=="T2",]) #P<0.0001
kruskal.test(hts~dataset, data=df.hts4[df.hts4$tree=="T3",]) #P<0.0001
kruskal.test(hts~dataset, data=df.hts4[df.hts4$tree=="T1",]) #P=0.004378
