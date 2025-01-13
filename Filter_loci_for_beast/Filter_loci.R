#This script is written by Lukas J Musher. Please cite:
#Musher, L.J., M. Ferreira, A. Auerbach, J. Mckay, and J. Cracraft. 
#Why is Amazonia a source of biodiversity? Climate mediated dispersal and synchronous 
#speciation across the Andes in an avian group (Tityrinae) (In revision) Proceedings of the Royal Society B
rm(list=ls())
par(mfrow=c(1,1))
library(gplots)
library(phytools)
library(ape)
library(phangorn)
library(ips)

#must have exactly same genes and trees--if there are 20 genes there must be 20 genes 
#and named appropriately so you know you have the right tree with each gene
#ideally you will have uces and gts named exactly the same e.g. uce-10.nexus uce-10.tre
#the easiest way to do this is make a list of the genes, and use a for loop to change
# the names of the RAxML outputs to just the locus name
#e.g. in linux/unix: for i in `cat names.genes`; do mv RAxML*$i* $i.tre; done
#runs slowly on datasets with lots of tips
#Start fresh R to maximize efficiency

setwd("/Users/lukasmusher/Documents/ANSDU/Tinamous/TimeTree/Autosomes-mafft-nexus-clean-trimmed-100p/")

files <- list.files(path="./", pattern="*.nexus", full.names=F, recursive=FALSE)
trees <- list.files(path="./",pattern="*.treefile",full.names=F,recursive=FALSE)
length(files)
length(trees)

outgroup<-c("rhea_penn")
spec_tree<- root(read.tree("../../Trees/uce1000.full.100p.iqtree.treefile"),outgroup, resolve.root = T)
RF.dist(spec_tree,read.tree(trees[10]))

t<-drop.tip(read.tree(trees[10]),)
plot(spec_tree)

remove.taxa<-c("rhea_penn","anodid","crybou","cssou586295","crysouzanf","cbcos53413","cuund34614","crytatino","notmacmac","eeele16893","nncad473927","rhyrufpal") #any taxa in the tree that you don't want involved in the likelihood estimations, for example taxa that might overinfluence the likelihoodds due to missing data etc

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

#lambda is a smoothing parameter for estimating divergence times using penalized likelihood
#lambda is correlated with "clock-likeness" higher values indicate closer to clock i.e. even/lined up tips

find_lambda_best<-function(x,tree,remove.taxa){
  x<-read.nexus.data(x)
  x<-as.DNAbin(x)
  X<-as.phyDat(x)
  t<-c()
  tree1<-drop.tip(tree,tree$tip.label[na.omit(match(remove.taxa,tree$tip.label))])
  fit<-pml(tree1, X, model="GTRG")   #do likelihood estimation of the distance tree l<-10^(-1:6)
  l<-c(0.001,1,100)
  cv<-sapply(l, function(x) sum(attr(chronopl(tree1, lambda=x, CV=T), "D2")))
  df<-as.data.frame(cbind(l,cv))
  l.min<-df$l[df$cv==min(df$cv)]
  returnValue(l.min)
}

#Use likelihood ratio tests to estimate likelihood ratio of clocklike to not clocklike model
is.clocklike<-function(dna,tree,lambda,outgroup,remove.taxa){
  x<-read.nexus.data(dna)
  x<-as.DNAbin(x)
  X<-as.phyDat(x)
  t<-c()
  counts=0
  tree1<-drop.tip(tree,tree$tip.label[na.omit(match(remove.taxa,tree$tip.label))])
  fit<-pml(tree1, X, model="GTRG")   #do likelihood estimation of tree 
  chr<-chronopl(phy=tree1, lambda=lambda, age.min=1, CV=T) #create chronogram
  fit2<-pml(chr, X, model="GTRG")   #do likelihood estimation of chronogram 
  LR=(2*(fit$logLik-fit2$logLik))
  return(LR)
}

#create the dataframe
Make_df_lambda<-function(files,trees,outgroup,remove.taxa,spec_tree,outfile){
  loc<-c()
  lambda<-c() #create empty vector
  log10lambda<-c()
  var.sites<-c()
  var.frac<-c()
  seq.len<-c()
  df<-c()
  LR<-c()
  p<-c()
  tst<-c()
  Tre<-c()
  RF<-c()
  for(i in 1:length(files)){ #for each locus in the directory
    loc[i]<-files[i]
    tree<-read.tree(trees[i])
    Tre[i]<-trees[i]
    tree<-drop.tip(tree,tree$tip.label[na.omit(match(remove.taxa,tree$tip.label))])
    print(paste("Locus",i,": ",loc[i],sep=""))
    try(lambda[i]<-find_lambda_best(files[i],tree,remove.taxa)) #the ith element is the result of finding the var.frac of the ith locus in the directory
    try(log10lambda[i]<-log10(lambda[i]))
    var.sites[i]<-pis(files[i])
    seq.len[i]<-site_count(files[i])
    var.frac[i]<-var.sites[i]/seq.len[i]
    try(LR[i]<-is.clocklike(files[i],tree,lambda[i],outgroup,remove.taxa))
    df[i]<-length(tree$tip.label)-2
    try(p[i]<-(pchisq(LR[i],df[i],lower.tail = FALSE)))
    spec_tree1<-drop.tip(spec_tree,spec_tree$tip.label[na.omit(match(remove.taxa,spec_tree$tip.label))])
    RF[i]<-RF.dist(tree,spec_tree1)
    try(data.f<-as.data.frame(cbind(loc,Tre,seq.len,var.sites,var.frac,lambda,log10lambda,LR,RF,p,df)))
    write.csv(data.f,file = outfile)
  }
  return(data.f)
}

#run code to make dataframe:

start <- Sys.time()

uce_table<-Make_df_lambda(files,trees,outgroup,remove.taxa,spec_tree,outfile="../tinamous-27Aug24.csv")

print( Sys.time() - start )

###############################
##########ANALYSES#############
###############################

par(mfrow=c(1,1))
tab<-na.omit(read.csv("../tinamous-27Aug24.csv"))
hist(tab$seq.len)
names(tab)
tab$LR<-as.numeric(tab$LR)
summary(lm(tab$LR~tab$RF))
plot(tab$RF,tab$LR)
abline(lm(tab$LR~tab$RF))
mod<-lm(LR~var.sites+var.frac+RF+log10lambda+seq.len, data = tab)

LR.opt<-mean((tab$LR))-2*sd(tab$LR)
RF.opt<-mean(tab$RF)#+2*sd(tab$RF)
tab1<-tab[tab$LR<=LR.opt,]
tab2<-tab1[tab1$RF<=RF.opt,]

length(tab2[,1])

pdf(file = "../../figs_30Oct24/loci_clock_plots.pdf", width = 8, height = 6)
par(mfrow=c(2,2))
plot(tab$var.sites,tab$var.frac,pch=19, xlab="# PIS", ylab="Prop. PIS")
points(tab1$var.sites,tab1$var.frac, col="red",pch=3,cex=2)
points(tab2$var.sites,tab2$var.frac,col="green",pch=3,cex=2)

plot(tab$RF,tab$LR,pch=19, ylab="Likelihood ratio", xlab="Robinson-Fould's Distance")
points(tab1$RF,tab1$LR, col="red",pch=3,cex=2)
points(tab2$RF,tab2$LR,col="green",pch=3,cex=2)

plot(tab$var.frac,tab$LR,pch=19, ylab="Likelihood ratio", xlab="# PIS")
points(tab1$var.frac,tab1$LR, col="red",pch=3,cex=2)
points(tab2$var.frac,tab2$LR,col="green",pch=3,cex=2)

plot(tab$var.sites,tab$LR,pch=19, ylab="Likelihood ratio", xlab="Prop. PIS")
points(tab1$var.sites,tab1$LR, col="red",pch=3,cex=2)
points(tab2$var.sites,tab2$LR,col="green",pch=3,cex=2)

dev.off()

hist(tab$LR, main="Historgram of likelihood ratios",xlab="2*(lnL clock - lnL nonclock)")
hist(tab2$LR, main="Historgram of likelihood ratios after filtering",xlab="2*(lnL clock - lnL nonclock)")
hist(tab$RF, main="Historgram of RF distances",xlab="RF Distance")
hist(tab2$RF, main="Historgram of RF distances after filtering",xlab="RF Distance")

write.csv(x = tab2,file = "../filtered_loci.csv")
