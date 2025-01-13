library(ape)
library(phytools)
library(phangorn)
library(devtools)
library(TESS)
#We can run this with two expected changes, who knows? We come up with some expected survival prob after an extinction. 
#Using the expected survival probability, we compute the α and β parameters of
#the beta distribution. We set the value of β to be large, which focuses the prior
#density more tightly around the expected survival probability. Then, we compute
#α based on the expected survival probability and the specified β value
par(mfrow=c(1,1))
citation('TESS')
setwd("~/Documents/ANSDU/Tinamous/TimeTree/")

t <- treeio::read.beast("~/Documents/ANSDU/Tinamous/TimeTree/Beuti_2.7_23_uces_6fossils+Moa_final.tre")
t<-drop.tip(as.phylo(t),c("anodid","rhea_penn","Eudromia_sp","E_olsoni","N_parvula","C_reai","Macn-sc-3610","Macn-sc-3613"))
plot(t)

expectedSurvivalProbability <- 0.05

numExpectedRateChanges<-2
numExpectedMassExtinctions<-0

# force.ultrametric<-function(tree,method=c("nnls","extend")){
#   method<-method[1]
#   if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
#                                      rooted=TRUE,trace=0)
#   else if(method=="extend"){
#     h<-diag(vcv(tree))
#     d<-max(h)-h
#     ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
#                y=tree$edge[,2])
#     tree$edge.length[ii]<-tree$edge.length[ii]+d
#   } else 
#     cat("method not recognized: returning input tree\n\n")
#   tree
# }
#t<-force.ultrametric(t)

tess.analysis(t,
              empiricalHyperPriors = TRUE,
              samplingProbability = 1.0,
              estimateNumberMassExtinctions = TRUE,
              MAX_ITERATIONS = 100000,
              dir = "comet_no_mass_extinctions_tyran_CC15")

output <- tess.process.output("comet_no_mass_extinctions_tyran_CC15",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)


par(mfrow=c(2,3))
tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "speciation Bayes factors",
                               "extinction rates",
                               "extinction shift times",
                               "extinction Bayes factors"),las=2)


