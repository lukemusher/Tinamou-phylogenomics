####First part looks at quartet frequencies, and is a little hard to use
library(ggplot2)
library(ape)
library(phytools)
library(MSCquartets)
citation("MSCquartets")
setwd("~/Documents/ANSDU/Tinamous/ASTRAL/")

#whole UCE1000 dataset
#read in file with genetrees. These should have complete taxon sampling (no missing taxa in any tree)
trees<-read.tree("autosomes-cladeA-100p.trees") #complete autosomes UCE1000 gene trees
plot(trees[[3]])
#define species to drop from tree
species<-c("cuund34614", "cvar27099","crybou")
counts=0
for(i in 1:length(trees)){
  tree<-drop.tip(trees[[i]],species)
  #print(length(trees[[i]]$tip.label))
  if(length(tree$tip.label)==8){
    counts=counts+1
    trees[[counts]]<-tree
  }
  
}
counts

#use mscquartets to evaluate all possible quartets for all gene trees. I suggest using a small number of species because the number of quartets gets excessive fast.

par(mfrow=c(1,1))
qtab<-quartetTable(trees[1:counts])

qtab.print<-quartetTablePrint(qtab)

#####################
#####Branch 1########
#####################

rqt<-quartetTableResolved(qtab)
rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("crdui419412","crker147057","crnoc59478","ceery21146","cratr320360","cstr9577"), taxaB = c("ccgol2213","ctra66583"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("ccgol2213ctra66583","crnoc59478","cratr320360"), taxaB = c("crdui419412","cstr9577","ceery21146","crker147057"))
#rqt.collapsed3<-quartetTableCollapse(rqt = rqt.collapsed2, taxaA = c("Ceery21146Crker147057","Ccgol2213CryCin_GCA_003342915Ctra66583","Cratr320360"), taxaB = c("Cstr9577","Crdui419412","CryUnd_GCA_013389825"))

tab<-quartetTablePrint(rqt.collapsed2)
rows<-length(tab[,1])

quartet<-c(rep(1:rows,3))
tax1<-rep(tab[,1],3)
tax2<-rep(tab[,2],3)
tax3<-rep(tab[,3],3)
tax4<-rep(tab[,4],3)
topology<-c(rep("12|34",rows),rep("13|24",rows),rep("14|23",rows))
no.genetrees<-as.numeric(c(tab[,5],tab[,6],tab[,7]))
prop.quartets<-no.genetrees/(sum(no.genetrees))

df<-data.frame(cbind(quartet,tax1,tax2,tax3,tax4,topology,no.genetrees,prop.quartets))

df$no.genetrees<-as.numeric(df$no.genetrees)

#prune tree to one member of each clade of interst
tree<-read.tree("../Trees/crypt.uce1000.full.75p.astral.tre")
tree<-drop.tip(tree,c("cuyap42937","cuver35355","cuund34614","cusim651488","crycin","ccgol2213","cryund","cstr9577","crybou","ceery21146","crdui419412"))
plot(tree)

exp2<-expectedCFs(speciestree = tree,plot = T, model = "T1") #calculate expected freqs given species tree

library(EMT)
citation("EMT")
#test if observed differs from expected
multinomial.test(observed=df$no.genetrees,prob = c(exp2[6],exp2[7],exp2[5]),useChisq = T, MonteCarlo = T, ntrial = 10000000)
as.numeric(df$prop.quartets)
c(exp2[6],exp2[7],exp2[5])

#test if observed differs from equal proportions
multinomial.test(observed=df$no.genetrees,prob = c(1/3,1/3,1/3),useChisq = T, MonteCarlo = T, ntrial = 10000000)

p1<-ggplot(data=df, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() + theme(legend.position = "none") +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  #ylim(0,0.5)+
  ylab("Relative frequency")

p1

#####################
#####Branch 2########
#####################

rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("crdui419412","crker147057","crnoc59478","ceery21146","cratr320360","cstr9577","crybou"), taxaB = c("ccgol2213","ctra66583"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("ccgol2213ctra66583","crnoc59478","cratr320360","ceery21146","crker147057","crybou"), taxaB = c("crdui419412","cstr9577"))
rqt.collapsed3<-quartetTableCollapse(rqt = rqt.collapsed2, taxaA = c("cratr320360","crdui419412cstr9577","ceery21146","crker147057","crybou"), taxaB = c("crnoc59478","ccgol2213ctra66583"))
rqt.collapsed4<-quartetTableCollapse(rqt = rqt.collapsed3, taxaA = c("ccgol2213ctra66583crnoc59478","cratr320360","crdui419412cstr9577"), taxaB = c("ceery21146","crker147057","crybou"))

tab<-quartetTablePrint(rqt.collapsed4)
rows<-length(tab[,1])

quartet<-c(rep(1:rows,3))
tax1<-rep(tab[,1],3)
tax2<-rep(tab[,2],3)
tax3<-rep(tab[,3],3)
tax4<-rep(tab[,4],3)
topology<-c(rep("12|34",rows),rep("13|24",rows),rep("14|23",rows))
no.genetrees<-as.numeric(c(tab[,5],tab[,6],tab[,7]))
prop.quartets<-no.genetrees/(sum(no.genetrees))

df2<-data.frame(cbind(quartet,tax1,tax2,tax3,tax4,topology,no.genetrees,prop.quartets))

df2
df2$no.genetrees<-as.numeric(df2$no.genetrees)

tree<-read.tree("../Trees/crypt.uce1000.full.75p.astral.tre")
tree<-drop.tip(tree,c("cuyap42937","cuver35355","cuund34614","cusim651488","crycin","ccgol2213","cryund","ctra66583","crybou","ceery21146","crdui419412"))
plot(tree)

exp3<-expectedCFs(speciestree = tree,plot = T, model = "T1")

#test if observed differs from expected freqs based on species tree
multinomial.test(observed=df2$no.genetrees,prob = c(exp3[7],exp3[6],exp3[5]),useChisq = T, MonteCarlo = T, ntrial = 10000000)
as.numeric(df2$prop.quartets)
c(exp3[7],exp3[6],exp3[5])

#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p2<-ggplot(data=df2, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() + theme(legend.position = "none") +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  #ylim(0,0.5)+
  ylab("Relative frequency")

p2


###CHR Z ONLY
#read in file with genetrees. These should have complete taxon sampling (no missing taxa in any tree)
#trees<-read.tree("bestML_UCE1000_collapsed.trees")
trees<-read.tree("chrz-cladeA-100p.trees") #comp Z-chrom linked crypturellus

plot(drop.tip(trees[[2]],species))
#drop unwanted tips from all genetrees

counts=0
for(i in 1:length(trees)){
  tree<-drop.tip(trees[[i]],species)
  #print(length(trees[[i]]$tip.label))
  if(length(tree$tip.label)==9){
    counts=counts+1
    trees[[counts]]<-tree
  }
  
}
trees[[counts]]

#use mscquartets to evaluate all possible quartets for all gene trees. I suggest using a small number of species because the number of quartets gets excessive fast.

par(mfrow=c(1,1))
qtab<-quartetTable(trees[1:counts])

qtab.print<-quartetTablePrint(qtab)

#####################
#####Branch 1########
#####################

rqt<-quartetTableResolved(qtab)
rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("crdui419412","crker147057","crnoc59478","ceery21146","cratr320360","cstr9577","crybou"), taxaB = c("ccgol2213","ctra66583"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("ccgol2213ctra66583","crnoc59478","cratr320360"), taxaB = c("crdui419412","cstr9577","ceery21146","crker147057","crybou"))
#rqt.collapsed3<-quartetTableCollapse(rqt = rqt.collapsed2, taxaA = c("Ceery21146Crker147057","Ccgol2213CryCin_GCA_003342915Ctra66583","Cratr320360"), taxaB = c("Cstr9577","Crdui419412","CryUnd_GCA_013389825"))

tab<-quartetTablePrint(rqt.collapsed2)
rows<-length(tab[,1])

quartet<-c(rep(1:rows,3))
tax1<-rep(tab[,1],3)
tax2<-rep(tab[,2],3)
tax3<-rep(tab[,3],3)
tax4<-rep(tab[,4],3)
topology<-c(rep("12|34",rows),rep("13|24",rows),rep("14|23",rows))
no.genetrees<-as.numeric(c(tab[,5],tab[,6],tab[,7]))
prop.quartets<-no.genetrees/(sum(no.genetrees))

df3<-data.frame(cbind(quartet,tax1,tax2,tax3,tax4,topology,no.genetrees,prop.quartets))

df3

df3$no.genetrees<-as.numeric(df3$no.genetrees)

multinomial.test(observed=df3$no.genetrees,prob = c(exp2[6],exp2[7],exp2[5]),useChisq = T, MonteCarlo = T, ntrial = 10000000)
as.numeric(df3$prop.quartets)
c(exp2[6],exp2[7],exp2[5])
#p<0.0001

multinomial.test(observed=df3$no.genetrees,prob = as.numeric(df$prop.quartets),useChisq = T, MonteCarlo = T, ntrial = 10000000)
as.numeric(df3$prop.quartets)
as.numeric(df$prop.quartets)
#p=0.00034

multinomial.test(observed=df3$no.genetrees,prob = c(1/3,1/3,1/3),useChisq = T, MonteCarlo = T, ntrial = 10000000)
#p=0.876

#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p3<-ggplot(data=df3, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() + theme(legend.position = "none") +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  #ylim(0,0.5)+
  ylab("Relative frequency")

p3

#####################
#####Branch 2########
#####################

rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("crdui419412","crker147057","crnoc59478","ceery21146","cratr320360","cstr9577","crybou"), taxaB = c("ccgol2213","ctra66583"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("ccgol2213ctra66583","crnoc59478","cratr320360","ceery21146","crker147057","crybou"), taxaB = c("crdui419412","cstr9577"))
rqt.collapsed3<-quartetTableCollapse(rqt = rqt.collapsed2, taxaA = c("cratr320360","crdui419412cstr9577","ceery21146","crker147057","crybou"), taxaB = c("crnoc59478","ccgol2213ctra66583"))
rqt.collapsed4<-quartetTableCollapse(rqt = rqt.collapsed3, taxaA = c("ccgol2213ctra66583crnoc59478","cratr320360","crdui419412cstr9577"), taxaB = c("ceery21146","crker147057","crybou"))

tab<-quartetTablePrint(rqt.collapsed4)
rows<-length(tab[,1])

quartet<-c(rep(1:rows,3))
tax1<-rep(tab[,1],3)
tax2<-rep(tab[,2],3)
tax3<-rep(tab[,3],3)
tax4<-rep(tab[,4],3)
topology<-c(rep("12|34",rows),rep("13|24",rows),rep("14|23",rows))
no.genetrees<-as.numeric(c(tab[,5],tab[,6],tab[,7]))
prop.quartets<-no.genetrees/(sum(no.genetrees))

df4<-data.frame(cbind(quartet,tax1,tax2,tax3,tax4,topology,no.genetrees,prop.quartets))

df4

df4$no.genetrees<-as.numeric(df4$no.genetrees)

multinomial.test(observed=df4$no.genetrees,prob = c(exp3[7],exp3[6],exp3[5]),useChisq = T, MonteCarlo = T, ntrial = 10000000)
as.numeric(df4$prop.quartets)
c(exp3[7],exp3[6],exp3[5])
#p<0.0001

multinomial.test(observed=df4$no.genetrees,prob = as.numeric(df2$prop.quartets),useChisq = T, MonteCarlo = T, ntrial = 10000000)
as.numeric(df4$prop.quartets)
as.numeric(df2$prop.quartets)
#p<0.0001

#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p4<-ggplot(data=df4, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() + theme(legend.position = "none") +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  #ylim(0,0.5)+
  ylab("Relative frequency")

p4

library(ggpubr)
#pdf(file = "../figs_30Oct24//MSC_Quartets_UCE1000.pdf", width = 4, height = 5)
ggarrange(p1,p2,p3,p4,
          ncol = 2, nrow = 2)
#dev.off()


