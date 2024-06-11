####First part looks at quartet frequencies, and is a little hard to use
library(ggplot2)
library(ape)
library(phytools)
library(MSCquartets)
citation("MSCquartets")
setwd("~/Documents/ANSDU/Tinamous/ASTRAL/")

#read in file with genetrees. These should have complete taxon sampling (no missing taxa in any tree)
trees<-read.tree("bestML_UCE1000_collapsed.trees")

#define species to drop from tree
species<-c("Cuund34614","Cusim651488","Cuyap42937","Cuver35355","AnoDid_GCA_006937325", "Njul32810", "NotJul_GCA_013398735", "Nbbon33138", "Nncad473927", "Not_nig_GCA_013398345", "Tmper40426", "Tmser456416", "Tmmaj20810", "Ttkle11280", "Tsol25977", "Tiosg311159", "TinGut_GCA_000705375", "Tgut85944", "Cttat35476", "CparMNT800", "Coobs25896", "Cocas11195", "Copun86899", "Cooch113498", "Cvar27099", "Cbar42396", "Crcas434024", "Cbre625102", "CsalbGUR081", "CrySou_GCA_013389845", "Csnig6048", "Cssou586295", "Cscap613347", "Cbcos53413", "Crber182233", "Ccin21636", "Cpta650391", "Nbor79764", "Tanan133116AB", "Ndaga68765", "Nomin344282", "Nmnig16884", "Nmann614501", "Nccin2861", "Rrruf13855", "Rhmac145784", "Npnie62872", "NotOrn_GCA_013398335", "NotPen_GCA_013398315", "Ntac27158", "NotPer_GCA_003342845", "Nppen6601", "Nocur425359", "Nobra61359", "Noros645439", "Noorn61422", "Effor12668", "EudEle_GCA_003342815", "Eepat56833", "Eeele16893", "Tpen7597", "Tiing186371", "Rhea_penn_GCA_003342835")

#drop unwanted tips from all genetrees
for(i in 1:length(trees)){
  for (j in species){
    trees[[i]]<-drop.tip(trees[[i]],j)
  }
  #print(length(trees[[i]]$tip.label))
}

#use mscquartets to evaluate all possible quartets for all gene trees. I suggest using a small number of species because the number of quartets gets excessive fast.

par(mfrow=c(1,1))
qtab<-quartetTable(trees)

qtab.print<-quartetTablePrint(qtab)

#####################
#####Branch 1########
#####################

rqt<-quartetTableResolved(qtab)
rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("CryUnd_GCA_013389825","Crdui419412","Crker147057","Crnoc59478","Ceery21146","Cratr320360","Cstr9577"), taxaB = c("CryCin_GCA_003342915","Ccgol2213","Ctra66583"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("Ccgol2213CryCin_GCA_003342915Ctra66583","Crnoc59478","Cratr320360"), taxaB = c("CryUnd_GCA_013389825","Crdui419412","Cstr9577","Ceery21146","Crker147057"))
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

df


#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p1<-ggplot(data=df, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  ylim(0,0.75)+
  ylab("Relative frequency")

p1

#####################
#####Branch 2########
#####################

rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("CryUnd_GCA_013389825","Crdui419412","Crker147057","Crnoc59478","Ceery21146","Cratr320360","Cstr9577"), taxaB = c("CryCin_GCA_003342915","Ccgol2213","Ctra66583"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("Ccgol2213CryCin_GCA_003342915Ctra66583","Crnoc59478","Cratr320360","Ceery21146","Crker147057"), taxaB = c("CryUnd_GCA_013389825","Crdui419412","Cstr9577"))
rqt.collapsed3<-quartetTableCollapse(rqt = rqt.collapsed2, taxaA = c("Ccgol2213CryCin_GCA_003342915Ctra66583","Crdui419412CryUnd_GCA_013389825Cstr9577","Ceery21146","Crker147057"), taxaB = c("Crnoc59478","Cratr320360"))
rqt.collapsed4<-quartetTableCollapse(rqt = rqt.collapsed3, taxaA = c("Ccgol2213CryCin_GCA_003342915Ctra66583","Cratr320360Crnoc59478","Crdui419412CryUnd_GCA_013389825Cstr9577"), taxaB = c("Ceery21146","Crker147057"))

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


#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p2<-ggplot(data=df2, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  ylim(0,0.75)+
  ylab("Relative frequency")

p2

#####################
######branch 4#########
#####################

rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("CryUnd_GCA_013389825","Crdui419412","Crker147057","Ceery21146","Cstr9577"), taxaB = c("CryCin_GCA_003342915","Ccgol2213","Ctra66583","Crnoc59478","Cratr320360"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("Ccgol2213Cratr320360Crnoc59478CryCin_GCA_003342915Ctra66583","Ceery21146","Crker147057"), taxaB = c("CryUnd_GCA_013389825","Crdui419412","Cstr9577"))

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


#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p3<-ggplot(data=df3, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  ylim(0,0.75)+
  ylab("Relative frequency")

p3

#####################
#####Branch 5########
#####################

rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("CryUnd_GCA_013389825","Crdui419412","Cstr9577"), taxaB = c("CryCin_GCA_003342915","Ccgol2213","Ctra66583","Crnoc59478","Cratr320360","Crker147057","Ceery21146"))


tab<-quartetTablePrint(rqt.collapsed1)
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


#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p4<-ggplot(data=df4, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  ylim(0,0.75)+
  ylab("Relative frequency")

p4

#####################
#####Branch 3########
#####################

rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("CryUnd_GCA_013389825","Crdui419412","Crker147057","Ceery21146","Cstr9577","Ctra66583","Crnoc59478","Cratr320360"), taxaB = c("CryCin_GCA_003342915","Ccgol2213"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("Ccgol2213CryCin_GCA_003342915","CryUnd_GCA_013389825","Crdui419412","Crker147057","Ceery21146","Cstr9577","Ctra66583"), taxaB = c("Crnoc59478","Cratr320360"))
rqt.collapsed3<-quartetTableCollapse(rqt = rqt.collapsed2, taxaA = c("Ccgol2213CryCin_GCA_003342915","Cratr320360Crnoc59478","Ctra66583"), taxaB = c("CryUnd_GCA_013389825","Crdui419412","Crker147057","Ceery21146","Cstr9577"))

tab<-quartetTablePrint(rqt.collapsed3)
rows<-length(tab[,1])

quartet<-c(rep(1:rows,3))
tax1<-rep(tab[,1],3)
tax2<-rep(tab[,2],3)
tax3<-rep(tab[,3],3)
tax4<-rep(tab[,4],3)
topology<-c(rep("12|34",rows),rep("13|24",rows),rep("14|23",rows))
no.genetrees<-as.numeric(c(tab[,5],tab[,6],tab[,7]))
prop.quartets<-no.genetrees/(sum(no.genetrees))

df5<-data.frame(cbind(quartet,tax1,tax2,tax3,tax4,topology,no.genetrees,prop.quartets))

df5


#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p5<-ggplot(data=df5, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  ylim(0,0.75)+
  ylab("Relative frequency")

p5

#####################
#####Branch 6########
#####################

rqt.collapsed1<-quartetTableCollapse(rqt = rqt, taxaA = c("CryUnd_GCA_013389825","Crdui419412","Crker147057","Ceery21146","Cstr9577"), taxaB = c("CryCin_GCA_003342915","Ccgol2213","Ctra66583","Crnoc59478","Cratr320360"))
rqt.collapsed2<-quartetTableCollapse(rqt = rqt.collapsed1, taxaA = c("Ccgol2213Cratr320360Crnoc59478CryCin_GCA_003342915Ctra66583","CryUnd_GCA_013389825","Crdui419412","Cstr9577"), taxaB = c("Crker147057","Ceery21146"))
rqt.collapsed3<-quartetTableCollapse(rqt = rqt.collapsed2, taxaA = c("Ccgol2213Cratr320360Crnoc59478CryCin_GCA_003342915Ctra66583","Ceery21146Crker147057","Crdui419412"), taxaB = c("CryUnd_GCA_013389825","Cstr9577"))

tab<-quartetTablePrint(rqt.collapsed3)
rows<-length(tab[,1])

quartet<-c(rep(1:rows,3))
tax1<-rep(tab[,1],3)
tax2<-rep(tab[,2],3)
tax3<-rep(tab[,3],3)
tax4<-rep(tab[,4],3)
topology<-c(rep("12|34",rows),rep("13|24",rows),rep("14|23",rows))
no.genetrees<-as.numeric(c(tab[,5],tab[,6],tab[,7]))
prop.quartets<-no.genetrees/(sum(no.genetrees))

df6<-data.frame(cbind(quartet,tax1,tax2,tax3,tax4,topology,no.genetrees,prop.quartets))

df6


#Plot relative quartet frequencies (again there may be more than 5 to visualize)

p6<-ggplot(data=df6, aes(x=topology, y=as.numeric(prop.quartets), fill=topology)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal() +
  geom_hline(yintercept = 0.33, col="red", lty=2)+
  xlab("Topology")+
  ylim(0,0.75)+
  ylab("Relative frequency")

p6

library(ggpubr)

pdf(file = "../figs/MSC_Quartets_UCE300.pdf", width = 14, height = 4)
ggarrange(p1,p2,p5,p3,p4,
          ncol = 5, nrow = 1)
dev.off()
