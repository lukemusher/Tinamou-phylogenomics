#plotting smilogram of UCE data
#first, run phyluce:

#phyluce_align_get_smilogram_from_alignments \
#--alignments mafft-nexus-clean-50p/ --output smilogram \
#--cores 16 --input-format nexus

#next take the output from the phyluce script and create a csv
#sqlite3 -header -csv smilogram.sqlite "select * from by_locus;" > smilogram.csv

#now, in R:

setwd("~/Documents/ANSDU/Tinamous/smilograms/")

library(ggplot2)

data<-read.csv("smilogram-1000flank")

s<-qplot(data$distance_from_center,data$freq, colour=data$substitutions, 
      xlab = "Distance from center of alignment (bp)", 
      ylab = "Frequency of variant bases", main="UCEs 1000 bp flanking region") + 
  scale_colour_gradient(high="tomato2", low ="steelblue2")+
  labs(colour='# substitutions')


data2<-read.csv("smilogram-100flank")

s1<-qplot(data2$distance_from_center,data2$freq, colour=data2$substitutions, 
      xlab = "Distance from center of alignment (bp)", 
      ylab = "Frequency of variant bases", main="UCEs 100 bp flanking region") + 
  scale_colour_gradient(high="tomato2", low ="steelblue2")+
  labs(colour='# substitutions')

data3<-read.csv("smilogram-300flank")

s2<-qplot(data3$distance_from_center,data3$freq, colour=data3$substitutions, 
      xlab = "Distance from center of alignment (bp)", 
      ylab = "Frequency of variant bases", main="UCEs 300 bp flanking region") + 
  scale_colour_gradient(high="tomato2", low ="steelblue2")+
  labs(colour='# substitutions')

data4<-read.csv("smilogram-CDS-95p")

s3<-qplot(data4$distance_from_center,data4$freq, colour=data4$substitutions, 
          xlim=c(-5000,5000), xlab = "Distance from center of alignment (bp)", 
          ylab = "Frequency of variant bases", main="CDS") + 
  scale_colour_gradient(high="tomato2", low ="steelblue2")+
  labs(colour='# substitutions')

library(ggpubr)
pdf(file = "./smilograms.pdf", width = 13, height = 15)
ggarrange(s3,s1,s2,s,p2,p3,p,p1, 
          ncol = 2, nrow = 4)

dev.off()

data<-read.csv("smilogram-missing-1000flank")

p<-qplot(data$distance_from_center,data$freq, colour=data$substitutions, 
         xlab = "Distance from center of alignment (bp)", 
         ylab = "Frequency of missing data", main="UCEs 1000 bp flanking region") + 
  scale_colour_gradient(high="tomato2", low ="steelblue2")+
  labs(colour='# substitutions')

data2<-read.csv("smilogram-missing-100flank")

p1<-qplot(data2$distance_from_center,data2$freq, colour=data2$substitutions, 
          xlab = "Distance from center of alignment (bp)", 
          ylab = "Frequency of missing data", main="UCEs 100 bp flanking region") + 
  scale_colour_gradient(high="tomato2", low ="steelblue2")+
  labs(colour='# substitutions')

data3<-read.csv("smilogram-missing-300flank")

p2<-qplot(data3$distance_from_center,data3$freq, colour=data3$substitutions, 
          xlab = "Distance from center of alignment (bp)", 
          ylab = "Frequency of missing data", main="UCEs 300 bp flanking region") + 
  scale_colour_gradient(high="tomato2", low ="steelblue2")+
  labs(colour='# substitutions')

data4<-read.csv("smilogram-missing-CDS-95p")

p3<-qplot(data4$distance_from_center,data4$freq, colour=data4$substitutions, 
          xlim=c(-10000,10000), xlab = "Distance from center of alignment (bp)", 
          ylab = "Frequency of missing data", main="CDS") + 
  scale_colour_gradient(high="tomato2", low ="steelblue2")+
  labs(colour='# substitutions')

library(ggpubr)
pdf(file = "./smilograms-missing.pdf", width = 9, height = 6)
ggarrange(p1,p2,p,p3, 
          ncol = 2, nrow = 2)

dev.off()
