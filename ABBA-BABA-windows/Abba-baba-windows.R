#######################################
#############INTROGRESSION#############
#######################################

#laptop
setwd("~/Documents/ANSDU/Tinamous/github/Tinamou-phylogenomics2/ABBA-BABA-windows/")

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(FSA)

#standard error function

se<-funtion(x){
  std.err<-sd(x)/sqrt(length(x))
  print(std.err)
}
# Read the CSV file
data1 <- read_csv("ABBABABAwindows.w100k.T1.csv")
names(data1)
data1$scaffold<-factor(data1$scaffold,levels = 
              c("1","2","3","4","5","6","7",
                "8","9","10","11","12","13","14",
                "15","16","17","18","19","20","21",
                "22","23","24","25","26","27","28",
                "29","30","31","32","33","34","35",
                "36","37","38","39","40"))

m<-ggplot(data1,aes(x=scaffold, y=fdM, group=scaffold))+
  ggtitle("Introgression")+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +xlab("") + ylab("fdM")
m

m1<-ggplot(data1,aes(x=micro.macro, y=fdM, group=micro.macro))+
  ggtitle("Introgression")+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +xlab("") + ylab("fdM")
m1

mean(na.omit(data1$fdM[data1$type=="autosome"]))
sd(na.omit(data1$fdM[data1$type=="autosome"]))
se(data1$fdM[data1$type=="autosome"])
mean(na.omit(data1$fdM[data1$micro.macro=="macro"]))
sd(na.omit(data1$fdM[data1$micro.macro=="macro"]))
se(data1$fdM[data1$micro.macro=="macro"])
mean(na.omit(data1$fdM[data1$micro.macro=="micro"]))
sd(na.omit(data1$fdM[data1$micro.macro=="micro"]))
se(data1$fdM[data1$micro.macro=="micro"])
mean(na.omit(data1$fdM[data1$micro.macro=="Z-chromosome"]))
sd(na.omit(data1$fdM[data1$micro.macro=="Z-chromosome"]))
se(data1$fdM[data1$micro.macro=="Z-chromosome"])

kruskal.test(fdM~type,data=data1)

kruskal.test(fdM~micro.macro,data=data1)

dunnTest(as.numeric(fdM)~micro.macro,data=data1)

data_clean <- data1 %>%
  filter(!is.na(fdM) & !is.na(mid)) %>%
  mutate(scaffold = as.character(scaffold))

# Function to extract chromosome number/letter for proper sorting
extract_chr <- function(scaffold) {
  # Extract the chromosome part (number or letter after "Chr")
  chr_match <- regexpr("Chr(\\d+|[A-Z])", scaffold, perl = TRUE)
  if (chr_match > 0) {
    chr_part <- regmatches(scaffold, chr_match)
    chr_value <- gsub("Chr", "", chr_part)
    return(chr_value)
  }
  return(scaffold)
}

# Add chromosome extraction and create sorting key
data_clean <- data_clean %>%
  mutate(
    chr_extract = sapply(scaffold, extract_chr),
    # Create numeric version for sorting (letters get high numbers)
    chr_numeric = ifelse(
      grepl("^\\d+$", chr_extract),
      as.numeric(chr_extract),
      1000 + utf8ToInt(chr_extract)  # Letters get values > 1000
    )
  ) %>%
  # Sort by chromosome then by mid position
  arrange(chr_numeric, mid) %>%
  # Create sequential x-position for plotting
  mutate(x_position = row_number())

# Get scaffold labels for x-axis (sample for readability)
n_scaffolds <- length(unique(data_clean$scaffold))
label_step <- max(1, floor(n_scaffolds / 20))  # Show ~20 labels max

scaffold_labels <- data_clean %>%
  group_by(scaffold) %>%
  summarise(x_pos = first(x_position), .groups = 'drop') %>%
  slice(seq(1, n(), by = label_step))

# Create alternating colors for chromosomes
data_clean <- data_clean %>%
  mutate(
    # Create unique chromosome identifiers
    chr_id = dense_rank(chr_numeric),
    # Alternate colors for every other chromosome
    chr_color = ifelse(chr_id %% 2 == 1, "Color1", "Color2")
  )

# Create the scatter plot with alternating colors
p <- ggplot(data_clean, aes(x = x_position, y = fdM, color = chr_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Color1" = "steelblue", "Color2" = "darkred"),
    guide = "none"  # Hide legend since colors just indicate chromosome alternation
  ) +
  labs(
    title = "Genome-wide introgression assuming T1",
    subtitle = paste("Showing", nrow(data_clean), "data points across", 
                     length(unique(data_clean$scaffold)), "scaffolds"),
    x = "Chromosome",
    y = "fdM"
  ) +
  scale_x_continuous(
    breaks = scaffold_labels$x_pos,
    labels = scaffold_labels$scaffold
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(p)

data2 <- read_csv("ABBABABAwindows.w100k.T2.csv")
data2$scaffold<-factor(data2$scaffold,levels = 
                        c("1","2","3","4","5","6","7",
                          "8","9","10","11","12","13","14",
                          "15","16","17","18","19","20","21",
                          "22","23","24","25","26","27","28",
                          "29","30","31","32","33","34","35",
                          "36","37","38","39","40"))

m2<-ggplot(data2,aes(x=scaffold, y=fdM, group=scaffold))+
  ggtitle("Introgression")+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +xlab("") + ylab("fdM p2<->p3")
m2

m3<-ggplot(data2,aes(x=micro.macro, y=fdM, group=micro.macro))+
  ggtitle("Introgression")+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +xlab("") + ylab("fdM")
m3

mean(na.omit(data2$fdM[data1$type=="autosome"]))
sd(na.omit(data2$fdM[data1$type=="autosome"]))
se(data2$fdM[data1$type=="autosome"])
mean(na.omit(data2$fdM[data1$micro.macro=="macro"]))
sd(na.omit(data2$fdM[data1$micro.macro=="macro"]))
se(data2$fdM[data1$micro.macro=="macro"])
mean(na.omit(data2$fdM[data1$micro.macro=="micro"]))
sd(na.omit(data2$fdM[data1$micro.macro=="micro"]))
se(data2$fdM[data1$micro.macro=="micro"])
mean(na.omit(data2$fdM[data1$micro.macro=="Z-chromosome"]))
sd(na.omit(data2$fdM[data1$micro.macro=="Z-chromosome"]))
se(data2$fdM[data1$micro.macro=="Z-chromosome"])

kruskal.test(fdM~type,data=data2)

kruskal.test(fdM~micro.macro,data=data2)

dunnTest(as.numeric(fdM)~micro.macro,data=data2)

# Filter out rows with missing fdM values
data_clean <- data2 %>%
  filter(!is.na(fdM) & !is.na(mid)) %>%
  mutate(scaffold = as.character(scaffold))

# Add chromosome extraction and create sorting key
data_clean <- data_clean %>%
  mutate(
    chr_extract = sapply(scaffold, extract_chr),
    # Create numeric version for sorting (letters get high numbers)
    chr_numeric = ifelse(
      grepl("^\\d+$", chr_extract),
      as.numeric(chr_extract),
      1000 + utf8ToInt(chr_extract)  # Letters get values > 1000
    )
  ) %>%
  # Sort by chromosome then by mid position
  arrange(chr_numeric, mid) %>%
  # Create sequential x-position for plotting
  mutate(x_position = row_number())

# Get scaffold labels for x-axis (sample for readability)
n_scaffolds <- length(unique(data_clean$scaffold))
label_step <- max(1, floor(n_scaffolds / 20))  # Show ~20 labels max

scaffold_labels <- data_clean %>%
  group_by(scaffold) %>%
  summarise(x_pos = first(x_position), .groups = 'drop') %>%
  slice(seq(1, n(), by = label_step))

# Create alternating colors for chromosomes
data_clean <- data_clean %>%
  mutate(
    # Create unique chromosome identifiers
    chr_id = dense_rank(chr_numeric),
    # Alternate colors for every other chromosome
    chr_color = ifelse(chr_id %% 2 == 1, "Color1", "Color2")
  )

# Create the scatter plot with alternating colors
p2 <- ggplot(data_clean, aes(x = x_position, y = fdM, color = chr_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Color1" = "steelblue", "Color2" = "darkred"),
    guide = "none"  # Hide legend since colors just indicate chromosome alternation
  ) +
  labs(
    title = "Genome-wide introgression assuming T2",
    subtitle = paste("Showing", nrow(data_clean), "data points across", 
                     length(unique(data_clean$scaffold)), "scaffolds"),
    x = "Chromosome",
    y = "fdM"
  ) +
  scale_x_continuous(
    breaks = scaffold_labels$x_pos,
    labels = scaffold_labels$scaffold
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(p2)

# Read the CSV file
data3 <- read_csv("ABBABABAwindows.w100k.T3.csv")
data3$scaffold<-factor(data3$scaffold,levels = 
                        c("1","2","3","4","5","6","7",
                          "8","9","10","11","12","13","14",
                          "15","16","17","18","19","20","21",
                          "22","23","24","25","26","27","28",
                          "29","30","31","32","33","34","35",
                          "36","37","38","39","40"))

m4<-ggplot(data3,aes(x=scaffold, y=fdM, group=scaffold))+
  ggtitle("Introgression")+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +xlab("") + ylab("fdM p2<->p3")
m4

m5<-ggplot(data3,aes(x=micro.macro, y=fdM, group=micro.macro))+
  ggtitle("Introgression")+
  geom_boxplot(notch=F, outlier.shape=8, fill="red", alpha=0.4)+
  theme_bw() +xlab("") + ylab("fdM")
m5

mean(na.omit(data3$fdM[data1$type=="autosome"]))
sd(na.omit(data3$fdM[data1$type=="autosome"]))
se(data3$fdM[data1$type=="autosome"])
mean(na.omit(data3$fdM[data1$micro.macro=="macro"]))
sd(na.omit(data3$fdM[data1$micro.macro=="macro"]))
se(data3$fdM[data1$micro.macro=="macro"])
mean(na.omit(data3$fdM[data1$micro.macro=="micro"]))
sd(na.omit(data3$fdM[data1$micro.macro=="micro"]))
se(data3$fdM[data1$micro.macro=="micro"])
mean(na.omit(data3$fdM[data1$micro.macro=="Z-chromosome"]))
sd(na.omit(data3$fdM[data1$micro.macro=="Z-chromosome"]))
se(data3$fdM[data1$micro.macro=="Z-chromosome"])

kruskal.test(fdM~type,data=data3)

kruskal.test(fdM~micro.macro,data=data3)

dunnTest(as.numeric(fdM)~micro.macro,data=data3)

# Filter out rows with missing fdM values
data_clean <- data3 %>%
  filter(!is.na(fdM) & !is.na(mid)) %>%
  mutate(scaffold = as.character(scaffold))

# Function to extract chromosome number/letter for proper sorting
extract_chr <- function(scaffold) {
  # Extract the chromosome part (number or letter after "Chr")
  chr_match <- regexpr("Chr(\\d+|[A-Z])", scaffold, perl = TRUE)
  if (chr_match > 0) {
    chr_part <- regmatches(scaffold, chr_match)
    chr_value <- gsub("Chr", "", chr_part)
    return(chr_value)
  }
  return(scaffold)
}

# Add chromosome extraction and create sorting key
data_clean <- data_clean %>%
  mutate(
    chr_extract = sapply(scaffold, extract_chr),
    # Create numeric version for sorting (letters get high numbers)
    chr_numeric = ifelse(
      grepl("^\\d+$", chr_extract),
      as.numeric(chr_extract),
      1000 + utf8ToInt(chr_extract)  # Letters get values > 1000
    )
  ) %>%
  # Sort by chromosome then by mid position
  arrange(chr_numeric, mid) %>%
  # Create sequential x-position for plotting
  mutate(x_position = row_number())

# Get scaffold labels for x-axis (sample for readability)
n_scaffolds <- length(unique(data_clean$scaffold))
label_step <- max(1, floor(n_scaffolds / 20))  # Show ~20 labels max

scaffold_labels <- data_clean %>%
  group_by(scaffold) %>%
  summarise(x_pos = first(x_position), .groups = 'drop') %>%
  slice(seq(1, n(), by = label_step))

# Create alternating colors for chromosomes
data_clean <- data_clean %>%
  mutate(
    # Create unique chromosome identifiers
    chr_id = dense_rank(chr_numeric),
    # Alternate colors for every other chromosome
    chr_color = ifelse(chr_id %% 2 == 1, "Color1", "Color2")
  )

# Create the scatter plot with alternating colors
p3 <- ggplot(data_clean, aes(x = x_position, y = fdM, color = chr_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Color1" = "steelblue", "Color2" = "darkred"),
    guide = "none"  # Hide legend since colors just indicate chromosome alternation
  ) +
  labs(
    title = "Genome-wide introgression assuming T3",
    subtitle = paste("Showing", nrow(data_clean), "data points across", 
                     length(unique(data_clean$scaffold)), "scaffolds"),
    x = "Chromosome",
    y = "fdM"
  ) +
  scale_x_continuous(
    breaks = scaffold_labels$x_pos,
    labels = scaffold_labels$scaffold
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(p3)

library(ggpubr)
pdf(file = "../../../Figs_August2025/introgressionXchromosome.pdf", width = 14, height = 15)
ggarrange(p, p2,p3,
          ncol = 1, nrow = 3)
dev.off()

pdf(file = "../../../Figs_August2025/introgression.boxplots.pdf", width = 16, height = 5)
ggarrange(m, m2,m3,
          ncol = 3, nrow = 1)
dev.off()

pdf(file = "../../../Figs_August2025/introgressionXchromosome.T1.T2.pdf", width = 14, height = 10)
ggarrange(p, p2,
          ncol = 1, nrow = 2)
dev.off()

pdf(file = "../../../Figs_August2025/introgression.boxplots.T1.T2.pdf", width = 14, height = 10)
ggarrange(m, m2,
          ncol = 1, nrow = 2)
dev.off()


###distance from center
library(ggplot2)
library(ggpubr)
library(ape)
library(phytools)
library(phangorn)

tab1<-read.csv("ABBABABAwindows.w100k.T2.csv")
names(tab1)

counts=0
fdm.val<-c()
D.val<-c()
abs.fdm.val<-c()
abs.D.val<-c()
chrom.log.length<-c()
chrom.length<-c()
chrom.number<-c()
abba.count<-c()
baba.count<-c()
abba.prop<-c()
baba.prop<-c()
sites.used<-c()
position<-c()
chr.ctr<-c()
per.dist.fr.ctr<-c()
mid<-c()
type<-c()
micro.macro<-c()

for (i in 1:length(tab1[,1])){
  counts=counts+1
  abba.count[counts]=tab1$ABBA[i]
  baba.count[counts]=tab1$BABA[i]
  abba.prop[counts]=tab1$ABBA[i]/tab1$sitesUsed[i]
  baba.prop[counts]=tab1$BABA[i]/tab1$sitesUsed[i]
  abs.fdm.val[counts]=abs(tab1$fdM[i])
  fdm.val[counts]=tab1$fdM[i]
  D.val[counts]=tab1$D[i]
  abs.D.val[counts]=abs(tab1$D[i])
  chrom.log.length[counts]<-log10(max(tab1$end[tab1$scaffold==tab1$scaffold[i]])) #log length
  chrom.length[counts]<-max(tab1$end[tab1$scaffold==tab1$scaffold[i]]) #length
  chrom.number[counts]<-tab1$scaffold[i]
  sites.used[counts]<-tab1$sitesUsed[i]
  position[counts]<-tab1$mid[i]
  chr.ctr[counts]<-chrom.length[counts]/2
  per.dist.fr.ctr[counts]<-abs(position[counts]-chr.ctr[counts])/chr.ctr[counts]
  type[counts]<-tab1$type[i]
  micro.macro[counts]<-tab1$micro.macro[i]
  mid[counts]<-tab1$mid[i]
}

tab<-na.omit(data.frame(cbind(chrom.length,chrom.log.length,mid,chrom.number,sites.used,abba.count,baba.count,abba.prop,baba.prop,D.val,abs.fdm.val,fdm.val,per.dist.fr.ctr, type, micro.macro)))
tab$fdm.val<-as.numeric(tab$fdm.val)
tab$per.dist.fr.ctr<-as.numeric(tab$per.dist.fr.ctr)
tab$baba.prop<-as.numeric(tab$baba.prop)
tab$abs.fdm.val<-as.numeric(tab$abs.fdm.val)
tab$chrom.log.length<-as.numeric(tab$chrom.log.length)
tab$mid<-as.numeric(tab$mid)
tab.pos<-tab#[tab$fdm.val<=0,]

####Linear model
for(i in c(1:29,40)){
  tab2<-tab.pos[tab.pos$chrom.number==i,]
  tab2$per.dist.fr.ctr<-as.numeric(tab2$per.dist.fr.ctr)
  tab2$fdm.val<-as.numeric(tab2$fdm.val)
  tab2$baba.prop<-as.numeric(tab2$baba.prop)
  model <- lm(fdm.val ~ per.dist.fr.ctr, data=tab2)
  #summary(model)
  pval<-coef(summary(model))[2,4]
  slope<-coef(summary(model))[2,1]
  nam<-paste("p",i,sep="")
  assign(nam,ggplot(tab2,aes(x=per.dist.fr.ctr,y=fdm.val)) + ggtitle(paste("Chr",i,sep=""),subtitle = paste("p=", round(pval,digits = 2),sep = ""))+
           geom_point(alpha=0.5) +
           stat_smooth(colour="red", method="glm",se=TRUE)+
           theme_bw() + xlab("% dist. from center") + ylab("fdm"))
}


# ###polynomial model 
for(i in c(1:29,40)){
  tab2<-tab.pos[tab.pos$chrom.number==i,]
  tab2$per.dist.fr.ctr<-as.numeric(tab2$per.dist.fr.ctr)
  tab2$fdm.val<-as.numeric(tab2$fdm.val)
  tab2$abs.fdm.val<-as.numeric(tab2$abs.fdm.val)

  # Polynomial model: y = a + b*x + c*x^2
  model <- lm(fdm.val ~ poly(mid, 2, raw=TRUE), data=tab2)

  # Extract model info
  pval <- coef(summary(model))[3,4]  # p-value for x^2 term
  r_squared <- summary(model)$r.squared

  nam<-paste("q",i,sep="")
  assign(nam, ggplot(tab2, aes(x=mid, y=fdm.val)) +
           ggtitle(paste("Chr",i,sep=""),
                   subtitle = paste("p=", round(pval, digits = 3),
                                    ", R²=", round(r_squared, digits = 3), sep = "")) +
           geom_point(alpha=0.5) +
           stat_smooth(method = "lm",
                       formula = y ~ poly(x, 2, raw=TRUE),
                       colour="red", se=TRUE) +
           theme_bw() + xlab("window position") + ylab("fdm"))
}

for(i in c(1:29,40)){
  tab2<-tab.pos[tab.pos$chrom.number==i,]
  tab2$per.dist.fr.ctr<-as.numeric(tab2$per.dist.fr.ctr)
  tab2$fdm.val<-as.numeric(tab2$fdm.val)
  tab2$baba.prop<-as.numeric(tab2$baba.prop)
  
  # Remove any NA values
  tab2 <- tab2[complete.cases(tab2[c("per.dist.fr.ctr", "baba.prop")]), ]
  
  # Skip if not enough data points
  if(nrow(tab2) < 10) {
    cat("Skipping chromosome", i, "- insufficient data\n")
    next
  }
  
  # Try multiple exponential approaches
  model_fitted <- FALSE
  
  # Quadratic (often captures exponential-like behavior)
  if(!model_fitted) {
    tryCatch({
      model <- lm(fdm.val ~ mid + I(per.dist.fr.ctr^2), data=tab2)
      pval <- coef(summary(model))[3,4]  # p-value for quadratic term
      r_squared <- summary(model)$r.squared
      quad_coef <- coef(model)[3]
      
      nam<-paste("r",i,sep="")
      assign(nam, ggplot(tab2, aes(x=mid, y=fdm.val)) + 
               ggtitle(paste("Chr",i,sep=""), 
                       subtitle = paste("coef=", round(quad_coef, digits = 2), 
                                        ", p=", round(pval, digits = 3), 
                                        ", R²=", round(r_squared, digits = 3), sep = "")) +
               geom_point(alpha=0.5) +
               stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                           colour="orange", se=TRUE) +
               theme_bw() + xlab("window position") + ylab("fdm"))
      
      model_fitted <- TRUE
      cat("Quadratic fitted for chromosome", i, "\n")
      
    }, error = function(e) {})
  }
  
  }
}
#pdf(file = "../../../Figs_August2025/fdmXdist.from.ctr.T2.pdf", width = 13, height = 15)
ggarrange(p1,p2, p3,p4, 
          p5, p6, p7, p8,
          p9,p10,p11,p12,
          p13,p14,p15,p16,p17,p18,p19,p20,p21,
          p22,p23,p24,p25,p26,p27,p28,p29,p40,
          ncol = 5, nrow = 6)

#dev.off()


pdf(file = "../../../Figs_August2025/polynomial_fdmXmid.T2.pdf", width = 13, height = 15)
ggarrange(q1,q2, q3,q4, 
          q5, q6, q7, q8,
          q9,q10,q11,q12,
          q13,q14,q15,q16,q17,q18,q19,q20,q21,
          q22,q23,q24,q25,q26,q27,q28,q29,q40,
          ncol = 5, nrow = 6)

dev.off()

pdf(file = "../../../Figs_August2025/quadratic_fdmXmid.T2.pdf", width = 13, height = 15)
ggarrange(r1,r2, r3,r4, 
          r5, r6, r7, r8,
          r9,r10,r11,r12,
          r13,r14,r15,r16,r17,r18,r19,r20,
          r22,r23,r24,r25,r26,r27,r28,r29,r40,
          ncol = 5, nrow = 6)

dev.off()

mean.fdm<-c()
chrom.length<-c()
chrom.log.length<-c()
baba.prop<-c()
for(i in 1:40){
  mean.fdm[i]<-mean(tab.pos$fdm.val[tab$chrom.number==i])
  chrom.length[i]<-max(tab.pos$chrom.length[tab$chrom.number==i])
  chrom.log.length[i]<-max(tab.pos$chrom.log.length[tab$chrom.number==i])
  baba.prop[i]<-max(tab.pos$baba.prop[tab$chrom.number==i])
}
tab3<-na.omit(data.frame(cbind(mean.fdm,chrom.length,chrom.log.length, baba.prop)))
tab3$mean.fdm<-as.numeric(tab3$mean.fdm)
tab3$chrom.log.length<-as.numeric(tab3$chrom.log.length)
tab3$baba.prop<-as.numeric(tab3$baba.prop)

tab3<-tab3[-c(12,16,25,29,39),]

model <- lm(mean.fdm ~ chrom.log.length, data=tab3)
#summary(model)
pval<-coef(summary(model))[2,4]
slope<-coef(summary(model))[2,1]

ggplot(tab3,aes(x=chrom.log.length,y=mean.fdm)) +
  geom_point(alpha=0.5) +
  ggtitle("Chromosome length",subtitle = paste("p=", round(pval,digits = 2),sep = ""))+
  stat_smooth(colour="red", method="glm",se=TRUE)+
  theme_bw() + xlab("Chromosome ln(length") + ylab("fdm")

ggplot(tab3,aes(x=chrom.log.length,y=baba.prop)) +
  geom_point(alpha=0.5) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  theme_bw() + xlab("Chromosome ln(length") + ylab("BABA Proportion")

ggplot(tab.pos,aes(x=per.dist.fr.ctr,y=fdm.val)) +
  geom_point(alpha=0.5) +
  stat_smooth(colour="red", method="glm",se=TRUE)+
  theme_bw() + xlab("Chromosome ln(length") + ylab("fdM")


