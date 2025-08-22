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


