##script adapted from Jorge Doña (University of Illinois Champaign-Urbana)

knitr::opts_chunk$set(echo = TRUE)
require("readxl")
require("ggplot2")
require("segmented")
require("knitr")

#likelihoods pulled from phylonet output log-files

m<-0:10
likelihood.wg<-c(-443513.9631,-443331.5423,-443287.7948,-443189.7785,-443156.1047,-443138.3756,-443155.4182,-443162.5741,-443157.1914,-443157.7499,-443158.3962)
likelihood.z<-c(-20446.28941,-20195.90265,-20168.981104789356,-20152.36069,-20100.39856,-20111.74738,-20103.41958,-20108.44324,-20106.14877,-20100.29237,-20110.99554)

data<-data.frame(cbind(m,likelihood.wg,likelihood.z))

par(mfrow=c(2,1))

par(mfrow=c(1,1)+0.3)  # Leave space for z axis
plot(m, likelihood.wg, type='l',lwd=2, main="Network likelihoods", xlab="# migration edges", ylab="likelihood (autosomes)") # first plot
par(new = TRUE)
plot(m, likelihood.z, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd=2, col='red')
axis(side=4, at = pretty(range(likelihood.z)), col="red")
mtext("likelihood (Z-chromosome)", side=4, line=3, col = "red")

# Fit a linear model
my.lm <- lm(likelihood.wg ~ m)
summary(my.lm)

# Fit a segmented linear model
my.seg <- segmented(my.lm, seg.Z = ~ m)

# Output the summary of the segmented model
summary(my.seg)

# Extract breakpoints and slopes
my.seg$psi

my.slopes <- coef(my.seg)

# Compute coefficients for the lines
b0 <- coef(my.seg)[[1]]
b1 <- coef(my.seg)[[2]]
c1 <- coef(my.seg)[[2]] + coef(my.seg)[[3]]
break1 <- round(my.seg$psi[, 2])
c0 <- b0 + b1 * break1 - c1 * break1
d1 <- coef(my.seg)[[4]] + c1
break_point <- round(my.seg$psi[, 2])

a<-ggplot(data, aes(x = m, y = likelihood.wg)) + 
  geom_point(size=2) + theme_classic() +
  scale_x_continuous(labels=as.character(data$number_of_max_reticulations), breaks=data$number_of_max_reticulations) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +
  geom_vline(xintercept=break_point, linetype='dashed') + ylab("likelihood (autosomes)")+
  geom_abline(intercept = b0, slope = b1, aes(colour = "first part"), show.legend = TRUE) +
  geom_abline(intercept = c0, slope = c1, aes(colour = "second part"), show.legend = TRUE)

#Chr Z
# Fit a linear model
my.lm <- lm(likelihood.z ~ m)
summary(my.lm)

# Fit a segmented linear model
my.seg <- segmented(my.lm, seg.Z = ~ m)

# Output the summary of the segmented model
summary(my.seg)

# Extract breakpoints and slopes
my.seg$psi

my.slopes <- coef(my.seg)

# Compute coefficients for the lines
b0 <- coef(my.seg)[[1]]
b1 <- coef(my.seg)[[2]]
c1 <- coef(my.seg)[[2]] + coef(my.seg)[[3]]
break1 <- round(my.seg$psi[, 2])
c0 <- b0 + b1 * break1 - c1 * break1
d1 <- coef(my.seg)[[4]] + c1
break_point <- round(my.seg$psi[, 2])

b<-ggplot(data, aes(x = m, y = likelihood.z)) + 
  geom_point(size=2) + theme_classic() +
  scale_x_continuous(labels=as.character(data$number_of_max_reticulations), breaks=data$number_of_max_reticulations) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +
  geom_vline(xintercept=break_point, linetype='dashed') + ylab("likelihood (Z-chromosome)")+
  geom_abline(intercept = b0, slope = b1, aes(colour = "first part"), show.legend = TRUE) +
  geom_abline(intercept = c0, slope = c1, aes(colour = "second part"), show.legend = TRUE)

a
b
library(ggpubr)
ggarrange(a,b,ncol=1,nrow=2)

