#script adapted from Jorge Doña (University of Illinois Champaign-Urbana)

knitr::opts_chunk$set(echo = TRUE)
require("readxl")
require("ggplot2")
require("segmented")
require("knitr")

#likelihoods pulled from phylonet output log-files

likelihood<-c(-612221.1008,-609081.5674,-591372.2571,-588720.9330807882,-587729.4264,-587873.2769,-580148.2359,-580771.36826, -582384.8279, -580750.6224990899,-578359.019804455)
m<-0:10

data<-data.frame(cbind(m,likelihood))
plot(m,likelihood, type = 'l')

# Fit a linear model
my.lm <- lm(likelihood ~ m, data = data)
summary(my.lm)
abline(my.lm)

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

ggplot(data, aes(x = m, y = likelihood)) + 
  geom_point(size=2) + theme_classic() +
  scale_x_continuous(labels=as.character(data$number_of_max_reticulations), breaks=data$number_of_max_reticulations) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +
  geom_vline(xintercept=break_point, linetype='dashed') +
  geom_abline(intercept = b0, slope = b1, aes(colour = "first part"), show.legend = TRUE) +
  geom_abline(intercept = c0, slope = c1, aes(colour = "second part"), show.legend = TRUE)
