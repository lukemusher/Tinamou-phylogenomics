##script adapted from Jorge Doña (University of Illinois Champaign-Urbana)

knitr::opts_chunk$set(echo = TRUE)
require("readxl")
require("ggplot2")
require("segmented")
require("knitr")

#likelihoods pulled from phylonet output log-files

likelihood<-c(-613205.7546, -612719.7746, -606942.4648, -604702.3631, -603158.9451158985, -601408.9604, -601310.8391, -599843.1398, -601772.1617, -598549.8285,-600205.1753)
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
