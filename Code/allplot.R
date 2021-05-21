# Code to plot output from varioGAM simulation study

# Created: May 19, 2021
# Last modified: May 21, 2021

# Set working directory
setwd(paste(mypath, "variogam", sep = ""))

# Load packages
library(mgcv)
library(cmdstanr)
library(plotly)
library(lattice)
library(loo)
library(raster)
library(rasterVis)
library(sp)
library(HDInterval)
library(fields)

# Source the necessary scripts
source("Code/omega.R")		# Calculates penalty matrix of penalized splines
source("Code/alldata.R")		# Recreates all data sets, both with and without spatial autocorrelation
source("Code/allplot_fn.R")		# Contains all functions necessary to create the following figures


# Contents (ctrl-f; run in order):
#	0a. Common values
#	0b. Load all models
#	I. Figure 1
#	II. Figure 2
#	III. Figure 3
#	IV. Figure 4
#	V. Figure 5
#	VI. Figure 6
#	VII. Figure 7
#	VIII. Figure 8
#	IX. Figure 9
#	X. Figure 10
#	XI. Figure 11
#	XII. Figure 12
#	XIII. Figure 13
#	XIV. Figure 14
#	XV. Figure 15
#	XVI. Figure 16


########## 0a. Common values ##########

# Nine shades of gray for plotting as a function of knots (from 2 to 10)
lcol_knots = c("gray10","gray20","gray30","gray40","gray50","gray60","gray70","gray80","gray90")

# Ten shades of gray for plotting as a function of peakedness (variance)
lcol_peak = c("gray10","gray19","gray28","gray37","gray46","gray55","gray64","gray73","gray82","gray91")

# Nine shades of gray for three-dimensional plotly as a function of knots (from 2 to 10)
lcol_plotly = c("rgb(25,25,25)","rgb(50,50,50)","rgb(75,75,75)","rgb(100,100,100)","rgb(125,125,125)","rgb(150,150,150)","rgb(175,175,175)","rgb(200,200,200)","rgb(225,225,225)")

# Vector of opacity values for three-dimensional plotly as a function of knots (from 2 to 10)
ovec = seq(0.95,0.3,length.out = 9)


########## 0b. Load all models ##########

# Read in GAM+OK table
vario = read.table("Output/GAM/vario.txt")

# Read in GAM output
load("Output/GAM/gamone.rda")
load("Output/GAM/gamtwo.rda")

# Load linear models for variogram results
load("Output/LM/sillone_lm.rda")
load("Output/LM/rangeone_lm.rda")
load("Output/LM/silltwo_lm.rda")
load("Output/LM/rangetwo_lm.rda")

# Load all two-dimensional model results
load("Output/Grouped/one1stan.rda")
load("Output/Grouped/one2stan.rda")
load("Output/Grouped/one3stan.rda")
load("Output/Grouped/one4stan.rda")
load("Output/Grouped/one5stan.rda")
load("Output/Grouped/one6stan.rda")
load("Output/Grouped/one7stan.rda")
load("Output/Grouped/one8stan.rda")
load("Output/Grouped/one9stan.rda")
load("Output/Grouped/one10stan.rda")

# Load all three-dimensional model results
load("Output/Grouped/two1stan.rda")
load("Output/Grouped/two2stan.rda")
load("Output/Grouped/two3stan.rda")
load("Output/Grouped/two4stan.rda")
load("Output/Grouped/two5stan.rda")
load("Output/Grouped/two6stan.rda")
load("Output/Grouped/two7stan.rda")
load("Output/Grouped/two8stan.rda")
load("Output/Grouped/two9stan.rda")
load("Output/Grouped/two10stan.rda")

# Load no autocorrelation model output
load("Output/No Autocorrelation/no_one1stan.rda")
load("Output/No Autocorrelation/no_two1stan.rda")

# Load no autocorrelation GAM model output
load("Output/GAM/gamnoone.rda")
load("Output/GAM/gamnotwo.rda")
novario = read.table("Output/GAM/novario.txt")


########## I. Figure 1 ##########

png("Plots/allmed_one.png", width = 30, height = 45, units = 'cm', res = 300)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
one1.med = medplot2(one1, 50, 10, one1.stan, gam.one[[1]], "2d")
one2.med = medplot2(one2, 60, 20, one2.stan, gam.one[[2]], "2d")
one3.med = medplot2(one3, 8, 2, one3.stan, gam.one[[3]], "2d")
one4.med = medplot2(one4, 40, 10, one4.stan, gam.one[[4]], "2d")
one5.med = medplot2(one5, 110, 25, one5.stan, gam.one[[5]], "2d")
one6.med = medplot2(one6, 15, 5, one6.stan, gam.one[[6]], "2d")
one7.med = medplot2(one7, 40, 10, one7.stan, gam.one[[7]], "2d")
one8.med = medplot2(one8, 40, 10, one8.stan, gam.one[[8]], "2d")
one9.med = medplot2(one9, 30, 10, one9.stan, gam.one[[9]], "2d")
one10.med = medplot2(one10, 20, 5, one10.stan, gam.one[[10]], "2d")
mtext("Scaled linear coordinate", side = 1, line = 4, outer = T, cex = 2)
mtext("Simulated density", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()


########## II. Figure 2 ##########

png("Plots/sillbox_one.png", width = 30, height = 45, units = 'cm', res = 300)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
sill.one1 = sillbox(one1.stan)
boxplot(sill ~ knots, sill.one1, outline = F, ylim = c(0-200*0.1,200+200*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,200,50), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 1, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one2 = sillbox(one2.stan)
boxplot(sill ~ knots, sill.one2, outline = F, ylim = c(0-200*0.1,200+200*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,200,50), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 2, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one3 = sillbox(one3.stan)
boxplot(sill ~ knots, sill.one3, outline = F, ylim = c(0-4*0.1,4+4*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,4,1), cex.axis = 2); box(); points(seq(9), rep(4,9), pch = 8)
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 3, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one4 = sillbox(one4.stan)
boxplot(sill ~ knots, sill.one4, outline = F, ylim = c(0-80*0.1,80+80*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,80,20), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 4, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one5 = sillbox(one5.stan)
boxplot(sill ~ knots, sill.one5, outline = F, ylim = c(0-600*0.1,600+600*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,600,200), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 5, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one6 = sillbox(one6.stan)
boxplot(sill ~ knots, sill.one6, outline = F, ylim = c(0-10*0.1,10+10*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,10,2), cex.axis = 2); box(); points(c(1,2), rep(10,2), pch = 8)
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 6, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one7 = sillbox(one7.stan)
boxplot(sill ~ knots, sill.one7, outline = F, ylim = c(0-60*0.1,60+60*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,60,20), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 7, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one8 = sillbox(one8.stan)
boxplot(sill ~ knots, sill.one8, outline = F, ylim = c(0-60*0.1,60+60*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,60,20), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 8, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one9 = sillbox(one9.stan)
boxplot(sill ~ knots, sill.one9, outline = F, ylim = c(0-80*0.1,80+80*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,80,20), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 9, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one10 = sillbox(one10.stan)
boxplot(sill ~ knots, sill.one10, outline = F, ylim = c(0-60*0.1,60+60*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,60,20), cex.axis = 2); box(); points(8, 60, pch = 8)
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 10, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
mtext("Knots", side = 1, line = 4, outer = T, cex = 2)
mtext("Sill posterior", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()


########## III. Figure 3 ##########

png("Plots/rangebox_one.png", width = 30, height = 45, units = 'cm', res = 300)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
range.one1 = rangebox(one1.stan)
boxplot(range ~ knots, range.one1, outline = F, ylim = c(0-0.3*0.1,0.3+0.3*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.3,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 1, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.3,"3.41")
range.one2 = rangebox(one2.stan)
boxplot(range ~ knots, range.one2, outline = F, ylim = c(0-0.15*0.1,0.15+0.15*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.15,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 2, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.15,"3.41")
range.one3 = rangebox(one3.stan)
boxplot(range ~ knots, range.one3, outline = F, ylim = c(0-0.8*0.1,0.8+0.8*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.8,0.2), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 3, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
range.one4 = rangebox(one4.stan)
boxplot(range ~ knots, range.one4, outline = F, ylim = c(0-0.15*0.1,0.15+0.15*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.15,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 4, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
range.one5 = rangebox(one5.stan)
boxplot(range ~ knots, range.one5, outline = F, ylim = c(0-0.15*0.1,0.15+0.15*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.15,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 5, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.15,"3.41")
range.one6 = rangebox(one6.stan)
boxplot(range ~ knots, range.one6, outline = F, ylim = c(0-0.08*0.1,0.08+0.08*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.08,0.02), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 6, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.08,"3.41")
range.one7 = rangebox(one7.stan)
boxplot(range ~ knots, range.one7, outline = F, ylim = c(0-0.1*0.1,0.1+0.1*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.1,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 7, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.1,"3.41")
range.one8 = rangebox(one8.stan)
boxplot(range ~ knots, range.one8, outline = F, ylim = c(0-0.1*0.1,0.1+0.1*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.1,0.02), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 8, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.1,"3.41")
range.one9 = rangebox(one9.stan)
boxplot(range ~ knots, range.one9, outline = F, ylim = c(0-0.2*0.1,0.2+0.2*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,0.2,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 9, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.2,"3.41")
range.one10 = rangebox(one10.stan)
boxplot(range ~ knots, range.one10, outline = F, ylim = c(0-0.8*0.1,0.8+0.8*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,0.8,0.2), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 10, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.8,"3.41")
mtext("Knots", side = 1, line = 4, outer = T, cex = 2)
mtext("Range posterior", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()


########## IV. Figure 4 ##########

# Get sill values
y2.sill = list()
y2.sill[[1]] = unlist(lapply(one1.stan, getsill))
y2.sill[[2]] = unlist(lapply(one2.stan, getsill))
y2.sill[[3]] = unlist(lapply(one3.stan, getsill))
y2.sill[[4]] = unlist(lapply(one4.stan, getsill))
y2.sill[[5]] = unlist(lapply(one5.stan, getsill))
y2.sill[[6]] = unlist(lapply(one6.stan, getsill))
y2.sill[[7]] = unlist(lapply(one7.stan, getsill))
y2.sill[[8]] = unlist(lapply(one8.stan, getsill))
y2.sill[[9]] = unlist(lapply(one9.stan, getsill))
y2.sill[[10]] = unlist(lapply(one10.stan, getsill))

# Get range values
y2.range = list()
y2.range[[1]] = unlist(lapply(one1.stan, getrange))
y2.range[[2]] = unlist(lapply(one2.stan, getrange))
y2.range[[3]] = unlist(lapply(one3.stan, getrange))
y2.range[[4]] = unlist(lapply(one4.stan, getrange))
y2.range[[5]] = unlist(lapply(one5.stan, getrange))
y2.range[[6]] = unlist(lapply(one6.stan, getrange))
y2.range[[7]] = unlist(lapply(one7.stan, getrange))
y2.range[[8]] = unlist(lapply(one8.stan, getrange))
y2.range[[9]] = unlist(lapply(one9.stan, getrange))
y2.range[[10]] = unlist(lapply(one10.stan, getrange))

# Get maximum values
mm1 = medmax(one1.med)
mm2 = medmax(one2.med)
mm3 = medmax(one3.med)
mm4 = medmax(one4.med)
mm5 = medmax(one5.med)
mm6 = medmax(one6.med)
mm7 = medmax(one7.med)
mm8 = medmax(one8.med)
mm9 = medmax(one9.med)
mm10 = medmax(one10.med)
mm = list(mm1, mm2, mm3, mm4, mm5, mm6, mm7, mm8, mm9, mm10)

png("Plots/msr_one.png", width = 30, height = 15, units = 'cm', res = 300)
par(mfrow = c(1,2))
par(mar = c(5,5,3,1), oma = c(2,2,0.5,0.5))
plot(1, type = 'n', xlim = c(0.3,1), ylim = c(0,1), xlab = "", ylab = "Sill (normalized)", cex.lab = 2, cex.axis = 1.5)
for(i in 1:10) {
	points(mm[[i]], y2.sill[[i]]/max(y2.sill[[i]]), pch = i, col = lcol_peak[[i]])
	lines(mm[[i]], predict(sillone[[i]], type = "response"), col = lcol_peak[[i]])
}
legend("bottomleft", legend = as.expression(sapply(round(sigma_one,2), function(x) {bquote(sigma^2 == .(round(x,2)))})), pch = seq(10), col = lcol_peak, bty = 'n')
par(mar = c(5,5,3,1))
plot(1, type = 'n', xlim = c(0.3,1), ylim = c(0,1), xlab = "", ylab = "Range (normalized)", cex.lab = 2, cex.axis = 1.5)
for(i in 1:10) {
	points(mm[[i]], y2.range[[i]]/max(y2.range[[i]]), pch = i, col = lcol_peak[[i]])
	lines(mm[[i]], predict(rangeone[[i]], type = "response"), col = lcol_peak[[i]])
}
mtext("Maximum density (normalized)", side = 1, outer = T, cex = 2)
legend("bottomleft", legend = as.expression(sapply(round(sigma_one,2), function(x) {bquote(sigma^2 == .(round(x,2)))})), pch = seq(10), col = lcol_peak, bty = 'n')
dev.off()


########## V. Figure 5 ##########

png("Plots/reltomax_one.png", width = 45, height = 15, units = 'cm', res = 300)
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0.5,0.5))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "Sill (normalized)", cex.lab = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2)
box()
sillplot(one1.stan, lcol_peak[1])
sillplot(one2.stan, lcol_peak[2])
sillplot(one3.stan, lcol_peak[3])
sillplot(one4.stan, lcol_peak[4])
sillplot(one5.stan, lcol_peak[5])
sillplot(one6.stan, lcol_peak[6])
sillplot(one7.stan, lcol_peak[7])
sillplot(one8.stan, lcol_peak[8])
sillplot(one9.stan, lcol_peak[9])
sillplot(one10.stan, lcol_peak[10])
par(mar = c(5,5,3,1))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "Range (normalized)", cex.lab = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2)
box()
rangeplot(one1.stan, lcol_peak[1])
rangeplot(one2.stan, lcol_peak[2])
rangeplot(one3.stan, lcol_peak[3])
rangeplot(one4.stan, lcol_peak[4])
rangeplot(one5.stan, lcol_peak[5])
rangeplot(one6.stan, lcol_peak[6])
rangeplot(one7.stan, lcol_peak[7])
rangeplot(one8.stan, lcol_peak[8])
rangeplot(one9.stan, lcol_peak[9])
rangeplot(one10.stan, lcol_peak[10])
par(mar = c(5,5,3,1))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "Maximum density (normalized)", cex.lab = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2)
box()
lines(seq(2,10), mm1, lwd = 2, col = lcol_peak[1])
lines(seq(2,10), mm2, lwd = 2, col = lcol_peak[2])
lines(seq(2,10), mm3, lwd = 2, col = lcol_peak[3])
lines(seq(2,10), mm4, lwd = 2, col = lcol_peak[4])
lines(seq(2,10), mm5, lwd = 2, col = lcol_peak[5])
lines(seq(2,10), mm6, lwd = 2, col = lcol_peak[6])
lines(seq(2,10), mm7, lwd = 2, col = lcol_peak[7])
lines(seq(2,10), mm8, lwd = 2, col = lcol_peak[8])
lines(seq(2,10), mm9, lwd = 2, col = lcol_peak[9])
lines(seq(2,10), mm10, lwd = 2, col = lcol_peak[10])
mtext("Knots", side = 1, outer = T, cex = 2)
dev.off()


########## VI. Figure 6 ##########

# Load all LOOIC and DIC output
load("Output/IC/oneic.rda")
load("Output/IC/twoic.rda")

# Load all LOOIC estimates
looic = read.table("Output/IC/looic.txt", header = T)

# Load loo_compare results
one_lb = read.table("Output/IC/one_lb.txt", header = T)
two_lb = read.table("Output/IC/two_lb.txt", header = T)

png("Plots/looic_all.png", width = 30, height = 15, units = 'cm', res = 300)
par(mfrow = c(1,2))
par(mar = c(5,5,3,1), oma = c(2,2,0.5,0.5))
plot(1, type = 'n', lwd = 2, col = lcol_peak[1], axes = F, xlab = "", ylab = "", xlim = c(2,10), ylim = c(-0.4,1), main = "2d")
axis(1, at = seq(2,10), cex.axis = 1.2)
axis(2, at = seq(0,1,by=0.2), cex.axis = 1.2)
box()
icplot(one.ic[[1]], lcol_peak[1], ic.knot$One_looic[1], 0)
icplot(one.ic[[2]], lcol_peak[2], ic.knot$One_looic[2], 0)
icplot(one.ic[[3]], lcol_peak[3], ic.knot$One_looic[3], 0)
icplot(one.ic[[4]], lcol_peak[4], ic.knot$One_looic[4], 0)
icplot(one.ic[[5]], lcol_peak[5], ic.knot$One_looic[5], 0.05)
icplot(one.ic[[6]], lcol_peak[6], ic.knot$One_looic[6], 0.05)
icplot(one.ic[[7]], lcol_peak[7], ic.knot$One_looic[7], 0)
icplot(one.ic[[8]], lcol_peak[8], ic.knot$One_looic[8], 0)
icplot(one.ic[[9]], lcol_peak[9], ic.knot$One_looic[9], 0.1)
icplot(one.ic[[10]], lcol_peak[10], ic.knot$One_looic[10], 0.05)
abline(0,0)
plot(1, type = 'n', lwd = 2, col = lcol_peak[1], axes = F, xlab = "", ylab = "", xlim = c(2,10), ylim = c(-0.4,1), main = "3d")
axis(1, at = seq(2,10), cex.axis = 1.2)
axis(2, at = seq(0,1,by=0.2), cex.axis = 1.2)
box()
icplot(two.ic[[1]], lcol_peak[1], ic.knot$Two_looic[1], 0)
icplot(two.ic[[2]], lcol_peak[2], ic.knot$Two_looic[2], 0)
icplot(two.ic[[3]], lcol_peak[3], ic.knot$Two_looic[3], 0.05)
icplot(two.ic[[4]], lcol_peak[4], ic.knot$Two_looic[4], 0.1)
icplot(two.ic[[5]], lcol_peak[5], ic.knot$Two_looic[5], 0.15)
icplot(two.ic[[6]], lcol_peak[6], ic.knot$Two_looic[6], 0)
icplot(two.ic[[7]], lcol_peak[7], ic.knot$Two_looic[7], 0.2)
icplot(two.ic[[8]], lcol_peak[8], ic.knot$Two_looic[8], 0.25)
icplot(two.ic[[9]], lcol_peak[9], ic.knot$Two_looic[9], 0.3)
icplot(two.ic[[10]], lcol_peak[10], ic.knot$Two_looic[10], 0.35)
abline(0,0)
mtext("Knots", side = 1, outer = T, cex = 1.5)
mtext("LOOIC (normalized)", side = 2, outer = T, cex = 1.5)
dev.off()


########## VII. Figure 7 ##########

z1 = medplot3(two1.stan, gam.two[[1]], two1, "3d")
z2 = medplot3(two2.stan, gam.two[[2]], two2, "3d")
z3 = medplot3(two3.stan, gam.two[[3]], two3, "3d")
z4 = medplot3(two4.stan, gam.two[[4]], two4, "3d")
z5 = medplot3(two5.stan, gam.two[[5]], two5, "3d")
z6 = medplot3(two6.stan, gam.two[[6]], two6, "3d")
z7 = medplot3(two7.stan, gam.two[[7]], two7, "3d")
z8 = medplot3(two8.stan, gam.two[[8]], two8, "3d")
z9 = medplot3(two9.stan, gam.two[[9]], two9, "3d")
z10 = medplot3(two10.stan, gam.two[[10]], two10, "3d")

# Create pdat data frame
z.pdat = data.frame(xs = rep(seq(min(two1$xs),max(two1$xs),length.out=20),20), ys = rep(seq(min(two1$ys),max(two1$ys),length.out=20),each=20))

# Plot all coordinates
png("Plots/allmed_two.png", width = 30, height = 45, units = 'cm', res = 300)
par(mar = c(0,0,0,0), oma = c(8,0,0,0))
set.panel(10,10)
simcd(list(two1.stan[[1]],two1.stan[[2]],two1.stan[[3]],two1.stan[[4]],two1.stan[[5]],two2.stan[[1]],two2.stan[[2]],two2.stan[[3]],two2.stan[[4]],two2.stan[[5]]), two1, c(2,3,4,5,6,2,3,4,5,6), rep(T,10), 90)
simcd(list(two1.stan[[6]],two1.stan[[7]],two1.stan[[8]],two1.stan[[9]],gam.two[[1]],two2.stan[[6]],two2.stan[[7]],two2.stan[[8]],two2.stan[[9]],gam.two[[2]]), two1, c(7,8,9,10,NA,7,8,9,10,NA), c(T,T,T,T,F,T,T,T,T,F), 90)
simcd(list(two3.stan[[1]],two3.stan[[2]],two3.stan[[3]],two3.stan[[4]],two3.stan[[5]],two4.stan[[1]],two4.stan[[2]],two4.stan[[3]],two4.stan[[4]],two4.stan[[5]]), two1, c(2,3,4,5,6,2,3,4,5,6), rep(T,10), 90)
simcd(list(two3.stan[[6]],two3.stan[[7]],two3.stan[[8]],two3.stan[[9]],gam.two[[3]],two4.stan[[6]],two4.stan[[7]],two4.stan[[8]],two4.stan[[9]],gam.two[[4]]), two1, c(7,8,9,10,NA,7,8,9,10,NA), c(T,T,T,T,F,T,T,T,T,F), 90)
simcd(list(two5.stan[[1]],two5.stan[[2]],two5.stan[[3]],two5.stan[[4]],two5.stan[[5]],two6.stan[[1]],two6.stan[[2]],two6.stan[[3]],two6.stan[[4]],two6.stan[[5]]), two1, c(2,3,4,5,6,2,3,4,5,6), rep(T,10), 90)
simcd(list(two5.stan[[6]],two5.stan[[7]],two5.stan[[8]],two5.stan[[9]],gam.two[[5]],two6.stan[[6]],two6.stan[[7]],two6.stan[[8]],two6.stan[[9]],gam.two[[6]]), two1, c(7,8,9,10,NA,7,8,9,10,NA), c(T,T,T,T,F,T,T,T,T,F), 90)
simcd(list(two7.stan[[1]],two7.stan[[2]],two7.stan[[3]],two7.stan[[4]],two7.stan[[5]],two8.stan[[1]],two8.stan[[2]],two8.stan[[3]],two8.stan[[4]],two8.stan[[5]]), two1, c(2,3,4,5,6,2,3,4,5,6), rep(T,10), 90)
simcd(list(two7.stan[[6]],two7.stan[[7]],two7.stan[[8]],two7.stan[[9]],gam.two[[7]],two8.stan[[6]],two8.stan[[7]],two8.stan[[8]],two8.stan[[9]],gam.two[[8]]), two1, c(7,8,9,10,NA,7,8,9,10,NA), c(T,T,T,T,F,T,T,T,T,F), 90)
simcd(list(two9.stan[[1]],two9.stan[[2]],two9.stan[[3]],two9.stan[[4]],two9.stan[[5]],two10.stan[[1]],two10.stan[[2]],two10.stan[[3]],two10.stan[[4]],two10.stan[[5]]), two1, c(2,3,4,5,6,2,3,4,5,6), rep(T,10), 90)
simcd(list(two9.stan[[6]],two9.stan[[7]],two9.stan[[8]],two9.stan[[9]],gam.two[[9]],two10.stan[[6]],two10.stan[[7]],two10.stan[[8]],two10.stan[[9]],gam.two[[10]]), two1, c(7,8,9,10,NA,7,8,9,10,NA), c(T,T,T,T,F,T,T,T,T,F), 90)
box(lty = "blank")
grid.lines(x = c(0.001,0.001), y = c(grconvertY(min(two1$ys), "user", "ndc")*0.95,1), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0.5,0.5), y = c(grconvertY(min(two1$ys), "user", "ndc")*0.95,1), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0.999,0.999), y = c(grconvertY(min(two1$ys), "user", "ndc")*0.95,1), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.93)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.93)*0.0, gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.95)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.95)*0.2, gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.97)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.97)*0.4, gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.97)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.97)*0.6, gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.97)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.97)*0.8, gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc"))+(0.998-grconvertY(min(two1$ys), "user", "ndc"))*1.0, gp = gpar(col = "black", lwd = 3))
set.panel()
par(oma = c(0.5,0,0,0))
image.plot(zlim = c(0,90), legend.only = T, horizontal = T, col = hcl.colors(1000,"Grays",rev = T), border = NULL, axis.args = list(cex.axis = 2))
dev.off()


########## VIII. Figure 8 ##########

png("Plots/sillbox_two.png", width = 30, height = 45, units = 'cm', res = 300)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
sill.two1 = sillbox(two1.stan)
boxplot(sill ~ knots, sill.two1, outline = F, ylim = c(0-100*0.1,100+100*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 1, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two2 = sillbox(two2.stan)
boxplot(sill ~ knots, sill.two2, outline = F, ylim = c(0-150*0.1,150+150*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 2, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two3 = sillbox(two3.stan)
boxplot(sill ~ knots, sill.two3, outline = F, ylim = c(0-100*0.1,100+100*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 3, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9, 100, "4801")
sill.two4 = sillbox(two4.stan)
boxplot(sill ~ knots, sill.two4, outline = F, ylim = c(0-100*0.1,100+100*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 4, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two5 = sillbox(two5.stan)
boxplot(sill ~ knots, sill.two5, outline = F, ylim = c(0-150*0.1,150+150*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 5, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two6 = sillbox(two6.stan)
boxplot(sill ~ knots, sill.two6, outline = F, ylim = c(0-100*0.1,100+100*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 6, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9, 100, "1081")
sill.two7 = sillbox(two7.stan)
boxplot(sill ~ knots, sill.two7, outline = F, ylim = c(0-200*0.1,200+200*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 7, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two8 = sillbox(two8.stan)
boxplot(sill ~ knots, sill.two8, outline = F, ylim = c(0-80*0.1,80+80*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 8, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two9 = sillbox(two9.stan)
boxplot(sill ~ knots, sill.two9, outline = F, ylim = c(0-80*0.1,80+80*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, cex.axis = 2); box(); points(seq(9), rep(80,9), pch = 8)
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 9, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two10 = sillbox(two10.stan)
boxplot(sill ~ knots, sill.two10, outline = F, ylim = c(0-15*0.1,15+15*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, cex.axis = 2); box(); points(seq(9), rep(15,9), pch = 8)
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 10, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
mtext("Knots", side = 1, line = 4, outer = T, cex = 2)
mtext("Sill posterior", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()

########## IX. Figure 9 ##########

png("Plots/rangebox_two.png", width = 30, height = 45, units = 'cm', res = 300)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
range.two1 = rangebox(two1.stan)
boxplot(range ~ knots, range.two1, outline = F, ylim = c(0-0.3*0.1,0.3+0.3*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.3,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 1, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.3,"4.41")
range.two2 = rangebox(two2.stan)
boxplot(range ~ knots, range.two2, outline = F, ylim = c(0-0.2*0.1,0.2+0.2*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.2,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 2, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.2,"4.41")
range.two3 = rangebox(two3.stan)
boxplot(range ~ knots, range.two3, outline = F, ylim = c(0-0.3*0.1,0.3+0.3*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.3,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 3, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.3,"2.73")
range.two4 = rangebox(two4.stan)
boxplot(range ~ knots, range.two4, outline = F, ylim = c(0-0.4*0.1,0.4+0.4*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.4,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 4, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.4,"4.41")
range.two5 = rangebox(two5.stan)
boxplot(range ~ knots, range.two5, outline = F, ylim = c(0-0.4*0.1,0.4+0.4*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.4,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 5, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.4,"4.41")
range.two6 = rangebox(two6.stan)
boxplot(range ~ knots, range.two6, outline = F, ylim = c(0-0.3*0.1,0.3+0.3*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.3,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 6, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.3,"2.50")
range.two7 = rangebox(two7.stan)
boxplot(range ~ knots, range.two7, outline = F, ylim = c(0-0.2*0.1,0.2+0.2*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.2,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 7, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
range.two8 = rangebox(two8.stan)
boxplot(range ~ knots, range.two8, outline = F, ylim = c(0-0.2*0.1,0.2+0.2*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.2,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 8, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
range.two9 = rangebox(two9.stan)
boxplot(range ~ knots, range.two9, outline = F, ylim = c(0-1.0*0.1,1.0+1.0*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,1.0,0.2), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 9, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
range.two10 = rangebox(two10.stan)
boxplot(range ~ knots, range.two10, outline = F, ylim = c(0-2.0*0.1,2.0+2.0*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,2.0,0.2), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 10, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
mtext("Knots", side = 1, line = 4, outer = T, cex = 2)
mtext("Range posterior", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()


########## X. Figure 10 ##########

# Get sill values
y3.sill = list()
y3.sill[[1]] = unlist(lapply(two1.stan, getsill))
y3.sill[[2]] = unlist(lapply(two2.stan, getsill))
y3.sill[[3]] = unlist(lapply(two3.stan, getsill))
y3.sill[[4]] = unlist(lapply(two4.stan, getsill))
y3.sill[[5]] = unlist(lapply(two5.stan, getsill))
y3.sill[[6]] = unlist(lapply(two6.stan, getsill))
y3.sill[[7]] = unlist(lapply(two7.stan, getsill))
y3.sill[[8]] = unlist(lapply(two8.stan, getsill))
y3.sill[[9]] = unlist(lapply(two9.stan, getsill))
y3.sill[[10]] = unlist(lapply(two10.stan, getsill))

# Get range values
y3.range = list()
y3.range[[1]] = unlist(lapply(two1.stan, getrange))
y3.range[[2]] = unlist(lapply(two2.stan, getrange))
y3.range[[3]] = unlist(lapply(two3.stan, getrange))
y3.range[[4]] = unlist(lapply(two4.stan, getrange))
y3.range[[5]] = unlist(lapply(two5.stan, getrange))
y3.range[[6]] = unlist(lapply(two6.stan, getrange))
y3.range[[7]] = unlist(lapply(two7.stan, getrange))
y3.range[[8]] = unlist(lapply(two8.stan, getrange))
y3.range[[9]] = unlist(lapply(two9.stan, getrange))
y3.range[[10]] = unlist(lapply(two10.stan, getrange))

zm1 = medmax(z1)
zm2 = medmax(z2)
zm3 = medmax(z3)
zm4 = medmax(z4)
zm5 = medmax(z5)
zm6 = medmax(z6)
zm7 = medmax(z7)
zm8 = medmax(z8)
zm9 = medmax(z9)
zm10 = medmax(z10)
zm = list(zm1, zm2, zm3, zm4, zm5, zm6, zm7, zm8, zm9, zm10)

# Plot normalized sill and range as a function of normalized maximum density
png("Plots/msr_two.png", width = 30, height = 15, units = 'cm', res = 300)
par(mfrow = c(1,2))
par(mar = c(5,5,3,1), oma = c(2,2,0.5,0.5))
plot(1, type = 'n', xlim = c(0.6,1), ylim = c(0,1), xlab = "", ylab = "Sill (normalized)", cex.lab = 1.5, cex.axis = 1.2)
for(i in 1:10) {
	points(zm[[i]], y3.sill[[i]]/max(y3.sill[[i]]), pch = i, col = lcol_peak[[i]])
	lines(zm[[i]], predict(silltwo[[i]], type = "response"), col = lcol_peak[[i]])
}
legend("bottomleft", legend = as.expression(sapply(round(sigma_two,2), function(x) {bquote(sigma^2 == .(round(x,2)))})), pch = seq(10), col = lcol_peak, bty = 'n')
par(mar = c(5,5,3,1))
plot(1, type = 'n', xlim = c(0.6,1), ylim = c(0,1), xlab = "", ylab = "Range (normalized)", cex.lab = 1.5, cex.axis = 1.2)
for(i in 1:10) {
	points(zm[[i]], y3.range[[i]]/max(y3.range[[i]]), pch = i, col = lcol_peak[[i]])
	lines(zm[[i]], predict(rangetwo[[i]], type = "response"), col = lcol_peak[[i]])
}
mtext("Maximum estimated density (normalized)", side = 1, outer = T, cex = 1.5)
legend("bottomleft", legend = as.expression(sapply(round(sigma_two,2), function(x) {bquote(sigma^2 == .(round(x,2)))})), pch = seq(10), col = lcol_peak, bty = 'n')
dev.off()


########## XI. Figure 11 ##########

png("Plots/reltomax_two.png", width = 45, height = 15, units = 'cm', res = 300)
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0.5,0.5))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "Sill (normalized)", cex.lab = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2)
box()
sillplot(two1.stan, lcol_peak[1])
sillplot(two2.stan, lcol_peak[2])
sillplot(two3.stan, lcol_peak[3])
sillplot(two4.stan, lcol_peak[4])
sillplot(two5.stan, lcol_peak[5])
sillplot(two6.stan, lcol_peak[6])
sillplot(two7.stan, lcol_peak[7])
sillplot(two8.stan, lcol_peak[8])
sillplot(two9.stan, lcol_peak[9])
sillplot(two10.stan, lcol_peak[10])
par(mar = c(5,5,3,1))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "Range (normalized)", cex.lab = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2)
box()
rangeplot(two1.stan, lcol_peak[1])
rangeplot(two2.stan, lcol_peak[2])
rangeplot(two3.stan, lcol_peak[3])
rangeplot(two4.stan, lcol_peak[4])
rangeplot(two5.stan, lcol_peak[5])
rangeplot(two6.stan, lcol_peak[6])
rangeplot(two7.stan, lcol_peak[7])
rangeplot(two8.stan, lcol_peak[8])
rangeplot(two9.stan, lcol_peak[9])
rangeplot(two10.stan, lcol_peak[10])
par(mar = c(5,5,3,1))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "Maximum density (normalized)", cex.lab = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2)
box()
lines(seq(2,10), zm1, lwd = 2, col = lcol_peak[1])
lines(seq(2,10), zm2, lwd = 2, col = lcol_peak[2])
lines(seq(2,10), zm3, lwd = 2, col = lcol_peak[3])
lines(seq(2,10), zm4, lwd = 2, col = lcol_peak[4])
lines(seq(2,10), zm5, lwd = 2, col = lcol_peak[5])
lines(seq(2,10), zm6, lwd = 2, col = lcol_peak[6])
lines(seq(2,10), zm7, lwd = 2, col = lcol_peak[7])
lines(seq(2,10), zm8, lwd = 2, col = lcol_peak[8])
lines(seq(2,10), zm9, lwd = 2, col = lcol_peak[9])
lines(seq(2,10), zm10, lwd = 2, col = lcol_peak[10])
mtext("Knots", side = 1, outer = T, cex = 2)
dev.off()


########## XII. Figure 12 ##########

# Initialize plot
png("Plots/nosp_all.png", width = 20, height = 30, units = 'cm', res = 300)
par(mfrow = c(4,2))

# Top left - all median lines for 1-d
par(mar = c(5,5,3,1))
plot(y~xs, no_one1, pch = 16, ylim = c(0-25*0.1,25+25*0.1), cex.axis = 1.2, cex.lab = 1.5, xlab = "Scaled linear coordinate", ylab = "Simulated density", main = "2d")
no_one1.med = lapply(no_one1.stan, getmed, type = "2d")
for(i in 1:length(no_one1.med)) {
	lines(no_one1$xs, no_one1.med[[i]], lwd = 2, col = lcol_peak[i])
}
lines(seq(-1.7,1.7,by=0.01), predict(nogam.one1, newdata = list(xs = seq(-1.7,1.7,by=0.01)), type = "response"), lwd = 2, lty = 2, col = "red")

# Top right - all median surfaces for 2-d (placeholder)
plot(1, type = 'n', axes = F, xlab = "", ylab = "", main = "3d")

# Second left - sill boxplots for 1-d
sill.noone1 = data.frame(sill = c(c(no_one1.stan[[1]]$ostan$draws("sill")), 
					c(no_one1.stan[[2]]$ostan$draws("sill")),
					c(no_one1.stan[[3]]$ostan$draws("sill")),
					c(no_one1.stan[[4]]$ostan$draws("sill")),
					c(no_one1.stan[[5]]$ostan$draws("sill")),
					c(no_one1.stan[[6]]$ostan$draws("sill")),
					c(no_one1.stan[[7]]$ostan$draws("sill")),
					c(no_one1.stan[[8]]$ostan$draws("sill")),
					c(no_one1.stan[[9]]$ostan$draws("sill"))))
sill.noone1$knots = c(rep(2, no_one1.stan[[1]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[1]]$ostan$metadata()$id)),
				rep(3, no_one1.stan[[2]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[2]]$ostan$metadata()$id)),
				rep(4, no_one1.stan[[3]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[3]]$ostan$metadata()$id)),
				rep(5, no_one1.stan[[4]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[4]]$ostan$metadata()$id)),
				rep(6, no_one1.stan[[5]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[5]]$ostan$metadata()$id)),
				rep(7, no_one1.stan[[6]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[6]]$ostan$metadata()$id)),
				rep(8, no_one1.stan[[7]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[7]]$ostan$metadata()$id)),
				rep(9, no_one1.stan[[8]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[8]]$ostan$metadata()$id)),
				rep(10, no_one1.stan[[9]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[9]]$ostan$metadata()$id)))
par(mar = c(5,5,3,1))
boxplot(sill ~ knots, sill.noone1, outline = F, ylim = c(0-10*0.1,10+10*0.1), col = "gray", xlab = "Knots", ylab = "Sill posterior", cex.axis = 1.2, cex.lab = 1.5)
lines(c(0,12), rep(novario[novario$Dim == 1 & novario$Set == 1, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")

# Second right - sill boxplots for 2-d
sill.notwo1 = data.frame(sill = c(c(no_two1.stan[[1]]$ostan$draws("sill")), 
					c(no_two1.stan[[2]]$ostan$draws("sill")),
					c(no_two1.stan[[3]]$ostan$draws("sill")),
					c(no_two1.stan[[4]]$ostan$draws("sill")),
					c(no_two1.stan[[5]]$ostan$draws("sill")),
					c(no_two1.stan[[6]]$ostan$draws("sill")),
					c(no_two1.stan[[7]]$ostan$draws("sill")),
					c(no_two1.stan[[8]]$ostan$draws("sill")),
					c(no_two1.stan[[9]]$ostan$draws("sill"))))
sill.notwo1$knots = c(rep(2, no_two1.stan[[1]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[1]]$ostan$metadata()$id)),
				rep(3, no_two1.stan[[2]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[2]]$ostan$metadata()$id)),
				rep(4, no_two1.stan[[3]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[3]]$ostan$metadata()$id)),
				rep(5, no_two1.stan[[4]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[4]]$ostan$metadata()$id)),
				rep(6, no_two1.stan[[5]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[5]]$ostan$metadata()$id)),
				rep(7, no_two1.stan[[6]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[6]]$ostan$metadata()$id)),
				rep(8, no_two1.stan[[7]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[7]]$ostan$metadata()$id)),
				rep(9, no_two1.stan[[8]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[8]]$ostan$metadata()$id)),
				rep(10, no_two1.stan[[9]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[9]]$ostan$metadata()$id)))
par(mar = c(5,5,3,1))
boxplot(sill ~ knots, sill.notwo1, outline = F, ylim = c(0-10*0.1,10+10*0.1), col = "gray", xlab = "Knots", ylab = "Sill posterior", cex.axis = 1.2, cex.lab = 1.5)
lines(c(0,12), rep(novario[novario$Dim == 2 & novario$Set == 1, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")

# Third left - range boxplots for 1-d
range.noone1 = data.frame(range = c(c(no_one1.stan[[1]]$ostan$draws("range")), 
					c(no_one1.stan[[2]]$ostan$draws("range")),
					c(no_one1.stan[[3]]$ostan$draws("range")),
					c(no_one1.stan[[4]]$ostan$draws("range")),
					c(no_one1.stan[[5]]$ostan$draws("range")),
					c(no_one1.stan[[6]]$ostan$draws("range")),
					c(no_one1.stan[[7]]$ostan$draws("range")),
					c(no_one1.stan[[8]]$ostan$draws("range")),
					c(no_one1.stan[[9]]$ostan$draws("range"))))
range.noone1$knots = c(rep(2, no_one1.stan[[1]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[1]]$ostan$metadata()$id)),
				rep(3, no_one1.stan[[2]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[2]]$ostan$metadata()$id)),
				rep(4, no_one1.stan[[3]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[3]]$ostan$metadata()$id)),
				rep(5, no_one1.stan[[4]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[4]]$ostan$metadata()$id)),
				rep(6, no_one1.stan[[5]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[5]]$ostan$metadata()$id)),
				rep(7, no_one1.stan[[6]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[6]]$ostan$metadata()$id)),
				rep(8, no_one1.stan[[7]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[7]]$ostan$metadata()$id)),
				rep(9, no_one1.stan[[8]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[8]]$ostan$metadata()$id)),
				rep(10, no_one1.stan[[9]]$ostan$metadata()$iter_sampling*length(no_one1.stan[[9]]$ostan$metadata()$id)))
par(mar = c(5,5,3,1))
boxplot(range ~ knots, range.noone1, outline = F, ylim = c(0-0.04*0.1,0.04+0.04*0.1), col = "gray", xlab = "Knots", ylab = "Range posterior", cex.axis = 1.2, cex.lab = 1.5)
lines(c(0,12), rep(novario[novario$Dim == 1 & novario$Set == 1, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9, 0.04, "3.41")

# Third right - range boxplots for 2-d
range.notwo1 = data.frame(range = c(c(no_two1.stan[[1]]$ostan$draws("range")), 
					c(no_two1.stan[[2]]$ostan$draws("range")),
					c(no_two1.stan[[3]]$ostan$draws("range")),
					c(no_two1.stan[[4]]$ostan$draws("range")),
					c(no_two1.stan[[5]]$ostan$draws("range")),
					c(no_two1.stan[[6]]$ostan$draws("range")),
					c(no_two1.stan[[7]]$ostan$draws("range")),
					c(no_two1.stan[[8]]$ostan$draws("range")),
					c(no_two1.stan[[9]]$ostan$draws("range"))))
range.notwo1$knots = c(rep(2, no_two1.stan[[1]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[1]]$ostan$metadata()$id)),
				rep(3, no_two1.stan[[2]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[2]]$ostan$metadata()$id)),
				rep(4, no_two1.stan[[3]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[3]]$ostan$metadata()$id)),
				rep(5, no_two1.stan[[4]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[4]]$ostan$metadata()$id)),
				rep(6, no_two1.stan[[5]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[5]]$ostan$metadata()$id)),
				rep(7, no_two1.stan[[6]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[6]]$ostan$metadata()$id)),
				rep(8, no_two1.stan[[7]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[7]]$ostan$metadata()$id)),
				rep(9, no_two1.stan[[8]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[8]]$ostan$metadata()$id)),
				rep(10, no_two1.stan[[9]]$ostan$metadata()$iter_sampling*length(no_two1.stan[[9]]$ostan$metadata()$id)))
par(mar = c(5,5,3,1))
boxplot(range ~ knots, range.notwo1, outline = F, ylim = c(0-0.6*0.1,0.6+0.6*0.1), col = "gray", xlab = "Knots", ylab = "Range posterior", cex.axis = 1.2, cex.lab = 1.5)
lines(c(0,25), rep(novario[novario$Dim == 2 & novario$Set == 1, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9, 0.6, "4.41")

# Fourth left - normalized estimates for 1-d
par(mar = c(5,5,3,1))
y2.nosill = unlist(lapply(no_one1.stan, getsill))
y2.norange = unlist(lapply(no_one1.stan, getrange))
y2.nomed = medmax(lapply(no_one1.stan, getmed, type = "2d"))
plot(seq(2,10), y2.nosill/max(y2.nosill), type = 'l', lwd = 2, lty = 1, ylim = c(0,1), xlab = "Knots", ylab = "Normalized estimates", cex.axis = 1.2, cex.lab = 1.5)
lines(seq(2,10), y2.norange/max(y2.norange), lwd = 2, lty = 2)
lines(seq(2,10), y2.nomed, lwd = 2, lty = 3)
legend("bottomleft", c("Sill", "Range", "Median max density"), lty = c(1,2,3), lwd = 2, bty = 'n')

# Fourth right - normalized estimates for 1-d
par(mar = c(5,5,3,1))
y3.nosill = unlist(lapply(no_two1.stan, getsill))
y3.norange = unlist(lapply(no_two1.stan, getrange))
y3.nomed = medmax(lapply(no_two1.stan, getmed, type = "3d"))
plot(seq(2,10), y3.nosill/max(y3.nosill), type = 'l', lwd = 2, lty = 1, ylim = c(0,1), xlab = "Knots", ylab = "Normalized estimates", cex.axis = 1.2, cex.lab = 1.5)
lines(seq(2,10), y3.norange/max(y3.norange), lwd = 2, lty = 2)
lines(seq(2,10), y3.nomed, lwd = 2, lty = 3)
legend("bottomleft", c("Sill", "Range", "Median max density"), lty = c(1,2,3), lwd = 2, bty = 'n')
dev.off()

png("Plots/nosp_twomed.png", width = 15, height = 10, units = 'cm', res = 300)
par(mar = c(0,0,0,0), oma = c(8,0,0,0))
set.panel(2,5)
simcd(list(no_two1.stan[[1]],no_two1.stan[[2]],no_two1.stan[[3]],no_two1.stan[[4]],no_two1.stan[[5]],no_two1.stan[[6]],no_two1.stan[[7]],no_two1.stan[[8]],no_two1.stan[[9]],nogam.two1), no_two1, seq(2,10), c(rep(T,9),F), 25)
box(lty = "blank")
grid.lines(x = c(0.001,0.001), y = c(grconvertY(min(no_two1$ys), "user", "ndc"),1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0.999,0.999), y = c(grconvertY(min(no_two1$ys), "user", "ndc"),1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = grconvertY(min(no_two1$ys), "user", "ndc"), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = 0.999, gp = gpar(col = "black", lwd = 2))
set.panel()
par(oma = c(0.5,0,0,0))
image.plot(zlim = c(0,25), legend.only = T, horizontal = T, col = hcl.colors(1000,"Grays",rev = T), border = NULL, axis.args = list(cex.axis = 1))
dev.off()


########## XIII. Figure 13 ##########

# Scallop depth plot


########## XIV. Figure 14 ##########

# Scallop coordinates plot


########## XVI. Figure 15 ##########

# Scallop color level plot


########## XVI. Figure 16 ##########

# Predict at 1000 new locations for two-dimensional data set with sigma = 1.0 and knots = 2
pdat.12 = data.frame(xs = seq(min(one1$xs), max(one1$xs),length.out=1000))
one12.pred = pred(G = gcalc(one1, pdat.12, "2d", K = 2), mod = one1.stan[[1]]$ostan, Gmod = one1.stan[[1]]$G, dat = one1, pdat = pdat.12, type = "2d")
one12.mean = apply(one12.pred$mu, 1, mean)
one12.var = apply(one12.pred$mu, 1, var)
one12.km = apply(one12.pred$kr, 1, mean)
one12.kv = apply(one12.pred$kr, 1, var)
all = one12.pred$mu + one12.pred$kr
one12.sm = apply(one12.pred$mu + one12.pred$kr, 1, mean)
one12.sv = apply(one12.pred$mu + one12.pred$kr, 1, var)

png("Plots/var_one12.png", width = 45, height = 15, units = 'cm', res = 300)
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0.5,0.5))
plot(pdat.12$xs, one12.mean, type = 'l', lwd = 2, xlab = "", ylab = "", cex.axis = 2, ylim = c(0,60))
points(pdat.12$xs, one12.sm, pch = 16, cex = 0.25)
points(one1$xs, one1$y, pch = 9, cex = 0.5)
legend("topleft", bty = 'n', lty = c(1,NA,NA), pch = c(NA,16,9), lwd = c(2,NA,NA), pt.cex = c(1,0.25,0.5), cex = 1.5, c("Spline mean", "Spline + krige mean", "Observations"))
plot(pdat.12$xs, one12.var, type = 'l', lwd = 2, xlab = "", ylab = "", cex.axis = 2, ylim = c(0,0.5))
lines(pdat.12$xs, one12.kv, type = 'l', lwd = 2, lty = 3)
legend("topleft", bty = 'n', lty = c(1,3), lwd = 2, pt.cex = c(1,1,1), cex = 1.5, c("Spline variance", "Krige variance"))
plot(pdat.12$xs, one12.sv, type = 'l', lwd = 2, xlab = "", ylab = "", cex.axis = 2, ylim = c(0,0.15))
points(one1$xs, rep(0,100), pch = 9, cex = 0.5)
legend("topleft", bty ='n', lty = c(1,NA), pch = c(NA,9), lwd = c(2,1), pt.cex = c(1,1,1), cex = 1.5, c("Spline + krige variance", "Observation locations"))
mtext("Scaled linear coordinate", outer = T, cex = 1.5, side = 1)
mtext("Variance", outer = T, cex = 1.5, side = 2)
dev.off()