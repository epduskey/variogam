# Code to plot output from varioGAM simulation study

# Created: February 22, 2020
# Last modified: August 11, 2021

# Set working directory
setwd("/Users/epdus/OneDrive/School/Research/Splines")

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
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(grid)
library(cowplot)
library(geometry)

# Source the necessary scripts
source("Code/omega.R")		# Calculates penalty matrix of penalized splines
source("Code/alldata.R")		# Recreates all data sets, both with and without spatial autocorrelation
source("Code/allplot_fn.R")		# Contains all functions necessary to create the following figures

# Load Atlantic herring data
herring_full = read.table("Data/herring.txt", header = T)
herring = read.table("Data/herring_sub.txt", header = T)

# Load three-spined stickleback data
stickle_full = read.table("Data/stickle.txt", header = T)
stickle = read.table("Data/stickle_sub.txt", header = T)

# Load scallop data
mab = read.table("Data/allmab5000_corrected.txt", header = T)
dmv = mab[mab$Zone == "dmv", ]
dmv$x = dmv$sXutm
dmv$y = dmv$sYutm
dmv$z = dmv$Count

# Load grid data
load("Data/HabcamGridAll.RData")
load("Data/InfoAll.RData")

# Contents (ctrl-f; run in order):
#	0a. Common values
#	0b. Load all models
#	I. Figure 1
#	II. Figure 2
#	III. Figure 3
#	IV. Figure 4
#	V. Figure 5
#	VI. Figure 6
# VII. Figure 7
#	VII. Figure 7
#	VIII. Figure 8
# IX. Figure 9
#	X Figure S1
#	XI. Figure S2
#	XII. Figure S3
#	XIII. Figure S4
#	XIV. Figure S5
#	XV. Figure S6
#	XVI. Figure S7
#	XVII. Figure S8
#	XVIII. Figure S9
#	XIX. Figure S10
#	XX. Figure S11
#	XX. Figure S12


########## 0a. Common values ##########

# Nine shades of gray for plotting as a function of knots (from 2 to 10)
lcol_knots = c("gray10","gray20","gray30","gray40","gray50","gray60","gray70","gray80","gray90")

# Ten shades of gray for plotting as a function of peakedness (variance)
lcol_peak = c("gray10","gray19","gray28","gray37","gray46","gray55","gray64","gray73","gray82","gray91")

# Nine shades of gray for three-dimensional plotly as a function of knots (from 2 to 10)
lcol_plotly = c("rgb(25,25,25)","rgb(50,50,50)","rgb(75,75,75)","rgb(100,100,100)","rgb(125,125,125)","rgb(150,150,150)","rgb(175,175,175)","rgb(200,200,200)","rgb(225,225,225)")

# Vector of opacity values for three-dimensional plotly as a function of knots (from 2 to 10)
ovec = seq(0.95,0.3,length.out = 9)

# Sigma values
sigma_one = seq(from = 1.0, to = 0.5, length.out = 10)
sigma_two = seq(from = 0.5, to = 0.1, length.out = 10)

# Simulation magnitude values
mag_one = 10
mag_two = 20


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

# Load herring model output
load("Example/Herring/herring_stan.rda")
load("Example/Herring/herring_ic.rda")
herring.df = read.table("Example/Herring/results.txt")

# Load stickleback model output
load("Example/Stickleback/stickle_stan.rda")
load("Example/Stickleback/stickle_ic.rda")
stickle.df = read.table("Example/Stickleback/results.txt")

# Load scallop model output
load("Example/dmvstan5c.rda")
load("Example/dmvstan6c.rda")
load("Example/dmvstan7c.rda")
load("Example/dmvstan8c.rda")
load("Example/dmvstan9c.rda")
load("Example/dmvstan10c.rda")
load("Example/dmvstan11c.rda")
load("Example/dmvstan12c.rda")
load("Example/dmvstan13c.rda")
load("Example/dmvstan14c.rda")


########## I. Figure 1 ##########

# Get world map
world = ne_countries(scale = "medium", returnclass = "sf")

# Draw world map and study region boxes
world_map = ggplot(data = world) + 
  geom_sf(fill = "gray", color = "darkgray") +
  geom_rect(xmin = -80, xmax = -68, ymin = 35, ymax = 42, fill = NA, colour = "black", linewidth = 1) +
  geom_rect(xmin = 8, xmax = 32, ymin = 52, ymax = 66, fill = NA, colour = "black", linewidth = 1) +
  theme(panel.background = element_rect(fill = "azure"), panel.border = element_rect(fill = NA))

# Convert data to sf
herring_sf = st_as_sf(herring, coords = c("LogLongitude", "LogLatitude"), crs = 4326, agr = "constant")
stickle_sf = st_as_sf(stickle, coords = c("LogLongitude", "LogLatitude"), crs = 4326, agr = "constant")
scallop_sf = st_as_sf(dmv, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

# Acoustic legend points
herring_leg = st_as_sf(data.frame(x = c(11), y = c(64)), coords = c("x","y"), crs = 4326, agr = "constant")
stickle_leg = st_as_sf(data.frame(x = c(11), y = c(65)), coords = c("x","y"), crs = 4326, agr = "constant")
scallop_leg = st_as_sf(data.frame(x = c(-79), y = c(41)), coords = c("x","y"), crs = 4326, agr = "constant")

# Draw Atlantic herring and three-spined stickleback map
acoustic_map = ggplot(data = world) +
  geom_sf(fill = "gray", color = "darkgray") +
  geom_sf(data = herring_sf, size = 1, shape = 21, fill = "darkred") +
  geom_sf(data = stickle_sf, size = 1, shape = 24, fill = "lightblue") +
  geom_sf(data = herring_leg, size = 1, shape = 21, fill = "darkred") +
  geom_sf(data = stickle_leg, size = 1, shape = 24, fill = "lightblue") +
  coord_sf(xlim = c(8,32), ylim = c(52,66), expand = FALSE) + 
  annotate("text", x = 12, y = 65, label= "Herring", size = 4, hjust = 0) + 
  annotate("text", x = 12, y = 64, label= "Stickleback", size = 4, hjust = 0) + 
  theme_void() + 
  theme(panel.background = element_rect(fill = "azure"), panel.border = element_rect(linewidth = 2, fill = NA))

# Draw Atlantic sea scallop map
scallop_map = ggplot(data = world) +
  geom_sf(fill = "gray", color = "darkgray") +
  geom_sf(data = scallop_sf, size = 1, shape = 22, fill = "pink") +
  geom_sf(data = scallop_leg, size = 1, shape = 24, fill = "pink") +
  coord_sf(xlim = c(-80,-68), ylim = c(35,42), expand = FALSE) + 
  annotate("text", x = -78.6, y = 41, label= "Scallops", size = 4, hjust = 0) + 
  theme_void() + 
  theme(panel.background = element_rect(fill = "azure"), panel.border = element_rect(linewidth = 2, fill = NA))

# Draw the final map
png("Plots/figure_1.png", width = 5000, height = 2500, units = 'px', res = 600)
ggdraw(world_map) +
  draw_plot(acoustic_map, width = 0.82, height = 0.52, x = 0.30, y = 0.24) +
  draw_plot(scallop_map, width = 0.39, height = 0.33, x = -0.01, y = 0.31)
dev.off()


########## II. Figure 2 ##########

png("Plots/figure_2.png", width = 6000, height = 9000, units = 'px', res = 600)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
one1.med = medplot2(one1, sigma_one[1], 50, 10, one1.stan, gam.one[[1]], "2d")
one2.med = medplot2(one2, sigma_one[2], 60, 20, one2.stan, gam.one[[2]], "2d")
one3.med = medplot2(one3, sigma_one[3], 12, 3, one3.stan, gam.one[[3]], "2d")
one4.med = medplot2(one4, sigma_one[4], 40, 10, one4.stan, gam.one[[4]], "2d")
one5.med = medplot2(one5, sigma_one[5], 110, 25, one5.stan, gam.one[[5]], "2d")
one6.med = medplot2(one6, sigma_one[6], 15, 5, one6.stan, gam.one[[6]], "2d")
one7.med = medplot2(one7, sigma_one[7], 40, 10, one7.stan, gam.one[[7]], "2d")
one8.med = medplot2(one8, sigma_one[8], 40, 10, one8.stan, gam.one[[8]], "2d")
one9.med = medplot2(one9, sigma_one[9], 30, 10, one9.stan, gam.one[[9]], "2d")
one10.med = medplot2(one10, sigma_one[10], 20, 5, one10.stan, gam.one[[10]], "2d")
mtext("Scaled linear coordinate", side = 1, line = 4, outer = T, cex = 2)
mtext("Simulated density", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()


########## III. Figure 3 ##########

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

png("Plots/figure_3.png", width = 9000, height = 3000, units = 'px', res = 600)
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


########## IV. Figure 4 ##########

# Get median estimates
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

# Get maximum values
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

png("Plots/figure_4.png", width = 9000, height = 3000, units = 'px', res = 600)
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


########## V. Figure 5 ##########

# Plot all herring and stickleback results
png("Plots/figure_5.png", width = 7500, height = 2500, units = 'px', res = 600)
par(mfrow = c(1,4))
par(mar = c(5,3,3,1), oma = c(2,4,0.5,0.5))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "", main = "Sill", cex.lab = 2, cex.main = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2)
box()
lines(seq(2,10), herring.df[herring.df$Model == "Bayes+OK", ]$SillNorm, lwd = 2, col = "#D81B60")
points(6, herring.df[herring.df$Model == "Bayes+OK", ]$SillNorm[6-1], pch = 9, cex = 1.5, col = "#D81B60")
lines(seq(2,10), stickle.df[stickle.df$Model == "Bayes+OK", ]$SillNorm, lwd = 2, col = "#004D40")
points(6, stickle.df[stickle.df$Model == "Bayes+OK", ]$SillNorm[6-1], pch = 9, cex = 1.5, col = "#004D40")
par(mar = c(5,3,3,1))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "", main = "Range", cex.lab = 2, cex.main = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2, labels = F)
box()
lines(seq(2,10), herring.df[herring.df$Model == "Bayes+OK", ]$RangeNorm, lwd = 2, col = "#D81B60")
points(6, herring.df[herring.df$Model == "Bayes+OK", ]$RangeNorm[6-1], pch = 9, cex = 1.5, col = "#D81B60")
lines(seq(2,10), stickle.df[stickle.df$Model == "Bayes+OK", ]$RangeNorm, lwd = 2, col = "#004D40")
points(6, stickle.df[stickle.df$Model == "Bayes+OK", ]$RangeNorm[6-1], pch = 9, cex = 1.5, col = "#004D40")
par(mar = c(5,3,3,1))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "", main = "Maximum density", cex.lab = 2, cex.main = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2, labels = F)
box()
herring.max = unlist(lapply(lapply(herring.stan, getmedex, dat = herring),max))
stickle.max = unlist(lapply(lapply(stickle.stan, getmedex, dat = stickle),max))
lines(seq(2,10), herring.max/max(herring.max), lwd = 2, col = "#D81B60")
points(6, (herring.max/max(herring.max))[6-1], pch = 9, cex = 1.5, col = "#D81B60")
lines(seq(2,10), stickle.max/max(stickle.max), lwd = 2, col = "#004D40")
points(6, (stickle.max/max(stickle.max))[6-1], pch = 9, cex = 1.5, col = "#004D40")
par(mar = c(5,3,3,1))
plot(1, ylim = c(0,1), xlim = c(2,10), type = 'n', xlab = "", ylab = "", main = "LOOIC", cex.lab = 2, cex.main = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 2)
axis(2, cex.axis = 2, labels = F)
box()
herring.loo = unlist(lapply(herring.ic, function(x) x$loo$estimates[3,1]))
stickle.loo = unlist(lapply(stickle.ic, function(x) x$loo$estimates[3,1]))
lines(seq(2,10), herring.loo/max(herring.loo), lwd = 2, col = "#D81B60")
points(6, (herring.loo/max(herring.loo))[6-1], pch = 9, cex = 1.5, col = "#D81B60")
lines(seq(2,10), stickle.loo/max(stickle.loo), lwd = 2, col = "#004D40")
points(6, (stickle.loo/max(stickle.loo))[6-1], pch = 9, cex = 1.5, col = "#004D40")
mtext("Knots", side = 1, outer = T, cex = 2)
mtext("Relative to max", side = 2, outer = T, cex = 2, line = 1)
dev.off()


########## VI. Figure 6 ##########

# Data frame
scallop = data.frame(Kd = rep(seq(5,14), 10), Kc = rep(seq(5,14), each = 10))

dmv.sill = c(unlist(lapply(dmv.stan5c, getsill)),
             unlist(lapply(dmv.stan6c, getsill)),
             unlist(lapply(dmv.stan7c, getsill)),
             unlist(lapply(dmv.stan8c, getsill)),
             unlist(lapply(dmv.stan9c, getsill)),
             unlist(lapply(dmv.stan10c, getsill)),
             unlist(lapply(dmv.stan11c, getsill)),
             unlist(lapply(dmv.stan12c, getsill)),
             unlist(lapply(dmv.stan13c, getsill)),
             unlist(lapply(dmv.stan14c, getsill)))
dmv.range = c(unlist(lapply(dmv.stan5c, getrange)),
              unlist(lapply(dmv.stan6c, getrange)),
              unlist(lapply(dmv.stan7c, getrange)),
              unlist(lapply(dmv.stan8c, getrange)),
              unlist(lapply(dmv.stan9c, getrange)),
              unlist(lapply(dmv.stan10c, getrange)),
              unlist(lapply(dmv.stan11c, getrange)),
              unlist(lapply(dmv.stan12c, getrange)),
              unlist(lapply(dmv.stan13c, getrange)),
              unlist(lapply(dmv.stan14c, getrange)))
dmv.max = c(unlist(lapply(lapply(dmv.stan5c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan6c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan7c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan8c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan9c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan10c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan11c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan12c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan13c, getmedex, dat = dmv),max)),
            unlist(lapply(lapply(dmv.stan14c, getmedex, dat = dmv),max)))
#dmv5c.ic = lapply(dmv.stan5c, getic, dat = dmv, type = "3d")
#dmv6c.ic = lapply(dmv.stan6c, getic, dat = dmv, type = "3d")
#dmv7c.ic = lapply(dmv.stan7c, getic, dat = dmv, type = "3d")
#dmv8c.ic = lapply(dmv.stan8c, getic, dat = dmv, type = "3d")
#dmv9c.ic = lapply(dmv.stan9c, getic, dat = dmv, type = "3d")
#dmv10c.ic = lapply(dmv.stan10c, getic, dat = dmv, type = "3d")
#dmv11c.ic = lapply(dmv.stan11c, getic, dat = dmv, type = "3d")
#dmv12c.ic = lapply(dmv.stan12c, getic, dat = dmv, type = "3d")
#dmv13c.ic = lapply(dmv.stan13c, getic, dat = dmv, type = "3d")
#dmv14c.ic = lapply(dmv.stan14c, getic, dat = dmv, type = "3d")
#dmv.ic = list(dmv5c.ic, dmv6c.ic, dmv7c.ic, dmv8c.ic, dmv9c.ic, dmv10c.ic, dmv11c.ic, dmv12c.ic, dmv13c.ic, dmv14c.ic)
#save(object = dmv.ic, file = "Output/IC/dmvic.rda")
load("Output/IC/dmvic.rda")
plooic = function(ic) {return(ic$loo$estimates["looic","Estimate"])}
dmv.pic = c(unlist(lapply(dmv.ic[[1]], plooic)),
            unlist(lapply(dmv.ic[[2]], plooic)),
            unlist(lapply(dmv.ic[[3]], plooic)),
            unlist(lapply(dmv.ic[[4]], plooic)),
            unlist(lapply(dmv.ic[[5]], plooic)),
            unlist(lapply(dmv.ic[[6]], plooic)),
            unlist(lapply(dmv.ic[[7]], plooic)),
            unlist(lapply(dmv.ic[[8]], plooic)),
            unlist(lapply(dmv.ic[[9]], plooic)),
            unlist(lapply(dmv.ic[[10]], plooic)))
scallop$Sill = dmv.sill
scallop$Range = dmv.range
scallop$Max = dmv.max
scallop$LOOIC = dmv.pic

# Get LOO objects for loo_compare
getloo = function(obj) {return(obj$loo)}
listloo = function(obj) {out = lapply(obj, getloo); return(out)}
dmv.loolist = lapply(dmv.ic, listloo)
dmv.looall = unlist(dmv.loolist, recursive = F)
loo_compare(dmv.looall)

scallop.sill = scallop[,c("Kd","Kc","Sill")]
scallop.sill$Norm = scallop.sill$Sill/max(scallop.sill$Sill)
scallop.sill$Group = "Sill"
colnames(scallop.sill)[3] = "Value"
scallop.range = scallop[,c("Kd","Kc","Range")]
scallop.range$Norm = scallop.range$Range/max(scallop.range$Range)
scallop.range$Group = "Range"
colnames(scallop.range)[3] = "Value"
scallop.max = scallop[,c("Kd","Kc","Max")]
scallop.max$Norm = scallop.max$Max/max(scallop.max$Max)
scallop.max$Group = "Max"
colnames(scallop.max)[3] = "Value"
scallop.ic = scallop[,c("Kd","Kc","LOOIC")]
scallop.ic$Norm = scallop.ic$LOOIC/max(scallop.ic$LOOIC)
scallop.ic$Group = "LOOIC"
colnames(scallop.ic)[3] = "Value"
allop = rbind(scallop.sill, scallop.range, scallop.max, scallop.ic)
opt = data.frame(x = c(8), y = c(5))
coordinates(opt) = c("x","y")
png("Plots/figure_8.png", width = 20, height = 20, units = 'cm', res = 300)
rasterVis::levelplot(Norm ~ Kd * Kc | Group, allop, col.regions = gray.colors(100,rev=T), at = seq(0,1,by=0.01), pretty = T, xlab = list(label = expression("K"["cube"]), cex = 1.5), ylab = list(label = expression("K"["tps"]), cex = 1.5), 
                     colorkey = list(labels = list(cex = 1.2)), main = "", cex.lab = 5, cex.axis = 1.2, par.strip.text=list(cex = 1.2), par.settings = list(strip.background = list(col = "white")), scales = list(x = list(at = seq(5,14), cex = 1.2), 
                                                                                                                                                                                                                   y = list(at = seq(5,14), cex = 1.2)), strip = strip.custom(factor.levels = c("LOOIC", "Maximum density", "Range", "Sill"))) + latticeExtra::layer(sp.points(opt, pch = 9, cex = 1.2, lwd = 2, col = lcol_knots[1:3]))
dev.off()


########## VII. Figure 7 ##########

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

png("Plots/figure_7.png", width = 45, height = 15, units = 'cm', res = 300)
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


########## VIII. Figure S1 ##########

png("Plots/figure_S1.png", width = 6000, height = 9000, units = 'px', res = 600)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
sill.one1 = sillbox(one1.stan)
boxplot(sill ~ knots, sill.one1, outline = F, ylim = c(0-200*0.1,200+200*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,200,100), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 1, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.one2 = sillbox(one2.stan)
boxplot(sill ~ knots, sill.one2, outline = F, ylim = c(0-200*0.1,200+200*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,200,100), cex.axis = 2); box()
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


########## IX. Figure S2 ##########

png("Plots/figure_S2.png", width = 6000, height = 9000, units = 'px', res = 600)
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
boxplot(range ~ knots, range.one6, outline = F, ylim = c(0-0.08*0.1,0.08+0.08*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.08,0.04), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 6, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.08,"3.41")
range.one7 = rangebox(one7.stan)
boxplot(range ~ knots, range.one7, outline = F, ylim = c(0-0.1*0.1,0.1+0.1*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.1,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 7, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.1,"3.41")
range.one8 = rangebox(one8.stan)
boxplot(range ~ knots, range.one8, outline = F, ylim = c(0-0.1*0.1,0.1+0.1*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.1,0.05), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 8, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.1,"3.41")
range.one9 = rangebox(one9.stan)
boxplot(range ~ knots, range.one9, outline = F, ylim = c(0-0.2*0.1,0.2+0.2*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,0.2,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 9, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.2,"3.41")
range.one10 = rangebox(one10.stan)
boxplot(range ~ knots, range.one10, outline = F, ylim = c(0-0.8*0.1,0.8+0.8*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,0.8,0.2), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 1 & vario$Set == 10, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9,0.8,"3.41")
mtext("Knots", side = 1, line = 4, outer = T, cex = 2)
mtext("Range posterior", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()


########## X. Figure S3 ##########

# Get sill values
y2.sill = vector(mode = "list", length = 10)
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
y2.range = vector(mode = "list", length = 10)
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

png("Plots/figure_S3.png", width = 8000, height = 4000, units = 'px', res = 600)
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


########## XI. Figure S4 ##########

# Get sill values
y3.sill = vector(mode = "list", length = 10)
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
y3.range = vector(mode = "list", length = 10)
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

# Plot normalized sill and range as a function of normalized maximum density
png("Plots/figure_S4.png", width = 8000, height = 4000, units = 'px', res = 600)
par(mfrow = c(1,2))
par(mar = c(5,5,3,1), oma = c(2,2,0.5,0.5))
plot(1, type = 'n', xlim = c(0.6,1), ylim = c(0,1), xlab = "", ylab = "Sill (normalized)", cex.lab = 2, cex.axis = 1.5)
for(i in 1:10) {
  points(zm[[i]], y3.sill[[i]]/max(y3.sill[[i]]), pch = i, col = lcol_peak[[i]])
  lines(zm[[i]], predict(silltwo[[i]], type = "response"), col = lcol_peak[[i]])
}
legend("bottomleft", legend = as.expression(sapply(round(sigma_two,2), function(x) {bquote(sigma^2 == .(round(x,2)))})), pch = seq(10), col = lcol_peak, bty = 'n')
par(mar = c(5,5,3,1))
plot(1, type = 'n', xlim = c(0.6,1), ylim = c(0,1), xlab = "", ylab = "Range (normalized)", cex.lab = 2, cex.axis = 1.5)
for(i in 1:10) {
  points(zm[[i]], y3.range[[i]]/max(y3.range[[i]]), pch = i, col = lcol_peak[[i]])
  lines(zm[[i]], predict(rangetwo[[i]], type = "response"), col = lcol_peak[[i]])
}
mtext("Maximum estimated density (normalized)", side = 1, outer = T, cex = 2)
legend("bottomleft", legend = as.expression(sapply(round(sigma_two,2), function(x) {bquote(sigma^2 == .(round(x,2)))})), pch = seq(10), col = lcol_peak, bty = 'n')
dev.off()


########## XII. Figure S5 ##########

# Load all LOOIC output
load("Output/IC/oneic.rda")
load("Output/IC/twoic.rda")

# Load all LOOIC estimates
looic = read.table("Output/IC/looic.txt", header = T)

# Load loo_compare results
one_lb = read.table("Output/IC/one_lb.txt", header = T)
two_lb = read.table("Output/IC/two_lb.txt", header = T)

png("Plots/figure_S5.png", width = 7000, height = 3500, units = 'px', res = 600)
par(mfrow = c(1,2))
par(mar = c(5,5,3,1), oma = c(2,2,0.5,0.5))
plot(1, type = 'n', lwd = 2, col = lcol_peak[1], axes = F, xlab = "", ylab = "", xlim = c(2,10), ylim = c(-0.4,1), main = "2d")
axis(1, at = seq(2,10), cex.axis = 1.2)
axis(2, at = seq(0,1,by=0.2), cex.axis = 1.2)
box()
icplot(one.ic[[1]], lcol_peak[1], one_lb$lb[1], 0)
icplot(one.ic[[2]], lcol_peak[2], one_lb$lb[2], 0)
icplot(one.ic[[3]], lcol_peak[3], one_lb$lb[3], 0)
icplot(one.ic[[4]], lcol_peak[4], one_lb$lb[4], 0)
icplot(one.ic[[5]], lcol_peak[5], one_lb$lb[5], 0.05)
icplot(one.ic[[6]], lcol_peak[6], one_lb$lb[6], 0.05)
icplot(one.ic[[7]], lcol_peak[7], one_lb$lb[7], 0)
icplot(one.ic[[8]], lcol_peak[8], one_lb$lb[8], 0)
icplot(one.ic[[9]], lcol_peak[9], one_lb$lb[9], 0.1)
icplot(one.ic[[10]], lcol_peak[10], one_lb$lb[10], 0.05)
abline(0,0)
plot(1, type = 'n', lwd = 2, col = lcol_peak[1], axes = F, xlab = "", ylab = "", xlim = c(2,10), ylim = c(-0.4,1), main = "3d")
axis(1, at = seq(2,10), cex.axis = 1.2)
axis(2, at = seq(0,1,by=0.2), cex.axis = 1.2)
box()
icplot(two.ic[[1]], lcol_peak[1], two_lb$lb[1], 0)
icplot(two.ic[[2]], lcol_peak[2], two_lb$lb[2], 0)
icplot(two.ic[[3]], lcol_peak[3], two_lb$lb[3], 0.05)
icplot(two.ic[[4]], lcol_peak[4], two_lb$lb[4], 0.1)
icplot(two.ic[[5]], lcol_peak[5], two_lb$lb[5], 0.15)
icplot(two.ic[[6]], lcol_peak[6], two_lb$lb[6], 0)
icplot(two.ic[[7]], lcol_peak[7], two_lb$lb[7], 0.2)
icplot(two.ic[[8]], lcol_peak[8], two_lb$lb[8], 0.25)
icplot(two.ic[[9]], lcol_peak[9], two_lb$lb[9], 0.3)
icplot(two.ic[[10]], lcol_peak[10], two_lb$lb[10], 0.35)
abline(0,0)
mtext("Knots", side = 1, outer = T, cex = 1.5)
mtext("LOOIC (normalized)", side = 2, outer = T, cex = 1.5)
dev.off()


########## XIII. Figure S6 ##########

# Create pdat data frame
z.pdat = data.frame(xs = rep(seq(min(two1$xs),max(two1$xs),length.out=20),20), ys = rep(seq(min(two1$ys),max(two1$ys),length.out=20),each=20))

# Plot all coordinates
png("Plots/figure_S6.png", width = 6000, height = 9000, units = 'px', res = 600)
par(mar = c(0,0,0,0), oma = c(8,0,0,0))
set.panel(11,10)
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
getzmu(two1,sigma_two[1],90)
getzmu(two2,sigma_two[2],90)
getzmu(two3,sigma_two[3],90)
getzmu(two4,sigma_two[4],90)
getzmu(two5,sigma_two[5],90)
getzmu(two6,sigma_two[6],90)
getzmu(two7,sigma_two[7],90)
getzmu(two8,sigma_two[8],90)
getzmu(two9,sigma_two[9],90)
getzmu(two10,sigma_two[10],90)
box(lty = "blank")
grid.lines(x = c(0.001,0.001), y = c(grconvertY(min(two1$ys), "user", "ndc")*0.95,1), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0.5,0.5), y = c(grconvertY(min(two1$ys), "user", "ndc")*2.1,1), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0.999,0.999), y = c(grconvertY(min(two1$ys), "user", "ndc")*0.95,1), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.93)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.93)*0*(1/11), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.95)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.95)*1*(1/11), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.95)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.95)*3*(1/11), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.97)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.97)*5*(1/11), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.97)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.97)*7*(1/11), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc")*0.97)+(0.998-grconvertY(min(two1$ys), "user", "ndc")*0.97)*9*(1/11), gp = gpar(col = "black", lwd = 3))
grid.lines(x = c(0,1), y = (grconvertY(min(two1$ys), "user", "ndc"))+(0.998-grconvertY(min(two1$ys), "user", "ndc"))*1.0, gp = gpar(col = "black", lwd = 3))
set.panel()
par(oma = c(0.5,0,0,0))
image.plot(zlim = c(0,90), legend.only = T, horizontal = T, col = hcl.colors(1000,"Grays",rev = T), border = NULL, axis.args = list(cex.axis = 2))
dev.off()


########## XIV. Figure S7 ##########

png("Plots/figure_S7.png", width = 6000, height = 9000, units = 'px', res = 600)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
sill.two1 = sillbox(two1.stan)
boxplot(sill ~ knots, sill.two1, outline = F, ylim = c(0-100*0.1,100+100*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,90,30), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 1, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two2 = sillbox(two2.stan)
boxplot(sill ~ knots, sill.two2, outline = F, ylim = c(0-150*0.1,150+150*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,150,50), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 2, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two3 = sillbox(two3.stan)
boxplot(sill ~ knots, sill.two3, outline = F, ylim = c(0-100*0.1,100+100*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,90,30), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 3, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9, 100, "4801")
sill.two4 = sillbox(two4.stan)
boxplot(sill ~ knots, sill.two4, outline = F, ylim = c(0-100*0.1,100+100*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,90,30), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 4, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two5 = sillbox(two5.stan)
boxplot(sill ~ knots, sill.two5, outline = F, ylim = c(0-150*0.1,150+150*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,150,50), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 5, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two6 = sillbox(two6.stan)
boxplot(sill ~ knots, sill.two6, outline = F, ylim = c(0-100*0.1,100+100*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,90,30), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 6, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9, 100, "1081")
sill.two7 = sillbox(two7.stan)
boxplot(sill ~ knots, sill.two7, outline = F, ylim = c(0-200*0.1,200+200*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,200,100), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 7, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two8 = sillbox(two8.stan)
boxplot(sill ~ knots, sill.two8, outline = F, ylim = c(0-80*0.1,80+80*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,80,20), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 8, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two9 = sillbox(two9.stan)
boxplot(sill ~ knots, sill.two9, outline = F, ylim = c(0-80*0.1,80+80*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,80,20), cex.axis = 2); box(); points(seq(9), rep(80,9), pch = 8)
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 9, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
sill.two10 = sillbox(two10.stan)
boxplot(sill ~ knots, sill.two10, outline = F, ylim = c(0-15*0.1,15+15*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,15,5), cex.axis = 2); box(); points(seq(9), rep(15,9), pch = 8)
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 10, ]$Sill, 2), lwd = 1.5, lty = 2, col = "red")
mtext("Knots", side = 1, line = 4, outer = T, cex = 2)
mtext("Sill posterior", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()


########## XV. Figure S8 ##########

png("Plots/figure_S8.png", width = 6000, height = 9000, units = 'px', res = 600)
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
boxplot(range ~ knots, range.two7, outline = F, ylim = c(0-0.2*0.1,0.2+0.2*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.2,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 7, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
range.two8 = rangebox(two8.stan)
boxplot(range ~ knots, range.two8, outline = F, ylim = c(0-0.2*0.1,0.2+0.2*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(2, at = seq(0,0.2,0.1), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 8, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
range.two9 = rangebox(two9.stan)
boxplot(range ~ knots, range.two9, outline = F, ylim = c(0-1.0*0.1,1.0+1.0*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,1.0,0.5), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 9, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
range.two10 = rangebox(two10.stan)
boxplot(range ~ knots, range.two10, outline = F, ylim = c(0-2.0*0.1,2.0+2.0*0.1), col = "gray", xlab = "", ylab = "", axes = F); axis(1,at=seq(9),labels=seq(2,10), cex.axis = 2); axis(2, at = seq(0,2.0,0.5), cex.axis = 2); box()
lines(c(0,12), rep(vario[vario$Dim == 2 & vario$Set == 10, ]$Range, 2), lwd = 1.5, lty = 2, col = "red")
mtext("Knots", side = 1, line = 4, outer = T, cex = 2)
mtext("Range posterior", side = 2, line = 1.8, outer = T, cex = 2)
dev.off()


########## XVI. Figure S9 ##########

# Initialize plot
png("Plots/figure_S9.png", width = 6000, height = 9000, units = 'px', res = 600)
par(mfrow = c(4,2))

# Top left - all median lines for 1-d
par(mar = c(5,5,3,1))
plot(y~xs, no_one1, pch = 16, ylim = c(0-25*0.1,25+25*0.1), cex.axis = 1.5, cex.lab = 2, cex.main = 2, xlab = "Scaled linear coordinate", ylab = "Simulated density", main = "2d")
no_one1.med = lapply(no_one1.stan, getmed, type = "2d")
for(i in 1:length(no_one1.med)) {
	lines(no_one1$xs, no_one1.med[[i]], lwd = 2, col = lcol_peak[i])
}
lines(seq(-1.7,1.7,by=0.01), predict(nogam.one1, newdata = list(xs = seq(-1.7,1.7,by=0.01)), type = "response"), lwd = 2, lty = 2, col = "red")

# Top right - all median surfaces for 2-d (placeholder)
plot(1, type = 'n', axes = F, cex.main = 2, xlab = "", ylab = "", main = "3d")

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
boxplot(sill ~ knots, sill.noone1, outline = F, ylim = c(0-10*0.1,10+10*0.1), col = "gray", xlab = "Knots", ylab = "Sill posterior", cex.axis = 1.5, cex.lab = 2)
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
boxplot(sill ~ knots, sill.notwo1, outline = F, ylim = c(0-10*0.1,10+10*0.1), col = "gray", xlab = "Knots", ylab = "Sill posterior", cex.axis = 1.5, cex.lab = 2)
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
boxplot(range ~ knots, range.noone1, outline = F, ylim = c(0-0.04*0.1,0.04+0.04*0.1), col = "gray", xlab = "Knots", ylab = "Range posterior", cex.axis = 1.5, cex.lab = 2)
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
boxplot(range ~ knots, range.notwo1, outline = F, ylim = c(0-0.6*0.1,0.6+0.6*0.1), col = "gray", xlab = "Knots", ylab = "Range posterior", cex.lab = 2, axes = F)
axis(1, at = seq(9), labels = seq(2,10), cex.axis = 1.5); axis(2, at = seq(0,0.6,0.2), cex.axis = 1.5); box()
lines(c(0,25), rep(novario[novario$Dim == 2 & novario$Set == 1, ]$Range, 2), lwd = 1.5, lty = 2, col = "red"); graphics::text(9, 0.6, "4.41")

# Fourth left - normalized estimates for 1-d
par(mar = c(5,5,3,1))
y2.nosill = unlist(lapply(no_one1.stan, getsill))
y2.norange = unlist(lapply(no_one1.stan, getrange))
y2.nomed = medmax(lapply(no_one1.stan, getmed, type = "2d"))
plot(seq(2,10), y2.nosill/max(y2.nosill), type = 'l', lwd = 2, lty = 1, ylim = c(0,1), xlab = "Knots", ylab = "Normalized estimates", cex.lab = 2)
axis(1, at = seq(2,10), cex.axis = 1.5); axis(2, cex.axis = 1.5); box()
lines(seq(2,10), y2.norange/max(y2.norange), lwd = 2, lty = 2)
lines(seq(2,10), y2.nomed, lwd = 2, lty = 3)
legend("bottomleft", c("Sill", "Range", "Median max density"), lty = c(1,2,3), lwd = 2, bty = 'n')

# Fourth right - normalized estimates for 1-d
par(mar = c(5,5,3,1))
y3.nosill = unlist(lapply(no_two1.stan, getsill))
y3.norange = unlist(lapply(no_two1.stan, getrange))
y3.nomed = medmax(lapply(no_two1.stan, getmed, type = "3d"))
plot(seq(2,10), y3.nosill/max(y3.nosill), type = 'l', lwd = 2, lty = 1, ylim = c(0,1), xlab = "Knots", ylab = "Normalized estimates", cex.lab = 2, axes = F)
axis(1, at = seq(2,10), cex.axis = 1.5); axis(2, cex.axis = 1.5); box()
lines(seq(2,10), y3.norange/max(y3.norange), lwd = 2, lty = 2)
lines(seq(2,10), y3.nomed, lwd = 2, lty = 3)
legend("bottomleft", c("Sill", "Range", "Median max density"), lty = c(1,2,3), lwd = 2, bty = 'n')
dev.off()

png("Plots/nosp_twomed.png", width = 2700, height = 1800, units = 'px', res = 600)
par(mar = c(0,0,0,0), oma = c(8,0,0,0))
set.panel(2,5)
simcd(list(no_two1.stan[[1]],no_two1.stan[[2]],no_two1.stan[[3]],no_two1.stan[[4]],no_two1.stan[[5]],no_two1.stan[[6]],no_two1.stan[[7]],no_two1.stan[[8]],no_two1.stan[[9]],nogam.two1), no_two1, seq(2,10), c(rep(T,9),F), 25)
box(lty = "blank")
grid.lines(x = c(0.001,0.001), y = c(grconvertY(min(no_two1$ys)*1.12, "user", "ndc"),1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0.999,0.999), y = c(grconvertY(min(no_two1$ys)*1.12, "user", "ndc"),1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = grconvertY(min(no_two1$ys)*1.12, "user", "ndc"), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = 1, gp = gpar(col = "black", lwd = 2))
set.panel()
par(oma = c(0.5,0,0,0))
image.plot(zlim = c(0,25), legend.only = T, horizontal = T, col = hcl.colors(1000,"Grays",rev = T), border = NULL, axis.args = list(cex.axis = 1))
dev.off()


########## XVII. Figure S10 ##########

# Get predictive grid for herring and stickleback
baltic.fullgrid = read.csv("Data/acoustic_grid.csv")

# Only keep points in convex hull of observations, remove edge points
acoustic.chull = convhulln(as.matrix(herring[,c("LogLongitude","LogLatitude")]))
inds.keep = inhulln(acoustic.chull, as.matrix(baltic.fullgrid[,c("Longitude","Latitude")]))
baltic.grid = baltic.fullgrid[inds.keep,]

# Convert lat/long to UTM coordinates
baltic_sp = SpatialPoints(cbind(baltic.grid$Longitude, baltic.grid$Latitude), 
                           proj4string = CRS("+proj=longlat +ellps=WGS84"))
baltic_utm = spTransform(baltic_sp, CRS("+init=epsg:3857"))
baltic.grid$sXutm = baltic_utm@coords[,1]
baltic.grid$sYutm = baltic_utm@coords[,2]

# Scale UTM coordinates
baltic.grid$xs = (baltic.grid$sXutm - mean(stickle_full$sXutm)) / sd(stickle_full$sXutm)
baltic.grid$ys = (baltic.grid$sYutm - mean(stickle_full$sYutm)) / sd(stickle_full$sYutm)

# Make predictions
mu.stickle = vector(mode = "list", length = 9)
mu.herring = vector(mode = "list", length = 9)
for(i in 2:10) {
  
  # Get design matrices
  G.stickle = cbind(1, gcalc(stickle, baltic.grid, "3d", i))
  G.herring = cbind(1, gcalc(herring, baltic.grid, "3d", i))
  
  # Get model predictions
  mu.stickle[[i-1]] = pred_ac(G.stickle, stickle.stan[[i-1]]$ostan, stickle.stan[[i-1]]$G, stickle, baltic.grid)
  mu.herring[[i-1]] = pred_ac(G.herring, herring.stan[[i-1]]$ostan, herring.stan[[i-1]]$G, herring, baltic.grid)
  
}

# Plot density as a function of the number of knots
png("Plots/figure_S10.png", width = 6000, height = 9000, units = 'px', res = 600)
par(mar = c(0,0,0,0), oma = c(8,1,5,1))
set.panel(9,2)
accd(baltic.grid, mu.herring[[1]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[1]], 24, 0, 80000)
accd(baltic.grid, mu.herring[[2]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[2]], 24, 0, 80000)
accd(baltic.grid, mu.herring[[3]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[3]], 24, 0, 80000)
accd(baltic.grid, mu.herring[[4]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[4]], 24, 0, 80000)
accd(baltic.grid, mu.herring[[5]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[5]], 24, 0, 80000)
accd(baltic.grid, mu.herring[[6]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[6]], 24, 0, 80000)
accd(baltic.grid, mu.herring[[7]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[7]], 24, 0, 80000)
accd(baltic.grid, mu.herring[[8]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[8]], 24, 0, 80000)
accd(baltic.grid, mu.herring[[9]], 24, 0, 80000)
accd(baltic.grid, mu.stickle[[9]], 24, 0, 80000)
box(lty = "blank")
grid.lines(x = c(0,0), y = c(0.07,1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0.5,0.5), y = c(0.07,1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(1,1), y = c(0.07,1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(1,1), y = c(0.07,1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = c(0.07,0.07), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = c(1,1), gp = gpar(col = "black", lwd = 2))
set.panel()
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = c(0,1), ylim = c(0,1))
text(x = 0.22, y = 1.02, labels = "Herring", cex = 2, adj = 0.5)
text(x = 0.77, y = 1.02, labels = "Stickleback", cex = 2, adj = 0.5)
image.plot(zlim = c(0,80000), legend.only = T, horizontal = T, col = hcl.colors(1000,"inferno",rev = T), border = NULL, axis.args = list(cex.axis = 1.5))
dev.off()


########## XVIII. Figure S11 ##########

# Get predictive grid for dmv
dmv.pdat = HabcamGridAll[[which(grepl("DMV", InfoAll$PlotTitle))[1]]]@data
dmv.pdat$xscale = (dmv.pdat$sXutm - mean(mab$sXutm)) / sd(mab$sXutm)
dmv.pdat$yscale = (dmv.pdat$sYutm - mean(mab$sYutm)) / sd(mab$sYutm)
dmv.pdat$dscale = (dmv.pdat$Depth - mean(mab$Depth)) / sd(mab$Depth)

# Plot all depths
png("Plots/figure_S10.png", width = 6000, height = 9000, units = 'px', res = 600)
par(mfrow = c(5,2))
par(mar = c(0,3,0,0), oma = c(6,4,0.5,0.5))
par(mgp = c(2,1,0))
scdp(dmv.pdat, mab, dmv, dmv.stan5c, 0.8, 0.2) 
scdp(dmv.pdat, mab, dmv, dmv.stan6c, 0.8, 0.2) 
scdp(dmv.pdat, mab, dmv, dmv.stan7c, 0.8, 0.2) 
scdp(dmv.pdat, mab, dmv, dmv.stan8c, 0.8, 0.2) 
scdp(dmv.pdat, mab, dmv, dmv.stan9c, 0.8, 0.2) 
scdp(dmv.pdat, mab, dmv, dmv.stan10c, 0.8, 0.2) 
scdp(dmv.pdat, mab, dmv, dmv.stan11c, 0.8, 0.2) 
scdp(dmv.pdat, mab, dmv, dmv.stan12c, 0.8, 0.2) 
scdp(dmv.pdat, mab, dmv, dmv.stan13c, 0.8, 0.2)
axis(1, at = ((seq(40,100,by=10) - mean(mab$Depth)) / sd(mab$Depth)), labels = seq(40,100,by=10), cex.axis = 2) 
scdp(dmv.pdat, mab, dmv, dmv.stan14c, 0.8, 0.3)
axis(1, at = ((seq(40,100,by=10) - mean(mab$Depth)) / sd(mab$Depth)), labels = seq(40,100,by=10), cex.axis = 2) 
mtext("Depth (m)", side = 1, line = 4, outer = T, cex = 2)
mtext(expression(paste("Median density (m"^"-2",")")), side = 2, line = 0.6, outer = T, cex = 2)
dev.off()


########## XIX. Figure S12 ##########

# Plot all coordinates
png("Plots/figure_S12.png", width = 6000, height = 9000, units = 'px', res = 600)
par(mar = c(0,0,0,0), oma = c(8,0,0,0))
set.panel(10,10)
sccd(list(dmv.stan5c[[1]],dmv.stan6c[[1]],dmv.stan7c[[1]],dmv.stan8c[[1]],dmv.stan9c[[1]],dmv.stan5c[[2]],dmv.stan6c[[2]],dmv.stan7c[[2]],dmv.stan8c[[2]],dmv.stan9c[[2]]), dmv, dmv.pdat, c(5,5,5,5,5,6,6,6,6,6), c(5,6,7,8,9,5,6,7,8,9))
sccd(list(dmv.stan10c[[1]],dmv.stan11c[[1]],dmv.stan12c[[1]],dmv.stan13c[[1]],dmv.stan14c[[1]],dmv.stan10c[[2]],dmv.stan11c[[2]],dmv.stan12c[[2]],dmv.stan13c[[2]],dmv.stan14c[[2]]), dmv, dmv.pdat, c(5,5,5,5,5,6,6,6,6,6), c(10,11,12,13,14,10,11,12,13,14))
sccd(list(dmv.stan5c[[3]],dmv.stan6c[[3]],dmv.stan7c[[3]],dmv.stan8c[[3]],dmv.stan9c[[3]],dmv.stan5c[[4]],dmv.stan6c[[4]],dmv.stan7c[[4]],dmv.stan8c[[4]],dmv.stan9c[[4]]), dmv, dmv.pdat, c(7,7,7,7,7,8,8,8,8,8), c(5,6,7,8,9,5,6,7,8,9))
sccd(list(dmv.stan10c[[3]],dmv.stan11c[[3]],dmv.stan12c[[3]],dmv.stan13c[[3]],dmv.stan14c[[3]],dmv.stan10c[[4]],dmv.stan11c[[4]],dmv.stan12c[[4]],dmv.stan13c[[4]],dmv.stan14c[[4]]), dmv, dmv.pdat, c(7,7,7,7,7,8,8,8,8,8), c(10,11,12,13,14,10,11,12,13,14))
sccd(list(dmv.stan5c[[5]],dmv.stan6c[[5]],dmv.stan7c[[5]],dmv.stan8c[[5]],dmv.stan9c[[5]],dmv.stan5c[[6]],dmv.stan6c[[6]],dmv.stan7c[[6]],dmv.stan8c[[6]],dmv.stan9c[[6]]), dmv, dmv.pdat, c(9,9,9,9,9,10,10,10,10,10), c(5,6,7,8,9,5,6,7,8,9))
sccd(list(dmv.stan10c[[5]],dmv.stan11c[[5]],dmv.stan12c[[5]],dmv.stan13c[[5]],dmv.stan14c[[5]],dmv.stan10c[[6]],dmv.stan11c[[6]],dmv.stan12c[[6]],dmv.stan13c[[6]],dmv.stan14c[[6]]), dmv, dmv.pdat, c(9,9,9,9,9,10,10,10,10,10), c(10,11,12,13,14,10,11,12,13,14))
sccd(list(dmv.stan5c[[7]],dmv.stan6c[[7]],dmv.stan7c[[7]],dmv.stan8c[[7]],dmv.stan9c[[7]],dmv.stan5c[[8]],dmv.stan6c[[8]],dmv.stan7c[[8]],dmv.stan8c[[8]],dmv.stan9c[[8]]), dmv, dmv.pdat, c(11,11,11,11,11,12,12,12,12,12), c(5,6,7,8,9,5,6,7,8,9))
sccd(list(dmv.stan10c[[7]],dmv.stan11c[[7]],dmv.stan12c[[7]],dmv.stan13c[[7]],dmv.stan14c[[7]],dmv.stan10c[[8]],dmv.stan11c[[8]],dmv.stan12c[[8]],dmv.stan13c[[8]],dmv.stan14c[[8]]), dmv, dmv.pdat, c(11,11,11,11,11,12,12,12,12,12), c(10,11,12,13,14,10,11,12,13,14))
sccd(list(dmv.stan5c[[9]],dmv.stan6c[[9]],dmv.stan7c[[9]],dmv.stan8c[[9]],dmv.stan9c[[9]],dmv.stan5c[[10]],dmv.stan6c[[10]],dmv.stan7c[[10]],dmv.stan8c[[10]],dmv.stan9c[[10]]), dmv, dmv.pdat, c(13,13,13,13,13,14,14,14,14,14), c(5,6,7,8,9,5,6,7,8,9))
sccd(list(dmv.stan10c[[9]],dmv.stan11c[[9]],dmv.stan12c[[9]],dmv.stan13c[[9]],dmv.stan14c[[9]],dmv.stan10c[[10]],dmv.stan11c[[10]],dmv.stan12c[[10]],dmv.stan13c[[10]],dmv.stan14c[[10]]), dmv, dmv.pdat, c(13,13,13,13,13,14,14,14,14,14), c(10,11,12,13,14,10,11,12,13,14))
box(lty = "blank")
grid.lines(x = c(0.001,0.001), y = c(grconvertY(min(dmv$sYutm), "user", "ndc"),1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0.5,0.5), y = c(grconvertY(min(dmv$sYutm), "user", "ndc"),1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0.999,0.999), y = c(grconvertY(min(dmv$sYutm), "user", "ndc"),1), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = rep(grconvertY(min(dmv$sYutm), "user", "ndc"),2), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = rep(((1-grconvertY(min(dmv$sYutm), "user", "ndc"))*0.2)+grconvertY(min(dmv$sYutm), "user", "ndc"),2), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = rep(((1-grconvertY(min(dmv$sYutm), "user", "ndc"))*0.4)+grconvertY(min(dmv$sYutm), "user", "ndc"),2), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = rep(((1-grconvertY(min(dmv$sYutm), "user", "ndc"))*0.6)+grconvertY(min(dmv$sYutm), "user", "ndc"),2), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = rep(((1-grconvertY(min(dmv$sYutm), "user", "ndc"))*0.8)+grconvertY(min(dmv$sYutm), "user", "ndc"),2), gp = gpar(col = "black", lwd = 2))
grid.lines(x = c(0,1), y = 0.999, gp = gpar(col = "black", lwd = 2))
set.panel()
par(oma = c(0.5,0,0,0))
image.plot(zlim = c(0,8), legend.only = T, horizontal = T, col = hcl.colors(1000,"Grays",rev = T), border = NULL, axis.args = list(cex.axis = 2))
dev.off()