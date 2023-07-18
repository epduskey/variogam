# Run cubetps model on just one area, with increasing knots

# Created: February 19, 2021
# Last modified: May 19, 2021

# Set working directory
setwd("/Users/Elizabeth/OneDrive/School/Research/Scallops")

# Load scripts
source("Code/omega.R")
source("Code/cubetps_stan.R") 

# Load packages
library(rgl)
library(rgdal)
library(mgcv)

# Contents (ctrl-f):
#	I. Load and organize data
#	II. Run Stan models
#	III. Run GAM and krige


########## I. Load and organize data ##########

# Load all data
load("Data/Corrected/HabcamDataAggregateAll.RData")

# Combine all Mid-Atlantic Bight 5000m resolution
dmv = HabcamDataAggregateAll[[16]]@data
dmv$Zone = rep("dmv", nrow(dmv))
dmv_sp = SpatialPoints(cbind(dmv$Longitude, dmv$Latitude), proj4string = CRS("+proj=longlat +ellps=WGS84"))
dmv_utm = spTransform(dmv_sp, CRS("+init=epsg:26918"))
dmv$sXutm = dmv_utm@coords[,1]
dmv$sYutm = dmv_utm@coords[,2]

# Add count data
dmv$Count = dmv$ImageDensity * dmv$sumFov

# Add a presence column to the data
dmv$Presence = as.numeric(dmv$ImageDensity > 0)

# Scale the xutm and yutm coordinates
dmv$xscale = scale(dmv$sXutm)[,1]
dmv$yscale = scale(dmv$sYutm)[,1]
head(dmv)

# Scale depth
dmv$dscale = scale(dmv$Depth)


########## II. Run Stan models ##########

# Knots
Kds = seq(5,14)
Kcs = seq(5,14)

# Run the combo model on DMV
dmv.stan5_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_5, file = "Model Output/Stan/dmv_stan55.rda")

dmv.stan6_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_5, file = "Model Output/Stan/dmv_stan65.rda")

dmv.stan7_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_5, file = "Model Output/Stan/dmv_stan75.rda")

dmv.stan8_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_5, file = "Model Output/Stan/dmv_stan85.rda")

dmv.stan9_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_5, file = "Model Output/Stan/dmv_stan95.rda")

dmv.stan10_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_5, file = "Model Output/Stan/dmv_stan105.rda")

dmv.stan11_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_5, file = "Model Output/Stan/dmv_stan115.rda")

dmv.stan12_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_5, file = "Model Output/Stan/dmv_stan125.rda")

dmv.stan13_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_5, file = "Model Output/Stan/dmv_stan135.rda")

dmv.stan14_5 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[1], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_5, file = "Model Output/Stan/dmv_stan145.rda")

dmv.stan5_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_6, file = "Model Output/Stan/dmv_stan56.rda")

dmv.stan6_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_6, file = "Model Output/Stan/dmv_stan66.rda")

dmv.stan7_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_6, file = "Model Output/Stan/dmv_stan76.rda")

dmv.stan8_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_6, file = "Model Output/Stan/dmv_stan86.rda")

dmv.stan9_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_6, file = "Model Output/Stan/dmv_stan96.rda")

dmv.stan10_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_6, file = "Model Output/Stan/dmv_stan106.rda")

dmv.stan11_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_6, file = "Model Output/Stan/dmv_stan116.rda")

dmv.stan12_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_6, file = "Model Output/Stan/dmv_stan126.rda")

dmv.stan13_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_6, file = "Model Output/Stan/dmv_stan136.rda")

dmv.stan14_6 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[2], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_6, file = "Model Output/Stan/dmv_stan146.rda")

dmv.stan5_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_7, file = "Model Output/Stan/dmv_stan57.rda")

dmv.stan6_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_7, file = "Model Output/Stan/dmv_stan67.rda")

dmv.stan7_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_7, file = "Model Output/Stan/dmv_stan77.rda")

dmv.stan8_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_7, file = "Model Output/Stan/dmv_stan87.rda")

dmv.stan9_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_7, file = "Model Output/Stan/dmv_stan97.rda")

dmv.stan10_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_7, file = "Model Output/Stan/dmv_stan107.rda")

dmv.stan11_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_7, file = "Model Output/Stan/dmv_stan117.rda")

dmv.stan12_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_7, file = "Model Output/Stan/dmv_stan127.rda")

dmv.stan13_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_7, file = "Model Output/Stan/dmv_stan137.rda")

dmv.stan14_7 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[3], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_7, file = "Model Output/Stan/dmv_stan147.rda")

dmv.stan5_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_8, file = "Model Output/Stan/dmv_stan58.rda")

dmv.stan6_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_8, file = "Model Output/Stan/dmv_stan68.rda")

dmv.stan7_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_8, file = "Model Output/Stan/dmv_stan78.rda")

dmv.stan8_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_8, file = "Model Output/Stan/dmv_stan88.rda")

dmv.stan9_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_8, file = "Model Output/Stan/dmv_stan98.rda")

dmv.stan10_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_8, file = "Model Output/Stan/dmv_stan108.rda")

dmv.stan11_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_8, file = "Model Output/Stan/dmv_stan118.rda")

dmv.stan12_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_8, file = "Model Output/Stan/dmv_stan128.rda")

dmv.stan13_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_8, file = "Model Output/Stan/dmv_stan138.rda")

dmv.stan14_8 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[4], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_8, file = "Model Output/Stan/dmv_stan148.rda")

dmv.stan5_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_9, file = "Model Output/Stan/dmv_stan59.rda")

dmv.stan6_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_9, file = "Model Output/Stan/dmv_stan69.rda")

dmv.stan7_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_9, file = "Model Output/Stan/dmv_stan79.rda")

dmv.stan8_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_9, file = "Model Output/Stan/dmv_stan89.rda")

dmv.stan9_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_9, file = "Model Output/Stan/dmv_stan99.rda")

dmv.stan10_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_9, file = "Model Output/Stan/dmv_stan109.rda")

dmv.stan11_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_9, file = "Model Output/Stan/dmv_stan119.rda")

dmv.stan12_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_9, file = "Model Output/Stan/dmv_stan129.rda")

dmv.stan13_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_9, file = "Model Output/Stan/dmv_stan139.rda")

dmv.stan14_9 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[5], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_9, file = "Model Output/Stan/dmv_stan149.rda")

dmv.stan5_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_10, file = "Model Output/Stan/dmv_stan510.rda")

dmv.stan6_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_10, file = "Model Output/Stan/dmv_stan610.rda")

dmv.stan7_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_10, file = "Model Output/Stan/dmv_stan710.rda")

dmv.stan8_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_10, file = "Model Output/Stan/dmv_stan810.rda")

dmv.stan9_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_10, file = "Model Output/Stan/dmv_stan910.rda")

dmv.stan10_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_10, file = "Model Output/Stan/dmv_stan1010.rda")

dmv.stan11_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_10, file = "Model Output/Stan/dmv_stan1110.rda")

dmv.stan12_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_10, file = "Model Output/Stan/dmv_stan1210.rda")

dmv.stan13_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_10, file = "Model Output/Stan/dmv_stan1310.rda")

dmv.stan14_10 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[6], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_10, file = "Model Output/Stan/dmv_stan1410.rda")

dmv.stan5_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_11, file = "Model Output/Stan/dmv_stan511.rda")

dmv.stan6_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_11, file = "Model Output/Stan/dmv_stan611.rda")

dmv.stan7_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_11, file = "Model Output/Stan/dmv_stan711.rda")

dmv.stan8_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_11, file = "Model Output/Stan/dmv_stan811.rda")

dmv.stan9_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_11, file = "Model Output/Stan/dmv_stan911.rda")

dmv.stan10_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_11, file = "Model Output/Stan/dmv_stan1011.rda")

dmv.stan11_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_11, file = "Model Output/Stan/dmv_stan1111.rda")

dmv.stan12_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_11, file = "Model Output/Stan/dmv_stan1211.rda")

dmv.stan13_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_11, file = "Model Output/Stan/dmv_stan1311.rda")

dmv.stan14_11 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[7], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_11, file = "Model Output/Stan/dmv_stan1411.rda")

dmv.stan5_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_12, file = "Model Output/Stan/dmv_stan512.rda")

dmv.stan6_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_12, file = "Model Output/Stan/dmv_stan612.rda")

dmv.stan7_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_12, file = "Model Output/Stan/dmv_stan712.rda")

dmv.stan8_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_12, file = "Model Output/Stan/dmv_stan812.rda")

dmv.stan9_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_12, file = "Model Output/Stan/dmv_stan912.rda")

dmv.stan10_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_12, file = "Model Output/Stan/dmv_stan1012.rda")

dmv.stan11_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_12, file = "Model Output/Stan/dmv_stan1112.rda")

dmv.stan12_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_12, file = "Model Output/Stan/dmv_stan1212.rda")

dmv.stan13_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_12, file = "Model Output/Stan/dmv_stan1312.rda")

dmv.stan14_12 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[8], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_12, file = "Model Output/Stan/dmv_stan1412.rda")

dmv.stan5_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_13, file = "Model Output/Stan/dmv_stan513.rda")

dmv.stan6_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_13, file = "Model Output/Stan/dmv_stan613.rda")

dmv.stan7_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_13, file = "Model Output/Stan/dmv_stan713.rda")

dmv.stan8_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_13, file = "Model Output/Stan/dmv_stan813.rda")

dmv.stan9_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_13, file = "Model Output/Stan/dmv_stan913.rda")

dmv.stan10_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_13, file = "Model Output/Stan/dmv_stan1013.rda")

dmv.stan11_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_13, file = "Model Output/Stan/dmv_stan1113.rda")

dmv.stan12_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_13, file = "Model Output/Stan/dmv_stan1213.rda")

dmv.stan13_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_13, file = "Model Output/Stan/dmv_stan1313.rda")

dmv.stan14_13 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[9], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_13, file = "Model Output/Stan/dmv_stan1413.rda")

dmv.stan5_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[1], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan5_14, file = "Model Output/Stan/dmv_stan514.rda")

dmv.stan6_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[2], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan6_14, file = "Model Output/Stan/dmv_stan614.rda")

dmv.stan7_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[3], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan7_14, file = "Model Output/Stan/dmv_stan714.rda")

dmv.stan8_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[4], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan8_14, file = "Model Output/Stan/dmv_stan814.rda")

dmv.stan9_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[5], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan9_14, file = "Model Output/Stan/dmv_stan914.rda")

dmv.stan10_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[6], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan10_14, file = "Model Output/Stan/dmv_stan1014.rda")

dmv.stan11_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[7], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan11_14, file = "Model Output/Stan/dmv_stan1114.rda")

dmv.stan12_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[8], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan12_14, file = "Model Output/Stan/dmv_stan1214.rda")

dmv.stan13_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[9], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan13_14, file = "Model Output/Stan/dmv_stan1314.rda")

dmv.stan14_14 = cubetps(h = dmv$dscale, px = dmv$xscale, xc = dmv$xscale, yc = dmv$yscale, x = dmv$xscale, y = dmv$yscale, z = dmv$Count, Kd = Kds[10], Kc = Kcs[10], nburn = 1000, niter = 1000, offset = log(dmv$sumFov), smush.low = 0.01, smush.hi = 0.01, plot.k = T, par = T)
save(object = dmv.stan14_14, file = "Model Output/Stan/dmv_stan1414.rda")


########## Run GAM and krige ##########

# Hurdle model
dmv.hurdle = glm(Zero ~ xscale, data = dmv, family = binomial)
stickle_hdl = stickle_sub[stickle_sub $Zero == 0, ]

# Run GAM
dmv.gam = gam(Count ~ s(Depth) + s(xscale,yscale), offset = log(sumFov), data = dmv, family = nb())
plot(dmv.gam, shade = T, shift = coef(dmv.gam)[1], trans = poisson()$linkinv)

# Krige
dat.dmv = data.frame(x = dmv$xscale, y = dmv$yscale)