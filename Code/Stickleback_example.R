# All cubic and thin plate spline estimations for acoustic data

# Created: March 23, 2023
# Last modified: March 24, 2023

# Set working directory
setwd("/Users/epdus/OneDrive/School/Research/Splines")

# Load scripts
source("Code/omega.R")
source("Code/tpsstan_examples.R")

# Load lattice
library(rgl)
library(rgdal)
library(mgcv)
library(gstat)
library(ggplot2)
library(ggpubr)

# File layout:		
#	I. Load and organize data
#	II. Run Stan models
# III. Choose best model
#	IV. Run GAM and krige
# V. Plot results


########## I. Load and organize data ##########

# Load acoustic density data
stickle = read.csv("Data/acoustic_density.csv")

# Convert lat/long to UTM coordinates
stickle_sp = SpatialPoints(cbind(stickle$LogLongitude, stickle$LogLatitude), 
				proj4string = CRS("+proj=longlat +ellps=WGS84"))
stickle_utm = spTransform(stickle_sp, CRS("+init=epsg:3857"))
stickle$sXutm = stickle_utm@coords[,1]
stickle$sYutm = stickle_utm@coords[,2]

# Add a Zero column to the data
stickle$Zero = as.numeric(stickle$StickleDensity == 0)

# Scale the xutm and yutm coordinates
stickle$xscale = scale(stickle$sXutm)[,1]
stickle$yscale = scale(stickle$sYutm)[,1]
head(stickle)

# Use only nighttime observations in SD 28
stickle_sub = subset(stickle, SubDivision == 27 & Cluster == "Night")

# Save data
write.table(stickle, "Data/stickle.txt")
write.table(stickle_sub, "Data/stickle_sub.txt")


########## II. Run Stan models ##########

# Number of knots
Ks = seq(2,10)

# Run the tps model on three-spined stickleback acoustic data
stickle.stan2 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[1], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan2, file = "Example/Stickleback/stickle_stan2.rda")

stickle.stan3 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[2], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan3, file = "Example/Stickleback/stickle_stan3.rda")

stickle.stan4 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[3], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan4, file = "Example/Stickleback/stickle_stan4.rda")

stickle.stan5 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[4], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan5, file = "Example/Stickleback/stickle_stan5.rda")

stickle.stan6 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[5], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan6, file = "Example/Stickleback/stickle_stan6.rda")

stickle.stan7 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[6], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan7, file = "Example/Stickleback/stickle_stan7.rda")

stickle.stan8 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[7], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan8, file = "Example/Stickleback/stickle_stan8.rda")

stickle.stan9 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[8], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan9, file = "Example/Stickleback/stickle_stan9.rda")

stickle.stan10 = tps(x = stickle_sub$xscale, y = stickle_sub$yscale, z = stickle_sub$StickleDensity, K = Ks[9], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = stickle.stan10, file = "Example/Stickleback/stickle_stan10.rda")

# Load model objects
load("Example/Stickleback/stickle_stan2.rda")
load("Example/Stickleback/stickle_stan3.rda")
load("Example/Stickleback/stickle_stan4.rda")
load("Example/Stickleback/stickle_stan5.rda")
load("Example/Stickleback/stickle_stan6.rda")
load("Example/Stickleback/stickle_stan7.rda")
load("Example/Stickleback/stickle_stan8.rda")
load("Example/Stickleback/stickle_stan9.rda")
load("Example/Stickleback/stickle_stan10.rda")

# Save all model objects
stickle.stan = list(stickle.stan2,
                    stickle.stan3,
                    stickle.stan4,
                    stickle.stan5,
                    stickle.stan6,
                    stickle.stan7,
                    stickle.stan8,
                    stickle.stan9,
                    stickle.stan10)
save(object = stickle.stan, file = "Example/Stickleback/stickle_stan.rda")

# Load object list
load("Example/Stickleback/stickle_stan.rda")


########## III. Choose best model ##########

# Calculate LOOIC and DIC for each model
stickle.ic = lapply(stickle.stan, getic, dat = stickle_sub, dcol = 39)

# Get information criteria
stickle.loo = lapply(stickle.ic, function(x) x$loo)
stickle.dic = lapply(stickle.ic, function(x) x$dic)

# Compare LOOIC
loo_compare(stickle.loo)

# Compare DIC
unlist(stickle.dic)

# Save output objects
save(object = stickle.ic, file = "Example/Stickleback/stickle_ic.rda")


########## III. Run GAM and krige ##########

# Data frame to store results
gam.df = data.frame(Knots = 4:10)
gam.df[,c("Sill","Range","Nugget")] = NA

# Do GAM and Krige for each knot value
for(i in 1:nrow(gam.df)) {
  
    # Temporary data set
    stickle_temp = stickle_sub
    
    # Add coordinates to data
    coordinates(stickle_temp) = c("xscale","yscale")
  
  	# Run GAM
  	stickle.gam = gam(log(StickleDensity) ~ s(xscale,yscale,k=gam.df$Knots[i]), data = stickle_temp)
  	#plot(stickle.gam, shade = T, shift = coef(stickle.gam)[1], trans = gaussian()$linkinv)
  
	  # Krige residuals
  	stickle_temp$resid = stickle_temp$StickleDensity - exp(stickle.gam$fitted)
  	stickle.vario = variogram(resid ~ 1, stickle_temp)
  	#plot(gamma ~ dist, stickle.vario, main = "vario")
  	fit = fit.variogram.gls(resid ~ 1, stickle_temp, vgm(nugget = 0, psill = 1e8, range = 0.2, model = "Exp"), maxiter = 0)
  
  	# Store results
  	gam.df[i,]$Sill = fit$psill[1]
  	gam.df[i,]$Range = fit$range[2]
  	gam.df[i,]$Nugget = fit$psill[2]
  
}

# Normalize
gam.df$SillNorm = gam.df$Sill / max(gam.df$Sill)
gam.df$RangeNorm = gam.df$Range / max(gam.df$Range)
gam.df$NuggetNorm = gam.df$Nugget / max(gam.df$Nugget)


########## IV. Plot results ##########

# Results data frame
stickle.df = data.frame(Knots = Ks)
stickle.df[,c("Sill","Range","Nugget")] = NA

# Get variogram parameters
stickle.df[1,2:4] = as.vector(stickle.stan2$ostan$summary(c("sill","range","nugget"))$median)
stickle.df[2,2:4] = as.vector(stickle.stan3$ostan$summary(c("sill","range","nugget"))$median)
stickle.df[3,2:4] = as.vector(stickle.stan4$ostan$summary(c("sill","range","nugget"))$median)
stickle.df[4,2:4] = as.vector(stickle.stan5$ostan$summary(c("sill","range","nugget"))$median)
stickle.df[5,2:4] = as.vector(stickle.stan6$ostan$summary(c("sill","range","nugget"))$median)
stickle.df[6,2:4] = as.vector(stickle.stan7$ostan$summary(c("sill","range","nugget"))$median)
stickle.df[7,2:4] = as.vector(stickle.stan8$ostan$summary(c("sill","range","nugget"))$median)
stickle.df[8,2:4] = as.vector(stickle.stan9$ostan$summary(c("sill","range","nugget"))$median)
stickle.df[9,2:4] = as.vector(stickle.stan10$ostan$summary(c("sill","range","nugget"))$median)

# Normalize
stickle.df$SillNorm = stickle.df$Sill / max(stickle.df$Sill)
stickle.df$RangeNorm = stickle.df$Range / max(stickle.df$Range)
stickle.df$NuggetNorm = stickle.df$Nugget / max(stickle.df$Nugget)

# Combine data
all.df = data.frame(Knots = c(stickle.df$Knots,gam.df$Knots))
all.df$Sill = c(stickle.df$Sill,gam.df$Sill)
all.df$Range = c(stickle.df$Range,gam.df$Range)
all.df$Nugget = c(stickle.df$Nugget,gam.df$Nugget)
all.df$SillNorm = c(stickle.df$SillNorm,gam.df$SillNorm)
all.df$RangeNorm = c(stickle.df$RangeNorm,gam.df$RangeNorm)
all.df$NuggetNorm = c(stickle.df$NuggetNorm,gam.df$NuggetNorm)
all.df$Model = c(rep("Bayes+OK",nrow(stickle.df)), rep("GAM+OK",nrow(gam.df)))

# Save results
write.table(all.df, "Example/Stickleback/results.txt")

# Plot three-spined stickleback results
p1 = ggplot(all.df, aes(x=Knots,y=SillNorm,color=Model)) + 
  	ylim(0,1) + xlab("") + ylab("Relative to max") +
  	geom_line() + 
  	geom_point() +
  	theme_bw() + theme(text = element_text(size = 20)) +
  	ggtitle("Partial Sill")
p2 = ggplot(all.df, aes(x=Knots,y=RangeNorm,color=Model)) + 
  	ylim(0,1) + xlab("Knots") + ylab("") +
  	geom_line() + 
  	geom_point() +
  	theme_bw() + theme(text = element_text(size = 20)) +
  	ggtitle("Range")
p3 = ggplot(all.df, aes(x=Knots,y=NuggetNorm,color=Model)) + 
  	ylim(0,1) + xlab("") + ylab("") +
  	geom_line() + 
  	geom_point() +
  	theme_bw() + theme(text = element_text(size = 20)) +
  	ggtitle("Nugget")
ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = "right")

# GAM and Bayes variogram values
b1 = ggplot(all.df, aes(x=Model,y=Sill)) +
  	geom_bar(stat="summary") + xlab("") + ylab("Value") + scale_y_log10() +
  	theme_bw() + theme(text = element_text(size = 20)) +
  	ggtitle("Sill")
b2 = ggplot(all.df, aes(x=Model,y=Range)) +
  	geom_bar(stat="summary") + xlab("Model") + ylab("") +
  	theme_bw() + theme(text = element_text(size = 20)) +
  	ggtitle("Range")
b3 = ggplot(all.df, aes(x=Model,y=Nugget)) +
  	geom_bar(stat="summary") + xlab("") + ylab("") +
  	theme_bw() + theme(text = element_text(size = 20)) +
  	ggtitle("Nugget")
ggarrange(b1, b2, b3, nrow = 1, common.legend = T, legend = "right")