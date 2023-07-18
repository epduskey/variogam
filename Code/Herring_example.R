# All cubic and thin plate spline estimations for acoustic data

# Created: March 23, 2023
# Last modified: June 30, 2023

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
herring = read.csv("Data/acoustic_density.csv")

# Convert lat/long to UTM coordinates
herring_sp = SpatialPoints(cbind(herring$LogLongitude, herring$LogLatitude), 
				proj4string = CRS("+proj=longlat +ellps=WGS84"))
herring_utm = spTransform(herring_sp, CRS("+init=epsg:3857"))
herring$sXutm = herring_utm@coords[,1]
herring$sYutm = herring_utm@coords[,2]

# Add a Zero column to the data
herring$Zero = as.numeric(herring$HerringDensity == 0)

# Scale the xutm and yutm coordinates
herring$xscale = scale(herring $sXutm)[,1]
herring$yscale = scale(herring $sYutm)[,1]
head(herring)

# Use only nighttime observations in SD 28
herring_sub = subset(herring, SubDivision == 27 & Cluster == "Night")

# Save data
write.table(herring, "Data/herring.txt")
write.table(herring_sub, "Data/herring_sub.txt")


########## II, Run Stan models ##########

# Number of knots
Ks = seq(2,10)

# Run the tps model on herring acoustic data
herring.stan2 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[1], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan2, file = "Example/Herring/herring_stan2.rda")

herring.stan3 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[2], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan3, file = "Example/Herring/herring_stan3.rda")

herring.stan4 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[3], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan4, file = "Example/Herring/herring_stan4.rda")

herring.stan5 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[4], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan5, file = "Example/Herring/herring_stan5.rda")

herring.stan6 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[5], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan6, file = "Example/Herring/herring_stan6.rda")

herring.stan7 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[6], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan7, file = "Example/Herring/herring_stan7.rda")

herring.stan8 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[7], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan8, file = "Example/Herring/herring_stan8.rda")

herring.stan9 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[8], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan9, file = "Example/Herring/herring_stan9.rda")

herring.stan10 = tps(x = herring_sub$xscale, y = herring_sub$yscale, z = herring_sub$HerringDensity, K = Ks[9], rorsd = 1, nchains = 6, nburn = 1000, niter = 1000, inc = 100, par = 6)
save(object = herring.stan10, file = "Example/Herring/herring_stan10.rda")

# Load model objects
load("Example/Herring/herring_stan2.rda")
load("Example/Herring/herring_stan3.rda")
load("Example/Herring/herring_stan4.rda")
load("Example/Herring/herring_stan5.rda")
load("Example/Herring/herring_stan6.rda")
load("Example/Herring/herring_stan7.rda")
load("Example/Herring/herring_stan8.rda")
load("Example/Herring/herring_stan9.rda")
load("Example/Herring/herring_stan10.rda")

# Save all objects as list
herring.stan = list(herring.stan2, 
                    herring.stan3, 
                    herring.stan4, 
                    herring.stan5, 
                    herring.stan6, 
                    herring.stan7, 
                    herring.stan8, 
                    herring.stan9, 
                    herring.stan10)
save(object = herring.stan, file = "Example/Herring/herring_stan.rda")

# Load object list
load("Example/Herring/herring_stan.rda")


########## III. Choose best model ##########

# Calculate LOOIC and DIC for each model
herring.ic = lapply(herring.stan, getic, dat = herring_sub, dcol = 38)

# Get information criteria
herring.loo = lapply(herring.ic, function(x) x$loo)
herring.dic = lapply(herring.ic, function(x) x$dic)

# Compare LOOIC
loo_compare(herring.loo)

# Compare DIC
unlist(herring.dic)

# Save output objects
save(object = herring.ic, file = "Example/Herring/herring_ic.rda")


########## IV. Run GAM and Krige ##########

# Data frame to store results
gam.df = data.frame(Knots = 4:10)
gam.df[,c("Sill","Range","Nugget")] = NA

# Do GAM and Krige for each knot value
for(i in 1:nrow(gam.df)) {
  
  # Temporary data set
  herring_temp = herring_sub
  
  # Add coordinates to data
  coordinates(herring_temp) = c("xscale","yscale")
  
  # Run GAM
  herring.gam = gam(log(HerringDensity) ~ s(xscale,yscale,k=gam.df$Knots[i]), data = herring_temp)
  #plot(herring.gam, shade = T, shift = coef(herring.gam)[1], trans = gaussian()$linkinv)
  
  # Krige residuals
  herring_temp$resid = herring_temp$HerringDensity - exp(herring.gam$fitted)
  herring.vario = variogram(resid ~ 1, herring_temp)
  #plot(gamma ~ dist, herring.vario, main = "vario")
  fit = fit.variogram.gls(resid ~ 1, herring_temp, vgm(nugget = 0, psill = 5e9, range = 0.2, model = "Exp"), maxiter = 0)
  
  # Store results
  gam.df[i,]$Sill = fit$psill[1]
  gam.df[i,]$Range = fit$range[2]
  gam.df[i,]$Nugget = fit$psill[2]
  
}

# Normalize
gam.df$SillNorm = gam.df$Sill / max(gam.df$Sill)
gam.df$RangeNorm = gam.df$Range / max(gam.df$Range)
gam.df$NuggetNorm = gam.df$Nugget / max(gam.df$Nugget)


########## V. Plot results ##########

# Results data frame
herring.df = data.frame(Knots = Ks)
herring.df[,c("Sill","Range","Nugget")] = NA

# Get variogram parameters
herring.df[1,2:4] = as.vector(herring.stan2$ostan$summary(c("sill","range","nugget"))$median)
herring.df[2,2:4] = as.vector(herring.stan3$ostan$summary(c("sill","range","nugget"))$median)
herring.df[3,2:4] = as.vector(herring.stan4$ostan$summary(c("sill","range","nugget"))$median)
herring.df[4,2:4] = as.vector(herring.stan5$ostan$summary(c("sill","range","nugget"))$median)
herring.df[5,2:4] = as.vector(herring.stan6$ostan$summary(c("sill","range","nugget"))$median)
herring.df[6,2:4] = as.vector(herring.stan7$ostan$summary(c("sill","range","nugget"))$median)
herring.df[7,2:4] = as.vector(herring.stan8$ostan$summary(c("sill","range","nugget"))$median)
herring.df[8,2:4] = as.vector(herring.stan9$ostan$summary(c("sill","range","nugget"))$median)
herring.df[9,2:4] = as.vector(herring.stan10$ostan$summary(c("sill","range","nugget"))$median)

# Normalize
herring.df$SillNorm = herring.df$Sill / max(herring.df$Sill)
herring.df$RangeNorm = herring.df$Range / max(herring.df$Range)
herring.df$NuggetNorm = herring.df$Nugget / max(herring.df$Nugget)

# Combine data
all.df = data.frame(Knots = c(herring.df$Knots,gam.df$Knots))
all.df$Sill = c(herring.df$Sill,gam.df$Sill)
all.df$Range = c(herring.df$Range,gam.df$Range)
all.df$Nugget = c(herring.df$Nugget,gam.df$Nugget)
all.df$SillNorm = c(herring.df$SillNorm,gam.df$SillNorm)
all.df$RangeNorm = c(herring.df$RangeNorm,gam.df$RangeNorm)
all.df$NuggetNorm = c(herring.df$NuggetNorm,gam.df$NuggetNorm)
all.df$Model = c(rep("Bayes+OK",nrow(herring.df)), rep("GAM+OK",nrow(gam.df)))

# Save results
write.table(all.df, "Example/Herring/results.txt")

# Plot herring results
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
  geom_bar(stat="summary") + xlab("") + ylab("Value") +
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