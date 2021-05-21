# Code to run standard GAM+OK method

# Created: May 19, 2021
# Last modified: May 21, 2021 by EPD

# Set working directory
setwd(paste(mypath, "variogam-main", sep = ""))

# Load packages
library(mgcv)
library(gstat)
library(sp)

# Source the necessary scripts
source("Code/alldata.R")

# Contents (ctrl-f):
#	0. Common functions
#	I. Run variogram for each two-dimensional data set
#	II. Run variogram for each three-dimensional data set
#	III. Create table of variogram parameters


########## 0. Common functions ##########

# Run and plot GAM+OK method on two-dimensional data sets
# 	dat: two-dimensional input data
# 	returns list with gam and variogram fits
gamok.one = function(dat) {
	
	# Initialize plot
	par(mfrow = c(1,2))
	
	# Run and plot GAM
	gam.dat = gam(y ~ s(xs,bs="cr"), data = dat, family = poisson())
	plot(gam.dat, shade = T, shift = coef(gam.dat)[1], trans = poisson()$linkinv, main = "gam")
	points(dat$xs, dat$y, pch = 16)
	
	# Run and plot variogram
	dat$ys = rep(0, length(dat$xs))
	dat$resid = gam.dat$residuals
	coordinates(dat) = c("xs","ys")
	dat.vario = variogram(resid ~ 1, dat)
	plot(gamma ~ dist, dat.vario, main = "vario")
	dat.fit = fit.variogram.gls(resid ~ 1, dat, vgm(nugget = 0, psill = 1, range = 0.001, model = "Exp"), maxiter = 0)
	lines(variogramLine(dat.fit, maxdist = 2), lwd = 2)
	
	# Return GAM and OK object
	return(list(gam = gam.dat, vario = dat.fit))
}

# Run and plot GAM+OK method on three-dimensional data sets
# 	dat: two-dimensional input data
# 	returns list with gam and variogram fits
gamok.two = function(dat) {
	
	# Initialize plot
	par(mfrow = c(1,2))
	
	# Run and plot GAM
	gam.dat = gam(z ~ ti(xs) + ti(ys) + ti(xs,ys), data = dat, family = poisson())
	dat.pred = data.frame(xs = rep(seq(min(dat$xs),max(dat$xs),length.out = 100),100), ys = rep(seq(min(dat$ys),max(dat$ys),length.out = 100),each=100))
	dat.pred$z = exp(predict(gam.dat, newdata = dat.pred))
	dat.sp = xyz2img(dat.pred, xcol = 1, ycol = 2, zcol = 3)
	image.plot(dat.sp, breaks = seq(0,max(dat.pred$z),length.out=1001), col = hcl.colors(1000,"Grays",rev=T), axes = F, main = "gam")

	# Run and plot variogram
	dat$resid = gam.dat$residuals
	coordinates(dat) = c("xs","ys")
	dat.vario = variogram(resid ~ 1, dat)
	plot(gamma ~ dist, dat.vario, main = "vario")
	dat.fit = fit.variogram.gls(resid ~ 1, dat, vgm(nugget = 0, psill = 1, range = 0.001, model = "Exp"), maxiter = 0)
	lines(variogramLine(dat.fit, maxdist = 2), lwd = 2)
	
	# Return GAM and OK objects
	return(list(gam = gam.dat, vario = dat.fit))
}

# Get variogram parameters from OLS fit
# 	fit: variogram fit with fit.variogram.gls
# 	model: variogram model type, usually "Exp"
# 	returns vector containing sill, range, and effective range
getvario = function(fit, model) {
	sill = sum(fit$psill)
	rn = sum(fit$range)
	ern = if(model == "Exp") {sum(fit$range)} else if(model == "Gau") {sqrt(sum(fit$range))}
	return(c(sill, rn, ern))
}

########## I. Run variogram for each two-dimensional data set ##########

# Run for all data sets
one1.fit = gamok.one(one1)
one2.fit = gamok.one(one2)
one3.fit = gamok.one(one3)
one4.fit = gamok.one(one4)
one5.fit = gamok.one(one5)
one6.fit = gamok.one(one6)
one7.fit = gamok.one(one7)
one8.fit = gamok.one(one8)
one9.fit = gamok.one(one9)
one10.fit = gamok.one(one10)

# Save grouped objects
gam.one = list(one1.fit[[1]], one2.fit[[1]], one3.fit[[1]], one4.fit[[1]], one5.fit[[1]], one6.fit[[1]], one7.fit[[1]], one8.fit[[1]], one9.fit[[1]], one10.fit[[1]])
save(object = gam.one, file = "Output/GAM/gamone.rda")


########## II. Run variogram for each three-dimensional data set ##########

# Run for all data sets
two1.fit = gamok.two(two1)
two2.fit = gamok.two(two2)
two3.fit = gamok.two(two3)
two4.fit = gamok.two(two4)
two5.fit = gamok.two(two5)
two6.fit = gamok.two(two6)
two7.fit = gamok.two(two7)
two8.fit = gamok.two(two8)
two9.fit = gamok.two(two9)
two10.fit = gamok.two(two10)

# Save grouped objects
gam.two = list(two1.fit[[1]], two2.fit[[1]], two3.fit[[1]], two4.fit[[1]], two5.fit[[1]], two6.fit[[1]], two7.fit[[1]], two8.fit[[1]], two9.fit[[1]], two10.fit[[1]])
save(object = gam.two, file = "Output/GAM/gamtwo.rda")


########## III. Create table of variogram parameters ##########

# Initialize table
vario = data.frame(Dim = c(rep(1,10), rep(2,10)), Set = rep(seq(10),2))
vario$Model = rep("Exp", 20)
vario$Sill = rep(NA, nrow(vario))
vario$Range = rep(NA, nrow(vario))
vario$Erange = rep(NA, nrow(vario))
vario$Nugget = rep(0, nrow(vario))

# Fill with fitted values and save
vario[1,4:6] = getvario(one1.fit[[2]], vario[vario$Dim == 1 & vario$Set == 1,]$Model)
vario[2,4:6] = getvario(one2.fit[[2]], vario[vario$Dim == 1 & vario$Set == 2,]$Model)
vario[3,4:6] = getvario(one3.fit[[2]], vario[vario$Dim == 1 & vario$Set == 3,]$Model)
vario[4,4:6] = getvario(one4.fit[[2]], vario[vario$Dim == 1 & vario$Set == 4,]$Model)
vario[5,4:6] = getvario(one5.fit[[2]], vario[vario$Dim == 1 & vario$Set == 5,]$Model)
vario[6,4:6] = getvario(one6.fit[[2]], vario[vario$Dim == 1 & vario$Set == 6,]$Model)
vario[7,4:6] = getvario(one7.fit[[2]], vario[vario$Dim == 1 & vario$Set == 7,]$Model)
vario[8,4:6] = getvario(one8.fit[[2]], vario[vario$Dim == 1 & vario$Set == 8,]$Model)
vario[9,4:6] = getvario(one9.fit[[2]], vario[vario$Dim == 1 & vario$Set == 9,]$Model)
vario[10,4:6] = getvario(one10.fit[[2]], vario[vario$Dim == 1 & vario$Set == 10,]$Model)
vario[11,4:6] = getvario(two1.fit[[2]], vario[vario$Dim == 2 & vario$Set == 1,]$Model)
vario[12,4:6] = getvario(two2.fit[[2]], vario[vario$Dim == 2 & vario$Set == 2,]$Model)
vario[13,4:6] = getvario(two3.fit[[2]], vario[vario$Dim == 2 & vario$Set == 3,]$Model)
vario[14,4:6] = getvario(two4.fit[[2]], vario[vario$Dim == 2 & vario$Set == 4,]$Model)
vario[15,4:6] = getvario(two5.fit[[2]], vario[vario$Dim == 2 & vario$Set == 5,]$Model)
vario[16,4:6] = getvario(two6.fit[[2]], vario[vario$Dim == 2 & vario$Set == 6,]$Model)
vario[17,4:6] = getvario(two7.fit[[2]], vario[vario$Dim == 2 & vario$Set == 7,]$Model)
vario[18,4:6] = getvario(two8.fit[[2]], vario[vario$Dim == 2 & vario$Set == 8,]$Model)
vario[19,4:6] = getvario(two9.fit[[2]], vario[vario$Dim == 2 & vario$Set == 9,]$Model)
vario[20,4:6] = getvario(two10.fit[[2]], vario[vario$Dim == 2 & vario$Set == 10,]$Model)

write.table(vario, "Output/GAM/vario.txt")
