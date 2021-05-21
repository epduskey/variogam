# Code to recreate all simulated data sets in Duskey et al 2021

# Created: May 19, 2021
# Last modified: May 21, 2021

# Set working directory
setwd(paste(mypath, "variogam-main", sep = ""))

# Source the necessary scripts

# To recreate simulated data sets
source("Code/datsim.R")

# Load packages
library(lattice)

# Contents (ctrl-f):
#	I. Two-dimensional data
#	II. Three-dimensional data
#	III. No autocorrelation two-dimensional data
#	IV. No autocorrelation three-dimensional data


########## I. Two-dimensional data ##########

# Number of points
n = 100	

# Data scale	
sc = 100	

# x-coordinates	
x = seq(1,100,length.out=n)*sc		
xscale = scale(x)

# Choose variogram parameters
psill = 100
rnge_one = abs(max(xscale)-min(xscale))/10
nugget = 5

# Choose sigma values for simulations
sigma_one = seq(from = 1.0, to = 0.5, length.out = 10)

# Set seed to replicate results in Duskey et al (2021)
set.seed(3194)

# Simulate data sets
one1 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[1], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one2 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[2], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one3 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[3], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one4 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[4], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one5 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[5], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one6 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[6], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one7 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[7], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one8 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[8], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one9 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[9], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
one10 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_one[10], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))

#plot(y~xs, one1)
#plot(y~xs, one2)
#plot(y~xs, one3)
#plot(y~xs, one4)
#plot(y~xs, one5)
#plot(y~xs, one6)
#plot(y~xs, one7)
#plot(y~xs, one8)
#plot(y~xs, one9)
#plot(y~xs, one10)


########## II. Three-dimensional data ##########

# Number of unique x-values
nx = 10

# Number of unique y-values	
ny = 10

# Data scale	
sc = 100

# x-coordinates	
x = rep(seq(1,100,length.out=nx), ny)*sc
xscale = scale(x)

# y-coordinates
y = rep(seq(1,100,length.out=ny), each = nx)*sc
yscale = scale(y)

# All points
xy = matrix(c(xscale,yscale), nrow = length(xscale), ncol = 2)
xc = matrix(c(xscale,yscale), nrow = length(x), ncol = 2)

# Choose variogram parameters
psill = 100
rnge_two = sqrt((max(xscale)-min(xscale))^2+(max(yscale)-min(yscale))^2)/10
nugget = 5

# Choose mean
mu = matrix(0, nrow = 2, ncol = 1)

# Choose sigma values for simulations
sigma_two = seq(from = 0.5, to = 0.1, length.out = 10)

# Set seed to replicate results in Duskey et al (2021)
set.seed(9128)

# Simulate data sets
two1 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[1], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two2 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[2], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two3 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[3], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two4 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[4], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two5 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[5], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two6 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[6], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two7 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[7], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two8 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[8], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two9 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[9], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))
two10 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_two[10], 20, xy, xc, psill, rnge_two, nugget, errvg = 0.1))

#wireframe(z ~ x + y, two1, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two2, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two3, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two4, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two5, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two6, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two7, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two8, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two9, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
#wireframe(z ~ x + y, two10, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)


########## III. No autocorrelation two-dimensional data ##########

# Number of points
n = 100

# Data scale	
sc = 100

# x-coordinates
x = seq(1,100,length.out=n)*sc
xscale = scale(x)

# Choose variogram parameters
psill = 1e-10
rnge_one = abs(max(xscale)-min(xscale))/10
nugget = 0

# Choose sigma value for simulation
sigma_noone = 1.0

# Set seed to replicate results in Duskey et al (2021)
set.seed(3194)

# Simulate data set
no_one1 = data.frame(x = x, xs = xscale, y = gauss(mean(xscale), sigma_noone[1], 10, xscale, xscale, psill, rnge_one, nugget, errvg = 0.1, type = "2d"))
#plot(y~xscale, no_one1)


########## IV. No autocorrelation three-dimensional data ##########

# Number of unique x-values
nx = 10		

# Number of unique y-values
ny = 10	

# Data scale	
sc = 100		

# x-coordinates
x = rep(seq(1,100,length.out=nx), ny)*sc	
xscale = scale(x)

# y-coordinates
y = rep(seq(1,100,length.out=ny), each = nx)*sc		
yscale = scale(y)

# All points
xy = matrix(c(xscale,yscale), nrow = length(xscale), ncol = 2)		
xc = matrix(c(xscale,yscale), nrow = length(x), ncol = 2)

# Choose mean value
mu = matrix(0, nrow = 2, ncol = 1)

# Choose variogram parameters
psill = 1e-10
rnge = sqrt((max(xscale)-min(xscale))^2+(max(yscale)-min(yscale))^2)/10
nugget = 0

# Choose sigma value for simulation
sigma_notwo = 0.5

# Set seed to replicate results in Duskey et al (2021)
set.seed(9128)

# Simulate data sets
no_two1 = data.frame(x = x, y = y, xs = xscale, ys = yscale, z = gauss(mu = mu, diag(2)*sigma_notwo[1], 20, xy, xc, psill, rnge, nugget, errvg = 0.1))
#wireframe(z ~ x + y, no_two1, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
