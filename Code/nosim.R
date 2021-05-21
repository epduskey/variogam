# Cubic and thin plate spline estimation for data sets sans autocorrelation

# Created: May 19, 2021
# Last modified: May 21, 2021

# Set working directory
setwd(paste(mypath, "variogam", sep = ""))

# Load packages
library(mgcv)
library(gstat)
library(sp)
library(lattice)

# Source the necessary scripts

# For simulating data from one- or two-dimensional Gaussian functions and adding correlated residuals
source("Code/datsim.R")

# To generate cubic and thin plate spline penalty matrices
source("Code/omega.R")	

# For cubic spline estimation	
source("Code/cubestan.R")

# For thin plate spline estimation
source("Code/tpsstan.R")

# Contents (ctrl-f):		
#	I. Simulate a one-dimensional data sets
#	II. Cubic spline estimation for increasing knot size
#	III. Simulate a two-dimensional data sets
#	IV. Thin plate spline estimation
#	V. Run GAM+OK on both data sets


########## I. Simulate a one-dimensional data sets ##########

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
plot(y~xscale, no_one1)


########## II. Cubic spline estimation ##########

# Choose range for number of knots
Ks = seq(2, 10)

# Apply Bayesian estimation for each
no_one1.stan2 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[1], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan2, file = "Output/No Autocorrelation/no_one1stan2.rda")
no_one1.stan3 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[2], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan3, file = "Output/No Autocorrelation/no_one1stan3.rda")
no_one1.stan4 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[3], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan4, file = "Output/No Autocorrelation/no_one1stan4.rda")
no_one1.stan5 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[4], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan5, file = "Output/No Autocorrelation/no_one1stan5.rda")
no_one1.stan6 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[5], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan6, file = "Output/No Autocorrelation/no_one1stan6.rda")
no_one1.stan7 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[6], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan7, file = "Output/No Autocorrelation/no_one1stan7.rda")
no_one1.stan8 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[7], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan8, file = "Output/No Autocorrelation/no_one1stan8.rda")
no_one1.stan9 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[8], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan9, file = "Output/No Autocorrelation/no_one1stan9.rda")
no_one1.stan10 = cubic(x = no_one1$xs, y = no_one1$y, K = Ks[9], offset = rep(0,n), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_one1.stan10, file = "Output/No Autocorrelation/no_one1stan10.rda")

# Save grouped object
no_one1.stan = list(no_one1.stan2, no_one1.stan3, no_one1.stan4, no_one1.stan5, no_one1.stan6, no_one1.stan7, no_one1.stan8, no_one1.stan9, no_one1.stan10)
save(object = no_one1.stan, file = "Output/No Autocorrelation/no_one1stan.rda")


########## III. Simulate ten two-dimensional data sets ##########

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
wireframe(z ~ x + y, no_two1, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)


########## IV. Thin plate spline estimation ##########

# Choose range for number of knots
Ks = seq(2, 10)

no_two1.stan2 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[1], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan2, file = "Output/No Autocorrelation/no_two1stan2.rda")
no_two1.stan3 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[2], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan3, file = "Output/No Autocorrelation/no_two1stan3.rda")
no_two1.stan4 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[3], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan4, file = "Output/No Autocorrelation/no_two1stan4.rda")
no_two1.stan5 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[4], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan5, file = "Output/No Autocorrelation/no_two1stan5.rda")
no_two1.stan6 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[5], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan6, file = "Output/No Autocorrelation/no_two1stan6.rda")
no_two1.stan7 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[6], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan7, file = "Output/No Autocorrelation/no_two1stan7.rda")
no_two1.stan8 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[7], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan8, file = "Output/No Autocorrelation/no_two1stan8.rda")
no_two1.stan9 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[8], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan9, file = "Output/No Autocorrelation/no_two1stan9.rda")
no_two1.stan10 = tps(x = no_two1$xs, y = no_two1$ys, z = no_two1$z, K = Ks[9], offset = rep(0, nx*ny), nchains = 6, nburn = 1000, niter = 5000, inc = 1000)
save(object = no_two1.stan10, file = "Output/No Autocorrelation/no_two1stan10.rda")

no_two1.stan = list(no_two1.stan2, no_two1.stan3, no_two1.stan4, no_two1.stan5, no_two1.stan6, no_two1.stan7, no_two1.stan8, no_two1.stan9, no_two1.stan10)
save(object = no_two1.stan, file = "Output/No Autocorrelation/no_two1stan.rda")


########## V. Run GAM+OK on both data sets ##########

# Run on 1-d
nogam.one1 = gam(y ~ s(xs,bs="cr"), data = no_one1, family = poisson())
plot(nogam.one1, shade = T, shift = coef(nogam.one1)[1], trans = poisson()$linkinv)
points(no_one1$xs, no_one1$y, pch = 16)
no_one1$ys = rep(0, length(no_one1$xs))
no_one1$resid = nogam.one1$residuals
coordinates(no_one1) = c("xs","ys")
no_one1.vario = variogram(resid ~ 1, no_one1)
plot(gamma ~ dist, no_one1.vario, xlim = c(0,1.2), ylim = c(0,1))
no_one1.fit = fit.variogram.gls(resid ~ 1, no_one1, vgm(nugget = 0.5, psill = 0.1, range = 0.001, model = "Exp"), maxiter = 0)
save(object = no_one1.fit, file = "Output/GAM/varionoone.rda")
lines(variogramLine(no_one1.fit, maxdist = 2), lwd = 2)

save(object = nogam.one1, file = "Output/GAM/gamnoone.rda")

# Run on 2-d
nogam.two1 = gam(z ~ ti(xs) + ti(ys) + ti(xs,ys), data = no_two1, family = poisson())
dat.two1 = data.frame(xs = rep(seq(min(no_two1$xs),max(no_two1$xs),length.out = 100),100), ys = rep(seq(min(no_two1$ys),max(no_two1$ys),length.out = 100),each=100))
dat.two1$z = exp(predict(nogam.two1, newdata = dat.two1))
wireframe(z ~ xs + ys, dat.two1, scales = list(x = list(arrows = T), y = list(arrows = T), z = list(arrows = F)), drape = T)
no_two1$resid = nogam.two1$residuals
coordinates(no_two1) = c("xs","ys")
two1.vario = variogram(resid ~ 1, no_two1)
plot(gamma ~ dist, two1.vario, xlim = c(0,1.5), ylim = c(0,1))
no_two1.fit = fit.variogram.gls(resid ~ 1, no_two1, vgm(nugget = 0.5, psill = 0.1, range = 0.2, model = "Exp"), maxiter = 0)
save(object = no_two1.fit, file = "Output/GAM/varionotwo.rda")
lines(variogramLine(no_two1.fit, maxdist = 2), lwd = 2)

save(object = nogam.two1, file = "Output/GAM/gamnotwo.rda")

novario = data.frame(Dim = c(1,2), Set = c(1,1), Model = c("Exp","Exp"))
novario$Sill = c(sum(no_one1.fit$psill),sum(no_two1.fit$psill))
novario$Range = c(sum(no_one1.fit$range),sum(no_two1.fit$range))
novario$Nugget = c(no_one1.fit[1,2],no_two1.fit[1,2])

write.table(novario, "Output/GAM/novario.txt")