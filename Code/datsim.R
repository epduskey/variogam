library(gstat)

# Simulate Gaussian data
# 	mu: mean of Gaussian function
# 	sigma: variance of Gaussian function
# 	mag: magnitude/multiplier of the function
# 	xy: points at which to evaluate the function
# 	xc: points at which to evaluate autocorrelation (may be scaled xy to reduce magnitude of distance values)
# 	psill: partial sill for auto function
# 	rnge: range for auto function
# 	nugget: nugget for auto function
# 	model: variogram model type for auto function
# 	errvg: variance of white noise for auto function
# 	type: "3d" by default i.e. two-variate Gaussian; "2d" for standard Gaussian
# 	returns simulated density values for each coordinate value or pair in xy
gauss = function(mu, sigma, mag, xy, xc, psill, rnge, nugget, model = "Exp", errvg, type = "3d") {
	
	# Return density
	if(type == "3d") {
		zmu = rep(NA, nrow(xy))
		for(i in 1:length(zmu)) {
			xi = t(t(xy[i,]))
			zmu[i] = mag*exp((-1/2)*(t(xi-mu)%*%solve(sigma)%*%(xi-mu)))
		}
		zres = exp(auto(xc[,1], xc[,2], psill, rnge, nugget, err = errvg))
		z = rpois(nrow(xy), zmu*zres)
		return(z)
	} else if (type == "2d") {
		zmu = mag*exp(-(xy - mu)^2/(2*sigma^2))
		zres = exp(auto(xc, rep(0, length(xc)), psill, rnge, nugget, err = errvg))
		z = rpois(length(zmu), zmu*zres)
		return(z)
	}
}

# Autocorrrelated residuals
# 	x: x-coordinate
# 	y: y-coordinate
# 	psill: partial sill of variogram
# 	rnge: range of variogram
# 	nugget: nugget of variogram
# 	model: variogram model, exponential by default
# 	err: variance of white noise added to residuals
# 	returns autocorrelated residual of length = length(x)
auto = function(x, y, psill, rnge, nugget, model = "Exp", err) {
	
	# Create variogram model
	vmod = vgm(psill = psill, range = rnge, nugget = nugget, model = model)
	
	# Calculate distance between points
	Mx = matrix(x, nrow = length(x), ncol = length(x))
	My = matrix(y, nrow = length(y), ncol = length(y))
	D = sqrt((Mx-t(Mx))^2 + (My-t(My))^2)
	
	# Create variance-covariance matrix from variogram
	vcov = matrix(variogramLine(vmod, maxdist = 2*max(D), dist_vector = c(D), covariance = T)$gamma, nrow = length(x), ncol = length(x))
	LL = t(chol(vcov))
	
	# Add autocorrelation to random normal draws
	res = (LL %*% (rnorm(length(x), 0, err)))
	return(res)
}