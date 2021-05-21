library(mvtnorm)

# A function to get estimated medians from Stan objects
#	obj: Stan model objects with output and design matrix
#	type: "2d" or "3d"
#	returns matrix of median density estimates
getmed = function(obj, type) {
	
	# Reconfigure draws object
	b = apply(obj$ostan$draws("b"), 3, cbind)
	
	# Generate posterior predictions
	out = matrix(NA, nrow = nrow(b), ncol = 100)
	for(i in 1:nrow(out)) {
		out[i,] = exp(obj$G %*% b[i,])
	}
	
	if(type == "2d") {return(apply(out, 2, median))}
	else if(type == "3d") {return(matrix(apply(out, 2, median), nrow = 10, ncol = 10, byrow = T))}
}

# Function to get covariance matrix from variogram parameters
#	params: a vector containing nugget, sill, range, mean residual, and residual, in that order
#	D: nxn distance matrix
#	returns covariance matrix
getdmvn = function(params, D) {
	nug = params[1]
	sill = params[2]
	rnge = params[3]
	mu_residual = params[4:(3+nrow(D))]
	residual = params[(3+nrow(D)+1):length(params)]
	dmvnorm(x = residual, mean = mu_residual, sigma = matrix(nug + sill*exp(-D/rnge), nrow = nrow(D), ncol = ncol(D)), log = T)
	#dmvnorm(x = params[104:203], mean = params[4:103], sigma = matrix(params[1] + params[2]*exp(-D/params[3]), nrow = 100, ncol = 100), log = T)
}

# Function to extract median nugget estimate from Stan model object
#	mod: Stan model object
#	returns median nugget estimate from a single model
getnugget = function(mod) {
	return(mod$ostan$summary("nugget")$median)
}

# Function to extract median sill estimate from Stan model object
#	mod: Stan model object
#	returns median sill estimate from a single model
getsill = function(mod) {
	return(mod$ostan$summary("sill")$median)
}

# Function to extract median range estimate from Stan model object
#	mod: JAGS model object
#	returns median range estimate from a single model
getrange = function(mod) {
	return(mod$ostan$summary("range")$median)
}

# Function to extract normalized maximum value from median estimates
#	medlist: a list containing vectors of all median estimates for a given model
#	returns a vector of normalized maximum density estimates
medmax = function(medlist) {
	return(unlist(lapply(medlist, max))/max(unlist(lapply(medlist, max))))
}

# Function to plot normalized sill estimates
#	modlist: a list of all models for a given data set, with knots ranging from 2 to 10
#	color: a color value for the line
#	type: norm = normalized (divided by max) OR raw = un-normalized
#	returns raw sill values if type = "raw"
sillplot = function(modlist, color, type = "norm") {
	if(type == "norm") {
		y.sill = unlist(lapply(modlist, getsill))
		sill.norm = y.sill/max(y.sill)
		lines(seq(2,10), y.sill/(max(y.sill)), lwd = 2, col = color)
		#return(y.sill)
	} else if(type == "raw") {
		y.sill = unlist(lapply(modlist, getsill))
		lines(seq(2,10), y.sill, lwd = 2, col = color)
		return(y.sill)
	}	
}

# Function to get sill posterior for boxplots
#	modlist: a list containing all models for one data set
#	returns data frame with sill and corresponding number of knots
sillbox = function(modlist) {
	sill = data.frame(sill = c(c(modlist[[1]]$ostan$draws("sill")), 
						c(modlist[[2]]$ostan$draws("sill")),
						c(modlist[[3]]$ostan$draws("sill")),
						c(modlist[[4]]$ostan$draws("sill")),
						c(modlist[[5]]$ostan$draws("sill")),
						c(modlist[[6]]$ostan$draws("sill")),
						c(modlist[[7]]$ostan$draws("sill")),
						c(modlist[[8]]$ostan$draws("sill")),
						c(modlist[[9]]$ostan$draws("sill"))))
	sill$knots = c(rep(2, modlist[[1]]$ostan$metadata()$iter_sampling*length(modlist[[1]]$ostan$metadata()$id)),
				rep(3, modlist[[2]]$ostan$metadata()$iter_sampling*length(modlist[[2]]$ostan$metadata()$id)),
				rep(4, modlist[[3]]$ostan$metadata()$iter_sampling*length(modlist[[3]]$ostan$metadata()$id)),
				rep(5, modlist[[4]]$ostan$metadata()$iter_sampling*length(modlist[[4]]$ostan$metadata()$id)),
				rep(6, modlist[[5]]$ostan$metadata()$iter_sampling*length(modlist[[5]]$ostan$metadata()$id)),
				rep(7, modlist[[6]]$ostan$metadata()$iter_sampling*length(modlist[[6]]$ostan$metadata()$id)),
				rep(8, modlist[[7]]$ostan$metadata()$iter_sampling*length(modlist[[7]]$ostan$metadata()$id)),
				rep(9, modlist[[8]]$ostan$metadata()$iter_sampling*length(modlist[[8]]$ostan$metadata()$id)),
				rep(10, modlist[[9]]$ostan$metadata()$iter_sampling*length(modlist[[9]]$ostan$metadata()$id)))
	return(sill)
}

# Function to plot normalized range estimates
#	modlist: a list of all models for a given data set, with knots ranging from 2 to 10
#	color: a color value for the line
#	type: norm = normalized (divided by max) OR raw = un-normalized
#	returns raw range values of type = "raw"
rangeplot = function(modlist, color, type = "norm") {
	if(type == "norm") {
		y.range = unlist(lapply(modlist, getrange))
		range.norm = y.range/max(y.range)
		lines(seq(2,10), y.range/(max(y.range)), lwd = 2, col = color)
		#return(y.range)
	} else if(type == "raw") {
		y.range = unlist(lapply(modlist, getrange))
		lines(seq(2,10), y.range, lwd = 2, col = color)
		return(y.range)
	}	
}

# Function to get range posterior for boxplots
#	modlist: a list containing all models for one data set
#	returns data frame with range and corresponding number of knots
rangebox = function(modlist) {
	range = data.frame(range = c(c(modlist[[1]]$ostan$draws("range")), 
						c(modlist[[2]]$ostan$draws("range")),
						c(modlist[[3]]$ostan$draws("range")),
						c(modlist[[4]]$ostan$draws("range")),
						c(modlist[[5]]$ostan$draws("range")),
						c(modlist[[6]]$ostan$draws("range")),
						c(modlist[[7]]$ostan$draws("range")),
						c(modlist[[8]]$ostan$draws("range")),
						c(modlist[[9]]$ostan$draws("range"))))
	range$knots = c(rep(2, modlist[[1]]$ostan$metadata()$iter_sampling*length(modlist[[1]]$ostan$metadata()$id)),
				rep(3, modlist[[2]]$ostan$metadata()$iter_sampling*length(modlist[[2]]$ostan$metadata()$id)),
				rep(4, modlist[[3]]$ostan$metadata()$iter_sampling*length(modlist[[3]]$ostan$metadata()$id)),
				rep(5, modlist[[4]]$ostan$metadata()$iter_sampling*length(modlist[[4]]$ostan$metadata()$id)),
				rep(6, modlist[[5]]$ostan$metadata()$iter_sampling*length(modlist[[5]]$ostan$metadata()$id)),
				rep(7, modlist[[6]]$ostan$metadata()$iter_sampling*length(modlist[[6]]$ostan$metadata()$id)),
				rep(8, modlist[[7]]$ostan$metadata()$iter_sampling*length(modlist[[7]]$ostan$metadata()$id)),
				rep(9, modlist[[8]]$ostan$metadata()$iter_sampling*length(modlist[[8]]$ostan$metadata()$id)),
				rep(10, modlist[[9]]$ostan$metadata()$iter_sampling*length(modlist[[9]]$ostan$metadata()$id)))
	return(range)
}

# Function to draw median predicted lines for two-dimensional models
#	dat: data set upon which the model was run
#	ymax: maximum value for the y-axis
#	yint: intervals for y-axis
#	modlist: a list containing all JAGS models for each knot value
#	modgam: the corresponding frequentist regression kriging model
#	type: "2d" or "3d"
#	returns a list of median density estimates for each model in modlist
medplot2 = function(dat, ymax, yint, modlist, modgam, type) {
	plot(y~xs, dat, pch = 16, ylim = c(0-ymax*0.1,ymax+ymax*0.1), cex = 1, axes = F, ylab = "")
	med = lapply(modlist, getmed, type = type)
	for(i in 1:length(med)) {
		lines(dat$xs, med[[i]], lwd = 2, col = lcol_knots[i])
	}
	lines(seq(-1.7,1.7,by=0.01), predict(modgam, newdata = list(xs = seq(-1.7,1.7,by=0.01)), type = "response"), lwd = 2, lty = 2, col = "red")
	axis(2, at = seq(0,ymax,by=yint), cex.axis = 2)
	box()
	return(med)
}

# Function to draw median predicted surfaces for three-dimensional models
#	modlist: a list containing all JAGS models for each knot value
#	modgam: the corresponding frequentist regression kriging model
#	dat: data set upon which the model was run
#	type: "2d" or "3d"
#	returns a list of median density estimates for each model in modlist
medplot3 = function(modlist, modgam, dat, type) {
	x = matrix(dat$xs, nrow=10, ncol = 10)[,1]
	y = matrix(dat$ys, nrow=10, ncol = 10)[1,]
	z = lapply(modlist, getmed, type = type)
	zgam = matrix(predict(modgam, type = "response"), nrow = 10, ncol = 10, byrow = T)
	fig = plot_ly(showscale = F)
	for(i in 1:length(z)) {
		fig = fig %>% add_surface(z = z[[i]], x = x, y = y, opacity = ovec[i], colorscale = list(c(0,1), c(lcol_plotly[i],lcol_plotly[i])))
	}
	fig = fig %>% add_surface(z = zgam, x = x, y = y, opacity = 0.3, colorscale = list(c(0,1), c("blue","red")))
	fig = fig %>% add_trace(showlegend = F, data = dat, x = dat$xs, y = dat$ys, z = dat$z, 
								mode = "markers", type = "scatter3d", marker = list(size = 3, color = "black"))
	print(fig)
	return(z)
}

# Function to plot normalized LOOIC values and their optimum
#	obj: a vector of all dic values for a given data set
#	color: color value for line and symbols
#	opt: the optimum knot choice according to LOOIC
#	adjust: placement of optimal knot symbol at bottom of graph
#	returns nothing, plots normalized LOOIC as line, adds point to optimal model
icplot = function (obj, color, opt, adjust) {
	ic = function(obj) {return(obj$loo$estimates["looic","Estimate"])}
	looic = unlist(lapply(obj, ic))
	lines(seq(2,10), looic/max(looic), lwd = 2, col = color)
	points(opt, looic[opt-1]/max(looic), pch = 9, cex = 1.2, lwd = 2, col = color)
	points(opt, adjust-0.4, pch = 9, cex = 1.2, lwd = 2, col = color)
}

# Function to calculate design matrix from data grid:
#	dat: data source that we used to run the model
#	pdat: data over which we are predicting
#	type: "2d" or "3d"
#	K: number of knots
#	returns design matrix
gcalc = function(dat, pdat, type, K) {
	
	if(type == "2d") {
		
		# Get knot values for depth
		knots = seq(min(dat$xs) + 0.01*diff(range(dat$xs)), max(dat$xs) - 0.01*diff(range(dat$xs)), length.out = K)
		
		# Set up the matrix Z_K where Z_K[i,j] = ((pdat$xs[i]-knots[j])_+)^3
		ZK = matrix(NA, nrow = length(pdat$xs), ncol = length(knots))
		for(i in 1:nrow(ZK)) {
			for(j in 1:ncol(ZK)) {
				ZK[i,j] = ifelse(pdat$xs[i] - knots[j] < 0, 0, (pdat$xs[i]-knots[j])^3) 
			}
		}
		
		# Calculate penalty matrix
		omega = pcube(knots, maxdat = max(pdat$xs))
		svd.omega = svd(omega)
		sqrt.omega = t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
		Z = t(solve(sqrt.omega, t(ZK)))
		X = matrix(c(rep(1, length(pdat$xs)), pdat$xs), nrow = length(pdat$xs), ncol = 2)
		G = cbind(X,Z)
	}
	
	if(type == "3d") {
		
		# Set up the knots (big K) for the coordinates, knots will concentrate in spatial data clusters
		K.df = data.frame(x = dat$xs, y = dat$ys)
		K.clust = hclust(dist(K.df))
		K.tree = cutree(K.clust, k = K)
		knots = matrix(NA, nrow = K, ncol = 2)
		for(i in 1:K) {
			knots[i, ] = apply(K.df[K.tree == i, ], 2, mean)
		} 
	
		# Set the matrix Z_K where Z_K[i,j] = r^2*log(r); r = sqrt((kx[i]-pdat$xscale[j])^2 + (ky[i]-pdat$yscale[j])^2)
		ZK = matrix(NA,  nrow = length(pdat$xs), ncol = nrow(knots))
		for(i in 1:nrow(ZK)) {
			for(j in 1:ncol(ZK)) {
				r = sqrt((pdat$xs[i]-knots[j,1])^2 + (pdat$ys[i]-knots[j,2])^2)
				ZK[i,j] = (r^2)*log(r)
			}
		}
	
		# Calculate the penalty matrix
		omega = pplate(knots)
		svd.omega = svd(omega)
		sqrt.omega = t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
		Z = t(solve(sqrt.omega, t(ZK)))
		X = matrix(c(pdat$xs, pdat$ys), nrow = length(pdat$xs), ncol = 2)
		G = cbind(X,Z)
	}
	
	return(G)
}

# Function to get expected response
#	G: design matrix
#	vec: vector of estimated coefficients
#	returns response vector
mmult = function(G, vec, offset) {
	return(exp((G %*% vec) + log(offset)))
}

# Function to calculate lagrange multiplier
#	D.tau: precision matrix
#	D.cov: covariance matrix for samples
#	D.eps: covariance of sample points with predictive grid points
#	nr: number of predictive grid points
#	returns lagrange multiplier
lagrange = function(D.tau, D.cov, D.eps, ns) {
	
	# Vector of ones
	ones = matrix(1, nrow = ns, ncol = 1)
	
	# Return lagrange multiplier
	return((t(ones) %*% D.tau %*% D.eps - 1)/(t(ones) %*% D.tau %*% ones))
	
}

# Function to calculate kriging weights
#	D.tau: data precision matrix
#	D.cov: data covariance matrix
#	nr: number of sample data points
# 	vars: covariance of predictive grid locations, and lagrange multiplier
#	returns kriging weights
weight = function(D.tau, D.cov, ns, vars) {
	
	# Get kriging weight and covariance
	D.eps = head(vars, -1)
	lambda = tail(vars, 1)
	
	# Calculate and return weights
	return(D.tau %*% (D.eps - lambda))
	
}

# Function to calculate data covariance matrices
#	D: data distance matrix
#	vars: sill and range
#	returns covariance matrix AND covariance of predicted with observed points
dcalc = function(D, vars) {
	
	# Get sill and range
	psill = vars[1]
	rnge = vars[2]
	
	# Return covariance and precision matrix
	D.cov = psill*exp(-D/rnge)
	D.tau = solve(D.tau)
	return(list(D.cov = D.cov, D.tau = D.tau))
		
}

# Function to krige residuals
#	D: distance vector within predictive grid
#	deps: distance vector between grid and samples
#	ns: number of sample points
#	np: number of predictive grid points
#	var: variogram parameters (psill, range)
#	returns kriged values
krg = function(D, deps, ns, np, vars) {
	
	# Get variogram values - calculate covariance matrix
	#nugget = vars[1]
	psill = vars[1]
	rnge = vars[2]
	resid = tail(vars, -2)
	
	# Calculate covariance matrix
	D.cov = psill*exp(-D/rnge) #+ nugget
	D.tau = solve(D.cov)
	
	# Calculate covariance of all observed values with all predicted values
	D.eps = psill*exp(-deps/rnge) #+ nugget
	
	# Get kriging weights
	lambda = apply(D.eps, 2, lagrange, D.tau = D.tau, D.cov = D.cov, ns = ns)
	w = apply(rbind(D.eps, lambda), 2, weight, D.tau = D.tau, D.cov = D.cov, ns = ns)
	
	# Krige and return
	k = (t(w) %*% resid)

	return(k)	
}

# Function to calculate distance of vector from point x
#	x: x coordinates of samples
#	y: y coordinates of samples
#	xy: coordinate pair of predictive grid point
#	returns distance vector between (x,y) and xy
dv = function(x, y, xy) {
	
	# Get coordinates
	px = xy[1]
	py = xy[2]
	
	# Return distance
	return(sqrt((x - px)^2 + (y - py)^2))
	
}

# Function to build posterior of predicted values from posterior distribution:
#	G: design matrix calculated by gcalc
#	mod: model object from which we are predicting
#	Gmod: design matrix arising from original model
#	grid: grid data over which we are predicting
#	dat: data used to run the original model
#	pdat: data over which we are predicting
#	type: coordinate data for probability of presence
#	returns a list with predicted median and kriged residuals
pred = function(G, mod, Gmod, dat, pdat, type) {
	
	# Get distance matrices
	if(type == "2d") {
		Mx = matrix(dat$xs, nrow = length(dat$xs), ncol = length(dat$xs), byrow = T)
		My = matrix(dat$xs, nrow = length(dat$xs), ncol = length(dat$xs), byrow = F)
		D = abs(Mx - My)
		deps = apply(matrix(c(pdat$xs, rep(0,length(pdat$xs))), ncol = 2), 1, dv, x = dat$xs, y = 0)	
	}
	if(type == "3d") {
		Mx = matrix(dat$xs, nrow = length(dat$xs), ncol = length(dat$xs))
		My = matrix(dat$ys, nrow = length(dat$ys), ncol = length(dat$ys))
		D = sqrt((Mx-t(Mx))^2 + (My-t(My))^2)
		deps = apply(matrix(c(pdat$xs, pdat$ys), ncol = 2), 1, dv, x = dat$xs, y = dat$ys)	
	}
	
	# Calculate mean response and residuals
	mumat = apply(apply(mod$draws("b"), 3, cbind), 1, mmult, G = G, offset = 1)
	resid = sweep(-apply(apply(mod$draws("b"), 3, cbind), 1, mmult, G = Gmod, offset = 1), 1, dat$y, FUN = "+")
	
	# Krige residuals
	krmat = apply(rbind(c(mod$draws("sill")), c(mod$draws("range")), resid), 2, krg, D = D, deps = deps, ns = nrow(dat), np = nrow(pdat))
	
	return(list(mu = mumat, kr = krmat))
}

# Function to plot heat maps of three-dimensional model output
#	modlist: list of stan model objects
#	dat: data upon which models were run
#	K: number of knots
#	stan: is the i^th object in modlist a stan model object; if false, assumes a GAM object
#	zmax: maximum z value for plots
#	returns nothing, plots heat maps from modlist
simcd = function(modlist, dat, K, stan, zmax) {
	for(i in 1:10) {
		if(stan[i]) {
			pc = dat
			pc$med = c(t(getmed(modlist[[i]],"3d")))
			dat.pc = xyz2img(pc, zcol = 6, xcol = 3, ycol = 4, tolerance = 1e-9)
			image(dat.pc, breaks = seq(0,zmax,length.out=1001), col = hcl.colors(1000,"Grays",rev=T), axes = F)
			box(lty=2)
		}
		else {
			pc = dat
			pc$mean = predict(modlist[[i]], newdata = list(xs = pc$xs, ys = pc$ys), type = "response")
			dat.pc = xyz2img(pc, zcol = 6, xcol = 3, ycol = 4, tolerance = 1e-9)
			image(dat.pc, breaks = seq(0,zmax,length.out=1001), col = hcl.colors(1000,"Grays",rev=T), axes = F)
			box(lty=2)
		}
	}
}
