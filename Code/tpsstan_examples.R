# install.packages("cmdstanr")
library(cmdstanr)
library(loo)

# Thin plate spline Bayesian estimation
# 	x: x-coordinate
# 	y: y-coordinate
# 	z: z-coordinate (i.e. surface height e.g. density)
# 	K: number of knots
#	rorsd: spline coefficient variance hyperparameter rate OR spline coefficient SD
# 	nchains: number of chains in Bayesian model run
# 	nburn: number of iterations to discard from each chain
# 	niter: number of iterations in each chain
# 	par: default true to run chains in parallel
#	returns stan output and design matrix
tps = function(x, y, z, K, rorsd, nchains, nburn, niter, inc, plot.k = F, par = 3) {
	
	# Calculate distances between points
	Mx = matrix(x, nrow = length(x), ncol = length(x))
	My = matrix(y, nrow = length(y), ncol = length(y))
	D = sqrt((Mx-t(Mx))^2 + (My-t(My))^2)
	
	# Set up the knots (big K) for the coordinates, knots will be clustered among observations in the grid
	K.df = data.frame(x = x, y = y)
	K.clust = hclust(dist(K.df))
	K.tree = cutree(K.clust, k = K)
	knots = matrix(NA, nrow = K, ncol = 2)
	for(i in 1:K) {
		knots[i, ] = apply(K.df[K.tree == i, ], 2, mean) + rnorm(2, sd = 0.01)
	} 
	if(plot.k) {plot(x, y, pch = 16); points(knots[,1], knots[,2], pch = 16, col = "red")}

	# Set the matrix Z_K where Z_K[i,j] = r^2*log(r); r = sqrt((kx[i]-x[j])^2 + (ky[i]-y[j])^2)
	ZK = matrix(NA,  nrow = length(x), ncol = nrow(knots))
	for(i in 1:nrow(ZK)) {
		for(j in 1:ncol(ZK)) {
			r = sqrt((x[i]-knots[j,1])^2 + (y[i]-knots[j,2])^2)
			ZK[i,j] = (r^2)*log(r)
		}
	}

	# Calculate the penalty matrix
	omega = pplate(knots)
	svd.omega = svd(omega)
	sqrt.omega = t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
	Z = t(solve(sqrt.omega, t(ZK)))
	X = matrix(c(rep(1, length(x)), x, y), nrow = length(x), ncol = 3)
	G = cbind(X,Z)
	
	# Write and save the Stan function
	cat("
	data
	{
		// Density observations
		int<lower=1> n;	
		vector<lower=0>[n] y;
		
		// Knot locations
		int<lower=1> k;
		
		// Design matrix
		array[n] vector[k+3] G;
		
		// Distance matrix
		matrix<lower=0>[n,n] D;
		
		// Spline coefficient variance rate
		real<lower=0> rorsd;
	}
	parameters
	{
		// Spline coefficients and variance
		vector[k+3] b;
		real<lower=0> btau;
		
		// Lognormal global variance
		real<lower=0> sigma;
		
		// Variogram parameters
		real<lower=0> nugget;
		real<lower=0> range;	
		real<lower=0> sill;
		
		// Local residual means and global variance
		vector[n] mu_residual;
	}
	transformed parameters
	{
		// Expected mean response and residual
		vector[n] mu;
		vector[n] residual;
		for(i in 1:n) {
			mu[i] = dot_product(b,G[i]);
			residual[i] = y[i] - exp(mu[i]);
		}
	}
	model
	{
		// Priors
		
		// Spline model coefficient and variance priors
		b[1:3] ~ normal(0,1);
		b[4:k+3] ~ normal(0,btau);
		btau ~ normal(0,rorsd)T[0,];
		
		// Variogram priors
		nugget ~ normal(0,1) T[0,];
		range ~ normal(0,1) T[0,];
		sill ~ normal(0,1000) T[0,];
		
		// Residual mean and variance priors
		mu_residual ~ normal(0,1);
		
		// Global lognormal variance prior
		sigma ~ normal(0,1)T[0,];
		
		// Likelihood
		
		// Mean
		y ~ lognormal(mu,sigma);
		
		// Residuals
		residual ~ multi_normal(mu_residual, sill*exp(-D/range)+nugget);
	}
		", file = "Code/tps_acoustic.stan")
	tpsstan = cmdstan_model("Code/tps_acoustic.stan")
	
	# Organize data for Stan run
	dat = list(y = as.vector(z),
				n = length(y),		# Number of total points
				k = K,		# Number of knots
				G = G,		# Restructured design matrix
				D = D,		# Matrix of distances between points
				rorsd = rorsd		# Spline coefficient variance rate OR spline coefficient SD
				)
	
	# Stan model run
	ostan = tpsstan$sample(
		data = dat,
		chains = nchains,
		parallel_chains = par,
		refresh = inc,
		adapt_delta = 0.8,
		max_treedepth = 15,
		iter_warmup = nburn,
		iter_sampling = niter,
		init = function() {list(range = 0.2, sill = 100, nugget = 1, btau = 1, b = rep(0,K+3), sigma = 1)}
	)
	ostan$draws()
	ostan$sampler_diagnostics()
	ostan$cmdstan_diagnose()
	ostan$summary()
	
	# Return model object
	return(list(ostan = ostan, G = G))
	
}

# Function to apply "tps" to a vector containing one or more rates
#	All arguments are identical to tps, save that rate may be a vector containing one or more values for 'rorsd'
#	returns a list of lists, each of which contains model results as applied with each value of 'rorsd'
tpsapply = function(x, y, z, K, rorsd, nchains, nburn, niter, inc, plot.k = F, par = 3, pv = T) {
	
	# Apply "tps" function to each value of rorsd
	out = sapply(rorsd, 
	tps, 
	x = x, 
	y = y, 
	z = z,
	K = K, 
	nchains = nchains, 
	nburn = nburn, 
	niter = niter, 
	inc = inc)
	
	# Return list of lists
	return(out)
}

# Function to extract information criteria (LOOIC and DIC) from Stan output
# obj: Stan model output
# dat: data input to Stan model
# dcol: column number in which density observations appear
# returns LOOIC and DIC
getic = function(obj, dat, dcol) {
  log_mu = obj$ostan$draws("mu")
  log_sd = obj$ostan$draws("sigma")
  y = dat[,dcol]
  y_arr = aperm(array(y, dim = dim(log_mu)[c(3,1,2)]), c(2,3,1))
  logp_y = array(mapply(dlnorm, meanlog = log_mu, sdlog = log_sd, x = y_arr, log = T), dim = dim(log_mu))
  dic = -4*mean(apply(logp_y,c(1,2),sum)) - -2*sum(dlnorm(y, obj$ostan$summary("mu")$mean, obj$ostan$summary("sigma")$mean, log = T))
  return(list(loo = loo(x = logp_y, r_eff = relative_eff(exp(logp_y))), dic = dic))
}