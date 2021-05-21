library(cmdstanr)

# Cubic spline Bayesian estimation
# 	x: x-coordinate
# 	y: observation e.g. organism density
# 	K: number of knots
# 	offset: area offset for density measurements
# 	nchains: number of chains in model run
# 	nburn: number of warmup iterations
# 	niter: number of post-warmup iterations in each chain
# 	nthin: thinning interval (i.e. keep every nthin^th iteration)
# 	par: number of chains to run in parallel at a given time
# 	pv: estimate btau if T, set btau to 0.1 if F; default is T
#	returns stan output and design matrix
cubic = function(x, y, K, offset, nchains, nburn, niter, nthin, inc, par = 3, pv = T) {
	
	# Calculate distance between points
	Mx = matrix(x, nrow = length(x), ncol = length(x), byrow = T)
	My = matrix(x, nrow = length(x), ncol = length(x), byrow = F)
	D = abs(Mx - My)
	
	# Get knot values
	knots = seq(min(x) + 0.01*sum(abs(range(x))), max(x) - 0.01*sum(abs(range(x))), length.out = K)
	
	# Set up the matrix Z_K where Z_K[i,j] = ((x[i]-knots[j])_+)^3
	ZK = matrix(NA, nrow = length(x), ncol = length(knots))
	for(i in 1:length(x)) {
		for(j in 1:length(knots)) {
			ZK[i,j] = ifelse(x[i] - knots[j] < 0, 0, (x[i]-knots[j])^3) 
		}
	}
	
	# Calculate penalty matrix
	omega = pcube(knots, maxdat = max(x))
	svd.omega = svd(omega)
	sqrt.omega = t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
	Z = t(solve(sqrt.omega, t(ZK)))
	X = matrix(c(rep(1, length(x)), x), nrow = length(x), ncol = 2)
	G = cbind(X,Z)
	
	# Write and save the Stan function: estimate btau if pv, set btau=0.1 otherwise
	if(pv) {
		cat("
		data
		{
			// Count observations and area offset
			int<lower=1> n;	
			int<lower=0> y[n];
			vector<lower=0>[n] offset;
			
			// Knots
			int<lower=1> k;
			
			// Design matrix
			vector[k+2] G[n];
			
			// Distance matrix
			matrix<lower=0>[n,n] D;
		}
		parameters
		{
			// Spline coefficients and variance
			vector[k+2] b;
			real<lower=0> btau;
			
			// Variogram parameters
			real<lower=0> nugget;
			real<lower=0> range;	
			real<lower=0> sill;
		
			// Residuals amd local residual means
			vector[n] mu_residual;
		}
		transformed parameters
		{	
			// Expected mean response and residual
			vector[n] log_mu;
			vector[n] residual;
			for(i in 1:n) {
				log_mu[i] = dot_product(b,G[i]) + offset[i];
				residual[i] = y[i] - exp(log_mu[i]);
			}
		}
		model
		{
			// Priors
			
			// Spline model coefficient and variance priors
			b[1:2] ~ normal(0,1);
			b[3:k+2] ~ normal(0,btau);
			btau ~ gamma(1,0.5);
		
			// Variogram priors
			nugget ~ normal(0,1) T[0,];
			range ~ normal(0,1) T[0,];
			sill ~ normal(0,1000) T[0,];
			
			// Residual mean and variance priors
			mu_residual ~ normal(0,1);

			// Likelihood
			
			// Count observations
			y ~ poisson(exp(log_mu));
			
			// Residuals
			residual ~ multi_normal(mu_residual,sill*exp(-D/range)+nugget);
		}
		", file = "Code/cube.stan")
		cubestan = cmdstan_model("Code/cube.stan")
	} else {
		cat("
		data
		{
			// Count observations and area offset
			int<lower=1> n;	
			int<lower=0> y[n];
			vector<lower=0>[n] offset;
			
			// Knots
			int<lower=1> k;
			
			// Design matrix
			vector[k+2] G[n];
			
			// Distance matrix
			matrix<lower=0>[n,n] D;
		}
		parameters
		{
			// Spline coefficients and variance
			vector[k+2] b;
			
			// Variogram parameters
			real<lower=0> nugget;
			real<lower=0> range;	
			real<lower=0> sill;
			
			// Residuals amd local residual means and variance
			vector[n] mu_residual;
		}
		transformed parameters
		{	
			// Expected mean response and residual
			vector[n] log_mu;
			vector[n] residual;
			for(i in 1:n) {
				log_mu[i] = dot_product(b,G[i]) + offset[i];
				residual[i] = y[i] - exp(log_mu[i]);
			}
		}
		model
		{
			// Priors
			
			// Spline model coefficient and variance priors
			b[1:2] ~ normal(0,1);
			b[3:k+2] ~ normal(0,0.1);
			
			// Variogram priors
			nugget ~ normal(0,1) T[0,];
			range ~ normal(0,1) T[0,];
			sill ~ normal(0,1000) T[0,];
			
			// Residual mean and variance priors
			mu_residual ~ normal(0,1);
			
			// Residuals
			residual ~ multi_normal(mu_residual,sill*exp(-D/range)+nugget);
			
			// Likelihood
			
			// Count observations
			y ~ poisson(exp(log_mu));
		}
		", file = "Code/cube.stan")
		cubestan = cmdstan_model("Code/cube.stan")
	}
	
	# Organize data for Stan run
	dat = list(y = as.vector(y),
				n = length(y),		# Number of total points
				k = K,		# Number of knots
				G = G,		# Restructured design matrix
				D = D,		# Matrix of distances between points
				offset = offset
				)
	
	# Stan model run
	ostan = cubestan$sample(
		data = dat,
		chains = nchains,
		parallel_chains = par,
		refresh = inc,
		adapt_delta = 0.8,
		max_treedepth = 15,
		iter_warmup = nburn,
		iter_sampling = niter,
		init = function() {list(range = 0.2, sill = 100, nugget = 1, btau = 1)}
	)
	ostan$draws()
	ostan$sampler_diagnostics()
	ostan$cmdstan_diagnose()
	ostan$summary()
	
	# Return model object
	return(list(ostan = ostan, G = G))
}