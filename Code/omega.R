library(Ryacas0)
x = Sym("x")
f1 = function(x) {return( x^2 )}
f2 = function(x) {return( x^3 )}
fu = function(x,kappa) {return( (x - kappa)^3 )}
f = c(f1, f2)

# Function to compute penalty matrix for cubic spline
#	knots: location of knots
#	maxdat: maximum of data
#	type: "K" to return just the non-parametric penalty, "all" to return parametric as well
#	returns penalty matrix
pcube = function(knots, maxdat, type = "K") {
	
	# Set x as a symbolic variable
	x = Sym("x")
	
	# Lower and upper integration range
	low = c(rep(min(knots), 2), knots)
	up = maxdat
	
	# Results matrix
	pmat = matrix(0, nrow = 2 + length(knots), ncol = 2 + length(knots))
	
	for(i in 3:nrow(pmat)) {
		for(j in 3:ncol(pmat)) {
			
			# Set the first function
			g1 = deriv(deriv(fu(x, knots[i-2]), x), x)
			
			# Set the second function
			g2 = deriv(deriv(fu(x, knots[j-2]), x), x)
			
			symbo = as.expression(Integrate(g1*g2, x))
			x = max(low[i], low[j])
			intlow = eval(symbo)
			x = up
			intup = eval(symbo)
			pmat[i,j] = intup - intlow
			x = Sym("x")
		}
	}
	if(type == "K") {return(pmat[3:dim(pmat)[1], 3:dim(pmat)[2]])}
	if(type == "all") {return(pmat)}
}

# Function to compute penalty matric for thin plate spline
#	knots: location of knots
#	returns penalty matrix
pplate = function(knots) {
	
	# Set up the matrix
	pmat = matrix(NA, nrow = nrow(knots), ncol = nrow(knots))
	
	# Calculate pairwise distances between knots
	for(i in 1:nrow(pmat)) {
		for(j in 1:ncol(pmat)) {
			pmat[i,j] = sqrt((knots[i,1]-knots[j,1])^2 + (knots[i,2]-knots[j,2])^2)
		}
	}
	
	# Return the penalty (note this matrix is only proportional to the penalty matrix)
	return(pmat)
}