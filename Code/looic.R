# Code to run and save LOOIC for all models

# Created: May 19, 2021
# Last modified: May 19, 2021

# Set working directory
setwd("/Users/Elizabeth/OneDrive/School/Research/Splines")

# Load packages
library(cmdstanr)
library(loo)

# Contents (ctrl-f):
#	0. Common functions
#	I. Load all model output
#	II. LOO for two-dimensional data sets
#	III. LOO for three-dimensional data sets
#	IV. Create LOOIC table
#	V. Run comparison
#	VI. Save knots of best models


########## 0. Common functions ##########

# Function to retrieve the LOOic and DIC from any Stan object
#	obj: Stan model object
#	dat: relevant simulated data
#	type: "2d" or "3d"
#	returns LOOic from a single model
getic = function(obj, dat, type) {
	lambda = exp(obj$ostan$draws("log_mu"))
	y = if(type == "2d") {dat$y} else{dat$z}
	if(type == "2d") {
		Mx = matrix(dat$x, nrow = length(dat$x), ncol = length(dat$x), byrow = T)
		My = matrix(dat$x, nrow = length(dat$x), ncol = length(dat$x), byrow = F)
		D = abs(Mx - My)
	} else if (type == "3d") {
		Mx = matrix(dat$x, nrow = length(dat$x), ncol = length(dat$x))
		My = matrix(dat$y, nrow = length(dat$y), ncol = length(dat$y))
		D = sqrt((Mx-t(Mx))^2 + (My-t(My))^2)
	}
	logp_y = aperm(apply(lambda, c(1,2), dpois, x = y, log = T), c(2,3,1))
	dic = -4*mean(apply(logp_y,c(1,2),sum)) - -2*sum(dpois(y, exp(obj$ostan$summary("log_mu")$mean), log = T))
	return(list(loo = loo(x = logp_y, r_eff = relative_eff(exp(logp_y))), dic = dic))
}

# Function to grab looic from loo object
#	obj: output from loo function
#	returns the looic estimate
plooic = function(obj) {
	return(obj$loo$estimates["looic","Estimate"])
}


########## I. Load all model output ##########

# Load all two-dimensional model results
load("Output/Grouped/one1stan.rda")
load("Output/Grouped/one2stan.rda")
load("Output/Grouped/one3stan.rda")
load("Output/Grouped/one4stan.rda")
load("Output/Grouped/one5stan.rda")
load("Output/Grouped/one6stan.rda")
load("Output/Grouped/one7stan.rda")
load("Output/Grouped/one8stan.rda")
load("Output/Grouped/one9stan.rda")
load("Output/Grouped/one10stan.rda")

# Load all three-dimensional model results
load("Output/Grouped/two1stan.rda")
load("Output/Grouped/two2stan.rda")
load("Output/Grouped/two3stan.rda")
load("Output/Grouped/two4stan.rda")
load("Output/Grouped/two5stan.rda")
load("Output/Grouped/two6stan.rda")
load("Output/Grouped/two7stan.rda")
load("Output/Grouped/two8stan.rda")
load("Output/Grouped/two9stan.rda")
load("Output/Grouped/two10stan.rda")


########## II. LOO for two-dimensional data sets ##########

one1.ic = lapply(one1.stan, getic, dat = one1, type = "2d")
one2.ic = lapply(one2.stan, getic, dat = one2, type = "2d")
one3.ic = lapply(one3.stan, getic, dat = one3, type = "2d")
one4.ic = lapply(one4.stan, getic, dat = one4, type = "2d")
one5.ic = lapply(one5.stan, getic, dat = one5, type = "2d")
one6.ic = lapply(one6.stan, getic, dat = one6, type = "2d")
one7.ic = lapply(one7.stan, getic, dat = one7, type = "2d")
one8.ic = lapply(one8.stan, getic, dat = one8, type = "2d")
one9.ic = lapply(one9.stan, getic, dat = one9, type = "2d")
one10.ic = lapply(one10.stan, getic, dat = one10, type = "2d")
one.ic = list(one1.ic, one2.ic, one3.ic, one4.ic, one5.ic, one6.ic, one7.ic, one8.ic, one9.ic, one10.ic)

save(object = one.ic, file = "Output/IC/oneic.rda")


########## III. LOO for three-dimensional data sets ##########

two1.ic = lapply(two1.stan, getic, dat = two1, type = "3d")
two2.ic = lapply(two2.stan, getic, dat = two1, type = "3d")
two3.ic = lapply(two3.stan, getic, dat = two1, type = "3d")
two4.ic = lapply(two4.stan, getic, dat = two1, type = "3d")
two5.ic = lapply(two5.stan, getic, dat = two1, type = "3d")
two6.ic = lapply(two6.stan, getic, dat = two1, type = "3d")
two7.ic = lapply(two7.stan, getic, dat = two1, type = "3d")
two8.ic = lapply(two8.stan, getic, dat = two1, type = "3d")
two9.ic = lapply(two9.stan, getic, dat = two1, type = "3d")
two10.ic = lapply(two10.stan, getic, dat = two1, type = "3d")
two.ic = list(two1.ic, two2.ic, two3.ic, two4.ic, two5.ic, two6.ic, two7.ic, two8.ic, two9.ic, two10.ic)

save(object = two.ic, file = "Output/IC/twoic.rda")


########## IV. Create LOOIC table ##########

looic = data.frame(Set = rep(seq(10),2), Dim = c(rep(1,10),rep(2,10)))
looic$k2 = rep(NA, nrow(looic))
looic$k3 = rep(NA, nrow(looic))
looic$k4 = rep(NA, nrow(looic))
looic$k5 = rep(NA, nrow(looic))
looic$k6 = rep(NA, nrow(looic))
looic$k7 = rep(NA, nrow(looic))
looic$k8 = rep(NA, nrow(looic))
looic$k9 = rep(NA, nrow(looic))
looic$k10 = rep(NA, nrow(looic))
looic[1,3:11] = unlist(lapply(one1.ic, plooic))
looic[2,3:11] = unlist(lapply(one2.ic, plooic))
looic[3,3:11] = unlist(lapply(one3.ic, plooic))
looic[4,3:11] = unlist(lapply(one4.ic, plooic))
looic[5,3:11] = unlist(lapply(one5.ic, plooic))
looic[6,3:11] = unlist(lapply(one6.ic, plooic))
looic[7,3:11] = unlist(lapply(one7.ic, plooic))
looic[8,3:11] = unlist(lapply(one8.ic, plooic))
looic[9,3:11] = unlist(lapply(one9.ic, plooic))
looic[10,3:11] = unlist(lapply(one10.ic, plooic))
looic[11,3:11] = unlist(lapply(two1.ic, plooic))
looic[12,3:11] = unlist(lapply(two2.ic, plooic))
looic[13,3:11] = unlist(lapply(two3.ic, plooic))
looic[14,3:11] = unlist(lapply(two4.ic, plooic))
looic[15,3:11] = unlist(lapply(two5.ic, plooic))
looic[16,3:11] = unlist(lapply(two6.ic, plooic))
looic[17,3:11] = unlist(lapply(two7.ic, plooic))
looic[18,3:11] = unlist(lapply(two8.ic, plooic))
looic[19,3:11] = unlist(lapply(two9.ic, plooic))
looic[20,3:11] = unlist(lapply(two10.ic, plooic))

write.table(looic, "Output/IC/looic.txt")


########## V. Run comparison ##########

# Extract LOO objects
one1.loo = list(one.ic[[1]][[1]]$loo,one.ic[[1]][[2]]$loo,one.ic[[1]][[3]]$loo,one.ic[[1]][[4]]$loo,one.ic[[1]][[5]]$loo,one.ic[[1]][[6]]$loo,one.ic[[1]][[7]]$loo,one.ic[[1]][[8]]$loo,one.ic[[1]][[9]]$loo)
one2.loo = list(one.ic[[2]][[1]]$loo,one.ic[[2]][[2]]$loo,one.ic[[2]][[3]]$loo,one.ic[[2]][[4]]$loo,one.ic[[2]][[5]]$loo,one.ic[[2]][[6]]$loo,one.ic[[2]][[7]]$loo,one.ic[[2]][[8]]$loo,one.ic[[2]][[9]]$loo)
one3.loo = list(one.ic[[3]][[1]]$loo,one.ic[[3]][[2]]$loo,one.ic[[3]][[3]]$loo,one.ic[[3]][[4]]$loo,one.ic[[3]][[5]]$loo,one.ic[[3]][[6]]$loo,one.ic[[3]][[7]]$loo,one.ic[[3]][[8]]$loo,one.ic[[3]][[9]]$loo)
one4.loo = list(one.ic[[4]][[1]]$loo,one.ic[[4]][[2]]$loo,one.ic[[4]][[3]]$loo,one.ic[[4]][[4]]$loo,one.ic[[4]][[5]]$loo,one.ic[[4]][[6]]$loo,one.ic[[4]][[7]]$loo,one.ic[[4]][[8]]$loo,one.ic[[4]][[9]]$loo)
one5.loo = list(one.ic[[5]][[1]]$loo,one.ic[[5]][[2]]$loo,one.ic[[5]][[3]]$loo,one.ic[[5]][[4]]$loo,one.ic[[5]][[5]]$loo,one.ic[[5]][[6]]$loo,one.ic[[5]][[7]]$loo,one.ic[[5]][[8]]$loo,one.ic[[5]][[9]]$loo)
one6.loo = list(one.ic[[6]][[1]]$loo,one.ic[[6]][[2]]$loo,one.ic[[6]][[3]]$loo,one.ic[[6]][[4]]$loo,one.ic[[6]][[5]]$loo,one.ic[[6]][[6]]$loo,one.ic[[6]][[7]]$loo,one.ic[[6]][[8]]$loo,one.ic[[6]][[9]]$loo)
one7.loo = list(one.ic[[7]][[1]]$loo,one.ic[[7]][[2]]$loo,one.ic[[7]][[3]]$loo,one.ic[[7]][[4]]$loo,one.ic[[7]][[5]]$loo,one.ic[[7]][[6]]$loo,one.ic[[7]][[7]]$loo,one.ic[[7]][[8]]$loo,one.ic[[7]][[9]]$loo)
one8.loo = list(one.ic[[8]][[1]]$loo,one.ic[[8]][[2]]$loo,one.ic[[8]][[3]]$loo,one.ic[[8]][[4]]$loo,one.ic[[8]][[5]]$loo,one.ic[[8]][[6]]$loo,one.ic[[8]][[7]]$loo,one.ic[[8]][[8]]$loo,one.ic[[8]][[9]]$loo)
one9.loo = list(one.ic[[9]][[1]]$loo,one.ic[[9]][[2]]$loo,one.ic[[9]][[3]]$loo,one.ic[[9]][[4]]$loo,one.ic[[9]][[5]]$loo,one.ic[[9]][[6]]$loo,one.ic[[9]][[7]]$loo,one.ic[[9]][[8]]$loo,one.ic[[9]][[9]]$loo)
one10.loo = list(one.ic[[10]][[1]]$loo,one.ic[[10]][[2]]$loo,one.ic[[10]][[3]]$loo,one.ic[[10]][[4]]$loo,one.ic[[10]][[5]]$loo,one.ic[[10]][[6]]$loo,one.ic[[10]][[7]]$loo,one.ic[[10]][[8]]$loo,one.ic[[10]][[9]]$loo)

two1.loo = list(two.ic[[1]][[1]]$loo,two.ic[[1]][[2]]$loo,two.ic[[1]][[3]]$loo,two.ic[[1]][[4]]$loo,two.ic[[1]][[5]]$loo,two.ic[[1]][[6]]$loo,two.ic[[1]][[7]]$loo,two.ic[[1]][[8]]$loo,two.ic[[1]][[9]]$loo)
two2.loo = list(two.ic[[2]][[1]]$loo,two.ic[[2]][[2]]$loo,two.ic[[2]][[3]]$loo,two.ic[[2]][[4]]$loo,two.ic[[2]][[5]]$loo,two.ic[[2]][[6]]$loo,two.ic[[2]][[7]]$loo,two.ic[[2]][[8]]$loo,two.ic[[2]][[9]]$loo)
two3.loo = list(two.ic[[3]][[1]]$loo,two.ic[[3]][[2]]$loo,two.ic[[3]][[3]]$loo,two.ic[[3]][[4]]$loo,two.ic[[3]][[5]]$loo,two.ic[[3]][[6]]$loo,two.ic[[3]][[7]]$loo,two.ic[[3]][[8]]$loo,two.ic[[3]][[9]]$loo)
two4.loo = list(two.ic[[4]][[1]]$loo,two.ic[[4]][[2]]$loo,two.ic[[4]][[3]]$loo,two.ic[[4]][[4]]$loo,two.ic[[4]][[5]]$loo,two.ic[[4]][[6]]$loo,two.ic[[4]][[7]]$loo,two.ic[[4]][[8]]$loo,two.ic[[4]][[9]]$loo)
two5.loo = list(two.ic[[5]][[1]]$loo,two.ic[[5]][[2]]$loo,two.ic[[5]][[3]]$loo,two.ic[[5]][[4]]$loo,two.ic[[5]][[5]]$loo,two.ic[[5]][[6]]$loo,two.ic[[5]][[7]]$loo,two.ic[[5]][[8]]$loo,two.ic[[5]][[9]]$loo)
two6.loo = list(two.ic[[6]][[1]]$loo,two.ic[[6]][[2]]$loo,two.ic[[6]][[3]]$loo,two.ic[[6]][[4]]$loo,two.ic[[6]][[5]]$loo,two.ic[[6]][[6]]$loo,two.ic[[6]][[7]]$loo,two.ic[[6]][[8]]$loo,two.ic[[6]][[9]]$loo)
two7.loo = list(two.ic[[7]][[1]]$loo,two.ic[[7]][[2]]$loo,two.ic[[7]][[3]]$loo,two.ic[[7]][[4]]$loo,two.ic[[7]][[5]]$loo,two.ic[[7]][[6]]$loo,two.ic[[7]][[7]]$loo,two.ic[[7]][[8]]$loo,two.ic[[7]][[9]]$loo)
two8.loo = list(two.ic[[8]][[1]]$loo,two.ic[[8]][[2]]$loo,two.ic[[8]][[3]]$loo,two.ic[[8]][[4]]$loo,two.ic[[8]][[5]]$loo,two.ic[[8]][[6]]$loo,two.ic[[8]][[7]]$loo,two.ic[[8]][[8]]$loo,two.ic[[8]][[9]]$loo)
two9.loo = list(two.ic[[9]][[1]]$loo,two.ic[[9]][[2]]$loo,two.ic[[9]][[3]]$loo,two.ic[[9]][[4]]$loo,two.ic[[9]][[5]]$loo,two.ic[[9]][[6]]$loo,two.ic[[9]][[7]]$loo,two.ic[[9]][[8]]$loo,two.ic[[9]][[9]]$loo)
two10.loo = list(two.ic[[10]][[1]]$loo,two.ic[[10]][[2]]$loo,two.ic[[10]][[3]]$loo,two.ic[[10]][[4]]$loo,two.ic[[10]][[5]]$loo,two.ic[[10]][[6]]$loo,two.ic[[10]][[7]]$loo,two.ic[[10]][[8]]$loo,two.ic[[10]][[9]]$loo)

# Compare objects
lc_one1 = loo_compare(one1.loo)
lc_one2 = loo_compare(one2.loo)
lc_one3 = loo_compare(one3.loo)
lc_one4 = loo_compare(one4.loo)
lc_one5 = loo_compare(one5.loo)
lc_one6 = loo_compare(one6.loo)
lc_one7 = loo_compare(one7.loo)
lc_one8 = loo_compare(one8.loo)
lc_one9 = loo_compare(one9.loo)
lc_one10 = loo_compare(one10.loo)

# Compare objects
lc_two1 = loo_compare(two1.loo)
lc_two2 = loo_compare(two2.loo)
lc_two3 = loo_compare(two3.loo)
lc_two4 = loo_compare(two4.loo)
lc_two5 = loo_compare(two5.loo)
lc_two6 = loo_compare(two6.loo)
lc_two7 = loo_compare(two7.loo)
lc_two8 = loo_compare(two8.loo)
lc_two9 = loo_compare(two9.loo)
lc_two10 = loo_compare(two10.loo)

loo_one = list(lc_one1, lc_one2, lc_one3, lc_one4, lc_one5, lc_one6, lc_one7, lc_one8, lc_one9, lc_one10)
loo_two = list(lc_two1, lc_two2, lc_two3, lc_two4, lc_two5, lc_two6, lc_two7, lc_two8, lc_two9, lc_two10)

save(object = loo_one, "Output/IC/loocompare_one.rda")
save(object = loo_two, "Output/IC/loocompare_two.rda")

# Get number of points with pareto k > 0.7
getct = function(obj) {
	return(pareto_k_table(obj)[1,1] + pareto_k_table(obj)[2,1])
}

pareto_one = matrix(NA, nrow = 10, ncol = 9)
rownames(pareto_one) = c("one1", "one2", "one3", "one4", "one5", "one6", "one7", "one8", "one9", "one10")
colnames(pareto_one) = seq(2,10)
pareto_one[1,] = unlist(lapply(one1.loo, getct))
pareto_one[2,] = unlist(lapply(one2.loo, getct))
pareto_one[3,] = unlist(lapply(one3.loo, getct))
pareto_one[4,] = unlist(lapply(one4.loo, getct))
pareto_one[5,] = unlist(lapply(one5.loo, getct))
pareto_one[6,] = unlist(lapply(one6.loo, getct))
pareto_one[7,] = unlist(lapply(one7.loo, getct))
pareto_one[8,] = unlist(lapply(one8.loo, getct))
pareto_one[9,] = unlist(lapply(one9.loo, getct))
pareto_one[10,] = unlist(lapply(one10.loo, getct))

pareto_two = matrix(NA, nrow = 10, ncol = 9)
rownames(pareto_two) = c("two1", "two2", "two3", "two4", "two5", "two6", "two7", "two8", "two9", "two10")
colnames(pareto_two) = seq(2,10)
pareto_two[1,] = unlist(lapply(two1.loo, getct))
pareto_two[2,] = unlist(lapply(two2.loo, getct))
pareto_two[3,] = unlist(lapply(two3.loo, getct))
pareto_two[4,] = unlist(lapply(two4.loo, getct))
pareto_two[5,] = unlist(lapply(two5.loo, getct))
pareto_two[6,] = unlist(lapply(two6.loo, getct))
pareto_two[7,] = unlist(lapply(two7.loo, getct))
pareto_two[8,] = unlist(lapply(two8.loo, getct))
pareto_two[9,] = unlist(lapply(two9.loo, getct))
pareto_two[10,] = unlist(lapply(two10.loo, getct))

write.table(pareto_one, file = "Output/IC/pareto_one.txt")
write.table(pareto_two, file = "Output/IC/pareto_two.txt")


########## VI. Save knots of best models ##########

one_lb = data.frame(set = seq(10), lb = c(9,8,7,5,8,9,10,6,9,10))
two_lb = data.frame(set = seq(10), lb = c(4,2,2,2,2,3,2,2,2,2))

write.table(one_lb, "Output/IC/one_lb.txt")
write.table(two_lb, "Output/IC/two_lb.txt")
