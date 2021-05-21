# Fit simple linear models to the maximum median estimates and variogram parameters

# Created: May 19, 2021
# Last modified: May 21, 2021 by EPD

# Set working directory
setwd(paste(mypath, "variogam-main", sep = ""))

# Load packages
library(mgcv)
library(xtable)

# Source the necessary scripts
source("Code/allplot_fn.R")
source("Code/alldata.R")

# Contents (ctrl-f):
#	0. Common functions
#	I. Load all models
#	II. Construct data sets
#	III. Run models


########## 0. Common functions ##########

# Return the p-value of the slope of linear model obj
#	obj: a linear model
#	returns the p-value of the slope
getp = function(obj) {
	return(summary(obj)$coefficients[2,4])
}


########## I. Load all models ##########

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


########## II. Construct data sets ##########

# Maximum
all.one1 = data.frame(max = unlist(lapply(lapply(one1.stan, getmed, type = "2d"), max)))
all.one2 = data.frame(max = unlist(lapply(lapply(one2.stan, getmed, type = "2d"), max)))
all.one3 = data.frame(max = unlist(lapply(lapply(one3.stan, getmed, type = "2d"), max)))
all.one4 = data.frame(max = unlist(lapply(lapply(one4.stan, getmed, type = "2d"), max)))
all.one5 = data.frame(max = unlist(lapply(lapply(one5.stan, getmed, type = "2d"), max)))
all.one6 = data.frame(max = unlist(lapply(lapply(one6.stan, getmed, type = "2d"), max)))
all.one7 = data.frame(max = unlist(lapply(lapply(one7.stan, getmed, type = "2d"), max)))
all.one8 = data.frame(max = unlist(lapply(lapply(one8.stan, getmed, type = "2d"), max)))
all.one9 = data.frame(max = unlist(lapply(lapply(one9.stan, getmed, type = "2d"), max)))
all.one10 = data.frame(max = unlist(lapply(lapply(one10.stan, getmed, type = "2d"), max)))

all.two1 = data.frame(max = unlist(lapply(lapply(two1.stan, getmed, type = "3d"), max)))
all.two2 = data.frame(max = unlist(lapply(lapply(two2.stan, getmed, type = "3d"), max)))
all.two3 = data.frame(max = unlist(lapply(lapply(two3.stan, getmed, type = "3d"), max)))
all.two4 = data.frame(max = unlist(lapply(lapply(two4.stan, getmed, type = "3d"), max)))
all.two5 = data.frame(max = unlist(lapply(lapply(two5.stan, getmed, type = "3d"), max)))
all.two6 = data.frame(max = unlist(lapply(lapply(two6.stan, getmed, type = "3d"), max)))
all.two7 = data.frame(max = unlist(lapply(lapply(two7.stan, getmed, type = "3d"), max)))
all.two8 = data.frame(max = unlist(lapply(lapply(two8.stan, getmed, type = "3d"), max)))
all.two9 = data.frame(max = unlist(lapply(lapply(two9.stan, getmed, type = "3d"), max)))
all.two10 = data.frame(max = unlist(lapply(lapply(two10.stan, getmed, type = "3d"), max)))

# Sill
all.one1$sill = unlist(lapply(one1.stan, getsill))
all.one2$sill = unlist(lapply(one2.stan, getsill))
all.one3$sill = unlist(lapply(one3.stan, getsill))
all.one4$sill = unlist(lapply(one4.stan, getsill))
all.one5$sill = unlist(lapply(one5.stan, getsill))
all.one6$sill = unlist(lapply(one6.stan, getsill))
all.one7$sill = unlist(lapply(one7.stan, getsill))
all.one8$sill = unlist(lapply(one8.stan, getsill))
all.one9$sill = unlist(lapply(one9.stan, getsill))
all.one10$sill = unlist(lapply(one10.stan, getsill))

all.two1$sill = unlist(lapply(two1.stan, getsill))
all.two2$sill = unlist(lapply(two2.stan, getsill))
all.two3$sill = unlist(lapply(two3.stan, getsill))
all.two4$sill = unlist(lapply(two4.stan, getsill))
all.two5$sill = unlist(lapply(two5.stan, getsill))
all.two6$sill = unlist(lapply(two6.stan, getsill))
all.two7$sill = unlist(lapply(two7.stan, getsill))
all.two8$sill = unlist(lapply(two8.stan, getsill))
all.two9$sill = unlist(lapply(two9.stan, getsill))
all.two10$sill = unlist(lapply(two10.stan, getsill))

# Range
all.one1$range = unlist(lapply(one1.stan, getrange))
all.one2$range = unlist(lapply(one2.stan, getrange))
all.one3$range = unlist(lapply(one3.stan, getrange))
all.one4$range = unlist(lapply(one4.stan, getrange))
all.one5$range = unlist(lapply(one5.stan, getrange))
all.one6$range = unlist(lapply(one6.stan, getrange))
all.one7$range = unlist(lapply(one7.stan, getrange))
all.one8$range = unlist(lapply(one8.stan, getrange))
all.one9$range = unlist(lapply(one9.stan, getrange))
all.one10$range = unlist(lapply(one10.stan, getrange))

all.two1$range = unlist(lapply(two1.stan, getrange))
all.two2$range = unlist(lapply(two2.stan, getrange))
all.two3$range = unlist(lapply(two3.stan, getrange))
all.two4$range = unlist(lapply(two4.stan, getrange))
all.two5$range = unlist(lapply(two5.stan, getrange))
all.two6$range = unlist(lapply(two6.stan, getrange))
all.two7$range = unlist(lapply(two7.stan, getrange))
all.two8$range = unlist(lapply(two8.stan, getrange))
all.two9$range = unlist(lapply(two9.stan, getrange))
all.two10$range = unlist(lapply(two10.stan, getrange))

# Add normalized values
all.one1$maxnorm = all.one1$max/max(all.one1$max)
all.one2$maxnorm = all.one2$max/max(all.one2$max)
all.one3$maxnorm = all.one3$max/max(all.one3$max)
all.one4$maxnorm = all.one4$max/max(all.one4$max)
all.one5$maxnorm = all.one5$max/max(all.one5$max)
all.one6$maxnorm = all.one6$max/max(all.one6$max)
all.one7$maxnorm = all.one7$max/max(all.one7$max)
all.one8$maxnorm = all.one8$max/max(all.one8$max)
all.one9$maxnorm = all.one9$max/max(all.one9$max)
all.one10$maxnorm = all.one10$max/max(all.one10$max)

all.two1$maxnorm = all.two1$max/max(all.two1$max)
all.two2$maxnorm = all.two2$max/max(all.two2$max)
all.two3$maxnorm = all.two3$max/max(all.two3$max)
all.two4$maxnorm = all.two4$max/max(all.two4$max)
all.two5$maxnorm = all.two5$max/max(all.two5$max)
all.two6$maxnorm = all.two6$max/max(all.two6$max)
all.two7$maxnorm = all.two7$max/max(all.two7$max)
all.two8$maxnorm = all.two8$max/max(all.two8$max)
all.two9$maxnorm = all.two9$max/max(all.two9$max)
all.two10$maxnorm = all.two10$max/max(all.two10$max)

all.one1$sillnorm = all.one1$sill/max(all.one1$sill)
all.one2$sillnorm = all.one2$sill/max(all.one2$sill)
all.one3$sillnorm = all.one3$sill/max(all.one3$sill)
all.one4$sillnorm = all.one4$sill/max(all.one4$sill)
all.one5$sillnorm = all.one5$sill/max(all.one5$sill)
all.one6$sillnorm = all.one6$sill/max(all.one6$sill)
all.one7$sillnorm = all.one7$sill/max(all.one7$sill)
all.one8$sillnorm = all.one8$sill/max(all.one8$sill)
all.one9$sillnorm = all.one9$sill/max(all.one9$sill)
all.one10$sillnorm = all.one10$sill/max(all.one10$sill)

all.two1$sillnorm = all.two1$sill/max(all.two1$sill)
all.two2$sillnorm = all.two2$sill/max(all.two2$sill)
all.two3$sillnorm = all.two3$sill/max(all.two3$sill)
all.two4$sillnorm = all.two4$sill/max(all.two4$sill)
all.two5$sillnorm = all.two5$sill/max(all.two5$sill)
all.two6$sillnorm = all.two6$sill/max(all.two6$sill)
all.two7$sillnorm = all.two7$sill/max(all.two7$sill)
all.two8$sillnorm = all.two8$sill/max(all.two8$sill)
all.two9$sillnorm = all.two9$sill/max(all.two9$sill)
all.two10$sillnorm = all.two10$sill/max(all.two10$sill)

all.one1$rangenorm = all.one1$range/max(all.one1$range)
all.one2$rangenorm = all.one2$range/max(all.one2$range)
all.one3$rangenorm = all.one3$range/max(all.one3$range)
all.one4$rangenorm = all.one4$range/max(all.one4$range)
all.one5$rangenorm = all.one5$range/max(all.one5$range)
all.one6$rangenorm = all.one6$range/max(all.one6$range)
all.one7$rangenorm = all.one7$range/max(all.one7$range)
all.one8$rangenorm = all.one8$range/max(all.one8$range)
all.one9$rangenorm = all.one9$range/max(all.one9$range)
all.one10$rangenorm = all.one10$range/max(all.one10$range)

all.two1$rangenorm = all.two1$range/max(all.two1$range)
all.two2$rangenorm = all.two2$range/max(all.two2$range)
all.two3$rangenorm = all.two3$range/max(all.two3$range)
all.two4$rangenorm = all.two4$range/max(all.two4$range)
all.two5$rangenorm = all.two5$range/max(all.two5$range)
all.two6$rangenorm = all.two6$range/max(all.two6$range)
all.two7$rangenorm = all.two7$range/max(all.two7$range)
all.two8$rangenorm = all.two8$range/max(all.two8$range)
all.two9$rangenorm = all.two9$range/max(all.two9$range)
all.two10$rangenorm = all.two10$range/max(all.two10$range)

# All data
all.one = rbind(all.one1, all.one2, all.one3, all.one4, all.one5, all.one6, all.one7, all.one8, all.one9, all.one10)
all.one$sigma = rep(sigma_one, each = 9)
all.two = rbind(all.two1, all.two2, all.two3, all.two4, all.two5, all.two6, all.two7, all.two8, all.two9, all.two10)
all.two$sigma = rep(sigma_two, each = 9)

# Save data
write.table(all.one1, "Output/LM/allone1.txt")
write.table(all.one2, "Output/LM/allone2.txt")
write.table(all.one3, "Output/LM/allone3.txt")
write.table(all.one4, "Output/LM/allone4.txt")
write.table(all.one5, "Output/LM/allone5.txt")
write.table(all.one6, "Output/LM/allone6.txt")
write.table(all.one7, "Output/LM/allone7.txt")
write.table(all.one8, "Output/LM/allone8.txt")
write.table(all.one9, "Output/LM/allone9.txt")
write.table(all.one10, "Output/LM/allone10.txt")

write.table(all.two1, "Output/LM/alltwo1.txt")
write.table(all.two2, "Output/LM/alltwo2.txt")
write.table(all.two3, "Output/LM/alltwo3.txt")
write.table(all.two4, "Output/LM/alltwo4.txt")
write.table(all.two5, "Output/LM/alltwo5.txt")
write.table(all.two6, "Output/LM/alltwo6.txt")
write.table(all.two7, "Output/LM/alltwo7.txt")
write.table(all.two8, "Output/LM/alltwo8.txt")
write.table(all.two9, "Output/LM/alltwo9.txt")
write.table(all.two10, "Output/LM/alltwo10.txt")

write.table(all.one, "Output/LM/allone.txt")
write.table(all.two, "Output/LM/alltwo.txt")

# Plots
plot(sillnorm ~ max, all.one1)
plot(sillnorm ~ max, all.one2)
plot(sillnorm ~ max, all.one3)
plot(sillnorm ~ max, all.one4)
plot(sillnorm ~ max, all.one5)
plot(sillnorm ~ max, all.one6)
plot(sillnorm ~ max, all.one7)
plot(sillnorm ~ max, all.one8)
plot(sillnorm ~ max, all.one9)
plot(sillnorm ~ max, all.one10)

plot(sillnorm ~ max, all.two1)
plot(sillnorm ~ max, all.two2)
plot(sillnorm ~ max, all.two3)
plot(sillnorm ~ max, all.two4)
plot(sillnorm ~ max, all.two5)
plot(sillnorm ~ max, all.two6)
plot(sillnorm ~ max, all.two7)
plot(sillnorm ~ max, all.two8)
plot(sillnorm ~ max, all.two9)
plot(sillnorm ~ max, all.two10)

plot(rangenorm ~ max, all.one1)
plot(rangenorm ~ max, all.one2)
plot(rangenorm ~ max, all.one3)
plot(rangenorm ~ max, all.one4)
plot(rangenorm ~ max, all.one5)
plot(rangenorm ~ max, all.one6)
plot(rangenorm ~ max, all.one7)
plot(rangenorm ~ max, all.one8)
plot(rangenorm ~ max, all.one9)
plot(rangenorm ~ max, all.one10)

plot(rangenorm ~ max, all.two1)
plot(rangenorm ~ max, all.two2)
plot(rangenorm ~ max, all.two3)
plot(rangenorm ~ max, all.two4)
plot(rangenorm ~ max, all.two5)
plot(rangenorm ~ max, all.two6)
plot(rangenorm ~ max, all.two7)
plot(rangenorm ~ max, all.two8)
plot(rangenorm ~ max, all.two9)
plot(rangenorm ~ max, all.two10)

plot(sillnorm ~ max, all.one)
plot(rangenorm ~ max, all.one)

plot(sillnorm ~ max, all.two)
plot(rangenorm ~ max, all.one)


########## III. Run models ##########

# Individual sets
sillone1.lm = lm(sillnorm ~ max, all.one1)
sillone2.lm = lm(sillnorm ~ max, all.one2)
sillone3.lm = lm(sillnorm ~ max, all.one3)
sillone4.lm = lm(sillnorm ~ max, all.one4)
sillone5.lm = lm(sillnorm ~ max, all.one5)
sillone6.lm = lm(sillnorm ~ max, all.one6)
sillone7.lm = lm(sillnorm ~ max, all.one7)
sillone8.lm = lm(sillnorm ~ max, all.one8)
sillone9.lm = lm(sillnorm ~ max, all.one9)
sillone10.lm = lm(sillnorm ~ max, all.one10)

silltwo1.lm = lm(sillnorm ~ max, all.two1)
silltwo2.lm = lm(sillnorm ~ max, all.two2)
silltwo3.lm = lm(sillnorm ~ max, all.two3)
silltwo4.lm = lm(sillnorm ~ max, all.two4)
silltwo5.lm = lm(sillnorm ~ max, all.two5)
silltwo6.lm = lm(sillnorm ~ max, all.two6)
silltwo7.lm = lm(sillnorm ~ max, all.two7)
silltwo8.lm = lm(sillnorm ~ max, all.two8)
silltwo9.lm = lm(sillnorm ~ max, all.two9)
silltwo10.lm = lm(sillnorm ~ max, all.two10)

rangeone1.lm = lm(rangenorm ~ max, all.one1)
rangeone2.lm = lm(rangenorm ~ max, all.one2)
rangeone3.lm = lm(rangenorm ~ max, all.one3)
rangeone4.lm = lm(rangenorm ~ max, all.one4)
rangeone5.lm = lm(rangenorm ~ max, all.one5)
rangeone6.lm = lm(rangenorm ~ max, all.one6)
rangeone7.lm = lm(rangenorm ~ max, all.one7)
rangeone8.lm = lm(rangenorm ~ max, all.one8)
rangeone9.lm = lm(rangenorm ~ max, all.one9)
rangeone10.lm = lm(rangenorm ~ max, all.one10)

rangetwo1.lm = lm(rangenorm ~ max, all.two1)
rangetwo2.lm = lm(rangenorm ~ max, all.two2)
rangetwo3.lm = lm(rangenorm ~ max, all.two3)
rangetwo4.lm = lm(rangenorm ~ max, all.two4)
rangetwo5.lm = lm(rangenorm ~ max, all.two5)
rangetwo6.lm = lm(rangenorm ~ max, all.two6)
rangetwo7.lm = lm(rangenorm ~ max, all.two7)
rangetwo8.lm = lm(rangenorm ~ max, all.two8)
rangetwo9.lm = lm(rangenorm ~ max, all.two9)
rangetwo10.lm = lm(rangenorm ~ max, all.two10)

# All sets
sillone.lm = lm(sillnorm ~ max, all.one)
silltwo.lm = lm(sillnorm ~ max, all.two)

rangeone.lm = lm(rangenorm ~ max, all.one)
rangetwo.lm = lm(rangenorm ~ max, all.two)

# Save models
sillone = list(sillone1.lm, sillone2.lm, sillone3.lm, sillone4.lm, sillone5.lm, sillone6.lm, sillone7.lm, sillone8.lm, sillone9.lm, sillone10.lm)
rangeone = list(rangeone1.lm, rangeone2.lm, rangeone3.lm, rangeone4.lm, rangeone5.lm, rangeone6.lm, rangeone7.lm, rangeone8.lm, rangeone9.lm, rangeone10.lm)

save(object = sillone, file = "Output/LM/sillone_lm.rda")
save(object = rangeone, file = "Output/LM/rangeone_lm.rda")
save(object = sillone.lm, file = "Output/LM/sillone_all_lm.rda")
save(object = rangeone.lm, file = "Output/LM/rangeone_all_lm.rda")

silltwo = list(silltwo1.lm, silltwo2.lm, silltwo3.lm, silltwo4.lm, silltwo5.lm, silltwo6.lm, silltwo7.lm, silltwo8.lm, silltwo9.lm, silltwo10.lm)
rangetwo = list(rangetwo1.lm, rangetwo2.lm, rangetwo3.lm, rangetwo4.lm, rangetwo5.lm, rangetwo6.lm, rangetwo7.lm, rangetwo8.lm, rangetwo9.lm, rangetwo10.lm)

save(object = silltwo, file = "Output/LM/silltwo_lm.rda")
save(object = rangetwo, file = "Output/LM/rangetwo_lm.rda")
save(object = silltwo.lm, file = "Output/LM/silltwo_all_lm.rda")
save(object = rangetwo.lm, file = "Output/LM/rangetwo_all_lm.rda")


########## IV. Results table ##########

# Create table
res.lm = data.frame(sigma = c(round(sigma_one,2), "All", round(sigma_two,2), "All"))
res.lm$sill = round(c(coef(sillone1.lm)[2], coef(sillone2.lm)[2], coef(sillone3.lm)[2], coef(sillone4.lm)[2], coef(sillone5.lm)[2], coef(sillone6.lm)[2], coef(sillone7.lm)[2], coef(sillone8.lm)[2], coef(sillone9.lm)[2], coef(sillone10.lm)[2], coef(sillone.lm)[2],
			coef(silltwo1.lm)[2], coef(silltwo2.lm)[2], coef(silltwo3.lm)[2], coef(silltwo4.lm)[2], coef(silltwo5.lm)[2], coef(silltwo6.lm)[2], coef(silltwo7.lm)[2], coef(silltwo8.lm)[2], coef(silltwo9.lm)[2], coef(silltwo10.lm)[2], coef(silltwo.lm)[2]), 2)
res.lm$psill = round(c(getp(sillone1.lm), getp(sillone2.lm), getp(sillone3.lm), getp(sillone4.lm), getp(sillone5.lm), getp(sillone6.lm), getp(sillone7.lm), getp(sillone8.lm), getp(sillone9.lm), getp(sillone10.lm), getp(sillone.lm),
			getp(silltwo1.lm), getp(silltwo2.lm), getp(silltwo3.lm), getp(silltwo4.lm), getp(silltwo5.lm), getp(silltwo6.lm), getp(silltwo7.lm), getp(silltwo8.lm), getp(silltwo9.lm), getp(silltwo10.lm), getp(silltwo.lm)), 4)
res.lm$range = round(c(coef(rangeone1.lm)[2], coef(rangeone2.lm)[2], coef(rangeone3.lm)[2], coef(rangeone4.lm)[2], coef(rangeone5.lm)[2], coef(rangeone6.lm)[2], coef(rangeone7.lm)[2], coef(rangeone8.lm)[2], coef(rangeone9.lm)[2], coef(rangeone10.lm)[2], coef(rangeone.lm)[2],
			coef(rangetwo1.lm)[2], coef(rangetwo2.lm)[2], coef(rangetwo3.lm)[2], coef(rangetwo4.lm)[2], coef(rangetwo5.lm)[2], coef(rangetwo6.lm)[2], coef(rangetwo7.lm)[2], coef(rangetwo8.lm)[2], coef(rangetwo9.lm)[2], coef(rangetwo10.lm)[2], coef(rangetwo.lm)[2]), 2)
res.lm$prange = round(c(getp(rangeone1.lm), getp(rangeone2.lm), getp(rangeone3.lm), getp(rangeone4.lm), getp(rangeone5.lm), getp(rangeone6.lm), getp(rangeone7.lm), getp(rangeone8.lm), getp(rangeone9.lm), getp(rangeone10.lm), getp(rangeone.lm),
			getp(rangetwo1.lm), getp(rangetwo2.lm), getp(rangetwo3.lm), getp(rangetwo4.lm), getp(rangetwo5.lm), getp(rangetwo6.lm), getp(rangetwo7.lm), getp(rangetwo8.lm), getp(rangetwo9.lm), getp(rangetwo10.lm), getp(rangetwo.lm)), 4)

# Save table
write.table(res.lm, "Output/LM/reslm.txt")

# Convert table to latex xode
xtable(res.lm)
