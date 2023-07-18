# Processing acoustic data

# Created: March 22, 2023
# Last modified: March 24, 2023

# Set working directory
setwd("/Users/epdus/OneDrive/School/Research/Splines")

# Load packages
library(lattice)
library(proxy)
library(ggplot2)

# File layout:		
#	I. Read in acoustic data
#	II. Assign haul to backscatter data
#	III. Calculate proportion of each species in catch for each haul



########## I. Read in acoustic data ##########

# Backscatter data
acoustic = read.csv("Data/acoustic_backscatter.csv")

# Catch data
catch = read.csv("Data/acoustic_catch.csv")

# Haul data
haul = read.csv("Data/acoustic_haul.csv")

# Species codes
species = read.csv("Data/acoustic_species.csv")

# Cross section coefficient
sigma = read.csv("Data/acoustic_sigma.csv")


########## II. Assign haul to backscatter data ##########

# Calculate mean position of the reference hauls
haul_xy = data.frame(x = rowMeans(haul[,c("HaulStartLongitude","HaulStopLatitude")]))
haul_xy$y = rowMeans(haul[,c("HaulStartLongitude","HaulStopLatitude")])
head(haul_xy)

# Split night and day hauls
haul_night = haul_xy[haul$Cluster == "Night", ]
haul_day = haul_xy[haul$Cluster == "Day", ]

# Get minimum distance between hauls and backscatter measurements
night.mindist = apply(dist(haul_night, acoustic[acoustic$Cluster=="Night",c("LogLongitude","LogLatitude")]), 2, which.min)
day.mindist = apply(dist(haul_day, acoustic[acoustic$Cluster=="Day",c("LogLongitude","LogLatitude")]), 2, which.min)

# Assign reference haul to backscatter data
acoustic$Haul = NA
acoustic[acoustic$Cluster=="Night", ]$Haul = as.numeric(rownames(haul_night)[night.mindist])
acoustic[acoustic$Cluster=="Day", ]$Haul = as.numeric(rownames(haul_day)[day.mindist])


########## III. Calculate proportion of each species in catch for each haul ##########

# Assign species
catch$Species = species$Description[match(catch$CatchSpeciesCode, species$Code)]

# Mean cross section from BIAS manual for herring and three-spined stickleback
catch$Sigma = 9.533e-07 * (catch$CatchLengthClass/2)^2

# Split catch by haul number 
catch_haul = split(catch, catch$HaulNumber)

# Add new columns for proportion of each fish species
haul[c("HerringEchoIntegral","HerringCrossSection","StickleEchoIntegral","StickleCrossSection")] = NA

# Get mean proportion weighted by length class for each haul
for(i in 1:nrow(haul)) {
	
	# Get weighted number at length
	tempdat = catch_haul[[i]]
	tempdat$WeightedNumberAtLength = (tempdat$CatchLengthClass * tempdat$CatchNumberAtLength) / sum(tempdat$CatchLengthClass)
	
	# Calculate proportion for each weighted number
	tempdat$Proportion = tempdat$WeightedNumberAtLength / sum(tempdat$WeightedNumberAtLength)
	
	# Calculate cross section
	tempdat$CrossSection = sigma$d[match(tempdat$Species, sigma$Species)] * (tempdat$CatchLengthClass/2)^2
	
	# Split each haul by species
	tempdat.sp = split(tempdat, tempdat$Species)
	
	# Calculate weighted number at length
	prop = unlist(lapply(tempdat.sp, function(x) sum(x$Proportion)))
	sigma_mean = unlist(lapply(tempdat.sp, function(x) mean(x$CrossSection)))
	
	# Get proportion of echo integral allocated to herring and three-spined stickleback
	echo_herring = (prop["Clupea harengus"] * sigma_mean["Clupea harengus"]) / sum(prop*sigma_mean)
	echo_stickle = (prop["Gasterosteus aculeatus"] * sigma_mean["Gasterosteus aculeatus"]) / sum(prop*sigma_mean)
	
	# Get mean cross section for herring and three-spined stickleback
	sigma_herring = sigma_mean["Clupea harengus"]
	sigma_stickle = sigma_mean["Gasterosteus aculeatus"]
	
	# Assign to haul data frame
	haul[i,"HerringEchoIntegral"] = ifelse(is.na(echo_herring), 0, echo_herring)
	haul[i,"HerringCrossSection"] = ifelse(is.na(sigma_herring), 0, sigma_herring)
	haul[i,"StickleEchoIntegral"] = ifelse(is.na(echo_stickle), 0, echo_stickle)
	haul[i,"StickleCrossSection"] = ifelse(is.na(sigma_stickle), 0, sigma_stickle)
}

# Add proportion herring and three-spined stickleback echo integral to acoustic data
acoustic$HerringEchoIntegral = haul$HerringEchoIntegral[match(acoustic$Haul, haul$HaulNumber)]
acoustic$StickleEchoIntegral = haul$StickleEchoIntegral[match(acoustic$Haul, haul$HaulNumber)]

# Add herring and three-spined stickleback cross section to acoustic data
acoustic$HerringCrossSection = haul$HerringCrossSection[match(acoustic$Haul, haul$HaulNumber)]
acoustic$StickleCrossSection = haul$StickleCrossSection[match(acoustic$Haul, haul$HaulNumber)]

# Add density of herring and three-spined stickleback to acoustic data
acoustic$HerringDensity = ifelse(acoustic$HerringEchoIntegral == 0, 0, (acoustic$DataValue * acoustic$HerringEchoIntegral) / acoustic$HerringCrossSection)
acoustic$StickleDensity  = ifelse(acoustic$StickleEchoIntegral == 0, 0, (acoustic$DataValue * acoustic$StickleEchoIntegral) / acoustic$StickleCrossSection)

# Make bubble plot based on density
ggplot(acoustic, aes(x=LogLongitude,y=LogLatitude,size=HerringDensity)) +
	geom_point(alpha = 0.5) +
	scale_size(name = "Density", range = c(0.1,5)) +
	theme_bw() +
	theme(text = element_text(size=15)) +
	ggtitle("Herring")
	
# Make bubble plot based on density
ggplot(acoustic, aes(x=LogLongitude,y=LogLatitude,size=StickleDensity)) +
	geom_point(alpha = 0.5) +
	scale_size(name = "Density", range = c(0.1,5)) +
	theme_bw() +
	theme(text = element_text(size=15)) +
	ggtitle("Three-spined stickleback")

# Write density data
write.csv(acoustic, "Data/acoustic_density.csv", row.names = F)