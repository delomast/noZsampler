### writing MCMC sampler for composition with tag rates that doesn't use latent variable Z's
### trying to overcome orthogonal moves of Gibbs sampler which can lead to estimates for
### groups getting stuck near zero
## codename: no sleep till Brooklyn
## codename: noZfish

# loading for data generation function and potential comparisons
library(fishCompTools)
library(devMCpbt)

library(MCMCpack)

# so, for one strata, no variables, just PBT and GSI assignments

# generate some data

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)

nStrata <- 1
multStratData <- data.frame()
for(i in 1:nStrata){
	tempDataAll <- generatePBTGSIdata(sampRate = .025, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .9, true_wild = .1, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)
	tempData <- tempDataAll[[1]]
	tempData$StrataVar <- i
	# tempData$GSI <- 1
	multStratData <- rbind(multStratData, tempData)
}
multStratData$GSI <- paste0("GSIgroup", multStratData$GSI)
tags <- tempDataAll[[2]]

# organize data needed
# repurposing function, not everything is used
input <- prepStrata(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", variableCols = c(), variableColsOth = c(), "AdClip",
								AI = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA, verbose = FALSE, symPrior = 0.01)
input <- input[[1]]

# define prior for piTot
prior_piTot <- function(piTot){
	#let's still use a dirichlet for now
	if(sum(piTot < 0 | piTot > 1) != 0) return(0)
	return(ddirichlet(piTot, rep(1/length(piTot), length(piTot))))
}
# define log likelihood for piTot
llh_piTot <- function(piTot, t, ohnc, piGSI, utGSI, nPBT, nGSI){
	#piGSI as matrix with rows being groups in same order as piTot,
	# and cols being gsi groups, in same order as utGSI
	#utGSI is sum of untagged fish belonging to that GSI category
	llh <- 0
	if(sum(piTot < 0 | piTot > 1) != 0) return(llh)
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc * log(piTot[1:nPBT] * t[1:nPBT]))
	# then utGSI part
	untag <- 1 - t
	for(j in 1:nGSI){
		llh <- llh + log(sum(piTot * untag * piGSI[,j])) * utGSI[j]
	}
	return(llh)
}

#proposal for piTot
proposal_piTot <- function(stdevs, bigstdev){
	# nGroups <- length(stdevs)
	# rec <- matrix(NA, nrow = nGroups, ncol = 2)
	# for(i in 1:nGroups){
	# 	rec[i,] <- rnorm(2, 0, stdevs[i])
	# }
	# rec[,1] <- rec[,1] / sum(rec[,1])
	# rec[,2] <- rec[,2] / sum(rec[,2])
	# rec <- rec[,1] - rec[,2]
	# rec <- rec * rnorm(1,0, bigstdev)
	rec <- rdirichlet(2, rep(1, length(stdevs)))
	rec <- rec[1,] - rec[2,]
	return(rec * rnorm(1,0,bigstdev))
}


reps <- 10000
r_piTot <- matrix(nrow = reps, ncol = length(input$groups))
r_accept <- rep(NA, reps)

#initial values
piTot <- input$piTotInitial
piGSI <- input$pi_gsiInitial

#pull values out of input
ohnc <- input$ohnc
nPBT <- input$nPBT
t <- input$t
ohnc_gsi <- input$ohnc_gsi
utGSI <- c()
for(g in input$GSI_values){
	utGSI <- c(utGSI, sum(input$gsiUT == g))
}

#define some values used in function
nGSI <- length(input$groups) - nPBT
stdevs <- rep(10, nPBT + nGSI) # will have to have some sort of optimization for these
bigstdev <- .03

set.seed(7)

for(r in 1:reps){
	#propose new value
	piTotNew <- piTot + proposal_piTot(stdevs, bigstdev)
	
	#calculate ratio
	a <- (log(prior_piTot(piTotNew)) + llh_piTot(piTotNew, t, ohnc, piGSI, utGSI, nPBT, nGSI)) - 
				(log(prior_piTot(piTot)) + llh_piTot(piTot, t, ohnc, piGSI, utGSI, nPBT, nGSI))
	a <- exp(a)
	if(runif(1) < a){
		piTot <- piTotNew
		r_accept[r] <- TRUE
	} else {
		r_accept[r] <- FALSE
	}
	
	# sample from piGSI - just performing Gibbs with the ohnc right now as a temporary measure
	#			with improper 0,0,...,0 prior
	if(nPBT > 0){
		for(i in 1:nPBT){
			piGSI[i,] <- as.vector(rdirichlet(1, ohnc_gsi[i,]))
		}
	}

	
	#record values
	r_piTot[r,] <- piTot
	
	
}

sum(r_accept)/reps

plot(r_piTot[,3])
plot(r_piTot[,5])

#calculate true piTot proportions
pbtrel <- c(1,2,3)
gsirel <- c(1,2,1)
trueHNC <- .9
trueW <- .1
pbtrel <- (pbtrel / sum(pbtrel)) * trueHNC
gsirel <- (gsirel / sum(gsirel)) * trueW
truePiTot <- c(pbtrel, gsirel)
truePiTot

start <- 4500

apply(r_piTot[start:reps,], 2, mean)
lower <- apply(r_piTot[start:reps,], 2, quantile, .05)
upper <- apply(r_piTot[start:reps,], 2, quantile, .95)
cbind(round(lower,3), truePiTot, round(upper,3))

plot(r_accept)


varMat <- list(
	matrix(c(rep(c(.1,.9), 3), rep(c(.9,.1), 3)), nrow = 6, ncol = 2, byrow = TRUE),
	matrix(c(rep(c(.4,.6), 3), rep(c(.6,.4), 3)), nrow = 6, ncol = 2, byrow = TRUE)
)

plot(r_piTot[,6])


prior_piTot(r_piTot[10000,])
prior_piTot(r_piTot[1000,])


### looking at different proposals

testPropose <- function(){
	# rec <- rdirichlet(2, c(1,1,1,1,1,1))
	rec <- rdirichlet(1, c(3,1,1,1,1,1))[1,] - rep(1/6,6)
	return(rec * rnorm(1,0,.03))
}

hist(replicate(10000, testPropose())[1,])

## note that as difference in alphas becomes large, can result in groups 
##  that only move in opposite directions of each other b/c 1 is always greater
##  than 1/N and the other is always lower than 1/N
## will that cause problems?
## maybe throw in a proposal from all alphas equal 1 every now and then?