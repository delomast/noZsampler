### run lots of tests with noZ

# loading for data generation function and potential comparisons
library(devMCpbt)

library(MCMCpack)

# define prior for piTot
prior_piTot <- function(piTot){
	
	# if(sum(piTot <= 0 | piTot >= 1) != 0) return(-Inf)
	# # try in logit
	# prior <- 0
	# piTot <- log(piTot / (1-piTot))
	# for(p in piTot){
	# 	prior <- prior + log(dnorm(p, 0, 1000000))
	# }
	# return(prior)
	#let's still use a dirichlet for now
	if(sum(piTot < 0 | piTot > 1) != 0) return(-Inf)
	return(log(ddirichlet(piTot, rep(.01, length(piTot)))))
}
# define log likelihood for piTot
llh_piTot <- function(piTot, t, ohnc, piGSI, utGSI, nPBT, nGSI){
	#piGSI as matrix with rows being groups in same order as piTot,
	# and cols being gsi groups, in same order as utGSI
	#utGSI is sum of untagged fish belonging to that GSI category
	llh <- 0
	if(sum(piTot < 0 | piTot > 1) != 0) return(-Inf)
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
	
	
reps <- 8000
start <- 4500

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)

trials <- 50

means <- matrix(NA, nrow = trials, ncol = 6)

for(tr in 1:trials){
	print(tr)
	nStrata <- 1
	multStratData <- data.frame()
	for(i in 1:nStrata){
		tempDataAll <- generatePBTGSIdata(sampRate = .1, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
											 obsTagRates = c(.8, .85,.9), physTagRates = 0,
					    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)
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
	
	if(length(input$groups) != 6) next
	
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
	

	for(r in 1:reps){
		#propose new value
		piTotNew <- piTot + proposal_piTot(stdevs, bigstdev)
		
		#calculate ratio
		# a <- (log(prior_piTot(piTotNew)) + llh_piTot(piTotNew, t, ohnc, piGSI, utGSI, nPBT, nGSI)) - 
		# 			(log(prior_piTot(piTot)) + llh_piTot(piTot, t, ohnc, piGSI, utGSI, nPBT, nGSI))
		a <- (prior_piTot(piTotNew) + llh_piTot(piTotNew, t, ohnc, piGSI, utGSI, nPBT, nGSI)) - 
					(prior_piTot(piTot) + llh_piTot(piTot, t, ohnc, piGSI, utGSI, nPBT, nGSI))
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
	
	means[tr,] <- apply(r_piTot[start:reps,], 2, mean)
	
	
}
	
sum(r_accept)/reps

plot(r_piTot[,3])
plot(r_piTot[,5])

#calculate true piTot proportions
pbtrel <- c(1,2,3)
gsirel <- c(1,2,1)
trueHNC <- .3
trueW <- .7
pbtrel <- (pbtrel / sum(pbtrel)) * trueHNC
gsirel <- (gsirel / sum(gsirel)) * trueW
truePiTot <- c(pbtrel, gsirel)
truePiTot

start <- 4500

boxplot(means)
abline(h=truePiTot)

apply(means,2, mean)

apply(r_piTot[start:reps,], 2, mean)
lower <- apply(r_piTot[start:reps,], 2, quantile, .05)
upper <- apply(r_piTot[start:reps,], 2, quantile, .95)
cbind(round(lower,3), truePiTot, round(upper,3))


