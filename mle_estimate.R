# maybe fit by maximum likelihood estimation????
# or use mle to define a weak prior

library(devMCpbt)
library(MCMCpack)
library(bbmle)

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)

multStratData <- data.frame()
tempDataAll <- generatePBTGSIdata(sampRate = .1, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
									 obsTagRates = c(.8, .85,.9), physTagRates = 0,
			    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)
tempData <- tempDataAll[[1]]
tempData$StrataVar <- 1
multStratData <- rbind(multStratData, tempData)

multStratData$GSI <- paste0("GSIgroup", multStratData$GSI)
tags <- tempDataAll[[2]]

# organize data needed
# repurposing function, not everything is used
input <- prepStrata(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", variableCols = c(), variableColsOth = c(), "AdClip",
								AI = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA, verbose = FALSE, symPrior = 0.01)
input <- input[[1]]

#initial values
input$piTotInitial
input$pi_gsiInitial

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

piGSI <- rbind(pbtGSImat, diag(3))

#trying to use R's built in mle function, but leading to lots of errors

negllh_piTot <- function(piTot){
	#piGSI as matrix with rows being groups in same order as piTot,
	# and cols being gsi groups, in same order as utGSI
	#utGSI is sum of untagged fish belonging to that GSI category
	# piTot <- c(a,b,c,d,e,f)
	llh <- 0
	if(sum(piTot < 0 | piTot > 1) != 0) return(Inf)
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc * log(piTot[1:nPBT] * t[1:nPBT]))
	# then utGSI part
	untag <- 1 - t
	for(j in 1:nGSI){
		llh <- llh + log(sum(piTot * untag * piGSI[,j])) * utGSI[j]
	}
	return(-llh)
}

negllh_piTot(rep(1/6,6))
negllh_piTot(truePiTot)

stats4::mle(negllh_piTot, start = list(piTot = input$piTotInitial), fixed = list(ohnc = ohnc, t = t, utGSI = utGSI, nPBT = nPBT, nGSI = nGSI, piGSI = rbind(pbtGSImat, diag(1,3))))

# stats4::mle(negllh_piTot,start = list(a = 1/6,b = 1/6,c = 1/6,d = 1/6,e = 1/6,f = 1/6))
stats4::mle(negllh_piTot, start = list(piTot = rep(1/6,6)))
stats4::mle(negllh_piTot)


fit <- mle2(negllh_piTot, start = list(piTot = rep(1/6,6)))

optim(input$piTotInitial, negllh_piTot, control = list(trace=1))

fitMle <- optim(input$piTotInitial, negllh_piTot)

negllh_piTot(fitMle$par)

optim(rep(1/6,6), negllh_piTot)
# optim works, now need to parameterize in order to have unconstrained parameters instead of proportions

#first try as relative sizes with one fixed
rel_negllh_piTot <- function(piTot){
	#piGSI as matrix with rows being groups in same order as piTot,
	# and cols being gsi groups, in same order as utGSI
	#utGSI is sum of untagged fish belonging to that GSI category
	# piTot <- c(a,b,c,d,e,f)
	piTot <- c(1, piTot) #first group is fixed at 1
	piTot <- piTot / sum(piTot) #transform into proportions
	llh <- 0
	if(sum(piTot < 0 | piTot > 1) != 0) return(Inf)
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc * log(piTot[1:nPBT] * t[1:nPBT]))
	# then utGSI part
	untag <- 1 - t
	for(j in 1:nGSI){
		llh <- llh + log(sum(piTot * untag * piGSI[,j])) * utGSI[j]
	}
	return(-llh)
}

fitRel <- optim(rep(1,5), rel_negllh_piTot, control = list(maxit = 5000, trace = 1), hessian = TRUE)
estimPiTot <- c(1, fitRel$par)
estimPiTot <- estimPiTot / sum(estimPiTot)

negllh_piTot(estimPiTot)
negllh_piTot(truePiTot)

#now try as logs of ratios
logRatio_negllh_piTot <- function(piTotRatios){
	#piGSI as matrix with rows being groups in same order as piTot,
	# and cols being gsi groups, in same order as utGSI
	#utGSI is sum of untagged fish belonging to that GSI category
	# piTot <- c(a,b,c,d,e,f)
	denom <- 1 / sum(exp(piTotRatios)+1)
	piTot <- c()
	for(r in piTotRatios){
		piTot <- c(piTot, denom * exp(r))
	}
	piTot <- c(denom, piTot)
	
	llh <- 0
	if(sum(piTot < 0 | piTot > 1) != 0) return(Inf)
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc * log(piTot[1:nPBT] * t[1:nPBT]))
	# then utGSI part
	untag <- 1 - t
	for(j in 1:nGSI){
		llh <- llh + log(sum(piTot * untag * piGSI[,j])) * utGSI[j]
	}
	return(-llh)
}

fitRatio <- optim(rep(0,5), logRatio_negllh_piTot, control = list(maxit = 5000, trace = 1), hessian = TRUE)

denomR <- 1 / sum(exp(fitRatio$par)+1)
estimPiTotRat <- c()
for(r in fitRatio$par){
	estimPiTotRat <- c(estimPiTotRat, denomR * exp(r))
}
estimPiTotRat <- c(denomR, estimPiTotRat)
	
negllh_piTot(estimPiTotRat)
negllh_piTot(estimPiTot)
negllh_piTot(truePiTot)

det(fitRatio$hessian)
det(fitRel$hessian)

eigen(fitRatio$hessian)
eigen(fitRel$hessian)
optim(fitRatio$par, logRatio_negllh_piTot, control = list(maxit = 5000, trace = 1), hessian = TRUE)

optim(rep(0,5), logRatio_negllh_piTot, control = list(maxit = 5000), hessian = FALSE)
optim(rep(0,5), logRatio_negllh_piTot, method = "BFGS", control = list(maxit = 5000), hessian = FALSE)

optim(rep(1,5), rel_negllh_piTot, control = list(maxit = 5000))
optim(rep(1,5), rel_negllh_piTot, method = "BFGS", control = list(maxit = 5000))

# relative sizes with one fixed seem to lead to better results in this scenario
# BFGS is faster than Nelder-Mead, but both lead to the same results

# now add in piGSI to the likelihood algorithm
# params is piTot and piGSI (for only the PBT groups) (rowwise) concatenated with the first GSI group (column) in each missing (it is set at 1)
rel_negllh_piTot_piGSI <- function(params){
	# first, unpack params
	#piTot
	piTot <- c(1, params[1:(nPBT + nGSI - 1)])
	piTot <- piTot / sum(piTot) #transform into proportions
	if(sum(piTot < 0 | piTot > 1) != 0) return(Inf)
	#piGSI
	piGSItemp <- matrix(params[(nPBT + nGSI):length(params)], nrow = (nPBT), ncol = (nGSI-1), byrow = TRUE)
	piGSItemp <- cbind(1, piGSItemp)
	for(i in 1:nrow(piGSItemp)){
		piGSItemp[i,] <- piGSItemp[i,] / sum(piGSItemp[i,])
		if(sum(piGSItemp[i,] < 0 | piGSItemp[i,] > 1) != 0) return(Inf)
	}
	piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
	llh <- 0
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc * log(piTot[1:nPBT] * t[1:nPBT]))
	# then utGSI part
	untag <- 1 - t
	for(j in 1:nGSI){
		llh <- llh + log(sum(piTot * untag * piGSItemp[,j])) * utGSI[j]
	}
	# then ohnc GSI part
	if(nPBT > 0){
		for(i in 1:nPBT){
			llh <- llh + sum(ohnc_gsi[i,] * log(piGSItemp[i,]))
		}
	}
	
	# returning negative log-likelihood for minimization
	return(-llh)
}

#programming this interatively will be tricky b/c the group that is fixed must be a positively estimated group
# for piTot, this can be an observed pbt group, or any group if no pbt are observed
# for piGSI, this will change from group to group

bfgsFit <- optim(rep(1,11), rel_negllh_piTot_piGSI, method = "BFGS", control = list(maxit = 5000))
nmFit <- optim(rep(1,11), rel_negllh_piTot_piGSI, control = list(maxit = 5000))

# unpack piTot and piGSI estimates
tempFit <- bfgsFit
ptestim <- c(1, tempFit$par[1:(nPBT + nGSI - 1)])
ptestim <- ptestim / sum(ptestim)
gsiEstim <- matrix(tempFit$par[(nPBT + nGSI):length(par)], nrow = nPBT, ncol = (nGSI-1), byrow = TRUE)
gsiEstim <- cbind(1, gsiEstim)
for(i in 1:nrow(gsiEstim)){
	gsiEstim[i,] <- gsiEstim[i,] / sum(gsiEstim[i,])
}
pbtGSImat
ohnc_gsi

tempFit <- nmFit
ptestim <- c(1, tempFit$par[1:(nPBT + nGSI - 1)])
ptestim <- ptestim / sum(ptestim)
gsiEstim <- matrix(tempFit$par[(nPBT + nGSI):length(par)], nrow = nPBT, ncol = (nGSI-1), byrow = TRUE)
gsiEstim <- cbind(1, gsiEstim)
for(i in 1:nrow(gsiEstim)){
	gsiEstim[i,] <- gsiEstim[i,] / sum(gsiEstim[i,])
}


rel_negllh_piTot_piGSI(bfgsFit$par)
rel_negllh_piTot_piGSI(c(bfgsFit$par[1:5], 8, 1, 1/8, 1/8, 1, 8))
rel_negllh_piTot_piGSI(c(truePiTot[2:6]/truePiTot[1], 8, 1, 1/8, 1/8, 1, 8))
# pi_gsi estimate is off, but the neg llh is lower at the estimated parameters
# than with the estimate piTot and the true pi_gsi values
# is something messed up in my llh function?
# can run multiple simulated datasets to see

# try with lots of data

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)

multStratData <- data.frame()
tempDataAll <- generatePBTGSIdata(sampRate = .5, censusSize = 13000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
									 obsTagRates = c(.8, .85,.9), physTagRates = 0,
			    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)
tempData <- tempDataAll[[1]]
tempData$StrataVar <- 1
multStratData <- rbind(multStratData, tempData)

multStratData$GSI <- paste0("GSIgroup", multStratData$GSI)
tags <- tempDataAll[[2]]

# organize data needed
# repurposing function, not everything is used
input <- prepStrata(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", variableCols = c(), variableColsOth = c(), "AdClip",
								AI = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA, verbose = FALSE, symPrior = 0.01)
input <- input[[1]]

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

bfgsFit <- optim(rep(1,11), rel_negllh_piTot_piGSI, method = "BFGS", control = list(maxit = 5000))
nmFit <- optim(rep(1,11), rel_negllh_piTot_piGSI, control = list(maxit = 5000))

# unpack piTot and piGSI estimates
tempFit <- bfgsFit
ptestim <- c(1, tempFit$par[1:(nPBT + nGSI - 1)])
ptestim <- ptestim / sum(ptestim)
gsiEstim <- matrix(tempFit$par[(nPBT + nGSI):length(tempFit$par)], nrow = nPBT, ncol = (nGSI-1), byrow = TRUE)
gsiEstim <- cbind(1, gsiEstim)
for(i in 1:nrow(gsiEstim)){
	gsiEstim[i,] <- gsiEstim[i,] / sum(gsiEstim[i,])
}
ptestim
truePiTot

gsiEstim
pbtGSImat
ohnc_gsi

# pi tot looks good, pi_gsi does not. maybe llh function is missing something
### missign piece was mistake in unpacking GSI estimates, fixed in a few places
### wrong: length(par) correct: length(tempFit$par)


## nm did not converge
tempFit <- nmFit
ptestim <- c(1, tempFit$par[1:(nPBT + nGSI - 1)])
ptestim <- ptestim / sum(ptestim)
gsiEstim <- matrix(tempFit$par[(nPBT + nGSI):length(tempFit$par)], nrow = nPBT, ncol = (nGSI-1), byrow = TRUE)
gsiEstim <- cbind(1, gsiEstim)
for(i in 1:nrow(gsiEstim)){
	gsiEstim[i,] <- gsiEstim[i,] / sum(gsiEstim[i,])
}

ptestim
truePiTot

gsiEstim
pbtGSImat
ohnc_gsi


### now, need to:
###  compare against SD
###  make more flexible
###  allow handling of cats with 0 in the ohncGSI data for BFGS
###  reasonable wrapper script(s)


# making the function that calculates the likelihood more flexible
#' @param params list of paramaters to optimize
#' @param nPBT number of pbt groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc list of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t list of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param utGSI list of number of un-PBT assigned fish in each GSI group
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#' @param pbtGSIkey list with a vector telling which GSI groups (as integers giving their position in the order)
#'  are nonzero for each pbt group

flex_negllh <- function(params, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi, pbtGSIkey){
	# first, unpack params
	#piTot
	piTot <- c(1, params[1:(nPBT + nGSI - 1)])
	piTot <- piTot / sum(piTot) #transform into proportions
	if(sum(piTot < 0 | piTot > 1) != 0) return(Inf)
	#piGSI
	# piGSItemp <- matrix(params[(nPBT + nGSI):length(params)], nrow = (nPBT), ncol = (nGSI-1), byrow = TRUE)
	gsiParams <- params[(nPBT + nGSI):length(params)]
	piGSItemp <- matrix(0, nrow = (nPBT), ncol = (nGSI)) #initiate with zeros
	if (nPBT > 0){
		for(i in 1:nPBT){
			key <- pbtGSIkey[[i]]
			piGSItemp[i,key[1]] <- 1 #first non-zero group is fixed
			lk <- length(key)
			if(lk > 1){
				piGSItemp[i,key[2:lk]] <- gsiParams[1:(lk-1)]
				gsiParams <- gsiParams[lk:length(gsiParams)] #bump entries forward
				piGSItemp[i,] <- piGSItemp[i,] / sum(piGSItemp[i,]) #normalize
				if(sum(piGSItemp[i,] < 0 | piGSItemp[i,] > 1) != 0) return(Inf) #make sure all entries are valid
			}
		}
	}
	piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
	
	# now, caluculate the log likelihood
	llh <- 0
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc * log(piTot[1:nPBT] * t[1:nPBT]))
	# then utGSI part
	untag <- 1 - t
	for(j in 1:nGSI){
		llh <- llh + log(sum(piTot * untag * piGSItemp[,j])) * utGSI[j]
	}
	# then ohnc GSI part
	if(nPBT > 0){
		for(i in 1:nPBT){
			llh <- llh + sum(ohnc_gsi[i,] * log(piGSItemp[i,]))
		}
	}
	
	# returning negative log-likelihood for minimization
	return(-llh)
}


pbtGSIkey <- list(c(1,2,3), c(1,2,3), c(1,2,3))

flex_negllh(rep(1,11), nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi, pbtGSIkey)

# now automatically create pbtGSIkey
pbtGSIkey <- list()
if(nPBT > 0){
	for(i in 1:nPBT){
		pbtGSIkey[[i]] <- which(input$ohnc_gsi[i,] > 0)
	}
}

#now create a wrapper for this in a new file (MLEwrapper.R)
# now need to allow optional use of variables

# now test
source("./flex_negllh.R")
source("./flex_negllh_var.R")
source("./MLEwrapper.R")

varMat <- list(
	matrix(c(rep(c(.1,.9), 3), rep(c(.9,.1), 3)), nrow = 6, ncol = 2, byrow = TRUE),
	matrix(c(rep(c(.4,.6), 3), rep(c(.6,.4), 3)), nrow = 6, ncol = 2, byrow = TRUE)
)
pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)

multStratData <- data.frame()
tempDataAll <- generatePBTGSIdata(sampRate = .8, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
									 obsTagRates = c(.8, .85,.9), physTagRates = 0,
			    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat)
tempData <- tempDataAll[[1]]
tempData$StrataVar <- 1
multStratData <- rbind(multStratData, tempData)

multStratData$GSI <- paste0("GSIgroup", multStratData$GSI)
tags <- tempDataAll[[2]]

# MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", 5000)
# 
# 
# 
# testInput <- prepStrata(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", variableCols = c("Var1"), variableColsOth = c(), "AdClip",
# 								AI = TRUE, GSIgroups = NA,
# 									 variableValues = NA, variableValuesOth = NA, verbose = FALSE, symPrior = 0.5)
# tI <- testInput[[1]]


res <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", variableCols = c("Var1"), control = list(maxit = 5000))
res <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "CG", variableCols = c("Var1"), control = list(maxit = 5000))
res <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "SANN", variableCols = c("Var1"), control = list(maxit = 10000))

res <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "Nelder-Mead", variableCols = c("Var1"), control = list(maxit = 5000))



MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "Nelder-Mead", variableCols = c("Var1"), control = list(maxit = 5000), old = TRUE)
MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", variableCols = c("Var1"), control = list(maxit = 5000), old = TRUE)

st <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", variableCols = c("Var1"), control = list(maxit = 5000))$starting
varKey <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", variableCols = c("Var1"), control = list(maxit = 5000))$varKey
nCat <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", variableCols = c("Var1"), control = list(maxit = 5000))$nCat
utVar <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", variableCols = c("Var1"), control = list(maxit = 5000))$utVar
ohnc_var <- MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", variableCols = c("Var1"), control = list(maxit = 5000))$ohnc_var


MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", variableCols = c("Var1"), control = list(maxit = 5000), old = TRUE)[[1]]$piTot


## adding the variables makes the estimates worse
### is there something wrong in the llh equation?

pt <- c(1, st[1:5])
pt <- pt / sum(pt)

# piGSI
subParams <- st[(nPBT + nGSI):length(st)]
piGSItemp <- matrix(0, nrow = (nPBT), ncol = (nGSI)) #initiate with zeros
if(nPBT > 0){
	for(i in 1:nPBT){
		key <- pbtGSIkey[[i]]
		piGSItemp[i,key[1]] <- 1 #first non-zero group is fixed
		lk <- length(key)
		if(lk > 1){
			piGSItemp[i,key[2:lk]] <- subParams[1:(lk-1)]
			subParams <- subParams[lk:length(subParams)] #bump entries forward
			piGSItemp[i,] <- piGSItemp[i,] / sum(piGSItemp[i,]) #normalize
		}
	}
}
piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
colnames(piGSItemp) <- input$GSIkey[,1]
rownames(piGSItemp) <- input$groupsKey[,1]

#piVar
piVarList <- list()
if(nVar > 0){
	for(v in 1:nVar){
		piVarTemp <- matrix(0, nrow = (nPBT + nGSI), ncol = nCat[[v]]) #initiate with zeros
		tempKey <- varKey[[v]]
		for(i in 1:(nPBT +  nGSI)){
			key <- tempKey[[i]]
			piVarTemp[i,key[1]] <- 1 #first non-zero group is fixed
			lk <- length(key)
			if(lk > 1){
				piVarTemp[i,key[2:lk]] <- subParams[1:(lk-1)]
				subParams <- subParams[lk:length(subParams)] #bump entries forward
				piVarTemp[i,] <- piVarTemp[i,] / sum(piVarTemp[i,]) #normalize
				if(sum(piVarTemp[i,] < 0 | piVarTemp[i,] > 1) != 0) return(Inf) #make sure all entries are valid
			}
		}
		piVarList[[v]] <- piVarTemp
	}
}


#calculate true piTot proportions
pbtrel <- c(1,2,3)
gsirel <- c(1,2,1)
trueHNC <- .3
trueW <- .7
pbtrel <- (pbtrel / sum(pbtrel)) * trueHNC
gsirel <- (gsirel / sum(gsirel)) * trueW
truePiTot <- c(pbtrel, gsirel)
truePiTot


