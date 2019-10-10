#' define negative log-likelihood function
#' this version estimates piTot and piGSI and one categorical variable
#' @param params list of paramaters to optimize
#' @param nPBT number of pbt groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t vector of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param utGSI vector of number of un-PBT assigned fish in each GSI group
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#' @param pbtGSIkey list with a vector telling which GSI groups (as integers giving their position in the order) are nonzero for
#'   each PBT group
#' @param utVar matrix of numbers of un-PBT assigned fish in each GSI group (row) x category of the variable (column)
#' @param ohnc_var matrix of counts of PBT-assigned fish in each category (rows are PBT groups)
#' @param varKey list of vectors telling which categories (as integers giving their position in the order) are nonzero for each PBT
#'   and GSI group
#' @param nCat the numbers of categories in the variable
#'  

flex_negllh_var <- function(params, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi, pbtGSIkey,
									 utVar, ohnc_var, varKey, nCat){
	# first, unpack params
	#piTot
	piTot <- c(1, params[1:(nPBT + nGSI - 1)])
	piTot <- piTot / sum(piTot) #transform into proportions
	if(sum(piTot < 0 | piTot > 1) != 0) return(Inf)
	#piGSI
	subParams <- params[(nPBT + nGSI):length(params)]
	piGSItemp <- matrix(0, nrow = (nPBT), ncol = (nGSI)) #initiate with zeros
	if(nPBT > 0){
		for(i in 1:nPBT){
			key <- pbtGSIkey[[i]]
			piGSItemp[i,key] <- 1 # group that is fixed
			lk <- length(key)
			tempPos <- 1:nGSI
			tempPos <- tempPos[tempPos != key]
			if(length(tempPos) > 0){
				piGSItemp[i,tempPos] <- subParams[1:length(tempPos)]
				subParams <- subParams[(length(tempPos) + 1):length(subParams)] #bump entries forward
			}
			piGSItemp[i,] <- piGSItemp[i,] / sum(piGSItemp[i,]) #normalize
			if(sum(piGSItemp[i,] < 0 | piGSItemp[i,] > 1) != 0) return(Inf) #make sure all entries are valid
		}
	}
	piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
	#piVar
	piVar <- matrix(0, nrow = (nPBT + nGSI), ncol = nCat) #initiate with zeros
	for(i in 1:(nPBT +  nGSI)){
		key <- varKey[[i]]
		piVar[i,key] <- 1 #first non-zero group is fixed
		pos <- 1:nCat
		pos <- pos[pos != key]
		if(length(pos) > 0){
			piVar[i,pos] <- subParams[1:length(pos)]
			subParams <- subParams[(length(pos) + 1):length(subParams)]
		}
		piVar[i,] <- piVar[i,] / sum(piVar[i,])
		if(sum(piVar[i,] < 0 | piVar[i,] > 1) != 0) return(Inf) #make sure all entries are valid
	}

	# now, calculate the log likelihood
	llh <- 0
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc * log(piTot[1:nPBT] * t[1:nPBT]))
	# then ohnc GSI part
	if(nPBT > 0){
		for(i in 1:nPBT){
			llh <- llh + sum(ohnc_gsi[i,] * log(piGSItemp[i,]))
		}
	}
	#then ohnc var part
	if(nPBT > 0){
		for(i in 1:nPBT){
			llh <- llh + sum(ohnc_var[i,] * log(piVar[i,]))
		}
	}
	# then ut part
	untag <- 1 - t
	for(j in 1:nCat){
		for(k in 1:nGSI){
			if(utVar[k,j] > 0) {
				tempSumProp <- sum(piTot * untag * piGSItemp[,k] * piVar[,j])
				if(tempSumProp <= 0) return(Inf)
				llh <- llh + log(tempSumProp) * utVar[k,j]
			}
		}
	}


	
	# returning negative log-likelihood for minimization
	return(-llh)
}


