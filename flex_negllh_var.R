#' define negative log-likelihood function
#' this version estimates piTot and piGSI and extra variables that are independent of each other
#' @param params list of paramaters to optimize
#' @param nPBT number of pbt groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t vector of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param utGSI vector of number of un-PBT assigned fish in each GSI group
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#' @param pbtGSIkey list with a vector telling which GSI groups (as integers giving their position in the order)
#' @param utVar list of vectors (one for each variable) of number of un-PBT assigned fish in each category
#' @param ohnc_var list of matrices (one for each variable) of counts of PBT-assigned fish in each category (rows are PBT groups)
#' @param varKey list of lists (one for each variable) with a vector telling which categories (as integers giving their position in the order)
#' @param nCat list of numbers of categories in each variable
#'  are nonzero for each PBT and GSI group
#'  

flex_negllh_var <- function(params, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi, pbtGSIkey,
									 utVar, ohnc_var, varKey, nCat){
	# first, unpack params
	#piTot
	piTot <- c(1, params[1:(nPBT + nGSI - 1)])
	piTot <- piTot / sum(piTot) #transform into proportions
	if(sum(piTot < 0 | piTot > 1) != 0) return(Inf)
	#piGSI
	# piGSItemp <- matrix(params[(nPBT + nGSI):length(params)], nrow = (nPBT), ncol = (nGSI-1), byrow = TRUE)
	subParams <- params[(nPBT + nGSI):length(params)]
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
				if(sum(piGSItemp[i,] < 0 | piGSItemp[i,] > 1) != 0) return(Inf) #make sure all entries are valid
			}
		}
	}
	piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
	#piVar
	piVarList <- list()
	nVar <- length(varKey)
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
			# selecting by pbtGSIkey[[i]] to avoid unobserved groups (which are estimated at 0)
			llh <- llh + sum(ohnc_gsi[i,pbtGSIkey[[i]]] * log(piGSItemp[i,pbtGSIkey[[i]]]))
		}
	}
	# then utVar part
	if(nVar > 0){
		for(v in 1:nVar){
			tempVarMat <- piVarList[[v]] #for code readability
			for(j in 1:nCat[[v]]){
				if(utVar[[v]][j] > 0) {
					tempSumProp <- sum(piTot * untag * tempVarMat[,j])
					if(tempSumProp <= 0) return(Inf)
					llh <- llh + log(tempSumProp) * utVar[[v]][j]
				}
			}
		}
	}
	#then ohnc var part
	if(nVar > 0 && nPBT > 0){
		for(v in 1:nVar){
			tempVarMat <- piVarList[[v]] #for code readability
			tempKey <- varKey[[v]]
			for(i in 1:nPBT){
				# selecting by tempKey[[i]] to avoid unobserved groups (which are estimated at 0)
				llh <- llh + sum(ohnc_var[[v]][i,tempKey[[i]]] * log(tempVarMat[i,tempKey[[i]]]))
			}
		}
	}
	
	# returning negative log-likelihood for minimization
	return(-llh)
}
