#' define negative log-likelihood function
#' this version estimates piTot and piGSI, but not other variables
#' this version fixes piGSI estiamtes for PBT groups that are not observed in ohnc_GSI fixed at zero
#' @param params list of paramaters to optimize
#' @param nPBT number of pbt groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc list of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t list of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param utGSI list of number of un-PBT assigned fish in each GSI group
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#' @param pbtGSIkey list with a vector telling which GSI groups (as integers giving their position in the order)
#'  are nonzero for each pbt group
#'  

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
	if(nPBT > 0){
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
			# selecting by pbtGSIkey[[i]] to avoid unobserved groups (which are estimated at 0)
			llh <- llh + sum(ohnc_gsi[i,pbtGSIkey[[i]]] * log(piGSItemp[i,pbtGSIkey[[i]]]))
		}
	}
	
	# returning negative log-likelihood for minimization
	return(-llh)
}
