#' wrapper to process data and determine mle estimates
#' calls different likelihood function depending on if there is a variable or not
#' @param trapData a dataframe with a rwo for each individual and columns for GSI assignment, PBT assignment, etc.
#' @param tags a dataframe with the first column containing names of PBT groups, and the second column containing
#'   tag rates
#' @param GSIcol
#' @param PBTcol
#' @param strataCol
#' @param adFinCol
#' @param AI
#' @param optimMethod the method to use first with \code{optim}. If this method fails, "Nelder-Mead" is attempted
#'   and a warning is issued.
#' @param variableCols
#' @param ... other arguments to pass to \code{optim}

MLEwrapper <- function(trapData, tags, GSIcol, PBTcol, strataCol, adFinCol, AI = TRUE, optimMethod = "Nelder-Mead",
							  variableCols = NULL, ...){
	
	# determine if variables used or not
	varBool <- !is.null(variableCols) #TRUE if variables are used
	
	# input checking
	if(varBool && length(variableCols) > 1){
		stop("variableCols must be either NULL or only one variable")
	}
	
	#don't need all the output from this, but it includes most things we need
	allInput <- prepStrata(trapData, tags, GSIcol, PBTcol, strataCol, variableCols = variableCols, variableColsOth = c(), adFinCol,
								AI = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA, verbose = FALSE, symPrior = .5)
	
	estimates <- list()
	#get estimates for each strata
	for(input in allInput){
		#pull values out of input
		ohnc <- input$ohnc
		nPBT <- input$nPBT
		t <- input$t
		ohnc_gsi <- input$ohnc_gsi
		utGSI <- c()
		for(g in input$GSI_values){
			utGSI <- c(utGSI, sum(input$gsiUT == g))
		}

		#define some values used in the function
		nGSI <- length(input$groups) - nPBT
		
		
		#pull variables values out of input
		ohnc_var <- data.frame() #define in case no PBT groups
		if(varBool){
			nCat <- length(input$values[[1]])
			if(nPBT > 0) ohnc_var <- input$pi_Vohnc[[1]][1:nPBT,]
			utVar <- matrix(0, nrow = nGSI, ncol = nCat)
			for(v in 1:nCat){
				val <- input$values[[1]][v]
				for(g in 1:nGSI){
					gVal <- input$GSI_values[g]
					utVar[g,v] <- sum(input$v_ut[,1] == val & input$gsiUT == gVal)
				}
			}
		}
		#number of parameters to optimize - this was before reasonable starting values were implemented and
		#  start was all equal
		# nParam <- nPBT + nGSI - 1
		
		# create pbtGSIkey
		pbtGSIkey <- list()
		if(nPBT > 0){
			for(i in 1:nPBT){
				pbtGSIkey[[i]] <- which(input$ohnc_gsi[i,] > 0)[1]
				# nParam <- nParam + length(pbtGSIkey[[i]]) - 1
			}
		}

		# create varKey
		if(varBool){
			varKey <- list()
			if(nPBT > 0){
				for(i in 1:nPBT){
					varKey[[i]] <- which(input$pi_Vohnc[[1]][i,] > 0)[1]
					# nParam <- nParam + length(varKey[[i]]) - 1
				}
			}
			for(i in 1:nGSI){
				#determine which categories have un-tagged fish with this GSI assignment
				tempCounts <- c()
				for(val in input$values[[1]]){
					tempCounts <- c(tempCounts, sum(input$gsiUT == input$GSI_values[i] & input$v_ut[,1] == val))
				}
				varKey[[nPBT + i]] <- which(tempCounts > 0)[1]
				# nParam <- nParam + length(varKey[[nPBT + i]]) - 1
			}
		}
		
		# determine reasonable starting values
		#for piTot
		start_piTot <- c()
		if(nPBT > 0) start_piTot <- ohnc[1:nPBT] / t[1:nPBT] #scale PBT by tag rates
		start_piTot <- c(start_piTot, utGSI) #just use observed GSI
		start_piTot <- start_piTot / sum(start_piTot) #normalize
		start_piTot <- start_piTot / start_piTot[1] # scale against first group
		start_piTot <- start_piTot[2:length(start_piTot)] #remove first group
		
		#for piGSI
		start_piGSI <- c()
		for(i in 1:nrow(ohnc_gsi)){
			temp <- ohnc_gsi[i,] #just use ohnc assignments
			temp[temp < 1] <- .1
			temp <- temp / sum(temp)
			temp <- temp / temp[pbtGSIkey[[i]]]
			pos <- 1:length(temp)
			pos <- pos[pos != pbtGSIkey[[i]]]
			temp <- temp[pos]
			start_piGSI <- c(start_piGSI, temp)
		}
		
		#for piVar
		start_piVar <- c()
		if(varBool){
			for(i in 1:nrow(ohnc_var)){
				temp <- ohnc_var[i,] #just use ohnc assignments
				temp[temp < 1] <- .1
				temp <- temp / sum(temp)
				temp <- temp / temp[varKey[[i]]]
				pos <- 1:length(temp)
				pos <- pos[pos != varKey[[i]]]
				temp <- temp[pos]
				start_piVar <- c(start_piVar, temp)
			}
			# then ut fish
			for(i in 1:nGSI){
				#determine which categories have un-tagged fish with this GSI assignment
				temp <- c()
				for(val in input$values[[1]]){
					temp <- c(temp, sum(input$gsiUT == input$GSI_values[i] & input$v_ut[,1] == val))
				}
				temp <- temp / sum(temp)
				temp <- temp[temp > 0]
				temp <- temp / temp[1]
				if(length(temp) > 1) start_piVar <- c(start_piVar, temp[2:length(temp)])
			}
		}


		# find mle
		if(varBool){
			tempFit <- optim(c(start_piTot, start_piGSI, start_piVar), flex_negllh_var, method = optimMethod, ...,
					 #arguments to pass to flex_negllh
					 nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey,
					 utVar = utVar, ohnc_var = ohnc_var, varKey = varKey, nCat = nCat)
		} else {
			tempFit <- optim(c(start_piTot, start_piGSI), flex_negllh_allGSI, method = optimMethod, ...,
								 #arguments to pass to flex_negllh
								 nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey)
		}
		if(tempFit$convergence != 0) cat("\nOptimizer gave convergence code of", tempFit$convergence, "in strata", input$strataName, "\n")
		
		# return(tempFit)


		#unpack values
		# unpack piTot
		ptestim <- c(1, tempFit$par[1:(nPBT + nGSI - 1)])
		ptestim <- ptestim / sum(ptestim)
		names(ptestim) <- input$groupsKey[,1]
		
		# piGSI
		subParams <- tempFit$par[(nPBT + nGSI):length(tempFit$par)]
		piGSItemp <- matrix(0, nrow = (nPBT), ncol = (nGSI)) #initiate with zeros
		if(nPBT > 0){
			for(i in 1:nPBT){
				key <- pbtGSIkey[[i]]
				piGSItemp[i,key] <- 1 #first non-zero group is fixed
				pos <- 1:nGSI
				pos <- pos[pos != key]
				if(length(pos) > 0){
					piGSItemp[i,pos] <- subParams[1:length(pos)]
					subParams <- subParams[(length(pos) + 1):length(subParams)]
				}
				piGSItemp[i,] <- piGSItemp[i,] / sum(piGSItemp[i,])
			}
		}
		piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
		colnames(piGSItemp) <- input$GSIkey[,1]
		rownames(piGSItemp) <- input$groupsKey[,1]
		
		#piVar
		piVar <- matrix(nrow = 0, ncol = 0)
		if(varBool){
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
			}
		}
		
		estimates[[input$strataName]] <- list(piTot = ptestim, piGSI = piGSItemp, piVar = piVar, par = tempFit$par, strataName = input$strataName)
	}
	
	return(estimates)
}
