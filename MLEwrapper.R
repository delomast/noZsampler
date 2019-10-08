#' wrapper to process data and determine mle estimates

MLEwrapper <- function(trapData, tags, GSIcol, PBTcol, strataCol, optimMethod, optimMaxIter){
	
	#don't need all the output from this, but it includes most things we need
	allInput <- prepStrata(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", variableCols = c(), variableColsOth = c(), "AdClip",
								AI = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA, verbose = FALSE, symPrior = .5)
	
	estimates <- list()
	#make estiamte for each strata
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
		
		#define some values used in function
		nGSI <- length(input$groups) - nPBT
		
		#number of parameters to optimize
		nParam <- nPBT + nGSI - 1
		
		# create pbtGSIkey
		pbtGSIkey <- list()
		if(nPBT > 0){
			for(i in 1:nPBT){
				pbtGSIkey[[i]] <- which(input$ohnc_gsi[i,] > 0)
				nParam <- nParam + length(pbtGSIkey[[i]]) - 1
			}
		}
		
		# return(list(nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey))
		
		# find mle
		tempFit <- optim(rep(1,nParam), flex_negllh, method = optimMethod, control = list(maxit = optimMaxIter),
							 #arguments to pass to flex_negllh
							 nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey)
		if(tempFit$convergence != 0) cat("\nOptimizer gave convergence code of", tempFit$convergence, "in strata", input$strataName, "\n")
		#unpack values
		# unpack piTot and piGSI estimates
		ptestim <- c(1, tempFit$par[1:(nPBT + nGSI - 1)])
		ptestim <- ptestim / sum(ptestim)
		names(ptestim) <- input$groupsKey[,1]
		
		gsiParams <- tempFit$par[(nPBT + nGSI):length(tempFit$par)]
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
				}
			}
		}
		piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
		colnames(piGSItemp) <- input$GSIkey[,1]
		rownames(piGSItemp) <- input$groupsKey[,1]
		
		estimates[[input$strataName]] <- list(piTot = ptestim, piGSI = piGSItemp, strataName = input$strataName)
	}
	
	return(estimates)
}
