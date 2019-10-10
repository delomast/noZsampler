#' wrapper to process data and determine mle estimates

MLEwrapper_allGSI <- function(trapData, tags, GSIcol, PBTcol, strataCol, optimMethod = "BFGS",
							  variableCols = c(), old = FALSE, ...){
	
	#don't need all the output from this, but it includes most things we need
	allInput <- prepStrata(multStratData, tags, GSIcol, PBTcol, strataCol, variableCols = variableCols, variableColsOth = c(), "AdClip",
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
				pbtGSIkey[[i]] <- which(input$ohnc_gsi[i,] > 0)[1]
				nParam <- nParam + length(pbtGSIkey[[i]]) - 1
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
		start_piGSIMAT <- ohnc_gsi #just use ohnc assignments
		for(i in 1:nrow(start_piGSIMAT)){
			temp <- start_piGSIMAT[i,]
			temp[temp < 1] <- .1
			temp <- temp / sum(temp)
			temp <- temp / temp[pbtGSIkey[[i]]]
			pos <- 1:length(temp)
			pos <- pos[pos != pbtGSIkey[[i]]]
			temp <- temp[pos]
			start_piGSI <- c(start_piGSI, temp)
		}
		

		
		# return(list(nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey,
		# 				utVar = utVar, ohnc_var = ohnc_var, varKey = varKey, nCat = nCat,
		# 				starting = c(start_piTot, start_piGSI, start_piVar)))
		
		# find mle
			print("start")
			print(flex_negllh_allGSI(c(start_piTot, start_piGSI), nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, 
									t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey))
			tempFit <- optim(c(start_piTot, start_piGSI), flex_negllh_allGSI, method = optimMethod, ...,
								 #arguments to pass to flex_negllh_allGSI
								 nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey)

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
		
		estimates[[input$strataName]] <- list(piTot = ptestim, piGSI = piGSItemp, par = tempFit$par, strataName = input$strataName)
	}
	
	return(estimates)
}
