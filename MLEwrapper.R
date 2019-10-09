#' wrapper to process data and determine mle estimates

MLEwrapper <- function(trapData, tags, GSIcol, PBTcol, strataCol, optimMethod = "BFGS",
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
		nVar <- length(input$values)
		
		#pull variables values out of input
		nCat <- list()
		ohnc_var <- list()
		utVar <- list()
		if(nVar > 0){
			for(i in 1:nVar){
				nCat[[i]] <- length(input$values[[i]])
				if(nPBT > 0) ohnc_var[[i]] <- input$pi_Vohnc[[i]][1:nPBT,]
				tempUtVar <- c()
				for(v in input$values[[i]]){
					tempUtVar <- c(tempUtVar, sum(input$v_ut[,i] == v))
				}
				utVar[[i]] <- tempUtVar
			}
		}
		
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

		# create varKey
		varKey <- list()
		if(nVar > 0){
			for(v in 1:nVar){
				tempKey <- list()
				if(nPBT > 0){
					for(i in 1:nPBT){
						tempKey[[i]] <- which(input$pi_Vohnc[[v]][i,] > 0)
						nParam <- nParam + length(tempKey[[i]]) - 1
					}
				}
				tempValues <- input$values[[v]]
				for(i in 1:nGSI){
					#determine which categories have un-tagged fish with this GSI assignment
					tempCounts <- c()
					for(val in tempValues){
						tempCounts <- c(tempCounts, sum(input$gsiUT == input$GSI_values[i] & input$v_ut[,v] == val))
					}
					tempKey[[nPBT + i]] <- which(tempCounts > 0)
					nParam <- nParam + length(tempKey[[nPBT + i]]) - 1
				}
				varKey[[v]] <- tempKey
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
			temp <- start_piGSIMAT[i,] / sum(start_piGSIMAT[i,])
			temp <- temp[temp > 0]
			temp <- temp / temp[1]
			if(length(temp) > 1) start_piGSI <- c(start_piGSI, temp[2:length(temp)])
		}
		
		#for piVar
		start_piVar <- c()
		if(nVar > 0){
			for(v in 1:nVar){
				start_piVarMAT <- ohnc_var[[v]] #just use ohnc assignments
				for(i in 1:nrow(start_piVarMAT)){
					temp <- start_piVarMAT[i,] / sum(start_piVarMAT[i,])
					temp <- temp[temp > 0]
					temp <- temp / temp[1]
					if(length(temp) > 1) start_piVar <- c(start_piVar, temp[2:length(temp)])
				}
				# then ut fish
				tempValues <- input$values[[v]]
				for(i in 1:nGSI){
					#determine which categories have un-tagged fish with this GSI assignment
					temp <- c()
					for(val in tempValues){
						temp <- c(temp, sum(input$gsiUT == input$GSI_values[i] & input$v_ut[,v] == val))
					}
					temp <- temp / sum(temp)
					temp <- temp[temp > 0]
					temp <- temp / temp[1]
					if(length(temp) > 1) start_piVar <- c(start_piVar, temp[2:length(temp)])
				}
			}
		}

		
		
		# return(list(nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey,
		# 				utVar = utVar, ohnc_var = ohnc_var, varKey = varKey, nCat = nCat,
		# 				starting = c(start_piTot, start_piGSI, start_piVar)))
		
		# find mle
		if(old){
			print("start")
			print(flex_negllh(c(start_piTot, start_piGSI), nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, 
									t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey))
			tempFit <- optim(c(start_piTot, start_piGSI), flex_negllh, method = optimMethod, ...,
								 #arguments to pass to flex_negllh
								 nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey)
			nVar <- 0
		} else {
						print("start")
			print(flex_negllh_var(c(start_piTot, start_piGSI, start_piVar), nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, 
									t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey,
									utVar = utVar, ohnc_var = ohnc_var, varKey = varKey, nCat = nCat))
					tempFit <- optim(c(start_piTot, start_piGSI, start_piVar), flex_negllh_var, method = optimMethod, ...,
							 #arguments to pass to flex_negllh
							 nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey,
							 utVar = utVar, ohnc_var = ohnc_var, varKey = varKey, nCat = nCat)
		}
		if(tempFit$convergence != 0) cat("\nOptimizer gave convergence code of", tempFit$convergence, "in strata", input$strataName, "\n")
		
		# return(tempFit)
		
					# print(flex_negllh_var(tempFit$par, nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, 
					# 				t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey,
					# 				utVar = utVar, ohnc_var = ohnc_var, varKey = varKey, nCat = nCat))
					# tryPar <- tempFit$par
					# tryPar[1:5] <- c(2,3,3.5,7,3.5)
					# print(flex_negllh_var(tryPar, nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, 
					# 				t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi, pbtGSIkey = pbtGSIkey,
					# 				utVar = utVar, ohnc_var = ohnc_var, varKey = varKey, nCat = nCat))
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
					}
				}
				piVarList[[v]] <- piVarTemp
			}
		}
		
		estimates[[input$strataName]] <- list(piTot = ptestim, piGSI = piGSItemp, piVarList = piVarList, par = tempFit$par, strataName = input$strataName)
	}
	
	return(estimates)
}
