# testing mle


source("./flex_negllh.R")
source("./MLEwrapper.R")

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)

ptMLE <- matrix(nrow = 1000, ncol = 6)

for(r in 1:1000){
	if(r %% 100 == 0) print(r) 
	multStratData <- data.frame()
	tempDataAll <- generatePBTGSIdata(sampRate = .1, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)
	tempData <- tempDataAll[[1]]
	tempData$StrataVar <- 1
	multStratData <- rbind(multStratData, tempData)
	
	multStratData$GSI <- paste0("GSIgroup", multStratData$GSI)
	tags <- tempDataAll[[2]]
	
	fit <- tryCatch(MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", 5000),
						 error = function(e){
						 	cat("\nNelder-Mead\n")
						 	return(MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "Nelder-Mead", 5000))
						 }
			)
	
	ptMLE[r,] <- fit[[1]]$piTot
	
}


truePiTot
boxplot(t(t(ptMLE) / truePiTot))
abline(h=1)
