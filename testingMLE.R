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


# compare estiamtes against SD
# they may be identical (or close to it)

reps <- 100

pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)

ptMLE <- matrix(nrow = reps, ncol = 6)
ptSD <- matrix(nrow = reps, ncol = 6)

all_names <- c("pbtGroup1",    "pbtGroup2",    "pbtGroup3",    "GSIgroup1",    "GSIgroup2", 
   "GSIgroup3")
colnames(ptMLE) <- all_names
colnames(ptSD) <- all_names

for(r in 1:reps){
	if(r %% floor(reps/10) == 0) print(r) 
	multStratData <- data.frame()
	tempDataAll <- generatePBTGSIdata(sampRate = .1, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.8, .85,.9), 
										 obsTagRates = c(.8, .85,.9), physTagRates = 0,
				    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = NA)
	tempData <- tempDataAll[[1]]
	tempData$StrataVar <- 1
	multStratData <- rbind(multStratData, tempData)
	
	multStratData$GSI <- paste0("GSIgroup", multStratData$GSI)
	tags <- tempDataAll[[2]]
	
	fit <- tryCatch(MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", control = list(maxit = 5000), old = TRUE),
						 error = function(e){
						 	cat("\nNelder-Mead\n")
						 	return(MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "Nelder-Mead", control = list(maxit = 5000), old = TRUE))
						 }
			)
	# if(length(fit[[1]]$piTot) != 6) next
	
	ptMLE[r,] <- fit[[1]]$piTot[all_names]
	
	## SD
	
	#create window count input
	window <- cbind(1, 3000, 1)

	#run to get PBT group compositions
	SCOBI_deux_fast(adultData = multStratData, windowData = window,
			 Run = "HNC_sim", RTYPE = "noclip_H", Hierarch_variables = c("GenParentHatchery"),
	                  SizeCut = NULL, alph = 0.1, B = 5, writeBoot = F, pbtRates = tags,
			 adClipVariable = "AdClip", physTagsVariable = "PhysTag", pbtGroupVariable = "GenParentHatchery",
			 screenOutput = "tempScreen.txt", dataGroupVariable = "StrataVar")

	#run to get wild group compositions
	SCOBI_deux_fast(adultData = multStratData, windowData = window,
			 Run = "W_sim", RTYPE = "wild", Hierarch_variables = c("GSI"),
	                  SizeCut = NULL, alph = 0.1, B = 5, writeBoot = F, pbtRates = tags,
			 adClipVariable = "AdClip", physTagsVariable = "PhysTag", pbtGroupVariable = "GenParentHatchery",
			 screenOutput = "tempScreen.txt", dataGroupVariable = "StrataVar")

	#record results
	pbt_res <- read.table("./HNC_sim_CI_Hier_GenParentHatchery.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	ptSD[r,1:3] <- pbt_res[match(all_names[1:3], pbt_res[,1]), 2]
	
	wild_res <- read.table("./W_sim_CI_Hier_GSI.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	ptSD[r,4:6] <- wild_res[match(all_names[4:6], wild_res[,1]), 2]
	
	ptSD[r,] <- ptSD[r,] / sum(window[,2])

}


ptMLE
ptSD

head(ptMLE)
head(ptSD)

all.equal(ptMLE, ptSD)
# [1] "Mean relative difference: 2.191572e-05"

truePiTot
boxplot(cbind(t(t(ptMLE) / truePiTot), t(t(ptSD) / truePiTot)))
abline(h=1)

par(mfrow=c(2,3))
for(i in 1:6) {
	plot(ptMLE[,i], ptSD[,i])
	abline(0,1)
}
# the estimates are essentially the same

system.time(
	replicate(10, MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "BFGS", control = list(maxit = 5000), old = TRUE))
)


system.time(
	replicate(10, {
			#run to get PBT group compositions
	SCOBI_deux_fast(adultData = multStratData, windowData = window,
			 Run = "HNC_sim", RTYPE = "noclip_H", Hierarch_variables = c("GenParentHatchery"),
	                  SizeCut = NULL, alph = 0.1, B = 0, writeBoot = F, pbtRates = tags,
			 adClipVariable = "AdClip", physTagsVariable = "PhysTag", pbtGroupVariable = "GenParentHatchery",
			 screenOutput = "tempScreen.txt", dataGroupVariable = "StrataVar")

	#run to get wild group compositions
	SCOBI_deux_fast(adultData = multStratData, windowData = window,
			 Run = "W_sim", RTYPE = "wild", Hierarch_variables = c("GSI"),
	                  SizeCut = NULL, alph = 0.1, B = 0, writeBoot = F, pbtRates = tags,
			 adClipVariable = "AdClip", physTagsVariable = "PhysTag", pbtGroupVariable = "GenParentHatchery",
			 screenOutput = "tempScreen.txt", dataGroupVariable = "StrataVar")
	})
)

# MLE might be a little faster

#################################################################################################################
### looking at all vs set unobserved to zero with very low tag rates
###### note that after this section was run, MLEwrapper was changed to always estimate all GSI categories

varMat <- list(
	matrix(c(rep(c(.1,.9), 3), rep(c(.9,.1), 3)), nrow = 6, ncol = 2, byrow = TRUE),
	matrix(c(rep(c(.4,.6), 3), rep(c(.6,.4), 3)), nrow = 6, ncol = 2, byrow = TRUE)
)
pbtGSImat <- matrix(c(.1, .8, .1,.8, .1, .1,.1, .1, .8), nrow = 3, ncol = 3, byrow = TRUE)
# pbtGSImat <- matrix(1/3, nrow = 3, ncol = 3)

reps <- 300

ptZero <- matrix(0, nrow = reps, ncol = 6)
ptAll <- matrix(0, nrow = reps, ncol = 6)

all_names <- c("pbtGroup1",    "pbtGroup2",    "pbtGroup3",    "GSIgroup1",    "GSIgroup2", 
   "GSIgroup3")
colnames(ptZero) <- all_names
colnames(ptAll) <- all_names


for(r in 1:reps){
	
	multStratData <- data.frame()
	tempDataAll <- generatePBTGSIdata(sampRate = .1, censusSize = 3000, relSizePBTgroups = c(1,2,3), tagRates = c(.08, .085,.09), 
								 obsTagRates = c(.08, .085,.09), physTagRates = 0,
		    true_clipped = 0, true_noclip_H = .3, true_wild = .7, relSizeGSIgroups = c(1,2,1), PBT_GSI_calls = pbtGSImat, varMatList = varMat)
	tempData <- tempDataAll[[1]]
	tempData$StrataVar <- 1
	multStratData <- rbind(multStratData, tempData)
	
	multStratData$GSI <- paste0("GSIgroup", multStratData$GSI)
	tags <- tempDataAll[[2]]
	if(sum(all_names[1:3] %in% multStratData$GenParentHatchery) != 3) next
	
	obs <- tryCatch(MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", adFinCol = "AdClip", AI = TRUE, 
			  optimMethod = "BFGS", variableCols = c(), control = list(maxit = 5000))[[1]]$piTot,
		error = function(e){
			MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", adFinCol = "AdClip", AI = TRUE, 
			  optimMethod = "Nelder-Mead", variableCols = c(), control = list(maxit = 5000))[[1]]$piTot
			}
		)
	all <- MLEwrapper_allGSI(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "Nelder-Mead", 
							variableCols = c(), control = list(maxit = 50000))[[1]]$piTot
	ptZero[r,names(obs)] <- obs
	ptAll[r,names(all)] <- all
	
}



ptZero <- ptZero[rowSums(ptZero) > 0,]
ptAll <- ptAll[rowSums(ptAll) > 0,]

boxplot(cbind(t(t(ptZero) / truePiTot), t(t(ptAll) / truePiTot)))
abline(h=1)

par(mfrow=c(2,3))
for(i in 1:6) {
	plot(ptZero[,i], ptAll[,i])
	abline(0,1)
}

par(mfrow=c(1,1))


compMatZero <- t(t(ptZero) / truePiTot)
compMatAll <- t(t(ptAll) / truePiTot)

apply(compMatZero,2, function(x) {
	x <- abs(1-x)
	mean(x)
	})
apply(compMatAll,2, function(x) {
	x <- abs(1-x)
	mean(x)
	})

apply(compMatZero,2, function(x) {
	x <- (1-x)^2
	mean(x)
	})
apply(compMatAll,2, function(x) {
	x <- (1-x)^2
	mean(x)
	})

MLEwrapper(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", adFinCol = "AdClip", AI = TRUE, 
			  optimMethod = "BFGS", variableCols = c(), control = list(maxit = 5000))[[1]]$piTot

MLEwrapper_allGSI(multStratData, tags, "GSI", "GenParentHatchery", "StrataVar", "Nelder-Mead", 
							variableCols = c(), control = list(maxit = 50000))[[1]]$piTot

#################################################################################################################


