lmm = function(z, dInfo, vInfo, mInfo)
{
    #####################################################################################################################
    #Make a linear mixed model.                                                                                         #
    #Arguments:                                                                                                         #
    # - "z" = dataframe of "nr" rows that contains the data. Missing data must be coded with NA.                        #
    # - "dInfo" = list of specifications for the dealing with the data.                                                 #
    # - "vInfo" = list of specifications for making and modelling the variogram.                                        #
    # - "mInfo" = list of specifications for the minimisation algorithm used to model the variogram.                    #
    #                                                                                                                   #
    #The possible components of list "dInfo" are:                                                                       #
    # - "y" = the column names of the response variables of "z".                                                        #
    # - "coords" = the columns names of the spatial coordinates of "z".                                                 #
    # - "scale" = logical to scale each response variable by it standard deviation.                                     #
    # - "trend" = a list with sub-components:                                                                           #
    #               - for each variable in "y", a string representing the formula needed; OR,                           #
    #               - for each variable in "y", an object of class "rpart"; OR,                                         #
    #               - "X" (a custom design matrix), and "Xv" (a corresponding matrix that indexes the variables in "X". #
    #               - "nesting" is the column numbers of the design matrix to zero, if nesting.                         #
    #                                                                                                                   #
    #The possible components of list "vInfo" are:                                                                       #
    # - "name" = string that names the covariance function to fit.                                                      #
    # - "effrMin" = numeric for the smallest value of the effective range that can be taken by the covariance function. #
    # - "effrMax" = numeric for the largest value of the effective range that can be taken by the covariance function.  #
    # - "phiLock" = logical to lock distance parameter "phi".                                                           #
    # - "kappaLock" = logical to lock curvature parameter "kappa".                                                      #
    # - "cvci" = logical to bootstrap a confidence interval about the theta statistics.                                 #
    # - "plotIt" = plot the experimental variogram(s) and the initial guess on the covariance function. Can be either a #
    #              logical (to save the plot into the default directory) or a string specifying the desired directory.  #
    # - "useSparse" = use sparse matrices to speed up computation with big datasets.                                    #
    #                                                                                                                   #
    #The possible components of list "mInfo" are:                                                                       #
    # - "method" = string to indicate the minimisation method ("lbfgsb", "sanneal" or "simplex").                       #
    # - ???                                                                                                             #
    # - "ncore" = numeric for number of cores to use in the analysis.                                                   #
    # - "silent" = logical to print output as you go.                                                                   #
    #
    
    ##UP TO HERE!
    
    # - "method" = minimisation method ("lbfgsb" or "sanneal").                                                         #
    # - "nStep" = used with "sanneal"; the number of cooling steps in the minimisation.                                 #
    # - "nMarkov" = used with "sanneal"; the number of Markov chains per cooling step.                                  #
    # - "heat" = used with "sanneal"; the initial temperature of the system.                                            #
    # - "alpha" = used with "sanneal"; the proportion by which the system is cooled at each step.                       #
    # - "nStop" = used with "sanneal"; if parameters aren"t changed, the number cooling steps needed before termination.#
    #                                                                                                                   #
    #For some list components, if it's not specified the script will stop with an error message; for others, it will be #
    #reset to a default value.                                                                                          #
	#####################################################################################################################
    #The operating system.
	os = Sys.info()["sysname"]
    
	#Bring in the required libraries and function.
	library(geoR)
	library(MASS)
	library(Matrix)
    library(parallel)
    library(rpart)
    #if(os == 'Windows')
    #{
    #    source('U:\\hg\\rscrepo\\rsc\\geostats\\tree\\tree.R')
    #}else{
    #    source('~/hg/rscrepo/rsc/geostats/tree/tree.R')
    #}
	
    #The maximum number of observations in each variable that can be modelled.
    nMax = 2000
    
    #Error check.
    tmp = lmm.errorCheck(z, dInfo, vInfo, mInfo, os)
    dInfo = tmp$dInfo
    vInfo = tmp$vInfo
    mInfo = tmp$mInfo
    nv = tmp$nv
    nesting = tmp$nesting
  print(nesting)
  print(dInfo)
  print("line 77")
    useSparse = tmp$useSparse
    
	#If required, compute the standard deviation of each variable, and use it to scale "z".
    zScl = z
	stdev = rep(NA, nv)
    if(dInfo$scale)
    {
    	for(i in 1:nv)
	    {
		    stdev[i] = sd(z[, dInfo$y[i]], na.rm = TRUE)
		    zScl[[dInfo$y[i]]] = z[[dInfo$y[i]]] / stdev[i]
	    }
    }
    
    #Matrix of correlations and their signs.
	rho = as.matrix(1)
    sgn = as.matrix(1)
	if(nv > 1)
	{
		rho = cor(zScl[, dInfo$y], use = "pairwise.complete.obs")
        sgn = ifelse(rho < 0, -1, 1)
        if(any(is.na(rho)))
        {
            #This will only be invoked when there is no colocations between
            #observations of a pair of variables.
            for(i in 1:(nv - 1))
            {
                for(j in i:nv)
                {
                    if(i != j)
                    {
                        rho[i, j] = 0.1
                        rho[j, i] = rho[i, j]
                    }
                }
            }
            sgn = matrix(0, nv, nv)
            diag(sgn) = 1
        }
	}
	
    #Stack "zScl".
    tmp = lmm.stack(zScl, dInfo, NULL)
    zStk = tmp$z
    nr = tmp$nr
    if(all(nr < 10))stop("lmm: There are not enough data to analyse under any circumstances!")
    nrTotal = tmp$nrTotal
    dThresh = NULL
    X = NULL
    Xv = NULL
    if(any(nr > nMax))
    {
        if(!dInfo$useLinear)stop("lmm: 19.2.2015 - I am not convinced subsampling can be done when using rpart!")
        tmp = lmm.subsample(nMax, nv, nr, zStk, dInfo$trend)
        dThresh = tmp$dThresh
        zStk = tmp$z
        nr = tmp$nr
        nrTotal = tmp$nrTotal
        if(any(names(dInfo$trend) == "X"))
        {
            X = tmp$X
            Xv = tmp$Xv
        }
    }
    
	#Distances between pairs of locations.
    d = lmm.distanz(os, nrTotal, 2, zStk[, 2:3], vInfo$effrMax, useSparse)
    
    #Make the design matrix for the fixed effects and the initial guesses
    #for the parameters of the random effects.
  print(nesting)
  print("line151")
    tmp = lmm.XandGuesses(zStk, dInfo, vInfo, mInfo, nesting, nv, nr, stdev, rho, X, Xv)
    vInfo = tmp$vInfo
    X = tmp$X
    Xv = tmp$Xv
    treeRules = tmp$treeRules
    
	#Fit. Note the use of multiple cores if using Linux.
	cat("#######################################\n") 
	cat("#FITTING THE LINEAR MIXED MODEL TO 'z'#\n") 
	cat("#######################################\n")
    model = NULL
    if(nv == 1 & mInfo$method == "lbfgsb")
	{
        if(mInfo$ncore == 1)
        {
            model = lmm.lbfgsb(zStk$value, os, nv, nr, nrTotal, X, Xv, d, vInfo, mInfo, stdev, sgn, nesting, useSparse)
        }else{
            tmp = list()
            for(i in 1:mInfo$ncore)
            {
                tmp[[i]] = zStk$value
            }
            modelList = mclapply(tmp, lmm.lbfgsb, os = os, nv = nv, nr = nr, nrTotal = nrTotal, Xmat = X, Xv = Xv,
                d = d, vInfo = vInfo, mInfo = mInfo, stdev = stdev, sgn = sgn,
                nesting = nesting, useSparse = useSparse, mc.cores = mInfo$ncore)
            l = 10 ^ 100
            for(i in 1:mInfo$ncore)
            {
                if(!is.character(modelList[[i]]) && modelList[[i]]$l[length(modelList[[i]]$l)] < l)
                {
                    l = modelList[[i]]$l[length(modelList[[i]]$l)]
                    model = modelList[[i]]
                }
            }
        }
	}else{
        if(mInfo$ncore == 1)
        {
            if(mInfo$method == "sanneal")
            {
                model = lmm.sanneal(zStk$value, os, nv, nr, nrTotal, X, Xv, d, vInfo, mInfo, stdev, sgn, nesting, useSparse)
            }else{
                model = lmm.simplex(zStk$value, os, nv, nr, nrTotal, X, Xv, d, vInfo, mInfo, stdev, sgn, nesting, useSparse)
            }
        }else{
            tmp = list()
            for(i in 1:mInfo$ncore)
            {
                tmp[[i]] = zStk$value
            }
            if(mInfo$method == "sanneal")
            {
                modelList = mclapply(tmp, lmm.sanneal, os = os, nv = nv, nr = nr, nrTotal = nrTotal, Xmat = X, Xv = Xv,
                    d = d, vInfo = vInfo, mInfo = mInfo, stdev = stdev, sgn = sgn,
                    nesting = nesting, useSparse = useSparse, mc.cores = mInfo$ncore)
            }else{
                modelList = mclapply(tmp, lmm.simplex, os = os, nv = nv, nr = nr, nrTotal = nrTotal, Xmat = X, Xv = Xv,
                    d = d, vInfo = vInfo, mInfo = mInfo, stdev = stdev, sgn = sgn,
                    nesting = nesting, useSparse = useSparse, mc.cores = mInfo$ncore)
            }
            l = 10 ^ 100
            for(i in 1:mInfo$ncore)
            {
                if(!is.character(modelList[[i]]) && modelList[[i]]$l[length(modelList[[i]]$l)] < l)
                {
                    l = modelList[[i]]$l[length(modelList[[i]]$l)]
                    model = modelList[[i]]
                }
            }
        }
    }
    cvSummary = NULL
    if(is.logical(nesting) && !is.null(model))
    {
        #Cross-validate.
        cvSummary = lmm.crossValidate(os, nv, nr, nrTotal, zStk, dInfo$scale, stdev, X, d,
            vInfo, model, useSparse, mInfo$ncore)
	}
    if(dInfo$scale)
    {
        #Rescale "sig" and "b" for output.
        #We need to stack the coregionalisation matrices first.
        sig = array(0, c(nv, nv, 3))
        sig[, , 1] = model$c0
        if(vInfo$nstr > 0)
        {
            sig[, , 2] = model$c1
            if(vInfo$nstr == 2)sig[, , 3] = model$c2
        }
        tmp = lmm.rescale(nv, stdev, Xv, sig, model$b, model$C)
        model$c0 = tmp$sig[, , 1]
        model$c1 = tmp$sig[, , 2]
        model$c2 = tmp$sig[, , 3]
	    model$b = tmp$b
        model$C = tmp$C
    }
    
    #Convert to log-likelihood.
    model$l = model$l * -1
    
    #If scaling has been used, or a subset been taken,
    #we want to return the original data.
    if(dInfo$scale | !is.null(dThresh))
    {
        #Stack "z".
        tmp = lmm.stack(z, dInfo, NULL)
	    zStk = tmp$z
        nr = tmp$nr
        if(!is.null(dThresh))
        {
            #Remake "d".
            if(!is.null(dThresh))d = lmm.distanz(os, sum(nr), 2, zStk[, 2:3], vInfo$effrMax, useSparse)
            
            #Remake "X", if necessary.
            if(any(names(dInfo$trend) == "X"))
            {
                X = dInfo$trend$X
            }else{
                tmp = lmm.XandGuesses(zStk, dInfo, vInfo, mInfo, nesting, nv, nr, stdev, NULL, NULL)
                X = tmp$X
            }
        }
    }
    
    #Output list.
    out = list("model" = model, "cvSummary" = cvSummary, "zStk" = zStk, "scale" = dInfo$scale,
        "nr" = nr, "d" = d, "treeRules" = NULL, "X" = X, "subsettingThreshold" = dThresh)
    if(length(names(treeRules)) > 0)out$treeRules = treeRules
	out
}


lmm.check = function(os, nv, name, nstr, phi, kappa, effrMin, effrMax, c0, c1, c2, sgn, doThis)
{
	############################################################
	#Check various things related to the covariance parameters.#
	############################################################
    #Deal with any zero correlations between variables.
    newSgn = sgn
    if(nv > 1 && any(sgn == 0))
    {
        ndx = which(sgn == 0, arr.ind = TRUE)
        for(i in 1:nrow(ndx))
        {
            if(ndx[i, 1] < ndx[i, 2] && sgn[ndx[i, 1], ndx[i, 2]] == 0)
            {
                newSgn[ndx[i, 1], ndx[i, 2]] = 1
                if(runif(1) < 0.5)newSgn[ndx[i, 1], ndx[i, 2]] = -1
                newSgn[ndx[i, 2], ndx[i, 1]] = newSgn[ndx[i, 1], ndx[i, 2]]
            }
        }
    }
    
    if(os == "Windows")
    {
    	#A very small number.
	    eps = 10^-300
        
	    #Assume that all parameters are valid to begin with.
	    out = TRUE
	    ##if(doThis == 1 | doThis == 4)
	    ##{
	    ##	stop("check boxcox not written yet")
	    ##	The Box-Cox parameters.
	    ##	do i=1,nv
	    ##		if(boxcox(i)<=-2.or.boxcox(i)>=2)then
	    ##			check=.false.
	    ##			return
	    ##		endif
	    ##	enddo
	    ##}
        if((doThis == 2 | doThis == 4) & name != "pureNugget")
	    {
		    #The curvature parameters.
            if((name == "dampedPeriodic" & phi > eps & kappa > (-1.569795) & kappa < (-0.001)) |
			    (name == "doubleSpherical" & phi[1] > eps & phi[2] > phi[1]) |
                ((name == "exponential" | name == "spherical") & phi > eps) |
			    (name == "matern" & phi > eps & kappa > 0.01 & kappa <= 5))
		    {
			    effr = lmm.effectiveRange(name, phi, kappa, effrMax)
                if(effr < effrMin | effr > effrMax)out = FALSE
		    }else{
			    out = FALSE
		    }
		    if(!out)return(out)
	    }
        if(doThis == 3 | doThis == 4)
	    {
            #The variance parameters.
		    for(i in 1:nv)
		    {
			    for(j in i:nv)
			    {
				    if(i == j)
				    {
                        #Does each structure have finite variance?
                        if(c0[i, j] <= 0)out = FALSE
                        if(nstr >= 1 & c1[i, j] <= 0)out = FALSE
                        if(nstr > 1 & c2[i, j] <= 0)out = FALSE
                        
                        #Does the variable have finite variance overall?
                        if((c0[i, j] + c1[i, j] + c2[i, j]) <= eps)out = FALSE
				    }else{
                        #Does the cross-covariance of each structure go the right way?
                        if((newSgn[i, j] < 0 & c0[i, j] >= -eps) | (newSgn[i, j] > 0 & c0[i, j] <= eps))out = FALSE
                        if(nstr > 0)
                        {
                            if((newSgn[i, j] < 0 & c1[i, j] >= -eps) | (newSgn[i, j] > 0 & c1[i, j] <= eps))out = FALSE
                            if(nstr == 2)
                            {
                                if((newSgn[i, j] < 0 & c2[i, j] >= -eps) | (newSgn[i, j] > 0 & c2[i, j] <= eps))out = FALSE
                            }
                        }
				    }
				    if(!out)return(out)
			    }
		    }
            if(nv > 1)
		    {
                #Is each coregionalisation matrix positive-definite?
                V = 0
                for(i in seq(0, nstr))
                {
                    if(i == 0)V = c0
                    if(i == 1)V = c1
                    if(i == 2)V = c2
                    Vchol = try(chol(V), silent = TRUE)
                    if(class(Vchol) == 'try-error')
                    {
                        out = FALSE
                        return(out)
                    }
                    
                    #Are the structural correlations valid?
                    for(j in c(1, (nv - 1)))
                    {
                        for(k in c((j + 1), nv))
                        {
                            scorr = V[j, k] / sqrt(V[j, j] * V[k, k])
                            if(abs(scorr) >= 0.999)
                            {
                                out = FALSE
                                return(out)
                            }
                        }
                    }
                }
		    }
	    }
    }else{
        #Use FORTRAN to do it quicker.
        dyn.load("~/hg/rscrepo/rsc/geostats/lmm/check.so")
        if(name == "dampedPeriodic")name = 1
        if(name == "doubleSpherical")name = 2
        if(name == "exponential")name = 3
        if(name == "matern")name = 4
        if(name == "pureNugget")name = 5
        if(name == "spherical")name = 6
        if(name != "doubleSpherical")phi = c(phi, 0)
        out = TRUE #Assume that all parameters are valid to begin with.
        tmp = .Fortran("check", as.integer(nv), as.integer(name), as.integer(nstr),
            as.double(phi), as.double(kappa),
            as.double(effrMin), as.double(effrMax),
            as.double(c0), as.double(c1), as.double(c2),
            as.integer(newSgn), as.integer(doThis), out)
        out = tmp[[13]]
    }
    out
}


lmm.covFun = function(d, name, phi, kappa, c0, c1, c2)
{
	######################
	#Covariance function.#
	######################
    sill = c0 + c1 + c2
	if(name == "dampedPeriodic")
	{
		tmp = (2 * pi * d / phi) - kappa
		out = ifelse(d == 0, sill, c1 * (sin(tmp) / tmp))
	}
    if(name == "doubleSpherical")
	{
		dOverPhi1 = d / phi[1]
        dOverPhi2 = d / phi[2]
        out = ifelse(d == 0, sill, 
            ifelse(d > 0 & d < phi[1], c1 * (1 - 1.5 * dOverPhi1 + 0.5 * (dOverPhi1 ^ 3)) + c2 * (1 - 1.5 * dOverPhi2 + 0.5 * (dOverPhi2 ^ 3)),
            ifelse(d >= phi[1] & d < phi[2], c2 * (1 - 1.5 * dOverPhi2 + 0.5 * (dOverPhi2 ^ 3)),
            0)))
	}
	if(name == "exponential")out = ifelse(d == 0, sill, c1 * exp(-d / phi))
	if(name == "matern")
	{
		ratio = d / phi
		out = ifelse(d == 0, sill, ifelse(ratio > 500, 0, 
			c1 * (1 / ((2 ^ (kappa - 1)) * gamma(kappa))) * (ratio ^ kappa) * besselK(ratio, kappa)))
	}
	if(name == "pureNugget")out = ifelse(d == 0, sill, 0)
	if(name == "spherical")
	{
		dOverPhi = d / phi
        out = ifelse(d == 0, sill, 
            ifelse(d > 0 & d < phi, c1 * (1 - 1.5 * dOverPhi + 0.5 * (dOverPhi ^ 3)),
            0))
	}
    out
}


lmm.crossValidate = function(os, nv, nr, nrTotal, z, scaling, stdev, X, d, vInfo, model, useSparse, ncore)
{
	#################################################
	#Cross-validate the random effects of the model.#
	#################################################
    #Empty output list.
    out = list()
    
    #Make covariance matrix of the random effects, "V".
    tmp = lmm.makeV(os, nv, nr, nrTotal, d, model$covFun, vInfo$nstr, model$phi, model$kappa,
        model$c0, model$c1, model$c2, useSparse, FALSE)
	if(is.null(tmp$V))
    {
        #"V" was not positive definite. Return "NULL".
        out = NULL
        return(out)
    }else{
        #"V" was positive definite, so unpack the list.
        V = tmp$V
        Vchol = tmp$Vchol
    }
    rm(tmp)
    
    #Invert "V", ensuring it's symmetric.
    iV = forceSymmetric(chol2inv(Vchol), "U")
    
    #Transpose "Vchol" - we want the lower triangle.
    Vchol = t(Vchol)
    
    #Setup the EBLUP list.
    eblupList = list()
    for(j in 1:nrTotal)
    {
        eblupList[[j]] = list("j" = j)
    }
    
    #Loop over the number of cross-validations to do.
    n2do = 1
    if(vInfo$cvci)
	{
		##n2do = 501
        n2do = 11 #FOR TESTING.
	}
    thetaMean = matrix(0, nrow = nv, ncol = n2do)
	thetaMedian = matrix(0, nrow = nv, ncol = n2do)
    for(i in 1:n2do)
    {
        if(i == 1)
        {
            #Use the actual data.
            zi = as.matrix(z$value)
            cat("Cross-validating the original data\n")
        }else{
            #Simulate a Gaussian random field based on the model.
            if(!useSparse)
            {
                zi = X %*% mvrnorm(n = 1, model$b, model$C) #Simulate some fixed effects.
                zi = zi + mvrnorm(n = 1, rep(0, nrTotal), V) #Add some simulated random effects.
            }else{
                #You first need to make the Cholesky decomposition of the covariance
                #matrix of the fixed effects, "C".
                if(i == 2)Cchol = chol(forceSymmetric(model$C, "U"))
                zi = X %*% (model$b + as.matrix(Cchol %*% as.matrix(rnorm(nrow(model$b))))) #Simulate some fixed effects.
                zi = zi + as.matrix(Vchol %*% as.matrix(rnorm(nrTotal))) #Add some simulated random effects.
            }
            cat("Cross-validating", i - 1, "of", n2do - 1, "Gaussian random fields\n")
        }
        if(ncore == 1)
        {
            #Do it sequentially.
            tmp = lapply(eblupList, lmm.eblupCV, i = i, nv = nv, nr = nr, nrTotal = nrTotal,
                z = zi, Xmat = X, d = d, model = model, V = V, iV = iV, Vchol = Vchol,
                useSparse = useSparse)
        }else{
            #Do it in parallel.
            tmp = mclapply(eblupList, lmm.eblupCV, i = i, nv = nv, nr = nr, nrTotal = nrTotal,
                z = zi, Xmat = X, model = model, d = d, V = V, iV = iV, Vchol = Vchol,
                useSparse = useSparse, mc.cores = ncore)
        }
        
        #Stitch the output back together.
        pred = Matrix(0, nrow = nrTotal, ncol = 2)
        for(j in 1:nrTotal)
        {
            ndx = which(pred[, 1] == 0 & !is.na(tmp[[j]][, 1]), arr.ind = TRUE)
            if(length(ndx) > 0)pred[ndx, ] = tmp[[j]][ndx, ]
        }
        if(i == 1)
        {
            out[["cvPred"]] = data.frame(cbind(z$value, as.matrix(pred)))
            names(out[["cvPred"]]) = c("observed", "prediction", "predictionVariance")
            if(scaling)
            {
                #Rescale the values.
                kount = rep(0, 2)
                for(j in 1:nv)
                {
                    kount[1] = kount[2] + 1
                    kount[2] = kount[2] + nr[j]
                    out[["cvPred"]][kount[1]:kount[2], 1:2] = out[["cvPred"]][kount[1]:kount[2], 1:2] * stdev[j]
                    out[["cvPred"]][kount[1]:kount[2], 3] = out[["cvPred"]][kount[1]:kount[2], 3] * stdev[j] * stdev[j]
                }
            }
        }
        
        #Standardised squared prediction error for each variable.
        theta = ((zi - pred[, 1]) ^ 2) / pred[, 2]
        kount = rep(0, 2)
	    for(j in 1:nv)
	    {
	        kount[1] = kount[2] + 1
	        kount[2] = kount[2] + nr[j]
            thetaMean[j, i] = mean(theta[kount[1]:kount[2]])
            thetaMedian[j, i] = median(theta[kount[1]:kount[2]])
        }
    }
    
	#Make the output.
	for(i in c("thetaMean", "thetaMedian"))
    {
        out[[i]] = data.frame(matrix(0, nrow = nv, ncol = 3))
        names(out[[i]]) = c("value", "q2.5", "q97.5")
        for(j in 1:nv)
	    {
            if(i == "thetaMean")
            {
                out[[i]]$value[j] = thetaMean[j, 1]
                if(vInfo$cvci)tmp = quantile(thetaMean[j, -1], c(0.025, 0.975))
            }else{
                out[[i]]$value[j] = thetaMedian[j, 1]
                if(vInfo$cvci)tmp = quantile(thetaMedian[j, -1], c(0.025, 0.975))
            }
            if(vInfo$cvci)
            {
    		    out[[i]]$"q2.5"[j] = tmp[1]
	    	    out[[i]]$"q97.5"[j] = tmp[2]
            }
	    }
    }
    out
}


lmm.distanz = function(os, nrTotal, ndim, xy, effrMax, useSparse)
{
	##################
	#Distance matrix.#
	##################
    if(!useSparse)
	{
		#Return the full matrix.
		if(nrTotal > 10000)stop("lmm.distanz: There are too many observations; either take a random subset or use sparse matrices!")
		out = matrix(0, nrow = nrTotal, ncol = nrTotal)
		for(i in 1:nrTotal)
		{
            if(ndim == 1)out[i, ] = abs(xy[, 1] - xy[i, 1])
            if(ndim == 2)out[i, ] = sqrt((xy[, 1] - xy[i, 1]) ^ 2 + (xy[, 2] - xy[i, 2]) ^ 2)
			if(ndim == 3)out[i, ] = sqrt((xy[, 1] - xy[i, 1]) ^ 2 + (xy[, 2] - xy[i, 2]) ^ 2 + (xy[, 3] - xy[i, 3]) ^ 2)
		}
	}else{
if(ndim != 2)stop('lmm.distanz - I think this only written for 2 dimensional data and needs to be checked!')
        #Return the lower triangle, in "coordinate"-sparse format.
        if((os == "Windows" & nrTotal > 10000) | (os == "Linux" & nrTotal > 25000))stop("lmm.distanz: There are too many observations, even for sparse matrices!")
		n = nrTotal + nrTotal * (nrTotal - 1) / 2
        if(os == "Windows")
        {
            out = data.frame(matrix(0, nrow = n, ncol = 3))
		    names(out) = c("row", "col", "d")
            kount = rep(0, 4)
    		for(i in 1:nrTotal)
	    	{
		    	kount[1] = kount[2] + 1
			    kount[2] = kount[1] + i - 1
                ndx = list("n" = seq(kount[1], kount[2]))
                ndx$row = rep(i, length(ndx$n))
                ndx$col = seq(length(ndx$n))
                d = sqrt((xy[ndx$row, 1] - xy[i, 1]) ^ 2 + (xy[ndx$col, 2] - xy[i, 2]) ^ 2)
			    tmp = d<= effrMax
                tmp = cbind(ndx$row[tmp], ndx$col[tmp], d[tmp])
                kount[3] = kount[4] + 1
                kount[4] = kount[3] + nrow(tmp) - 1
                ndx = kount[3]:kount[4]
                out$row[ndx] = tmp[, 1]
			    out$col[ndx] = tmp[, 2]
                out$d[ndx] = tmp[, 3]
		    }
            out = out[1:kount[4], ]
        }else{
            #Use FORTRAN to do it quicker.
            dyn.load("~/hg/rscrepo/rsc/geostats/lmm/distanz.so")
            out = matrix(0, nrow = n, ncol = 3)
            tmp = .Fortran("distanz", as.integer(nrTotal), as.integer(ndim), as.double(as.matrix(xy)), as.double(effrMax), as.integer(n), as.double(out))
            out = matrix(tmp[[6]], nrow = n, ncol = 3)
            out = out[1:tmp[[5]], ]
            if(tmp[[5]] == 1)out = t(as.matrix(out))
            out = data.frame(out)
            names(out) = c("row", "col", "d")
        }
	}
	out
}


lmm.eblup = function(os, nv, nr, nrTotal, ndim, z, model, X, n0, xy0, X0, useSparse, iVoo)
{
    ##############################################
    #EBLUP for prediction at unsampled locations.#
    ##############################################
    #Load the required library.
    library(Matrix)
    
    #Structure information.
    if(nv == 1)
    {
        model$c0 = as.matrix(model$c0)
        if(any(names(model) == "c1"))model$c1 = as.matrix(model$c1)
        if(any(names(model) == "c2"))model$c2 = as.matrix(model$c2)
    }
    model$phi = c(model$phi)
    structInfo = list("c1" = matrix(0, nv, nv), "c2" = matrix(0, nv, nv), "nstr" = 0)
    if(model$covFun != "pureNugget")
    {
        structInfo$c1 = model$c1
        structInfo$nstr = 1
        if(model$covFun == "doubleSpherical")
        {
            structInfo$c2 = model$c2
            structInfo$nstr = 2
        }
    }
    
    #Inverse of the covariance matrix between the observed locations, if required.
    if(is.null(iVoo))
    {
        #The covariance matrix between the observed locations.
        d = lmm.distanz(os, nrTotal, ndim, z[, 2:(2 + ndim - 1)], max(model$phi) + 1, useSparse)
        tmp = lmm.makeV(os, nv, nr, nrTotal, d, model$covFun, structInfo$nstr, model$phi, model$kappa,
            model$c0, structInfo$c1, structInfo$c2, useSparse, FALSE)
        if(is.null(tmp$V))
        {
            #"V" was not positive definite. Return an infinite likelihood.
            stop("lmm.eblup: Voo is not positive definite!")
        }else{
            #"V" was positive definite, so unpack the list.
            Vchol = tmp$Vchol
        }
        iVoo = chol2inv(tmp$Vchol)
        rm(tmp)
    }
    
    #Covariance matrix between the prediction locations and the observation locations.
    Vpo = Matrix(0, nrow = n0, ncol = nrTotal, sparse = TRUE)
    if(os == "Windows")
    {
        kount = rep(0, 5)
        for(i in 1:nv)
        {
            kount[1] = kount[2] + 1
            kount[2] = kount[2] + n0 / nv
            kount[3:5] = 0
            for(j in 1:nv)
            {
                kount[3] = kount[4] + 1
                kount[4] = kount[4] + nr[j]
                kount[5] = 0
                d = matrix(0, nrow = n0 / nv, ncol = nr[j])
                for(k in kount[3]:kount[4])
                {
                    kount[5] = kount[5] + 1
                    if(ndim == 1)d[, kount[5]] = abs(xy0[kount[1]:kount[2], 1] - z[k, 2])
                    if(ndim == 2)
                    {
                        d[, kount[5]] = sqrt((xy0[kount[1]:kount[2], 1] - z[k, 2]) ^ 2 + (xy0[kount[1]:kount[2], 2] - z[k, 3]) ^ 2)
                    }
                    if(ndim == 3)
                    {
                        d[, kount[5]] = sqrt((xy0[kount[1]:kount[2], 1] - z[k, 2]) ^ 2 + (xy0[kount[1]:kount[2], 2] - z[k, 3]) ^ 2 + (xy0[kount[1]:kount[2], 3] - z[k, 4]) ^ 2)
                    }
                }
                Vpo[kount[1]:kount[2], kount[3]:kount[4]] = lmm.covFun(d, model$covFun, model$phi, model$kappa, model$c0[i, j], structInfo$c1[i, j], structInfo$c2[i, j])
            }
        }
    }else{
if(ndim != 2)stop('lmm.eblup - I think this only written for 2 dimensional data and needs to be checked!')
        ###Use FORTRAN to do it quicker.
        ##phi = rep(0, 2)
        ##phi[1:structInfo$nstr] = model$phi
        ##dyn.load("~/hg/rscrepo/rsc/geostats/lmm/makeVpo.so")
        ##tmp = .Fortran("makeVpo", as.integer(nv), as.integer(nr), as.integer(nrTotal), as.integer(max(nr)),
        ##    as.integer(n0), as.integer(n0 / nv), as.double(as.matrix(z[, 2:3])), as.double(as.matrix(xy0)),
        ##    as.integer(structInfo$nstr), as.double(phi), as.double(model$c0), as.double(structInfo$c1),
        ##    as.double(structInfo$c2), as.double(as.matrix(Vpo)))
        ##Vpo = Matrix(tmp[[14]], nrow = n0, ncol = nrTotal, sparse = TRUE)
    }
    
    #The covariance matrix between the prediction locations.
    d = lmm.distanz(os, n0, ndim, xy0, max(model$phi) + 1, useSparse)
    tmp = lmm.makeV(os, nv, rep(n0 / nv, nv), n0, d, model$covFun, structInfo$nstr, model$phi, model$kappa, model$c0, structInfo$c1, structInfo$c2, useSparse, FALSE)
    if(is.null(tmp$V))
    {
        #"V" was not positive definite. Return an infinite likelihood.
        stop("lmm.eblup: Vpp is not positive definite!")
    }else{
        #"V" was positive definite, so unpack the list.
        Vpp = tmp$V
    }
    rm(d, tmp)
    
    #Put the EBLUP together.
    wgt = Vpo %*% iVoo
    mat1 = X0 - wgt %*% X
    tmp = forceSymmetric(t(X) %*% iVoo %*% X, "L")
    mat2 = chol2inv(chol(tmp))
    mat3 = Vpp  - wgt %*% t(Vpo)
    Cmat = mat1 %*% mat2 %*% t(mat1) + mat3
    Cmat = forceSymmetric(Cmat, "L")
    
    #(Co)kriging estimates.
    pred0 = mat1 %*% model$b + wgt %*% z[, ncol(z)]
    
    #Output list.
    out = list("pred0" = pred0, "Cmat" = Cmat)
    out
}


lmm.eblupCV = function(info, i, nv, nr, nrTotal, z, Xmat, model, d, V, iV, Vchol, useSparse)
{
    #############################
    #EBLUP for cross-validation.#
    #############################
    #Print progress.
    if(info$j %% 20 == 0)cat("j = ",info$j, "of", nrTotal, "\n")
    
    #Empty array for storage of the predictions.
    out = Matrix(0, nrow = nrTotal, ncol = 2)
    
    #Indices for the locations colocated with the "j"th location.
    if(!useSparse)
    {
        ndx = which(d[info$j, ] == 0, arr.ind = TRUE)
    }else{
        ndx = d$row[which(d$col == info$j & d$d == 0, arr.ind = TRUE)]
    }
    n2do = length(ndx)
    
    #Split the variables.
    Xo = Xmat[-ndx, ]
    tXo = t(Xo)
    Xp = Xmat[ndx, ]
    iiVpp = iV[ndx, ndx]
    iVpo = iV[ndx, -ndx]
    Vpo = V[ndx, -ndx]
    Vpp = V[ndx, ndx]
    if(n2do == 1)
    {
        Xp = Matrix(Xp, nrow = 1, ncol = ncol(Xmat))
        iiVpp = Matrix(iiVpp, nrow = 1, ncol = 1)
        iVpo = Matrix(iVpo, nrow = 1, ncol = nrTotal - 1)
        Vpo = Matrix(Vpo, nrow = 1, ncol = nrTotal - 1)
        Vpp = Matrix(Vpp, nrow = 1, ncol = 1)
    }
    tVpo = t(Vpo)
    iVoo = iV[-ndx, -ndx] - t(iVpo) %*% chol2inv(chol(forceSymmetric(iiVpp, "U"))) %*% iVpo
    
    #Simple (co)kriging weights.
    skWgt = Vpo %*% iVoo
    
    #EBLUP prediction and prediction variance.
    mat1 = Xp - skWgt %*% Xo
    mat2 = chol2inv(chol(forceSymmetric(tXo %*% iVoo %*% Xo, "U")))
    mat3 = mat1 %*% mat2
    ukWgt = mat3 %*% tXo %*% iVoo + skWgt
    out[ndx, 1] = ukWgt %*% z[-ndx]
    out[ndx, 2] = diag(mat3 %*% t(mat1) + Vpp - skWgt %*% tVpo)
    out[setdiff(seq(nrTotal), ndx), ] = NA
    out
}

                
lmm.effectiveRange = function(name, phi, kappa, effrMax)
{
	##################################################
	#Effective range of the autocorrelation function.#
	##################################################
	#Parameter.
	kountMax = 500
	
	#For the Matern and exponential functions, take it as 95% of the sill.
	if(name == "dampedPeriodic")
	{
        inc = effrMax / kountMax
		d = seq(inc, kountMax * inc, length.out = kountMax)
		tmp = which(lmm.covFun(d, "dampedPeriodic", phi[1], kappa, 0, 1, 0) < 0, arr.ind = TRUE)[1]
		out = d[tmp] - inc * 0.5 #Effective range is half-way back to the last value.
    }
    if(name == "doubleSpherical")out = phi[2]
    if(name == "exponential")out = 3 * phi[1]
	if(name == "matern")
    {
			inc = effrMax / kountMax
			d = seq(inc, kountMax * inc, length.out = kountMax)
			tmp = which(lmm.covFun(d, "matern", phi, kappa, 0, 1, 0) <= 0.05, arr.ind = TRUE)[1]
			if(is.na(tmp))
			{
				#It's outside the permissible range, so return with an unacceptably large value.
				out = effrMax + 1
				return(out)
			}
			out = d[tmp] - inc * 0.5 #Effective range is half-way back to the last value.
		}
	if(name == "spherical")out = phi[1]
	out
}


lmm.errorCheck = function(z, dInfo, vInfo, mInfo, os)
{
    #########################################
    #Check the input information for errors.#
    #########################################
    #On "dInfo".
    for(i in dInfo$y)
	{
		if(!any(names(z) == i))stop(paste("lmm.errorCheck: Data variable '", i, "' does not exist in 'z'!", sep = ""))
	}
	nv = length(dInfo$y) #The number of variables.
    if(nv > 9)stop("lmm.errorCheck: The maximum number of variables permitted is 9!")
    if(!any(names(dInfo) != "coords"))stop("lmm: Coordinates are not specified!")
    if(!any(names(dInfo) == "scale") || !is.logical(dInfo$scale))dInfo$scale = FALSE
 
  print("line 189, test")
  nesting = FALSE
  print(nesting)
    if(any(names(dInfo) == "trend") && any(names(dInfo$trend) == "nesting"))nesting = dInfo$trend$nesting
   # if(!is.logical(nesting))stop("lmm.errorCheck: 19.2.2015 - check that the nesting even works!")
    if(any(names(dInfo) == "trend") && any(names(dInfo$trend) == "X") && !any(names(dInfo$trend) == "Xv"))
    {
        stop("lmm: If specifying the design matrix, you must include a matrix of corresponding size, which labels each element with a variable number!")
    }
    dInfo$useLinear = TRUE
    if(any(names(dInfo) == "trend"))
    {
        kount = 0
        for(i in dInfo$y)
	    {
            kount = kount + 1
            if(class(dInfo$trend[[i]]) == 'rpart')dInfo$useLinear = FALSE
            if(kount > 1 & !dInfo$useLinear & class(dInfo$trend[[i]]) != 'rpart')stop("lmm.errorCheck: All data variables must use 'rpart'!")
	    }
    }
    if(!is.logical(nesting) & !dInfo$useLinear)stop("lmm.errorCheck: 19.2.2015 - I am not convinced that 'rpart' will work with nesting!")
    
    #On "vInfo".
    if(vInfo$name != "dampedPeriodic" & vInfo$name != "doubleSpherical" & vInfo$name != "exponential" &
        vInfo$name != "matern" & vInfo$name != "pureNugget" & vInfo$name != "spherical")stop("lmm.erroCheck: Invalid covariance function!")
    if(nv > 1 & vInfo$name == "dampedPeriodic")stop("lmm.errorCheck: 'dampedPeriodic' covariance cannot be used with >1 variable at this stage!")
    if(vInfo$name == "dampedPeriodic" & length(dInfo$coords) > 1)stop("lmm.errorCheck: 'dampedPeriodic' covariance should probably only be used with 1-dimensional coordinates!")
    vInfo$nstr = 1
    if(vInfo$name == "dampedPeriodic" | vInfo$name == "matern")
	{
        if(!any(names(vInfo) == "kappaLock"))vInfo$kappaLock = FALSE
        if(!any(names(vInfo) == "kappa"))
        {
            if(vInfo$kappaLock)stop("lmm.errorCheck: The optimum value for 'kappa' must be specified if locking it!")
            vInfo$kappa = rep(-pi / 4, nv)
            if(vInfo$name == "matern")vInfo$kappa = rep(0.5, nv)
        }
	}else{
        if(vInfo$name == "doubleSpherical")vInfo$nstr = 2
		if(vInfo$name ==  "pureNugget")
        {
            vInfo$nstr = 0
            vInfo$phi = -9999
            vInfo$phiLock = TRUE
        }
		vInfo$kappa = -9999
		vInfo$kappaLock = TRUE
		vInfo$kappaMax = -9999
	}
    useSparse = FALSE
    if(any(names(vInfo) == "useSparse"))useSparse = vInfo$useSparse
    if(!is.logical(useSparse))stop("lmm.errorCheck: 'useSparse' must be logical!")
	if(vInfo$name != "doubleSpherical" & vInfo$name != "pureNugget" & vInfo$name != "spherical")useSparse = FALSE
    if(!any(names(vInfo) == "phiLock"))vInfo$phiLock = FALSE
    if(!any(names(vInfo) == "effrMin"))stop("lmm.errorCheck: Provide 'effrMin' in vInfo (the smallest lag of the variogram)!")
    if(min(vInfo$effrMin) <= 0)stop("lmm.errorCheck: 'effrMin' must be >0!")
	if(!any(names(vInfo) == "effrMax"))stop("lmm.errorCheck: Provide 'effrMax' in vInfo (the largest lag of the variogram)!")
    if((vInfo$name == "doubleSpherical" | vInfo$name == "exponential" | vInfo$name == "spherical") & vInfo$phiLock)
    {
        if(!any(names(vInfo) == "effrMax"))stop("lmm.errorCheck: 'effrMax' must be provided in vInfo when choosing to lock 'phi'!")
        if(vInfo$name == "doubleSpherical")vInfo$phi = c(vInfo$effrMin, vInfo$effrMax)
        if(vInfo$name == "exponential")vInfo$phi = vInfo$effrMax / 3
        if(vInfo$name == "dampedPeriodic" | vInfo$name == "spherical")vInfo$phi = vInfo$effrMax
    }
    if((vInfo$name == "dampedPeriodic" | vInfo$name == "matern") & (vInfo$phiLock | vInfo$kappaLock))
    {
        if(vInfo$name == "dampedPeriodic")stop("lmm.errorCheck: 19.2.2015 - sort out the locking on dampedPeriodic!")
        tmp = cbind(seq(vInfo$effrMin, vInfo$effrMax, length.out = 500), NA, NA)
        for(i in 1:nrow(tmp))
        {
            tmp[i, 2] = lmm.effectiveRange(vInfo$name, tmp[i, 1], vInfo$kappa, vInfo$effrMax)
        }
        tmp[, 3] = abs(tmp[, 2] - vInfo$effrMax)
        ndx = which(tmp[, 3] == min(tmp[, 3]), arr.ind = TRUE)[1]
        vInfo$phi = tmp[ndx, 1]
    }
    if(!any(names(vInfo) == "cvci"))vInfo$cvci = FALSE
    if(!is.logical(nesting))vInfo$cvci = FALSE
	if(!any(names(vInfo) == "plotIt"))vInfo$plotIt = TRUE
    
    #On "mInfo".
    if(mInfo$method != "lbfgsb" & mInfo$method != "sanneal" & mInfo$method != "simplex")stop("lmm.errorCheck: Invalid minimisation method!")
	if(nv > 1 & mInfo$method == "lbfgsb")stop("lmm.errorCheck: Minimisation by 'lbfgsb' cannot be used when nv > 1!")
    if(mInfo$method == "sanneal")
    {
        if(!any(names(mInfo) == "method"))mInfo$nStep = 1000
        if(!any(names(mInfo) == "nMarkov"))mInfo$nMarkov = 30
        if(!any(names(mInfo) == "heat"))mInfo$heat = 100
        if(!any(names(mInfo) == "alpha"))mInfo$alpha = 0.95
        if(!any(names(mInfo) == "nStop"))mInfo$nStop = 10
    }
    if(os == "Windows")
    {
        mInfo$ncore = 1
    }else{
        if(any(names(mInfo) == "ncore"))
        {
            mInfo$ncore = min(mInfo$ncore, 8)
        }else{
            mInfo$ncore = 8
        }
    }
    if(!any(names(mInfo) == "silent"))mInfo$silent = FALSE
    
    #Output list.
    out = list("dInfo" = dInfo, "vInfo" = vInfo, "mInfo" = mInfo,
        "nv" = nv, "nesting" = nesting, "useSparse" = useSparse)
    out
}


lmm.findNA = function(x)
{
    ###########################################################
    #Find the rows with missing data in a matrix or dataframe.#
    ###########################################################
    any(is.na(x))
}


lmm.insideCI = function(nv, medianTheta)
{
    ############################################################
    #Interrogate the "medianTheta" statistic(s) to see if the  #
    #residuals of the model are plausibly normally distributed.#
    ############################################################
    out = rep(FALSE, nv)
    for(i in 1:nv)
    {
        if(!is.na(medianTheta$"q2.5"[i]) & !is.na(medianTheta$"q97.5"[i]))
        {
            out[i] = medianTheta$value[i] > medianTheta$"q2.5"[i] & medianTheta$value[i] < medianTheta$"q97.5"[i]
        }
    }
    out
}

lmm.lbfgsb = function(z, os, nv, nr, nrTotal, Xmat, Xv, d, vInfo, mInfo, stdev, sgn, nesting, useSparse)
{
	#######################################################################
	#L-BFGS-B algorithm to minimise the negative log. likelihood function.#
	#######################################################################
    #The total number of parameters.
	n = nv
	if(vInfo$name != "pureNugget" & vInfo$name != "dampedPeriodic")n = n * 2
    if(!vInfo$phiLock)n = n + 1
	if(!vInfo$kappaLock)n = n + 1
	
    #Randomly perturb the initial covariance parameter values.
    tmp = lmm.perturbInitial(nv, vInfo)
    phi = tmp$phi
    kappa = tmp$kappa
    sig = array(0, c(nv, nv, 3))
    sig[, , 1] = tmp$c0
    if(vInfo$nstr > 0)
    {
        sig[, , 2] = tmp$c1
        if(vInfo$nstr == 2)sig[, , 3] = tmp$c2
    }
    
    #Box constraints for the parameters.
    vInfoMin = vInfo
    vInfoMax = vInfo
    vInfoMin$sig = sig * 0 + 10 ^ -12
    vInfoMax$sig = sig * 3
    if(vInfo$name != "pureNugget" & !vInfo$phiLock)
    {
        if(vInfo$name == "dampedPeriodic")
		{
            vInfoMin$phi = 90 #Note the hard-coded override of "vInfo$effrMin".
            vInfoMax$phi = vInfo$effrMax
		}
        if(vInfo$name == "doubleSpherical")
		{
			vInfoMin$phi = c(vInfo$effrMin, vInfo$effrMin + vInfo$effrMin * 0.1)
			vInfoMax$phi = c(vInfo$effrMax - vInfo$effrMax * 0.1, vInfo$effrMax)
		}
		if(vInfo$name == "exponential")
		{
			vInfoMin$phi = vInfo$effrMin / 3
			vInfoMax$phi = vInfo$effrMax / 3
		}
		if(vInfo$name == "matern" | vInfo$name == "spherical")
		{
			vInfoMin$phi = vInfo$effrMin
			vInfoMax$phi = vInfo$effrMax
		}
	}
	if(!vInfo$kappaLock)
	{
		if(vInfo$name == "dampedPeriodic")
		{
            vInfoMin$kappa = -1.570795
            vInfoMax$kappa = -0.001
		}
		if(vInfo$name == "matern")
		{
			vInfoMin$kappa = 0.01
            vInfoMax$kappa = 5
		}
	}
    
    #Map the initial parameters (and constraints) into a single variable.
    p = rep(0, n)
	par0 = lmm.map(nv, vInfo, sig, phi, kappa, n, p, FALSE)$p
    parLower = lmm.map(nv, vInfoMin, vInfoMin$sig, vInfoMin$phi, vInfoMin$kappa, n, p, FALSE)$p
    parUpper = lmm.map(nv, vInfoMax, vInfoMax$sig, vInfoMax$phi, vInfoMax$kappa, n, p, FALSE)$p
    
    #List where output will be stored.
    out = list("l" = rep(NA, 2))
    
    #Initial likelihood.
    out$l[1] = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, phi, kappa,
        sig[, , 1], sig[, , 2], sig[, , 3], nesting, useSparse)$l
    if(is.finite(out$l[1]))
    {
        #The initial parameter values are valid, so proceed.
        if(!mInfo$silent)cat("Initial value of the objective function =", out$l[1], "\n")
        tmp = optim(par0, lmm.reml_lbfgsb, gr = NULL,
            os = os, nv = nv, nr = nr, nrTotal = nrTotal, y = z, X = Xmat, d = d,
            name = vInfo$name, nstr = vInfo$nstr, phi = vInfo$phi, phiLock = vInfo$phiLock,
            kappa = vInfo$kappa, kappaLock = vInfo$kappaLock, nesting = nesting, 
            useSparse = useSparse, lOnly = TRUE,
            method = "L-BFGS-B", lower = parLower, upper = parUpper)
        parList = lmm.map(nv, vInfo, sig, phi, kappa, n, tmp$par, TRUE)
        out$l[2] = tmp$value
        if(!mInfo$silent)
        {
            cat("Solution =\n")
            cat("   l =", out$l[2], "\n")
        }
        effr = NULL
        if(vInfo$name != "pureNugget" & (!vInfo$phiLock | !vInfo$kappaLock))
        {
            effr = lmm.effectiveRange(vInfo$name, parList$phi, parList$kappa, vInfo$effrMax)
            if(!mInfo$silent)
            {
                if(!vInfo$phiLock & !vInfo$kappaLock)cat("   phi =", parList$phi, "kappa =", parList$kappa, "effr =", effr, "\n")
		        if(!vInfo$phiLock & vInfo$kappaLock)cat("   phi =", parList$phi, "effr =", effr, "\n")
				if(vInfo$phiLock & !vInfo$kappaLock)cat("   kappa =", parList$kappa, "effr =", effr, "\n")
			}
        }
        
        #Put in one last call to "lmm.reml" with the final parameter values.
        finalCall = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, parList$phi, parList$kappa,
            parList$sig[, , 1], parList$sig[, , 2], parList$sig[, , 3], nesting, useSparse)
        if(!mInfo$silent)
		{
            #Print.
            tmp = lmm.rescale(nv, stdev, Xv, parList$sig, finalCall$b, NA) #Rescale "sig" and "b" for printing purposes.
		    cat("   sig =\n")
			print(tmp$sig) #Note that it is printed on the original scale.
			cat("   b =\n")
			print(tmp$b) #Note that it is printed on the original scale.
			cat("\n")
		}
    }else{
        stop("lmm.lbfgsb: Infinite initial likelihood!")
    }
    
    #Fill the output list.
    out$covFun = vInfo$name
	out$phi = parList$phi
	out$kappa = parList$kappa
    out$effr = effr
    out$c0 = parList$sig[, , 1]
    out$c1 = parList$sig[, , 2]
	out$c2 = parList$sig[, , 3]
	out$b = finalCall$b
    out$b_percentVarianceExplained = finalCall$b_percentVarianceExplained
    out$C = finalCall$C
    out$l = finalCall$l
    out
}


lmm.makeV = function(os, nv, nr, nrTotal, d, name, nstr, phi, kappa, c0, c1, c2, useSparse, returnCoordFormat)
{
	###############################################
	#Make the covariance matrix between locations.#
	###############################################
    #Adjust the inputs class if necessary.
    if(nv == 1)
    {
        c0 = as.matrix(c0)
        c1 = as.matrix(c1)
        c2 = as.matrix(c2)
    }
    phi = c(phi)
    
	#Make a point-to-point covariance matrix.
    if(!useSparse)
	{
        V = matrix(0, nrow = nrTotal, ncol = nrTotal)
        if(nv == 1)
		{
            V = lmm.covFun(d, name, phi, kappa, c0[1, 1], c1[1, 1], c2[1, 1])
		}else{
			kount = rep(0, 4)
			for(i in 1:nv)
			{
				kount[1] = kount[2] + 1
				kount[2] = kount[2] + nr[i]
				for(j in i:nv)
				{
					if(i == j & i > 1)kount[4] = kount[2] - nr[j]
					kount[3] = kount[4] + 1
					kount[4] = kount[4] + nr[j]
                    V[kount[1]:kount[2], kount[3]:kount[4]] = lmm.covFun(d[kount[1]:kount[2], kount[3]:kount[4]], name, phi, kappa, c0[i, j], c1[i, j], c2[i, j])
					if(i != j)V[kount[3]:kount[4], kount[1]:kount[2]] = t(V[kount[1]:kount[2], kount[3]:kount[4]])
				}
			}
		}
	}else{
        V = Matrix(0, nrow = nrTotal, ncol = nrTotal, sparse = TRUE)
        if(os == "Windows")
        {
    		for(i in 1:nrow(d))
	    	{
		    	#Get the relevant variable numbers.
			    variable = list()
			    tmp = cumsum(nr)
			    for(j in c("row", "col"))
			    {
                    ndx = which(tmp >= d[[j]][i])
                    if(length(ndx) == 0)ndx = 1
				    variable[[j]] = ndx[1]
			    }
			    
			    #Write the covariance into the relevant element.
                V[d$row[i], d$col[i]] = lmm.covFun(d$d[i], name, phi, kappa, c0[variable$row, variable$col], c1[variable$row, variable$col], c2[variable$row, variable$col])
		    }
        }else{
            #Use FORTRAN to do it quicker.
            dyn.load("~/hg/rscrepo/rsc/geostats/lmm/makeV.so")
            if(name == "spherical")phi = c(phi, 0)
            n = nrow(d)
            if(n == 0)stop("lmm.makeV: There are no data in the sparse distance matrix!")
            tmp = .Fortran("makeV", as.integer(nrow(d)), as.double(as.matrix(d)), as.integer(nv), as.integer(nr),
                as.integer(nstr), as.double(phi), as.double(c0), as.double(c1), as.double(c2),
                as.double(as.matrix(d * 0)), as.integer(n))
            n = tmp[[11]] #The number of non-zero covariances.
            if(n == 0)stop("lmm.makeV: There are no data in the sparse covariance matrix!")
            tmp = matrix(tmp[[10]], nrow = nrow(d), ncol = 3)
            tmp = tmp[1:n, ]
            if(n == 1)tmp = t(as.matrix(tmp))
            if(returnCoordFormat)
            {
                #Return the lower triangle in coordinate-sparse format.
                V = data.frame(tmp)
                names(V) = c("row", "column", "covariance")
                return(V)
            }
            covariances = tmp[, 3]
            ndx = sum(nr) * (tmp[, 2] - 1) + tmp[ ,1]
            V[ndx] = covariances #Write into the sparse matrix.
        }
		V = forceSymmetric(V, "L")
	}
    
    #Check if "V" is positive definite.
    tmp = try(chol(V), silent = TRUE)
    if(is.character(tmp))
    {
        #"V" is not positive definite. Return an empty list.
        out = list()
        return(out)
    }else{
        #"V" is positive definite.
        #While we're here find the Cholesky decomposition and the log-determinant of "V".
        out = list("V" = V, "Vchol" = tmp)
        out$Vdet = 2 * sum(log(diag(out$Vchol)))
    }
	out
}


lmm.map = function(nv, vInfo, sig, phi, kappa, n, p, invert)
{
	##########################################################
	#Maps the covariance parameters from one form to another.#
	##########################################################
	if(!invert)
    {
        #Map the disparate parameters to a single variable.
		kount = 0
        for(i in 1:(vInfo$nstr + 1))
        {
    		for(j in 1:nv)
            {
                for(k in j:nv)
                {
                    kount = kount + 1
                    p[kount] = sig[j, k, i]
                }
            }
        }
		if(!vInfo$phiLock)
        {
            for(i in 1:vInfo$nstr)
            {
                kount = kount + 1
			    p[kount] = phi[i]
            }
		}
        if(!vInfo$kappaLock)
        {
            kount = kount + 1
			p[kount] = kappa
        }
	}else{
		#Map the single variable back to the disparate parameters.
		kount = 0
        for(i in 1:(vInfo$nstr + 1))
        {
    		for(j in 1:nv)
            {
                for(k in j:nv)
                {
                    kount = kount + 1
                    sig[j, k, i] = p[kount]
                    if(j != k)sig[k, j, i] = sig[j, k, i]
                }
            }
        }
		if(!vInfo$phiLock)
        {
            for(i in 1:vInfo$nstr)
            {
                kount = kount + 1
                phi[i] = p[kount]
            }
        }
        if(!vInfo$kappaLock)
        {
            kount = kount + 1
            kappa = p[kount]
        }
    }
    out = list("sig" = sig, "phi" = phi, "kappa" = kappa, "p" = p)
    out
}


lmm.metrop = function(l, l0, heat)
{
	#######################
	#Metropolis criterion.#
	#######################
	out = FALSE
	if(runif(1) <= exp((l0 - l) / heat))out = TRUE
	out
}


lmm.perturbInitial = function(nv, vInfo)
{
    ####################################################################
    #Perturb some of the initial covariance parameters except "kappa").#
    ####################################################################
    phi = vInfo$phi
    if(vInfo$name != "pureNugget" & !vInfo$phiLock)phi = phi * runif(length(phi), 0.8, 1.2)
    if((vInfo$nstr == 2 && phi[1] >= phi[2]) | min(phi) < vInfo$effrMin | max(phi) > vInfo$effrMax)phi = vInfo$phi
    for(i in c("c0", "c1", "c2"))
    {
        if(any(vInfo[[i]] != 0))
        {
            diag(vInfo[[i]]) = diag(vInfo[[i]]) + runif(nv) * diag(vInfo[[i]]) * 0.5
		}
	}
    out = list("phi" = phi, "kappa" = vInfo$kappa, "c0" = vInfo$c0, "c1" = vInfo$c1, "c2" = vInfo$c2)
    out
}


lmm.reml = function(os, nv, nr, nrTotal, y, X, d, name, nstr, phi, kappa, c0, c1, c2, nesting, useSparse)
{
	##############################
	#Residual maximum likelihood.#
	##############################
	#Make covariance matrix "V", it's Cholesky decomposition, and it's determinant.
	tmp = lmm.makeV(os, nv, nr, nrTotal, d, name, nstr, phi, kappa, c0, c1, c2, useSparse, FALSE)
    if(is.null(tmp$V))
    {
        #"V" was not positive definite. Return an infinite liklihood.
        out = list("l" = Inf)
        return(out)
    }else{
        #"V" was positive definite, so unpack the list.
        Vchol = tmp$Vchol
        logVdet = tmp$Vdet
    }
    rm(tmp)
    
    #Inverse of "V".
    Vi = chol2inv(Vchol)
    
    #Matrix multiply the transpose of "X" by the inverse of "V".
    XtVi = t(X) %*% Vi
    
    #Matrix multiply "XtVi" by "X".
	W = forceSymmetric(XtVi %*% X, "U")
    
    #Calculate the log determinant of "W".
	logWdet = 2 * sum(log(diag(chol(W))))
	
    #Generalised least-squares estimates of the trend parameter(s).
    if(is.logical(nesting))
    {
        #FULL MODEL.
        b = chol2inv(chol(W)) %*% XtVi %*% y
    }else{
        #NESTED MODEL.
        #Subset.
        X = X[, -nesting]
        if(length(X) == nrTotal)X = Matrix(X)
        
        #Matrix multiply the transpose of "X" by the inverse of "V".
        XtVi = t(X) %*% chol2inv(Vchol)
        
	    #Matrix multiply "XtVi" by "X".
	    W = forceSymmetric(XtVi %*% X, "U")
        
	    #Generalised least-squares estimates of the trend parameter(s).
        b = chol2inv(chol(W)) %*% XtVi %*% y
    }
    
    #Weighted squared error.
    pred = X %*% b
    e = y - pred
	wse = t(e) %*% Vi %*% e
	
    #Percent variance explained by the trend model.
    if(nv == 1)
    {
        b_percentVarianceExplained = var(as.matrix(pred)) / var(as.matrix(y))
    }else{
        obs = as.matrix(y)
        pred = as.matrix(pred)
        b_percentVarianceExplained = rep(0, nv)
        kount = rep(0, 2)
        for(i in 1:nv)
        {
            kount[1] = kount[2] + 1
            kount[2] = kount[2] + nr[i]
            b_percentVarianceExplained[i] = (var(pred[kount[1]:kount[2]]) / var(obs[kount[1]:kount[2]]))
        }
    }
    b_percentVarianceExplained = b_percentVarianceExplained * 100
    
	#Negative log-likelihood.
	bcPart = 0
	##if(doBoxcox)
	##{
	##	stop("reml boxcox stuff not written yet!")
	##	do i=1,nv
	##		if(lock_boxcox(i)<0)bcpart=bcpart+(boxcox(i)-1.0)*sumlogy(i)
	##	enddo
	##}
    out = list("b" = as.matrix(b),
        "b_percentVarianceExplained" = b_percentVarianceExplained,
        "l" = -1 * c(as.matrix((bcPart - 0.5 * logVdet - 0.5 * logWdet - 0.5 * wse))))
    if(is.logical(nesting))out$C = chol2inv(chol(W)) #Covariance matrix for the fixed effects.
    out
}


lmm.reml_lbfgsb = function(par, os, nv, nr, nrTotal, y, X, d, name, nstr, phi, phiLock, kappa, kappaLock, nesting, useSparse, lOnly)
{
	#######################################################################
	#Sets up the input from "lmm.lbfgsb" so that "lmm.reml" can be called.#
	#######################################################################
    #Invert the parameter mapping.
    vInfo = list("name" = name, "nstr" = nstr, "phiLock" = phiLock, "phi" = phi, "kappaLock" = kappaLock, "kappa" = kappa)
    sig = array(0, c(nv, nv, 3))
    phi = vInfo$phi
    kappa = vInfo$kappa
    parList = lmm.map(nv, vInfo, sig, phi, kappa, length(par), par, TRUE)
    
    #Calculate the likelihood.
    tmp = lmm.reml(os, nv, nr, nrTotal, y, X, d, name, nstr, parList$phi, parList$kappa,
        parList$sig[, ,1], parList$sig[, ,2], parList$sig[, ,3], nesting, useSparse)
    if(lOnly)
	{
		out = tmp$l
	}else{
		out = tmp
	}
	out
}

lmm.reml_simplex = function(par, os, nv, nr, nrTotal, y, X, d, name, nstr, effrMin, effrMax, phiLock, phi, kappaLock, kappa, sgn, nesting, useSparse, lOnly)
{
	########################################################################
	#Sets up the input from "lmm.simplex" so that "lmm.reml" can be called.#
	########################################################################
    #Invert the parameter mapping.
    vInfo = list("name" = name, "nstr" = nstr, "effrMin" = effrMin, "effrMax" = effrMax,
        "phiLock" = phiLock, "phi" = phi, "kappaLock" = kappaLock, "kappa" = kappa)
    sig = array(0, c(nv, nv, 3))
    phi = vInfo$phi
    kappa = vInfo$kappa
    parList = lmm.map(nv, vInfo, sig, phi, kappa, length(par), par, TRUE)
    
    #Check the parameters.
    chk = lmm.check(os, nv, name, nstr, parList$phi, parList$kappa, effrMin, effrMax,
		parList$sig[, , 1], parList$sig[, , 2], parList$sig[, , 3], sgn, 4)
    if(!chk)
    {
        #There is at least one invalid parameter. Return an infinite likelihood.
        out = Inf
    }else{
        #The parameters are valid. Compute the likelihood.
        tmp = lmm.reml(os, nv, nr, nrTotal, y, X, d, name, nstr, parList$phi, parList$kappa,
            parList$sig[, , 1], parList$sig[, , 2], parList$sig[, , 3], nesting, useSparse)
        if(lOnly)
	    {
            out = tmp$l
        }else{
		    out = tmp
	    }
    }
    out
}


lmm.rescale = function(nv, stdev, Xv, sig, b, Cmat)
{
	################################################################################
	#Rescale "sig", "b" and "Cmat" by the relevant standard deviation, if required.#
	################################################################################
	out = list("sig" = sig * 0, "b" = b * 0)
    if(all(is.finite(stdev)))
    {
        u = apply(Xv, 2, max)
	    for(i in 1:nv)
	    {
		    for(j in i:nv)
		    {
                out$sig[i, j, ] = sig[i, j, ] * stdev[i] * stdev[j]
			    if(i != j)out$sig[j, i, ] = out$sig[i, j, ]
		    }
            ndxi = which(u == i, arr.ind = TRUE)
		    out$b[ndxi, 1] = b[ndxi, 1] * stdev[i]
	    }
        if(any(!is.na(Cmat)))
        {
            out$C = Cmat * 0
            for(i in 1:nv)
	        {
                ndxi = which(u == i, arr.ind = TRUE)
		        for(j in 1:nv)
                {
                    ndxj = which(u == j, arr.ind = TRUE)
                    out$C[ndxi, ndxj] = Cmat[ndxi, ndxj] * stdev[i] * stdev[j]
                }
	        }
        }
    }else{
        #No rescaling is needed.
        out$sig = sig
        out$b = b
        out$C = Cmat
    }
    out
}


lmm.sanneal = function(z, os, nv, nr, nrTotal, Xmat, Xv, d, vInfo, mInfo, stdev, sgn, nesting, useSparse)
{
	########################################################################
	#Simulated annealing to minimise the negative log. likelihood function.#
	########################################################################
    #Randomly perturb the initial covariance parameter values.
    #Then stack the coregionalisation matrices.
    tmp = lmm.perturbInitial(nv, vInfo)
    vInfo$phi = tmp$phi
    vInfo$kappa = tmp$kappa
    vInfo$sig = array(0, c(nv, nv, 3))
    vInfo$sig[, , 1] = tmp$c0
    if(vInfo$nstr > 0)
    {
        vInfo$sig[, , 2] = tmp$c1
        if(vInfo$nstr == 2)vInfo$sig[, , 3] = tmp$c2
    }
    
    ###Apply the Box-Cox transform to the scaled values.
	##do i=1,nv
	##	call transform(nr,boxcox,lock_boxcox,i,y,ybc,-1)
	##enddo
    
	#Check that the initial variogram parameters are valid.
    chk = lmm.check(os, nv, vInfo$name, vInfo$nstr, vInfo$phi, vInfo$kappa,
		vInfo$effrMin, vInfo$effrMax,
		vInfo$sig[, , 1], vInfo$sig[, , 2], vInfo$sig[, , 3], sgn, 4)
    if(!chk)stop("lmm.sanneal: At least at one initial covariance parameter is invalid!")
	
	#Initialise the negative log. likelihood.
    out = list("l" = rep(NA, mInfo$nStep))
    tmp = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, vInfo$phi, vInfo$kappa,
        vInfo$sig[, , 1], vInfo$sig[, , 2], vInfo$sig[, , 3], nesting, useSparse)
    b = tmp$b
	out$l[1] = tmp$l
	if(!mInfo$silent)cat("Initial value of the objective function =", out$l[1], "\n")
	
	#Parameters that constrain how much the variogram parameters can change with each iteration.
	delta = list()
	##delta$boxcox = 0.1
	delta$phi = max(d) * 0.1
	delta$kappa = 0.1
	delta$sig = matrix(0, nv, nv)
    for(i in 1:nv)
	{
		for(j in i:nv)
		{
			if(i == j)
			{
				delta$sig[i, j] = (vInfo$c0[i, j] + vInfo$c1[i, j] + vInfo$c2[i, j]) * 0.05
			}else{
				delta$sig[i, j] = 0.05 * sqrt(vInfo$c0[i, i] * vInfo$c0[j, j] + vInfo$c1[i, i] * vInfo$c1[j, j] + vInfo$c2[i, i] * vInfo$c2[j, j])
			}
		}
	}
    
    #The initial temperature of the system.
	heat = mInfo$heat
	
    #These are used to count the number of Markov chains for which the log. likelihood has not been changed.
	llast = out$l[1]
	nUnch = 0
	
	#Loop over the cooling steps.
	iStep = 1
	reset = FALSE
    kountMax = 500
	while(iStep <= mInfo$nStep & nUnch <= mInfo$nStop)
	{
		if(iStep > mInfo$nStep | nUnch > mInfo$nStop)break
		iStep = iStep + 1
		
		#Reset to the initial values if necessary.
		if(iStep == 2)
		{
			##boxcox = variogramInfo$boxcox
			phi = vInfo$phi
			kappa = vInfo$kappa
			sig = vInfo$sig
		}
		
		#Iterate within the step.
		out$l[iStep] = out$l[iStep - 1]
		rej = 0
		acc = 0
		for(imc in 1:mInfo$nMarkov)
		{
###print(imc)
##			#Adjust "boxcox", if required.
##			if(doBoxcox)
##			{
##				kount = rep(0, 2)
##				for(i in 1:nv)
##				{
##					kount[1] = kount [2] + 1
##					kount[2] = kount[2] + nr[i]
##					if(!boxcoxLock[i])
##					{
##						l0 = out$l[iStep]
##						boxcox0 = boxcox[i]
##						ybc0[1:nrTotal, 1] = ybc[1:nrTotal, 1]
##					    kount(3)=0
##					    do
##						    kount(3)=kount(3)+1
##						    r=(ran2(seed)-0.5)*2*boxcox_delta
##						    boxcox(i)=boxcox0+r
##						    chk=check(nv,nr_total,np,nsvcl,d,covfun,phi,kappa,covariate,lambda,boxcox,meffr,effr,sig,prop,sgn,1)
##						    if(chk.or.kount(3)==mkount)exit
##					    enddo
##					    if(chk)then
##						    call transform(nr,boxcox,lock_boxcox,i,y,ybc,-1)
##						    call reml(nv,nr,nr_total,np,nsvcl,ybc,sumlogy,X_isnested,X,Xt,X_ndx,d,covfun,phi,kappa,covariate,lambda,boxcox,do_boxcox,lock_boxcox,sig,prop,iV,detV,0,b,C,l(is))
##						    if(l(is)<=l0)then
##							    acc=acc+1
##						    else
##							    if(metrop(seed,l(is),l0,heat)<0)then
##								    rej=rej+1
##								    boxcox(i)=boxcox0
##								    l(is)=l0
##							    else
##								    acc=acc+1
##							    endif
##						    endif
##					    else
##						    rej=rej+1
##						    boxcox(i)=boxcox0
##						    ybc(1:nr_total,1)=ybc0(1:nr_total,1)
##						    l(is)=l0
##					}
##				}
##			}
			
			#Adjust "phi", if required.
			if(!vInfo$phiLock)
			{
                for(istr in 1:vInfo$nstr)
                {
                    l0 = out$l[iStep]
	    			phi0 = phi[istr]
		    		kount = 0
			    	while(kount < kountMax)
				    {
					    kount = kount + 1
					    r = (runif(1) - 0.5) * 2 * delta$phi
					    phi[istr] = phi0 + r
                        chk = lmm.check(os, nv, vInfo$name, vInfo$nstr, phi, kappa,
						    vInfo$effrMin, vInfo$effrMax,
						    sig[, , 1], sig[, , 2], sig[, , 3], sgn, 2)
                        if(chk)break
				    }
                    if(chk)
				    {
					    tmp = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, phi, kappa,
                            sig[, , 1], sig[, , 2], sig[, , 3], nesting, useSparse)
					    b = tmp$b
                        Cmat = tmp$C
					    out$l[iStep] = tmp$l
					    if(out$l[iStep] <= l0)
					    {
						    acc = acc + 1
					    }else{
						    if(!lmm.metrop(out$l[iStep], l0, heat))
						    {
							    rej = rej + 1
							    phi[istr] = phi0
							    out$l[iStep] = l0
						    }else{
							    acc = acc + 1
						    }
					    }
				    }else{
					    rej = rej + 1
					    phi[istr] = phi0
					    out$l[iStep] = l0
				    }
                }
			}
            
			#Adjust "kappa", if required.
			if(!vInfo$kappaLock)
			{
				l0 = out$l[iStep]
				kappa0 = kappa
				kount = 0
				while(kount < kountMax)
				{
					kount = kount + 1
					r = (runif(1) - 0.5) * 2 * delta$kappa
					kappa = kappa0 + r
                    chk = lmm.check(os, nv, vInfo$name, vInfo$nstr, phi, kappa,
						vInfo$effrMin, vInfo$effrMax,
						sig[, , 1], sig[, , 2], sig[, , 3], sgn, 2)
					if(chk)break
				}
                if(chk)
				{
					tmp = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, phi, kappa,
                        sig[, , 1], sig[, , 2], sig[, , 3], nesting, useSparse)
                    b = tmp$b
                    Cmat = tmp$C
					out$l[iStep] = tmp$l
					if(out$l[iStep] <= l0)
					{
						acc = acc + 1
					}else{
						if(!lmm.metrop(out$l[iStep], l0, heat))
						{
							rej = rej + 1
							kappa = kappa0
							out$l[iStep] = l0
						}else{
							acc = acc + 1
						}
					}
				}else{
					rej = rej + 1
					kappa = kappa0
					out$l[iStep] = l0
				}
			}
			
			#Loop over the coregionalisation matrices.
            for(istr in 1:(vInfo$nstr + 1))
            {
                for(i in 1:nv)
			    {
				    for(j in i:nv)
				    {
					    #Adjust "sig".
					    l0 = out$l[iStep]
                        sig0 = sig[i, j, istr]
					    kount = 0
					    while(kount < kountMax)
					    {
						    kount = kount + 1
						    r = (runif(1) - 0.5) * 2 * delta$sig[i, j]
						    sig[i, j, istr] = sig0 + r
						    if(i != j)sig[j, i, istr] = sig[i, j, istr]
                            chk = lmm.check(os, nv, vInfo$name, vInfo$nstr, phi, kappa,
	    					    vInfo$effrMin, vInfo$effrMax,
		    				    sig[, , 1], sig[, , 2], sig[, , 3], sgn, 3)
                            if(chk)break
					    }
					    if(chk)
					    {
						    tmp = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, phi, kappa,
                                sig[, , 1], sig[, , 2], sig[, , 3], nesting, useSparse)
                            b = tmp$b
                            Cmat = tmp$C
						    out$l[iStep] = tmp$l
						    if(out$l[iStep] <= l0)
						    {
							    acc = acc + 1
						    }else{
							    if(!lmm.metrop(out$l[iStep], l0, heat))
							    {
								    rej = rej + 1
								    sig[i, j, istr] = sig0
								    if(i != j)sig[j, i, istr] = sig[i, j, istr]
								    out$l[iStep] = l0
							    }else{
								    acc = acc + 1
							    }
						    }
					    }else{
						    rej = rej + 1
						    sig[i, j, istr] = sig0
						    if(i != j)sig[j, i, istr] = sig[i, j, istr]
						    out$l[iStep] = l0
					    }
                    }
                }
			}
		}
		pacc = acc / (acc + rej)
		if(!mInfo$silent)
		{
			cat("Step =", iStep - 1, "of a possible", mInfo$nStep, "\n")
			cat("   heat =", heat, "l =", out$l[iStep], "pacc =", pacc, "\n")
		}
		##if(do_boxcox)print*," boxcox=",real(boxcox(1:nv))
		effr = NULL
		if(vInfo$name != "pureNugget" & (!vInfo$phiLock | !vInfo$kappaLock))
		{
			effr = lmm.effectiveRange(vInfo$name, phi, kappa, vInfo$effrMax)
			if(!mInfo$silent)
			{
				if(!vInfo$phiLock & !vInfo$kappaLock)cat("   phi =", phi, "kappa =", kappa, "effr =", effr, "\n")
				if(!vInfo$phiLock & vInfo$kappaLock)cat("   phi =", phi, "effr =", effr, "\n")
				if(vInfo$phiLock & !vInfo$kappaLock)cat("   kappa =", kappa, "effr =", effr, "\n")
			}
		}
		if(!mInfo$silent)
		{
			tmp = lmm.rescale(nv, stdev, Xv, sig, b, NA) #Rescale "sig" and "b" for printing purposes.
			cat("   sig =\n")
			print(tmp$sig)
			cat("   b =\n")
			print(tmp$b)
			cat("\n")
		}
		
		#Check after completion of the first Markov chain. Return if the algorithm requires different "heat".
		if(iStep == 2)
		{
			reset = FALSE
			if(pacc < 0.8 | pacc > 0.999)
			{
				#Try with a new "heat" on the next iteration.
				reset = TRUE
				inc = heat * 0.1
                if(heat <= 1000)inc = 50
				if(heat <= 100)inc = 5
				if(heat <= 10)inc = 0.5
				if(heat <= 1)inc = 0.05
				if(heat <= 0.1)inc = 0.005
				if(heat <= 0.01)inc = 0.0005
				if(heat <= 0.001)inc = 0.00005
				if(heat <= 0.0001)inc = heat * 0.5
				if(pacc < 0.9)heat = heat + inc
				if(pacc > 0.999)heat = heat - inc
                if(heat < 0)stop("lmm.sanneal: Negative heat! How has this happened?")
				iStep = 1
			}
		}
		if(!reset)
		{
			#Update "llast" and "nUnch".
			if(llast == out$l[iStep])
			{
				nUnch = nUnch + 1
			}else{
				llast = out$l[iStep]
				nUnch = 0
			}
			
			#Cool.
			heat = heat * mInfo$alpha
		}
	}
    
    #Output.
	out$covFun = vInfo$name
	out$phi = phi
	out$kappa = kappa
	out$effr = effr
	out$c0 = sig[, , 1]
    out$c1 = sig[, , 2]
	out$c2 = sig[, , 3]
    finalCall = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, phi, kappa, out$c0, out$c1, out$c2, nesting, useSparse)
	out$b = finalCall$b
    out$b_percentVarianceExplained = finalCall$b_percentVarianceExplained
    out$C = finalCall$C
	out$l = out$l[1:iStep]
	out
}


lmm.simplex = function(z, os, nv, nr, nrTotal, Xmat, Xv, d, vInfo, mInfo, stdev, sgn, nesting, useSparse)
{
	########################################################################
	#Nelder-Mead simplex to minimise the negative log. likelihood function.#
    ########################################################################
    #The total number of parameters.
    n = nv * (nv + 1) / 2
	if(vInfo$name != "pureNugget" & vInfo$name != "dampedPeriodic")
    {
        if(vInfo$nstr == 1)
        {
            n = n * 2
        }else{
            n = n * 3
        }
    }
    if(!vInfo$phiLock)n = n + length(vInfo$phi)
	if(!vInfo$kappaLock)n = n + 1
    
    #Randomly perturb the initial covariance parameter values.
    #Then stack the coregionalisation matrices.
    tmp = lmm.perturbInitial(nv, vInfo)
    phi = tmp$phi
    kappa = tmp$kappa
    sig = array(0, c(nv, nv, 3))
    sig[, , 1] = tmp$c0
    if(vInfo$nstr > 0)
    {
        sig[, , 2] = tmp$c1
        if(vInfo$nstr == 2)sig[, , 3] = tmp$c2
    }
    
    #Map the initial parameters into a single variable.
    p = rep(0, n)
	par0 = lmm.map(nv, vInfo, sig, phi, kappa, n, p, FALSE)$p
    
    #List where output will be stored.
    out = list("l" = rep(NA, 2))
    
    #Initial likelihood.
    out$l[1] = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, phi, kappa,
        sig[, , 1], sig[, , 2], sig[, , 3], nesting, useSparse)$l
    if(is.finite(out$l[1]))
    {
        #The initial parameter values are valid, so proceed.
        if(!mInfo$silent)cat("Initial value of the objective function =", out$l[1], "\n")
        tmp = optim(par0, lmm.reml_simplex,
            os = os, nv = nv, nr = nr, nrTotal = nrTotal, y = z, X = Xmat, d = d,
            name = vInfo$name, nstr = vInfo$nstr,
            effrMin = vInfo$effrMin, effrMax = vInfo$effrMax,
            phiLock = vInfo$phiLock, phi = vInfo$phi, 
            kappaLock = vInfo$kappaLock, kappa = vInfo$kappa, sgn = sgn,
            nesting = nesting, useSparse = useSparse, lOnly = TRUE,
            method = "Nelder-Mead", control = list("maxit" = mInfo$nStep))
        parList = lmm.map(nv, vInfo, sig, phi, kappa, n, tmp$par, TRUE)
        out$l[2] = tmp$value
        if(!mInfo$silent)
        {
            cat("Solution\n")
            cat("   l =", out$l[2], "\n")
        }
        effr = NULL
        if(vInfo$name != "pureNugget" & (!vInfo$phiLock | !vInfo$kappaLock))
        {
            effr = lmm.effectiveRange(vInfo$name, parList$phi, parList$kappa, vInfo$effrMax)
            if(!mInfo$silent)
            {
                if(!vInfo$phiLock & !vInfo$kappaLock)cat("   phi =", parList$phi, "kappa =", parList$kappa, "effr =", effr, "\n")
                if(!vInfo$phiLock & vInfo$kappaLock)cat("   phi =", parList$phi, "effr =", effr, "\n")
				if(vInfo$phiLock & !vInfo$kappaLock)cat("   kappa =", parList$kappa, "effr =", effr, "\n")
			}
        }
        
        #Put in one last call to "lmm.reml" with the final parameter values.
        finalCall = lmm.reml(os, nv, nr, nrTotal, z, Xmat, d, vInfo$name, vInfo$nstr, parList$phi, parList$kappa,
            parList$sig[, , 1], parList$sig[, , 2], parList$sig[, , 3], nesting, useSparse)
        if(!mInfo$silent)
		{
            #Print.
            tmp = lmm.rescale(nv, stdev, Xv, parList$sig, finalCall$b, NA) #Rescale "sig" and "b" for printing purposes.
		    cat("   sig =\n")
			print(tmp$sig)
			cat("   b =\n")
			print(tmp$b)
			cat("\n")
		}
    }else{
        stop("lmm.simplex: Infinite initial likelihood!")
    }
    
    #Fill the output list.
    out$covFun = vInfo$name
	out$phi = parList$phi
	out$kappa = parList$kappa
    out$effr = effr
	out$c0 = parList$sig[, , 1]
	out$c1 = parList$sig[, , 2]
    out$c2 = parList$sig[, , 3]
	out$b = finalCall$b
    out$b_percentVarianceExplained = finalCall$b_percentVarianceExplained
    out$C = finalCall$C
    out$l = finalCall$l
    out
}


lmm.simulater = function(nsim, n0, pred0, Cmat)
{
    ###########################################
    #Simulations on the parameters of the LMM.#
    ###########################################
    #Eigen decomposition of the covariance matrix of (co)kriging errors.
    #(THIS MAY CAUSE STORAGE PROBLEMS WITH LARGE DATATSETS, IN WHICH CASE USE CHOLESKY!)
    eCmat = eigen(Cmat, symmetric = TRUE)
    if(any(eCmat$values <= 0))
    {
        #"Cmat" is not positive definite, so drastic action is needed.
        #Try adjusting it to the nearest positive definite matrix...
        tmp = nearPD(Cmat, keepDiag = TRUE, maxit = 20)
        if(tmp$converged)
        {
            Cmat = tmp$mat
        }else{
            #...that didn't work, so set the off-diagonal elements to zero.
            #This return maximum uncertainty.
            tmp = Cmat * 0
            diag(tmp) = diag(Cmat)
            Cmat = tmp
        }
        eCmat = eigen(Cmat, symmetric = TRUE)
    }
    tmp = Matrix(0, n0, n0, sparse = TRUE)
    diag(tmp) = sqrt(eCmat$values)
    eCmat = eCmat$vectors %*% tmp
    
    #Adjust "pred0" with simulated residuals.
    gaussDev = Matrix(rnorm(n0 * nsim), nrow = n0, ncol = nsim)
    out = as.matrix(eCmat %*% gaussDev + pred0)
    out
}


lmm.stack = function(z, dInfo, trendNames)
{
	############
	#Stack "z".#
	############
	#Split the trend list into unique columns, if required.
	if(any(names(dInfo) == "trend"))
	{
        trendNames = NULL
        if(dInfo$useLinear & !any(names(dInfo$trend) == "X"))
        {
            for(i in dInfo$y)
	    	{
		    	if(any(names(dInfo$trend) == i))
			    {
				    tmp = gsub(" ", "", dInfo$trend[[i]])
                    tmp = gsub("+",  "$", tmp, fixed = TRUE)
				    tmp = gsub("*",  "$", tmp, fixed = TRUE)
				    tmp = strsplit(tmp, "$", fixed = TRUE)[[1]]
				    trendNames = c(trendNames, tmp)
			    }
		    }
        }
        if(!dInfo$useLinear)
        {
            trendNames = NULL
    		for(i in dInfo$y)
	    	{
                trendNames = c(trendNames, attr(dInfo$trend[[i]]$terms, "term.labels"))
            }
        }
        trendNames = unique(trendNames)
        for(i in trendNames)
		{
		    #Error check.
			if(!any(names(z) == i))stop(paste("lmm.stack: Trend variable '", i, "' does not exist in 'z'!", sep = ""))
		}
	}
	
	#Put it all together
	out = list()
	kount = 0
	for(i in dInfo$y)
	{
		kount = kount + 1
		ndx = !is.na(z[[i]])
		if(sum(ndx) < 2)stop(paste("lmm.stack: Data variable '", i, "' has fewer than 2 observations!", sep = ""))
		if(is.null(trendNames))
		{
			tmp = cbind(kount, z[ndx, c(dInfo$coords, i)])
		}else{
			tmp = cbind(kount, z[ndx, c(dInfo$coords, i, trendNames)])
		}
        tmp = data.frame(tmp)
        names(tmp)[c(1, 2 + length(dInfo$coords))] = c("variable", "value")
        if(kount == 1)
		{
			out$z = tmp
			out$nr = sum(ndx)
		}else{
			out$z = rbind(out$z, tmp)
			out$nr = c(out$nr, sum(ndx))
		}
	}
    out$z$variable = factor(out$z$variable, levels = sort(unique(out$z$variable)))
	out$nrTotal = sum(out$nr)
	out
}


lmm.subsample = function(nMax, nv, nr, z, XInfo)
{
    #######################################################################
    #At least one variable has too many observations, so take a subsample.#
    #######################################################################
    #Loop over the variables.
    for(i in 1:nv)
    {
        if(nr[i] > nMax)
        {
            #Exclude one of each pair of points located within a
            #threshold distance of each other.
            ndx1 = which(z$variable == i, arr.ind = TRUE)
            zi = z[ndx1, ]
            d = c(dist(zi[, 2:3]))
            dThresh = min(d)
            ndx2 = rep(TRUE, nr[i])
            while(sum(ndx2) > nMax)
            {
                ndx2 = rep(TRUE, nr[i])
                for(j in 1:(nr[i] - 1))
                {
                    if(ndx2[j])
                    {
                        d = sqrt((zi[(j + 1):nr[i], 2] - zi[j, 2]) ^ 2 + (zi[(j + 1):nr[i], 3] - zi[j, 3]) ^ 2)
                        tmp = j + which(d < dThresh, arr.ind = TRUE)
                        ndx2[tmp] = FALSE
                    }
                }
                if(sum(ndx2) > nMax)dThresh = dThresh + 1.5 * dThresh
            }
            
            #There will be < "nMax" now, so make it up with a random sample.
            set.seed(1234567890) #Ensure the same sample is taken every time.
            tmp = which(!ndx2, arr.ind = TRUE)
            ndx3 = sample(tmp, nMax - sum(ndx2))
            set.seed(NULL) #Reset the seed.
            ndx2[ndx3] = TRUE
            nr[i] = nMax
            zSub = zi[ndx2, ]
            z[ndx1, ] = NA
            z[ndx1[1]:(ndx1[1] + nMax - 1), ] = zSub
            if(any(names(XInfo) == "X"))
            {
                XSub = XInfo$X[ndx2, ]
                X = XInfo$X
                X[ndx1, ] = NA
                X[ndx1[1]:(ndx1[1] + nMax - 1), ] = XSub
                XvSub = XInfo$Xv[ndx2, ]
                Xv = XInfo$Xv
                Xv[ndx1, ] = NA
                Xv[ndx1[1]:(ndx1[1] + nMax - 1), ] = XvSub
            }
        }
    }
    
    #Take out the rows with the missing values.
    ndx = !apply(z, 1, lmm.findNA)
    z = z[ndx, ]
    out = list("dThresh" = dThresh, "z" = z, "nr" = nr, "nrTotal" = sum(nr))
    if(any(names(XInfo) == "X"))
    {
        X = X[ndx, ]
        Xv = Xv[ndx, ]
        out$X = X
        out$Xv = Xv
    }
    out
}


lmm.XandGuesses = function(z, dInfo, vInfo, mInfo, nesting, nv, nr, stdev, rho, X, Xv)
{
    ###############################################
    #Make the design matrix for the fixed effects,#
    #and the initial guesses for the parameters of#
    #the random effects.                          # 
    ###############################################
    #Loop over the variables.
    makePlot = FALSE
    if((is.logical(vInfo$plotIt) && vInfo$plotIt) | !is.logical(vInfo$plotIt))
    {
        makePlot = TRUE
        path = ""
        if(!is.logical(vInfo$plotIt))path = vInfo$plotIt
        if(nv == 1)
        {
            pdf(file = paste(path, "residualVariogram_", vInfo$name, ".pdf", sep = ""), width = 7, height = 7)
        }else{
            pdf(file = paste(path, "residualVariograms_", vInfo$name, ".pdf", sep = ""), width = 7, height = 7)
        }
        par(mfrow = c(nv, 1))
        if(nv > 2)par(mfrow = c(2, 2))
        if(nv > 4)par(mfrow = c(3, 3))
    }
    kount = rep(0, 2)
    vInfo$c0 = matrix(0, nv, nv) #Coregionalisation matrix for the nugget structure.
    vInfo$c1 = matrix(0, nv, nv) #Coregionalisation matrix for the 1st autocorrelated structure.
    vInfo$c2 = matrix(0, nv, nv) #Coregionalisation matrix for the 2nd autocorrelated structure. 
    if(vInfo$phiLock)tmp = vInfo$phi
    vInfo$phi = matrix(0, nrow = nv, ncol = max(1, vInfo$nstr))
	if(vInfo$nstr > 0)
	{
        if(vInfo$phiLock)
        {
            vInfo$phi[, 1] = tmp[1]
            if(vInfo$nstr == 2)vInfo$phi[, 2] = tmp[2]
        }
	}
    treeRules = list()
    for(i in 1:nv)
	{
        if(is.null(dInfo$trend$X))
        {
            #Ordinary least squares model. Retain the design matrix.
            fmla = list("full" = formula("value ~ 1"), "nested" = formula("value ~ 1"))
			if(any(names(dInfo$trend) == dInfo$y[i]))fmla$full = formula(paste("value ~", dInfo$trend[[dInfo$y[i]]]))
		    ndx = z$variable == i
            if(dInfo$useLinear)
            {
                #Trend is described by strictly linear parameters.
                ols = list("full" = lm(fmla$full, data = z, subset = ndx, x = TRUE),
                    "nested" = lm(fmla$nested, data = z, subset = ndx))
            }
            if(i == 1)
			{
                if(dInfo$useLinear)
                {
    			    X = ols$full$x
                }else{
                    tmp = treeSplits(dInfo$trend[[dInfo$y[i]]])
                    treeRules[[dInfo$y[i]]] = treeSplits2rules(tmp)
                    X = treeRules2X(z[ndx, ], treeRules[[dInfo$y[i]]])
                }
				Xv = X * 0 + i
			}else{
                if(dInfo$useLinear)
                {
			        tmp1 = matrix(0, nrow = nrow(ols$full$x), ncol = ncol(X))
				    tmp2 = matrix(0, nrow = nrow(X), ncol = ncol(ols$full$x))
				    X = cbind(rbind(X, tmp1), rbind(tmp2, ols$full$x))
                    Xv = cbind(rbind(Xv, tmp1), rbind(tmp2, ols$full$x * 0 + i))
                }else{
                    tmp = treeSplits(dInfo$trend[[dInfo$y[i]]])
                    treeRules[[dInfo$y[i]]] = treeSplits2rules(tmp)
                    Xi = treeRules2X(z[ndx, ], treeRules[[dInfo$y[i]]])
                    tmp1 = matrix(0, nrow = nrow(Xi), ncol = ncol(X))
				    tmp2 = matrix(0, nrow = nrow(X), ncol = ncol(Xi))
                    X = cbind(rbind(X, tmp1), rbind(tmp2, Xi))
                    Xv = cbind(rbind(Xv, tmp1), rbind(tmp2, Xi * 0 + i))
                }
			}
            if(dInfo$useLinear)
            {
                e = residuals(ols$full)
            }else{
                e = residuals(dInfo$trend[[dInfo$y[i]]])
            }
			if(!is.logical(nesting))e = residuals(ols$nested)
        }else{
            kount[1] = kount[2] + 1
            kount[2] = kount[2] + nr[i]
            Xi = dInfo$trend$X[kount[1]:kount[2], ]
            Xi = Xi[, colSums(Xi) != 0]
            if(!is.matrix(Xi))Xi = matrix(Xi, nrow = length(Xi), ncol = 1)
            ndx = which(colSums(Xi) == nrow(Xi), arr.ind = TRUE)
            if(length(ndx) > 1)
            {
                if(any(Xi != 1))
                {
                    Xi = Xi[, c(1, which(colSums(Xi) != nrow(Xi), arr.ind = TRUE))]
                }else{
                    Xi = matrix(1, nrow = nrow(Xi), ncol = 1)
                }
            }
            ndx = z$variable == i
            b = solve(t(Xi) %*% Xi) %*% t(Xi) %*% z$value[ndx]
            e = z$value[ndx] - Xi %*% b
        }
        if(length(e) >= 20)
        {
            #Plot the variogram of the residuals.
	        #Note that plotting is the original scale.
            tmp = cbind(z[ndx, dInfo$coords], e)
            lags = seq(vInfo$effrMin, vInfo$effrMax, length.out = 20)
            if(lags[1] == 0)lags = lags[-1]
            zv = variog(coords = z[ndx, dInfo$coords], data = e,
                uvec = lags, messages = FALSE)
            tmp = zv$v
            if(all(is.finite(stdev)))tmp = zv$v * stdev[i] ^ 2
            if(makePlot)plot(zv$u, tmp, main = dInfo$y[i], xlab = "Lag", ylab = "Semivariance",
                pch = 19, ylim = c(0, max(tmp)), xlim=c(0, vInfo$effrMax))
		    
            #Superimpose the initial guess.
		    if(vInfo$nstr > 0 & !vInfo$phiLock)
            {
                if(vInfo$nstr == 1)
                {
                    vInfo$phi[i, 1] = zv$u[zv$v == max(zv$v)][1]
                    if(vInfo$name != "spherical")vInfo$phi[i, 1] = zv$u[zv$v == max(zv$v)][1] / 3
                }else{
                    vInfo$phi[i, 2] = zv$u[zv$v == max(zv$v)][1]
                    vInfo$phi[i, 1] = 0.5 * vInfo$phi[i, 2]
                }
            }
            tmp = sort.list(zv$v, decreasing = TRUE)
            tmp = zv$v[tmp]
            if(vInfo$nstr == 0)
            {
                vInfo$c0[i, i] = mean(tmp)
            }else{
                if(vInfo$nstr == 1)
                {
                    vInfo$c1[i, i] = mean(tmp[1:min(length(tmp), 5)]) * 0.5
                    if(vInfo$name != "dampedPeriodic")vInfo$c0[i, i] = vInfo$c1[i, i]
                }else{
                    vInfo$c2[i, i] = mean(tmp[1:min(length(tmp), 5)]) * 0.333
                    vInfo$c1[i, i] = vInfo$c2[i, i]
                    vInfo$c0[i, i] = vInfo$c2[i, i]
                }
            }
            c2 = 0
            if(all(is.finite(stdev)))
            {
                c0 = vInfo$c0[i, i] * stdev[i] ^ 2
		        c1 = vInfo$c1[i, i] * stdev[i] ^ 2
                if(vInfo$nstr == 2)c2 = vInfo$c2[i, i] * stdev[i] ^ 2
            }else{
                c0 = vInfo$c0[i, i]
		        c1 = vInfo$c1[i, i]
                if(vInfo$nstr == 2)c2 = vInfo$c2[i, i]
            }
		    tmp = vInfo$phi[i, ]
            if(vInfo$nstr > 0)
            {
                tmp[1] = max(vInfo$effrMin + 0.001, min(vInfo$phi[i, 1], vInfo$effrMax - 0.001))
                if(vInfo$nstr == 2)tmp[2] = max(vInfo$effrMin + 0.001, min(vInfo$phi[i, 2], vInfo$effrMax - 0.001))
            }
            pred = (c0 + c1 + c2) - lmm.covFun(zv$u, vInfo$name, tmp, vInfo$kappa, c0, c1, c2)
            if(makePlot)lines(zv$u, pred, col = "red")
        }else{
            if(vInfo$name != "pureNugget")stop("lmm: Not enough points for a variogram; use the pureNugget covariance function!")
            vInfo$c0[i, i] = var(e)
            if(makePlot)
            {
                plot(c(0, vInfo$effrMax), rep(vInfo$c0[i, i], 2), main = dInfo$y[i], xlab = "Lag", ylab = "Semivariance", pch = 19, ylim = c(0, 2 * vInfo$c0[i, i]), xlim=c(0, vInfo$effrMax), type = "n")
		        lines(c(0, vInfo$effrMax), rep(vInfo$c0[i, i], 2), col = "red")
            }
        }
	}
    if(makePlot)dev.off()
    if(nv > 1)
	{
		#Take the mean of "vInfo$phi" and "vInfo$kappa".
		if(vInfo$name != "doubleSpherical")
        {
            vInfo$phi = mean(vInfo$phi)
        }else{
            vInfo$phi = colMeans(vInfo$phi)
        }
        vInfo$kappa = mean(vInfo$kappa)
		
        #Put some numbers on the cross-covariances.
        for(i in 1: (nv - 1))
		{
			for(j in (i + 1) :nv)
			{
                for(k in c("c0", "c1", "c2"))
                {
                    vInfo[[k]][i, j] = (sqrt(vInfo[[k]][i, i] * vInfo[[k]][j, j]) * rho[i, j]) * 0.001
                    vInfo[[k]][j, i] = vInfo[[k]][i, j]
                }
			}
		}
	}
    
    #Output list.
    out = list("vInfo" = vInfo)
    if(is.null(dInfo$trend$X))
    {
        out$X = X
        out$Xv = Xv
    }else{
        out$X = dInfo$trend$X
        out$Xv = dInfo$trend$Xv
    }
    out$treeRules = treeRules
    out
}
