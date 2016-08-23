#############################
### Code for fitting LMCR ###
#############################

# we have simplified landuse to irrigated and non-irrigated
# one landuse only 
# - if only 2015 sampled - 2015 LU used
# - if only 2002 sampled - 2002 LU used
# - if both years sampled - 2015 LU used

#Read functions
source('lmm_fun.R')

#z=read.csv("restacked_merged.csv")
z=read.csv("C:/Users/bishop/Dropbox/Patrick/Data/final.csv")
str(z)
names(z)

# select 0-10 cm only
z=z[which(as.character(z$depth)=="0-10 cm"),]

names(z)

#The possible components of list "dInfo" are:                                                                       #
# - "y" = the column names of the response variables of "z".                                                        #
# - "coords" = the columns names of the spatial coordinates of "z".                                                 #
# - "scale" = logical to scale each response variable by it standard deviation.                                     #
# - "trend" = a list with sub-components:                                                                           #
#               - for each variable in "y", a string representing the formula needed; OR,                           #
#               - for each variable in "y", an object of class "rpart"; OR,                                         #
#               - "X" (a custom design matrix), and "Xv" (a corresponding matrix that indexes the variables in "X". #
#               - "nesting" is the column numbers of the design matrix to zero, if nesting.   

dInfo = list ("y" = c("pH_02","pH_15"),"coords" = c("Eastings", "Northings"),
              "trend"=  list("pH_02" = "Landuse+PC1.x+PC5.x",
                             "pH_15" = "Landuse+Geo_3+PC1.x+PC2.x+PC5.x"))

#The possible components of list "vInfo" are:                                                                       #
# - "name" = string that names the covariance function to fit.                                                      #
# - "effrMin" = numeric for the smallest value of the effective range that can be taken by the covariance function. #
# - "effrMax" = numeric for the largest value of the effective range that can be taken by the covariance function.  #
# - "phiLock" = logical to lock distance parameter "phi".                                                           #
# - "kappaLock" = logical to lock curvature parameter "kappa".                                                      #
# - "cvci" = logical to bootstrap a confidence interval about the theta statistics.                                 #
# - "plotIt" = plot the experimental variogram(s) and the initial guess on the covariance function. Can be either a #
#              logical (to save the plot into the default directory) or a string specifying the desired directory.  #
#                                             

##vInfo = list("name" = "doubleSpherical", "effrMin" = 0.2, "effrMax" = 8, "phiLock" = TRUE)
##vInfo = list("name" = "exponential", "effrMin" = 0.2, "effrMax" = 8, "phiLock" = TRUE)
##vInfo = list("name" = "matern", "effrMin" = 1000, "effrMax" = 80000)
##vInfo = list("name" = "matern", "effrMin" = 0.2, "effrMax" = 8, "phiLock" = TRUE)
##vInfo = list("name" = "matern", "effrMin" = 0.2, "effrMax" = 8, "kappaLock" = TRUE, "kappa" = 2)

vInfo = list("name" = "exponential",
             "effrMin" = 30,
             "effrMax" = 80000,
             "phiLock" = FALSE)

##vInfo = list("name" = "pureNugget", "effrMin" = 0.2, "effrMax" = 8)
##vInfo = list("name" = "spherical", "effrMin" = 0.2, "effrMax" = 8, "phiLock" = TRUE)


#The possible components of list "mInfo" are:                                                                       #
# - "method" = string to indicate the minimisation method ("lbfgsb", "sanneal" or "simplex").                       #
# - ???                                                                                                             #
# - "ncore" = numeric for number of cores to use in the analysis.                                                   #
# - "silent" = logical to print output as you go.                                                                   #
# - "method" = minimisation method ("lbfgsb" or "sanneal").                                                         #
# - "nStep" = used with "sanneal"; the number of cooling steps in the minimisation.                                 #
# - "nMarkov" = used with "sanneal"; the number of Markov chains per cooling step.                                  #
# - "heat" = used with "sanneal"; the initial temperature of the system.                                            #
# - "alpha" = used with "sanneal"; the proportion by which the system is cooled at each step.                       #
# - "nStop" = used with "sanneal"; if parameters aren"t changed, the number cooling steps needed before termination.#

mInfo = list("method" = "sanneal",
             "silent" = FALSE,
             "nStep" = 2000,
             "nMarkov" = 50,
             "heat" = 2000,
             "alpha" = 0.98,
             "nStop"=20)

##################
##Pat's comments## 
##################
#"mInfo"#                                                                                  
# - "method" = use "sanneal" always 
# - "nStep" = the number of iterations, the number of cooling steps in the minimisation 
# - "nMarkov" = the number of Markov chains per cooling step                                  
# - "heat" = the initial temperature of the system   
#            increasing the temperature of the system makes the accpetance of a particular transition more likely 
# - "alpha" = the proportion by which the system is cooled at each step, 
#             the proportion of legal changes accepted
#             the rate of acceptance of changes               
#             metropolis criterion
# - "nStop" = if parameters aren"t changed, the number cooling steps needed before termination

# Initial guess at the model parameters & cooling schedule: 
# c1-the initial system temperature, 
# ??c, 
# Nm-the number of perturbations of each model parameter in a single Markov chain 
# Nt-a threshold number of Markov chains so that, if there is no change in WSS over Nt successive chains, the algorithm will stop.
# examples from lark and papritz paper:
# ??c=0.975 
# Nm=25 
# Nt=15
# ????,max=1.0
# c1 = was selected by trial and error so that 99% of the perturbations were accepted in the initial Markov chain

# nStep = Nm
# nMarkov = Nt
# heat = c1
# alpha = ??c
# nStop = NA 

model = lmm(z, dInfo, vInfo, mInfo)

#pH_11=model

# FINAL MODEL
##############################################################################################
pH_11$zStk
pH_11$cvSummary
pH_11$model
# plot obs vs pred
pH_11$cvSummary$cvPred
plot(pH_11$cvSummary$cvPred$observed,pH_11$cvSummary$cvPred$prediction)
abline(0,1)
# matrix 
pH_11_matrix=pH_11$model$C
pH_11_matrix
#############################   Extracting p values from model ###############################
datas=z
coefficients <- pH_11$model$b 
se_error <- sqrt(diag(pH_11$model$C))  # variance-covariance matrix
t_value <- coefficients/se_error## get t values
t_prob <- 2 * pt(-abs(t_value), df = (nrow(datas) -length (coefficients)))## and probabilities
## make pretty
coef_mat <- cbind(coefficients, se_error, t_value, t_prob)
colnames(coef_mat) <- c("Estimate", "Std.Err","t value", "Pr(>|t|)")
# order of the trend
trend1=~Landuse+PC1.x+PC5.x #need to populate with names from trend
trend2=~Landuse+Geo_3+PC1.x+PC2.x+PC5.x
trend=c(colnames(model.matrix(trend1, z)),colnames(model.matrix(trend2, z)))#model.matrix does n-1 for categories, not sure how it works
trend #bind two lots of trend, because bivariate model^^
rownames(coef_mat)=trend
printCoefmat(coef_mat)
##############################################################################################
###### CI's to test for overlap: ---------------
tcrit<-abs(qt(0.025,lower.tail=F,df=(nrow(datas) -length (coefficients))))
U95<-coefficients+(tcrit*se_error)
L95<-coefficients-(tcrit*se_error)
#
coef_mat <- cbind(coefficients,L95, U95, se_error, t_value, t_prob)
colnames(coef_mat) <- c("Estimate", "Lower 95CI", "Upper 95CI","Std.Err","t value", "Pr(>|t|)")
# order of the trend
trend1=~Landuse+PC1.x+PC5.x #need to populate with names from trend
trend2=~Landuse+Geo_3+PC1.x+PC2.x+PC5.x
trend=c(colnames(model.matrix(trend1, z)),colnames(model.matrix(trend2, z)))#model.matrix does n-1 for categories, not sure how it works
trend #bind two lots of trend, because bivariate model^^
rownames(coef_mat)=trend
printCoefmat(coef_mat, signif.stars = T)
##############################################################################################

# plot this to see if the same covariates overlap each other in different trends
# if they overlap, we can have it in a single trend
plot(coef_mat[c(2:4,6:11),3],xaxt="n",ylim=c(-1,1))
points(coef_mat[c(2:4,6:11),2],add=T)
labels=c("LanduseNon-irrigated","PC1.x","PC5.x","LanduseNon-irrigated","Geo_3Other","Geo_3Red","PC1.x","PC2.x","PC5.x")
axis(side=1,at=c(1:9),labels=labels,cex.axis=0.7)
grid()
# PC1 overlaps
# PC5 overlaps

save(pH_11,file="pH_finalmodel.RData")