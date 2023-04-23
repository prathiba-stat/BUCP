#Only intercepts are estimated
#R code to fit Bayesian unknown change-point model with autocorrelation
#Program begins 
#lines that begin#F with '#' are comments and will not be executed
#install.packages(runjags)
#install runjags if it is not already installed
#this extension of SCD-part2 includes loops

setwd("~/Desktop/Research/Bayes talk for small is beautiful")
#change the working directory

source("plots.R")
source("plot_SSD.R")
source("plotPost.R")
source("HDIofMCMC.R") 
source("openGraphSaveGraph.R") 
source("plotPost_upper.R")

datafl <- "SCD1-part2.csv"
postfl <- "SCD1-part2-Posterior-ES"
outfl  <- "SCD1-par2-results.csv"

library(runjags)

massage.data <- function(data){
  #data are fed in two columns, the first column is the DV and 
  #the second column in the phase number 1 or 2
  first <- data[data[, 2]==1, 1]
  second <- data[data[, 2]==2, 1]
  if (length(first)>length(second)){
    second <- c(second, rep(NA, length(first) - length(second))) 
    } else {
    first <- c(first, rep(NA, length(second) - length(first)))
    } 
  return(rbind(first, second))
  
}
#If there are 2 datasets all in sheets 1-2 in an xlsx spreadsheet
#ndata <- 1
#change the number of rows below depending on the number of datasets
#sum.results <- matrix(NA, ndata, 20) 
sum.results <- matrix(NA, 1, 20) 
colnames(sum.results) <- c("mean1", "mean2", "beta[1,1]", "beta[2,1]", "ES",
                           "sigma.epsilon", "rho", "cp",
                           "l.95b1", "u.95b1",
                           "l.95b2", "u.95b2", 
                           "l.95es", "u.95es", 
                           "l.95sigma", "u.95sigma",
                           "l.95rho", "u.95rho",
                           "l.95cp", "u.95cp")
#Read data from a file
#for (i in 1:ndata){
data <- read.csv("SCD1-part2.csv")
#data <- read.csv(file=datafl)
plot_SSD(massage.data(data), "data.jpg")
y <- data[,1]
mean1 <- mean(data[data[,2]==1, 1])
mean2 <- mean(data[data[,2]==2, 1])

T <- length(y) 
#Total number of time-points in the data 
P <- 2 
#two phases - baseline and treatment phase
pi <- c(rep(0, 2), rep(1/(T-4), (T-4)), rep(0, 2)) 
#probabilities for the categorical prior for change-point
#i.e. pi will equal (0, 0, 1/16, ., 1/16, 0, 0) because the change-point
#can only be between time point 3 and T-2
b1 <- mean(y[1:(T/2)]) 
#mean of the distribution. This is used to generate starting values for #intercept 1 when the model is fitted
b2 <- mean(y[((T/2)+1):T])  
#mean of the distribution. This is used to generate starting values for #intercept 2 when the model is fitted
#model definition begins here, the model is denoted using variable called 
#UCP.model
UCP.model <- "model {
                x[1] <- 0 
                #x is a vector of means of distributions from which the intercepts are
                #drawn
                yhat[1] <- beta[1, 1] 
                #beta[1, 1] is the mean of the distribution from which intercept 1 is
                #drawn
                y[1] ~ dnorm(yhat[1], tau.epsilon) 
                #y[1] is generated from normal distribution with mean yhat[1] and 
                #precision tau.epsilon, tau.epsilon is the inverse of the variance in 
                #equation 1
                  for (i in 2:T) {
                      dummy[i] <- step(cp - i) 
                # dummy = 1 if i in the baseline phase,  dummy = 0 if i in the
                #intervention phase
                      x[i] <- dummy[i] * beta[1, 1] + (1 - dummy[i]) * beta[2, 1] 
                #identifies if x for the ith time-point should be centered 
                #around intercept 1 or intercept 2
                      yhat[i] <- ifelse(cp == i + 1, x[i], x[i] + rho * (y[i - 1] - x[i 
                                  - 1])) 
                #if i is the first time-point of phase 2 no autocorrelation
                      y[i] ~ dnorm(yhat[i], tau.epsilon)
                #y[i] is generated from normal distribution with mean yhat[i] and
                #precision tau.epsilon
                  }
                #Assigning priors 
                for (i in 1:P){
                  beta[i, 1] ~ dnorm(mu[i], prec[i])
                #intercept[i] is drawn from a normal distribution with mean mu[i]
                #and precision prec[i]
                  mu[i] ~ dnorm(0, .0001)
                # mu is drawn from a normal distribution with mean 0 and 
                #precision .0001 (i.e. variance 10000)
                  prec[i] ~ dgamma(1, 1)
                #precision is drawn from a gamma distribution with shape 1
                }
  sigma.epsilon ~ dunif(0.1, 5)
                # standard deviation of y within each phase is drawn from a uniform
                #distribution ranging from 0.1 to 5
  tau.epsilon <- pow(sigma.epsilon, -2)
                #tau.epsilon is the precision of y obtained from sigma.epsilon
  es <- (beta[1, 1] - beta[2, 1])/sigma.epsilon
  
  cp ~ dcat(pi)
                #prior of the change-point is categorical with probabilities as defined
                #in the pi vector on line 6
  rho ~ dunif(-1, 1)
                #autocorrelation is drawn from a uniform distribution ranging from -1
                #to 1
  }"
#end of model definition
#run jags to estimate the parameters
#autorun.jags calls runjags to run the model in JAGS until convergence
results <- autorun.jags(
  model = UCP.model,
  #the model defined in the previous block is assigned as the model to be
  #fitted
  data = list(y = y, T = T, P = P, pi = pi),
  #data input: observations, #time-points, #phases, #priors for change-point
  monitor = c("beta", "es", "sigma.epsilon", "rho", "cp"),
  #parameters to monitor
  n.chains = 4,
  #run 4 chains
  startsample = 50000,
  #burn-in the first 50000
  inits = function() { #initial values are generated within this function
    list(
      beta = rbind(rnorm(1, b1, 1), rnorm(1, b2, 1)),
      sigma.epsilon = runif(1, 0.5, 2),
      rho = runif(1, -1, 1),
      cp = sample(3:(T-2), 1)
    )
  }, 
  method = "rjparallel"
  #run parallel chains
)
#end of running jags  
# combine all chains into a single chain for convenience
results$draws <- round(combine.mcmc(results$mcmc), digits = 4)
results #object that contains all posteriors
results.tab <- summary(results)
##############################
#####OVERLAP
length(which(results$draws[,"beta[2,1]"] > min(results$draws[,"beta[1,1]"])))
######################################################
#######EFFECT SIZE#####################
######################################################
lims <- c(min(c(results$draws[,"beta[1,1]"], results$draws[,"beta[2,1]"])),
          max(c(results$draws[,"beta[1,1]"], results$draws[,"beta[2,1]"])))
#pdf(file = paste0("Posterior-plots-es", i, ".pdf"))
pdf(file = paste0(postfl,".pdf"))
par(mar=c(3.5,3.5,2,1),mgp=c(2,0.7,0), mfcol = c(3, 2))
plot(density(results$draws[,"beta[1,1]"]), main = "Intercept Phase 1",
     xlab = expression(beta[1]), ylab =" ", xlim = lims)
abline(v=results.tab["beta[1,1]", c(1, 3)], lty = "dashed")
plot(density(results$draws[,"beta[2,1]"]), main = "Intercept Phase 2",
     xlab = expression(beta[2]), ylab =" ", xlim = lims)
abline(v=results.tab["beta[2,1]", c(1, 3)], lty = "dashed")

###############################################################
#If you want to change the "significant ES do it in this block
###############################################################

##############Small ROPE with 1 SD ABOVE 0
#in the plot below the mid-point of comparison for ROPE is 0.5
#the interval ranges from 0.5-0.5 to 0.5 +0.5
#change this value based on the hypothesis for ROPE
lims1 = c(min(results$draws[,"es"]), max(results$draws[,"es"]))
plots(results$draws[,"es"], compVal = 2, ropeRad = 0.2, "Effect size",
      lims1 = lims1)
hist(results$draws[, "cp"], main = "Change-point", xlab = "change-point", 
     ylab = " ", yaxt = 'n')
plot(density(results$draws[,"rho"]), main = "Autocorrelation",
     xlab = "rho", ylab =" ")
plot(density(results$draws[,"sigma.epsilon"]), main = "Standard deviation",
     xlab = "SD", ylab =" ")
dev.off()

#plot(results)
#all estimates in one row
sum.results[1,] <- c(mean1, mean2, results.tab[-6, "Mean"], results.tab["cp", "Mode"],
                     c(t(results.tab[, c("Lower95", "Upper95")])))
#}
#sum.results

#write.csv(sum.results, "results-es.csv")
write.csv(t(sum.results), outfl)

