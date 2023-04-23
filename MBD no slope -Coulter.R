##2023
setwd("~/Desktop/Research/Bayes talk for small is beautiful") #change your working directory here

library(runjags)
library(ggplot2)
library(rjags)

source('plots.R')
source('plot_MBD.R')
source('plot_SSD.R')
source("plotPost.R")
source("plotPost_upper.R")
source("HDIofMCMC.R") 
source("openGraphSaveGraph.R") 
#############FUNCTIONS###########################
pad <- function(a, b){
  if (length(a) > length(b)){
    temp = length(a) - length(b)
    b <- c(b, rep(NA, temp))
  } else {
    temp <- length(b) - length(a)
    a <- c(a, rep(NA, temp))
  }
  return(rbind(a, b))
}

##############################################
#data <- read.csv("Coulter_MBD.csv")
#data <- read.csv("C:/Users/Ratna/Dropbox/AERA2018Study/Data/Barber-Initiations-MBD.csv")
#data <- read.csv("C:/Users/Ratna/Dropbox/AERA2018Study/Data/Barber-Responses-MBD.csv")
data <- read.csv("Coulter_MBD.csv")

names <- c("Child1", "Child2", "Child3") #change the names here
dv <- c("Peer-Initiations") #change the DV here
n.cases <- 3
plotname <- "Initiations-ssd-plot.jpg"
############PLOT##############################################

yplot.1 <- pad(data[data[,1]==1, 2], 
               data[data[,1]==2, 2])
yplot.2 <- pad(data[data[,3]==1, 4], 
               data[data[,3]==2, 4])
yplot.3 <- pad(data[data[,5]==1, 6], 
               data[data[,5]==2, 6])
m <- max(ncol(yplot.1), ncol(yplot.2), ncol(yplot.3))
yplot <- rbind(cbind(yplot.1, matrix(NA, 2, m - ncol(yplot.1))), 
               cbind(yplot.2, matrix(NA, 2, m - ncol(yplot.2))), 
               cbind(yplot.3, matrix(NA, 2, m - ncol(yplot.3))))
lims <- nrow(data)
plot_MBD(yplot, names, dv, lims, plotname)
#######################END PLOT##################

sum.results <- data.frame(t(rep(NA, 28))) 
colnames(sum.results) <- c("data", "mean1", "mean2", "beta[1,1]", "beta[2,1]", "ES",
                           "sigma.epsilon", "rho", "cp", "p.calc",
                           "l.95b1","median.b1" , "u.95b1",
                           "l.95b2", "median.b2", "u.95b2", 
                           "l.95es", "median.es", "u.95es", 
                           "l.95sigma", "median", "u.95sigma",
                           "l.95rho", "median", "u.95rho",
                           "l.95cp", "median", "u.95cp")

###################MODELS###########
#suggest not using the UCP model with slope & autocorrelation
UCP.model.slope <- "model{
  x[1] <- 0 
  yhat[1] <- beta[1, 1] 
  y[1] ~ dnorm(yhat[1], tau.epsilon) 
  
  for (i in 2:T) {
    dummy[i] <- step(cp - i) 
    
    x[i] <- (dummy[i] * (beta[1,1] + (beta[1,2] * (i - 1)))) + 
            ((1 - dummy[i]) * (beta[2,1] + (beta[2,2]*(cp-i-1))))
    
    yhat[i] <- ifelse(cp == i + 1, x[i], 
                      x[i] + rho * (y[i - 1] - x[i - 1])) 

    y[i] ~ dnorm(yhat[i], tau.epsilon)
    
  }
  for (i in 1:P){
    beta[i,1] ~ dnorm(mu[i], prec[i])
    mu[i] ~ dnorm(30, 0.1)
    prec[i] ~ dgamma(1, 1)
    beta[i,2] ~ dnorm(mu.s[i], prec.s[i])
    mu.s[i] ~ dnorm(30, 0.1)
    prec.s[i] ~ dgamma(1, 1)

  }
  sigma.epsilon ~ dunif(0.1, 5)
  tau.epsilon <- pow(sigma.epsilon, -2)
  es <- (beta[1, 1] - beta[2, 1])/sigma.epsilon
  
  cp ~ dcat(pi)
  
  rho ~ dunif(-1, 1)
  
}"
########################################
##no slopes model
########################################

UCP.no.slope.model <- "model{
  x[1] <- 0 
  yhat[1] <- beta[1, 1] 
  y[1] ~ dnorm(yhat[1], tau.epsilon) 

  for (i in 2:T) {
    dummy[i] <- step(cp - i) 

    x[i] <- dummy[i] * beta[1, 1] + (1 - dummy[i]) * beta[2, 1] 

    yhat[i] <- ifelse(cp == i + 1, x[i], 
    x[i] + rho * (y[i - 1] - x[i - 1])) 

    y[i] ~ dnorm(yhat[i], tau.epsilon)

  }
  for (i in 1:P){
    beta[i, 1] ~ dnorm(mu[i], prec[i])
    mu[i] ~ dnorm(30, 0.1)
    prec[i] ~ dgamma(1, 1)
  }
  sigma.epsilon ~ dunif(0.1, 5)
  tau.epsilon <- pow(sigma.epsilon, -2)
  es <- abs((beta[1, 1] - beta[2, 1]))/sigma.epsilon

  cp ~ dcat(pi)

  rho ~ dunif(-1, 1)

}"
###############END MODELS############

results.tab <- list()
results <- list()
p.calc.int <- list()
p.calc.slope <- list()

for (i in 1:n.cases){
  y <- data[complete.cases(data[,i*2]), i*2]
  mean1 <- mean(data[data[,(2*i - 1)]==1, 2*i], na.rm = T)
  mean2 <- mean(data[data[,(2*i - 1)]==2, 2*i], na.rm = T)
  T <- length(y) 
  P <- 2 
  pi <- c(rep(0, 2), rep(1/(T-4), (T-4)), rep(0, 2)) 
  b1 <- mean(y[1:(T/2)]) 
  b2 <- mean(y[((T/2)+1):T])  
  results[[i]] <- autorun.jags(
    model = UCP.no.slope.model,
    data = list(y = y, T = T, P = P, pi = pi),
    monitor = c("beta", "es", "sigma.epsilon", "rho", "cp"),
    n.chains = 4,
    startsample = 50000,
    inits = function() { 
      list(
        beta = rbind(rnorm(1, b1, 1), rnorm(1, b2, 1)),
        sigma.epsilon = runif(1, 0.5, 2),
        rho = runif(1, -1, 1),
        cp = sample(3:(T-2), 1)
      )
    }, 
    method = "rjparallel"
  )

results[[i]]$draws <- round(combine.mcmc(results[[i]]$mcmc), digits = 4)
results.tab <- rbind(results.tab, summary(results[[i]]))
p.calc.int[i] <- length(intersect(results[[i]]$draws[,"beta[1,1]"], results[[i]]$draws[,"beta[2,1]"]))/
  length(results[[i]]$draws[,"beta[1,1]"])

pdf(file = paste0(i,"-MBD-posteriors.pdf"))
par(mar=c(3.5,3.5,2,1),mgp=c(2,0.7,0), mfcol = c(3, 2))
lims <- c(min(c(results[[i]]$draws[,"beta[1,1]"], results[[i]]$draws[,"beta[2,1]"])),
          max(c(results[[i]]$draws[,"beta[1,1]"], results[[i]]$draws[,"beta[2,1]"])))
plot(density(results[[i]]$draws[,"beta[1,1]"]), main = "Intercept Phase 1",
     xlab = expression(beta[1]), ylab =" ", xlim = lims)
abline(v=summary(results[[i]])["beta[1,1]",c(1, 3)], lty = "dashed")
plot(density(results[[i]]$draws[,"beta[2,1]"]), main = "Intercept Phase 2",
     xlab = expression(beta[2]), ylab =" ", xlim = lims)
abline(v=summary(results[[i]])["beta[2,1]",c(1, 3)], lty = "dashed")

###############################################################
#If you want to change the "significant ES do it in this block
###############################################################

##############Small ROPE with 1 SD ABOVE 0
#in the plot below the mid-point of comparison for ROPE is 0.5
#the interval ranges from 0.5-0.5 to 0.5 +0.5
#change this value based on the hypothesis for ROPE
lims1 = c(min(results[[i]]$draws[,"es"]), max(results[[i]]$draws[,"es"]))
plots(results[[i]]$draws[,"es"], compVal = 2, ropeRad = 0.2, "Effect size",
      lims1 = lims1)
hist(results[[i]]$draws[, "cp"], main = "Change-point", xlab = "change-point", 
     ylab = " ", yaxt = 'n')
plot(density(results[[i]]$draws[,"rho"]), main = "Autocorrelation",
     xlab = "rho", ylab =" ")
plot(density(results[[i]]$draws[,"sigma.epsilon"]), main = "Standard deviation",
     xlab = "SD", ylab =" ")
dev.off()

}

write.csv(results.tab, "Barber-initiations.csv")

