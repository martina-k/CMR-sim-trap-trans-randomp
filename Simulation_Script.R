# CMR-sim-trap-trans-randomp
#Simulation of CMR data with transcience, trap-dep and random recapture rates
#Simulation using multi-state approach, to be analyzed using single-state-models

rm(list = ls())
setwd("/Users/Martina/Documents/Synchrony Project/Shag_analysis/Simulation/")

library(jagsUI)
library(gdata)
library(R2ucare)

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
  # Unobservable: number of state that is unobservable
  n.occasions <- dim(PSI.STATE)[4] + 1
  CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
  g <- colSums(marked)
  for (s in 1:dim(PSI.STATE)[1]){
    if (g[s]==0) next  # To avoid error message if nothing to replace
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
  } #s
  for (i in 1:sum(marked)){
    for (s in 1:dim(PSI.STATE)[1]){
      if (mark.occ[i,s]==0) next
      first <- mark.occ[i,s]
      CH[i,first] <- s
      CH.TRUE[i,first] <- s
    } #s
    for (t in (first+1):n.occasions){
      # Multinomial trials for state transitions
      if (first==n.occasions) next
      state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
      CH.TRUE[i,t] <- state
      # Multinomial trials for observation process
      event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
      CH[i,t] <- event
    } #t
  } #i
  # Replace the NA and the highest state number (dead) in the file by 0
  CH[is.na(CH)] <- 0
  CH[CH==dim(PSI.STATE)[1]] <- 0
  CH[CH==unobservable] <- 0
  id <- numeric(0)
  for (i in 1:dim(CH)[1]){
    z <- min(which(CH[i,]!=0))
    ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
  }
  return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
  # CH: capture histories to be used
  # CH.TRUE: capture histories with perfect observation
}


# Starting from Kery & Schaub BPA book, example 9.4.2. 
# Generation of simulated data
# Define mean survival, recapture, as well as number of occasions, states, observations and released individuals
phi.j <- 0.4
phi.ad <- 0.8
p.seenbef <- 0.8
p.nseenbef <- 0.3
n.occasions <- 7  
n.states <- 4 #(Juvenile, Adult-seen year before, Adult-not seen year before, Dead)
n.obs <- 3 #(Seen as juv, Seen as ad, Not seen)
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)	# Releases only as juveniles

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      0, phi.j,   0,        1-phi.j,
      0, phi.ad,  phi.ad,   1-(2*phi.ad),
      0, phi.ad,  phi.ad,   1-(2*phi.ad),
      0, 0,       0,        1), nrow = n.states, byrow  = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      0, 0,         1,
      0, p.seenbef, 1-p.seenbef,
      0, p.nseenbef,1-p.nseenbef,
      0, 0,         1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed!
# 1 = seen as juv, 2 = seen, 3 = not seen
rCH <- CH
rCH[rCH==0] <- 3

# Analysis of the model
# Specify model in BUGS language
sink("trans_trap.jags")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi.juv: first year survival probability
    # phi.ad: adult survival probability
    # p.NO: recapture probability of non-seen-year-before
    # p.O: recapture probability of seen-year-before
    # -------------------------------------------------
    # States (S):
    # 1 juvenile
    # 2 ad, seen the year before
    # 3 ad, not seen the year before
    # 4 dead
    # Observations (O):
    # 1 seen as juvenile
    # 2 seen as adult
    # 3 not seen
    # -------------------------------------------------
    
    # Priors and constraints	
    for (t in 1:(n.occasions-1)){
    phi.juv[t] <- mean.phij
    phi.ad[t] <- mean.phiad
    p.NO[t] <- mean.pnonobs
    p.O[t] <- mean.pobs
    }
    mean.phij ~ dunif(0, 1)     # Prior for mean juv survival
    mean.phiad ~ dunif(0, 1)    # Prior for mean ad survival
    mean.pnonobs ~ dunif(0, 1)  # Prior for mean recapture non-obs
    mean.pobs ~ dunif(0, 1)     # Prior for mean recapture obs year before
    
    # Define state-transition and observation matrices 	
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- 0
    ps[1,i,t,2] <- phi.j[t]
    ps[1,i,t,3] <- 0
    ps[1,i,t,4] <- 1-phi.j[t]
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- phi.ad[t] 
    ps[2,i,t,3] <- phi.ad[t]
    ps[2,i,t,4] <- 1-2*phi.ad[t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- phi.ad[t]
    ps[3,i,t,3] <- phi.ad[t]
    ps[3,i,t,4] <- 1-2*phi.ad[t]
    ps[4,i,t,1] <- 0
    ps[4,i,t,2] <- 0
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- 0
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 1
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p.O[t]
    po[2,i,t,3] <- 1-p.O[t]
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- p.NO[t]
    po[3,i,t,3] <- 1-p.NO[t]
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 1
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()


# Initial values
# We have to provide initial values for the true latent state z. Since the observed state numbers do not always match the true state numbers, we cannot use simple the observed multistate capture history as initial value for z. The following function creates initial values for z from our observation under the model we want to analyse.

agefirst.init <- function(ch, f){
  age <- array(NA, dim=dim(ch))
  for (i in 1:nrow(ch)){
    for (t in f[i]:ncol(ch)){
      age[i,t] <- min(c(t-f[i]+1, 4))
    }
  }
  ini <- array(NA, dim=dim(ch))
  for (i in 1:nrow(ch)){
    for (t in f[i]:ncol(ch)){
      if(ch[i,t]==1) ini[i,t] <- 1
      if(ch[i,t]==2&age[i,t]==2) ini[i,t] <- 2
      if(ch[i,t]==3&age[i,t]==2) ini[i,t] <- 4
      if(ch[i,t]==3&age[i,t]==3) ini[i,t] <- 4
      if(ch[i,t]==4&age[i,t]==4) ini[i,t] <- 4
      if(ch[i,t]==2&age[i,t]==3) ini[i,t] <- 3
    }
  }
  ini[which(is.na(ini))] <- age[which(is.na(ini))]
  for (i in 1:nrow(ch)){
    ini[i,f[i]] <- NA
    for (t in f[i]:(ncol(ch)-1)){
      if(ini[i,t]==4 & ini[i,t+1]==3) ini[i,t+1] <- 4
    }
  }
  return(ini)
}

# Add-in for JAGS
# Since we have created a matrix with initial values for the true state z, we can use part of this information as data (see chapter 7.3.1) which can help with convergence and computing time). Here we give those initial values that are based on an actual observation. Since the first observation is deterministic, it must be excluded. The following code constructs the data matrix:

ch <- rCH
ch[ch==3] <- NA
z.known <- agefirst.init(rCH, f)
z.known[is.na(ch)] <- NA
for (i in 1:nrow(ch)){
  z.known[i,f[i]] <- NA
}
z <- agefirst.init(rCH, f)
z[!is.na(z.known)] <- NA

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = z.known)

# Initial values
inits <- function(){list(mean.phij = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), 
                         mean.pNO = runif(1, 0, 1), mean.pO = runif(1, 0, 1), 
                         z = z)}  

# Parameters monitored
parameters <- c("mean.phij", "mean.phiad",  "mean.pNO", "mean.pO")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 2 min)
transtrap  <- jags(jags.data, inits, parameters, "trans_trap.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)

#Change code below to match model name and JAGS
par(mfrow = c(3, 3), las = 1)
hist(agefirst$BUGSoutput$sims.list$mean.phi1, col = "gray", main = "",  xlab = expression(phi[1]))
abline(v = phi.1, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.phi2, col = "gray", main = "",  xlab = expression(phi[2]), ylab = "")
abline(v = phi.2, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.phiad, col = "gray", main = "",  xlab = expression(phi[ad]) , ylab = "")
abline(v = phi.ad, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.alpha1, col = "gray", main = "",  xlab = expression(alpha[1]))
abline(v = alpha.1, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.alpha2, col = "gray", main = "",  xlab = expression(alpha[2]) , ylab = "")
abline(v = alpha.2, col = "red", lwd = 2)
plot(0, type = "n", axes = F, ylab = "", xlab = "")
hist(agefirst$BUGSoutput$sims.list$mean.pNB, col = "gray", main = "",  xlab = expression(p[NB]))
abline(v = p.NB, col = "red", lwd = 2)
hist(agefirst$BUGSoutput$sims.list$mean.pB, col = "gray", main = "",  xlab = expression(p[B]) , ylab = "")
abline(v = p.B, col = "red", lwd = 2)





####################
