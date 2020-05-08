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


# Starting from Kery & Schaub BPA book, example 9.3.3 and 9.4.2. 
# Generation of simulated data
# Define mean survival, recapture, as well as number of occasions, states, observations and released individuals
phi.j <- 0.4
phi.ad <- 0.8 #Different survival rates to incorporate transcience
hide <- 0.2
hideagain <- 0.6
p.seenbef <- 0.8 # = 1 - hide
p.nseenbef <- 0.4 # = 1 - hideagain
n.occasions <- 7  
n.states <- 4 #(Juvenile, Adult-observable, Adult-hiding (unobservable), Dead)
n.obs <- 3 #(Seen as juv, Seen as ad (Not hiding), Not seen)
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
      0, phi.j,                 0,                  1-phi.j,
      0, phi.ad*(1-hide),       phi.ad*hide,        1-phi.ad,
      0, phi.ad*(1-hideagain),  phi.ad*hideagain,   1-phi.ad,
      0, 0,                     0,                  1), nrow = n.states, byrow  = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      0, 0,         1,
      0, 1,         0,
      0, 0,         1, #Unobservable state
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
# Specify model - with constant phi:s and p:s
sink("trans_trap_c-c.jags")
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
    ps[1,i,t,2] <- phi.juv[t]
    ps[1,i,t,3] <- 0
    ps[1,i,t,4] <- 1-phi.juv[t]
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- phi.ad[t] * p.O[t]
    ps[2,i,t,3] <- phi.ad[t] * (1-p.O[t])
    ps[2,i,t,4] <- 1-phi.ad[t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- phi.ad[t] * p.NO[t]
    ps[3,i,t,3] <- phi.ad[t] * (1-p.NO[t])
    ps[3,i,t,4] <- 1-phi.ad[t]
    ps[4,i,t,1] <- 0
    ps[4,i,t,2] <- 0
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- 0
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 1
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- 1
    po[2,i,t,3] <- 0
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- 1
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


# Function from BPA book for creating known latent states z
known.state.ms <- function(ms, notseen){
  # notseen: label for not seen
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}

known.state.cjs <- function(ch){ state <- ch
for (i in 1:dim(ch)[1]){
  n1 <- min(which(ch[i,]==1)) 
  n2 <- max(which(ch[i,]==1)) 
  state[i,n1:n2] <- 1 
  state[i,n1] <- NA
}
state[state==0] <- NA 
return(state)
}

#Trying to create function for known latent states
known.state.trap <- function(ch, notseen){
  state <- ch
  state[state==notseen] <- NA
  for (i in 1:dim(ch)[1]){
    m1 <- min(which(!is.na(state[i,])))
    m2 <- max(which(!is.na(state[i,])))
    if(m1<m2){
    for (j in (m1+1):m2){
      if(is.na(j)==T){
      if((j-1)<notseen){state[i,j] <- 2} #If event = 2 (thus, seen), latent state is 2 (seen on previous occasion)
      if((j-1)==notseen){state[i,j] <- 3} #If not seen on previous occ, latent state is 3 (not seen on previous occasion)
    }}}
    state[i,m1] <- NA #Because we condition on first capture?
  }
  return(state)
}

# Function to create initial values for unknown z
inits.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  notseen <- max(ch, na.rm = TRUE)
  v <- which(ch==notseen)
  ch[-v] <- NA
  ch[v] <- 2 #Since state 2 and 3 are possible on recap occasions but only 2 is observed?
  return(ch)
}


# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.trap(rCH, 3))

# Initial values
inits <- function(){list(mean.phij = runif(1, 0, 1), mean.phiad = runif(1, 0, 1), 
                         mean.pnonobs = runif(1, 0, 1), mean.pobs = runif(1, 0, 1), 
                         z = inits.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phij", "mean.phiad",  "mean.pnonobs", "mean.pobs")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R
transtrap.cc  <- jags(jags.data, inits, parameters, "trans_trap_c-c.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)
#Error in checkForRemoteErrors(val) : 
#3 nodes produced errors; first error: Error in node z[495,6]
#Cannot normalize density

#Something wrong with initial values??? Or known latent states?
#Nodes producing the first errors are those with 1 3 3 as last three events (only seen on marking occasion)

#Plots
par(mfrow = c(3, 3), las = 1)
hist(transtrap.cc$sims.list$mean.phij, col = "gray", main = "",  xlab = expression(phi[juv]))
abline(v = phi.j, col = "red", lwd = 2)
hist(transtrap.cc$sims.list$mean.phiad, col = "gray", main = "",  xlab = expression(phi[ad]) , ylab = "")
abline(v = phi.ad, col = "red", lwd = 2)
plot(0, type = "n", axes = F, ylab = "", xlab = "")
hist(transtrap.cc$sims.list$mean.pnonobs, col = "gray", main = "",  xlab = expression(p["Not observed year before"]))
abline(v = p.nseenbef, col = "red", lwd = 2)
hist(transtrap.cc$sims.list$mean.pobs, col = "gray", main = "",  xlab = expression(p["Observed year before"]) , ylab = "")
abline(v = p.seenbef, col = "red", lwd = 2)
####################




####################


