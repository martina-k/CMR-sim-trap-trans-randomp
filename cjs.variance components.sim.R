#Simulate a CJS data set for 3 sites 
#Evaluate identifiability of the ICC 
#Code written by Sarah Converse, sconver@uw.edu 
#CJS simulation code and m-array analysis code adapted from Kery and Schaub

#############################################################################

#############################################################################

#Simulate CJS data 
n.occ <- 15 #number of occasions
sites <- 3 #number of sites 
marked <- rep(20,n.occ) #number of inds marked per occasion
#matrix of values for phi and p (occasion - 1)*inds

b0 <- 3
sigma.s1 <- sigma.s2 <- sigma.s3 <- 0.5
sigma.nrw <- 0.3
sigma.all <- 0.3 

sims <- 50 
results <- matrix(NA,nrow=sims,ncol=7)
colnames(results) <- c("b0","mean.p","sigma[1]","sigma[2]","sigma[3]","sigma.nrw","sigma.all")
for(s in 1:sims){

delta.all <- rnorm((n.occ-1),0,sigma.all) 
delta.nrw <- rnorm((n.occ-1),0,sigma.nrw)
eps.s1 <- rnorm((n.occ-1),0,sigma.s1)
eps.s2 <- rnorm((n.occ-1),0,sigma.s1)
eps.s3 <- rnorm((n.occ-1),0,sigma.s1)

Phi.s1 <- 1/(1+exp(-(b0 + eps.s1 + delta.all))) #Isle of May 
Phi.s2 <- 1/(1+exp(-(b0 + eps.s2 + delta.nrw + delta.all))) #Sklinna
Phi.s3 <- 1/(1+exp(-(b0 + eps.s3 + delta.nrw + delta.all))) #Rost 

PHI.s1 <- matrix(Phi.s1, ncol = n.occ-1, nrow = sum(marked), byrow=TRUE)
PHI.s2 <- matrix(Phi.s2, ncol = n.occ-1, nrow = sum(marked), byrow=TRUE)
PHI.s3 <- matrix(Phi.s3, ncol = n.occ-1, nrow = sum(marked), byrow=TRUE)
P.s1 <- matrix(rep(0.7,n.occ-1), ncol = n.occ-1, nrow = sum(marked),byrow=TRUE)
P.s2 <- matrix(rep(0.7,n.occ-1), ncol = n.occ-1, nrow = sum(marked),byrow=TRUE)
P.s3 <- matrix(rep(0.7,n.occ-1), ncol = n.occ-1, nrow = sum(marked),byrow=TRUE)

#simulate CJS data 
CJS.sim <- function(marked,P,PHI){
  n.occ <- dim(PHI)[2]+1 #number of occasions
  EH <- matrix(0, ncol = n.occ, nrow = sum(marked)) #assign matrix
  first <- rep(1:length(marked),marked[1:length(marked)]) #get first encounter occasion for each ind
  for(i in 1:sum(marked)){ #loop over inds
    EH[i,first[i]] <- 1 #condition on first capture
    if(first[i]==n.occ) next #if first capture occasion for ind is last occasion, next
    for(t in (first[i]+1):n.occ){ #loop over time 
      surv <- rbinom(1,1,PHI[i,t-1]) #Bernoulli survival trial
      if(surv == 0) break #if didn't survive, done
      cap <- rbinom(1, 1, P[i,t-1]) #if survived, Bernoulli capture trial
      if (cap==1) EH[i,t] <- 1 #fill in the EH
    }
  }
  return(EH)
}
#run the simulation to get encounter histories 
EH.s1 <- CJS.sim(marked,P.s1,PHI.s1)
EH.s2 <- CJS.sim(marked,P.s2,PHI.s2)
EH.s3 <- CJS.sim(marked,P.s3,PHI.s3)
EH <- rbind(EH.s1,EH.s2,EH.s3)
site <- c(rep(1,nrow(EH.s1)), rep(2,nrow(EH.s2)), rep(3,nrow(EH.s3)))
nrw <- c(rep(0,nrow(EH.s1)), rep(1,nrow(EH.s2)), rep(1,nrow(EH.s3)))

#############################################################################

#############################################################################

#Fit in JAGS using a state-space model 

#JAGS MODEL 
cat("
model{
    
#state-space likelihood 
for(i in 1:n.ind){
  #  z[i,first[i]] <- 1  #condition on first capture 
  for(j in (first[i]+1):n.occ){
    z[i,j] ~ dbern(z[i,j-1]*phi[i,j-1]) #state process 
    
    y[i,j] ~ dbern(z[i,j]*p[i,j-1]) #observation process 
  }
}
    
#constraints and priors 
for(i in 1:n.ind){
  for(j in 1:(n.occ-1)){
    logit(phi[i,j]) <- b0 + delta.all[j] + delta.nrw[j]*nrw[i] + eps[site[i],j]  
    
    p[i,j] <- mean.p
  }
}

for(i in 1:3){
  tau[i] <- pow(sigma[i],-2)
  sigma[i] ~ dunif(0,10)

  for(j in 1:(n.occ-1)){
    eps[i,j] ~ dnorm(0,tau[i])
  }
}

for(j in 1:(n.occ-1)){
  delta.nrw[j] ~ dnorm(0,tau.nrw)
  delta.all[j] ~ dnorm(0,tau.all)
}



tau.nrw <- pow(sigma.nrw,-2)
sigma.nrw ~ dunif(0,10)
tau.all <- pow(sigma.all,-2)
sigma.all ~ dunif(0,10)

mean.p ~ dunif(0,1)
b0 ~ dunif(-5,5)

}
",file="CJS.txt")

## JAGS input
n.ind <- nrow(EH)
n.occ <- ncol(EH)
first <- rep(NA,n.ind)
for(i in 1:n.ind){ first[i] <- min(which(EH[i,]==1))}

#Provide known state as data  
state.data <- function(EH){
  state <- EH
  for (i in 1:dim(EH)[1]){
    first <- min(which(EH[i,]==1))
    last <- max(which(EH[i,]==1))
    state[i,first:last] <- 1
    #   state[i,first] <- NA
  }
  state[state==0] <- NA
  return(state)
}

z.state <- state.data(EH)
data <- list(n.ind=n.ind,n.occ=n.occ,y=EH,first=first,z=z.state,site=site,nrw=nrw)
params <- c("b0","mean.p","sigma","sigma.nrw","sigma.all") 
inits =  function() {list(b0=runif(1),mean.p=runif(1))} 

library(jagsUI)
## RUN
out.jags = jags(data, inits, params, model.file="CJS.txt",
                n.chains=3, n.iter=20000, n.burnin=5000, n.thin=1)

results[s,1] <- mean(out.jags$sims.list$b0)
results[s,2] <- mean(out.jags$sims.list$mean.p)
results[s,3] <- mean(out.jags$sims.list$sigma[,1])
results[s,4] <- mean(out.jags$sims.list$sigma[,2])
results[s,5] <- mean(out.jags$sims.list$sigma[,3])
results[s,6] <- mean(out.jags$sims.list$sigma.nrw)
results[s,7] <- mean(out.jags$sims.list$sigma.all)

print(s)

}

