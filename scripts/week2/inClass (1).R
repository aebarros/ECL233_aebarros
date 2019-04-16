# practicing loops

# function for Beverton-Holt model
# inputs:
# nt = population size at time t
# alpha = growth at low population sizes
# beta = determines saturation point
# output: population size at time t+1
bh = function(nt, alpha=1.2, beta=0.001) { # note: setting defaults
  nt1 = alpha*nt/(1+beta*nt) # Beverton-Holt over one time step
  return(nt1) # return the population size in the next time step
}

# set our parameters
n0 = 100 # initial population size
tmax = 10000 # maximum time allowed
minDiff = 0.01#5 # minimum difference to singals reach equilibrium
n = rep(NA,2) # start an empty vector
n[1] = n0 # fill in initial popula size
n[2] = bh(n[1]) # population size at the second time step
t=2 # starting at time two
while(abs(n[t]-n[t-1])>minDiff & t<tmax) { # run until equilbrium: difference in population size between two time steps is small - use the tmax to avoid an infinite loop
  n[t+1] = bh(n[t]) # updates population size
  t = t+1 # update time step
  #print(n)
  #print(t)
}

rm(n) # clear out population size
tf = 100 # final time (instead of having our minimum difference)
n = rep(NA, tf) # set aside a vector of population sizes to fill in: good for efficient memory use
n[1] = n0 # set initial population size
for(t in 1:(tf-1)) { # loops through each value of t from 1 to tf-1
  n[t+1] = bh(n[t]) # run Beverton-Holt model for one time step
  #print(t)
  #print(n)
}
plot(1:tf,n) # plot time series

# try Beverton-Holt through time with different initial conditions
vn0 = seq(10,300,by=10) # vector of initial conditions
Nn0 = length(vn0) # how many initial conditions we're using
n = matrix(NA, nrow=tf, ncol=Nn0) # setting up empty matrix to fill in with population sizes
n[1,] = vn0 # fill in initial conditions
for(t in 1:(tf-1)) { # loop through time
  n[t+1,] = bh(n[t,]) # fill in the next time step -- next row, all columns for all initial conditions
}
matplot(1:tf,n,type='l', xlab="Time", ylab="poplation size")

# how to do subplots
mf = par(mfrow=c(2,1)) # steps up 2 subplots
for(stp in 1:2) plot(1:tf,n[,stp],type='l', xlab="Time",ylab="N",main=paste("n0=",vn0[stp])) # plot two initial conditions in different subplots, have the title say what the initial condition was
par(mfrow=mf) # wraps up the suplots