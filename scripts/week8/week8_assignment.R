rm(list=ls()) #clear all
library(deSolve) #loads the ODE solving package
# Compare simulations with and without demographic stochasticity for the following model:
# 
# 1. SIR model with loss of immunity at rate alpha:
# dS/dt = -beta*S*I/N + alpha*R
# dI/dt = beta*S*I/N - gamma*I
# dR/dt = gamma*I-alpha*R
# which non-dimensionalizes in time to:
# dS/dtau = -R0*S*I/N + rho*R
# dI/dtau = R0*S*I/N - I
# dR/dtau = I - rho*R
# where R0=beta/gamma, rho=alpha/gamma, and tau=gamma*t
# note that, in the deterministic case, dS/dt+dI/dt+dR/dt=0, so total population size N=S+I+R is constant and you can replace R=N-S-I to reduce the number of equations for the deterministic simulations
# Suggested parameter values*: R0=4, rho=0.2, N=1000, tf=100

#For this I will maintain the same assumptions from class that of the initial starting population, ~40% are susceptible to the disease, and .1 % are infected

#--------------
#functions
#deterministic model
SIRdet = function(t,x,parms){
  with(as.list(c(parms,x)),{
    dS = -R0*S*(I/N) + rho*(N-S-I) #susceptibles
    dI = R0*S*(I/N) - I #infecteds
    list(c(dS,dI)) #return dS and dI in a vector
  })
}

#stochastic model
SIRstoch = function(R0,rho,N,I0,S0,tf=100){
  #error checking
  #making sure pop sizes are whole numbers
  N = as.integer(N)
  I0 = as.integer(I0)
  S0 = as.integer(S0)
  if((I0+S0)>N) stop("Total population sizes S0+I0 are greater than total population size N")
  if(rho>1) stop("Cannot have loss of immunity rate rho > 1")
  #assign remaining parameters and intial conditions
  alpha = 1-rho #recovery rate
  Rec0 = N-S0-I0 #initial population of recovered individuals
  vt = 0 #vector of time each event occured
  SIR = as.matrix(c(S0,I0,Rec0)) #build initial population sizes matrix to add to over time
  tstep = 1 #counter of events that have occured
  #outcome matrix
  dN = cbind(c(-1,1,0),c(0,-1,1),c(1,0,-1)) #the three possible outcomes of an event:  infection, recovery, and loss of immunity for recovered individuals
  nOut = ncol(dN) #number of possible outcomes
  while(vt[tstep]<tf) {#run until we hit tf
    #find rates of each event and their sum
    Nt = sum(SIR[,tstep])
    wt = c(R0*SIR[1,tstep]*SIR[2,tstep]/Nt, SIR[2,tstep],rho*SIR[3,tstep]) #rate for each event -  make sure in same order as outcome matrix above
    sum_wt = sum(wt) #sum of rate wait times
    #calculate wait time to next event
    dt = rexp(n=1,rate=sum_wt) #exponential distribution draw
    #calculate which event occurs
    wtn = wt/sum_wt #normalized rates
    outIdx = sample(1:nOut, size=1, prob=wtn) #determine which event occurs
    #update S,I, and/or R populations
    SIR = cbind(SIR, SIR[,tstep]+dN[,outIdx])
    #update our time
    vt[tstep+1] = vt[tstep]+dt
    #update counter
    tstep = tstep+1
  }
  return(list(t=vt, N=SIR)) #returns the SIR matrix and time vector together as a list
}


#--------------
#parameters and initial population sizes
# Suggested parameter values*: R0=4, rho=0.2, N=1000, tf=100
R0= 4 #infection rate
rho= 0.2 #loss of immunity rate
N= 1000 #population size
tf=100
S0= floor(0.4*N) #initial number of susceptibles as a whole number
I0= floor(0.01*N) #floor cuts off everything after the decimal, rounding down and guranteeing a whole number initial infected


#--------------
#run everything
#run the stochastic model
StochOut = SIRstoch(R0=R0,rho=rho,N=N, I0=I0, S0=S0)

#run the deterministic model
parms = c(R0=R0,rho=rho,N=N) #parameters vector
x0= c(S=S0, I=I0) #initial population sizes
detOut=as.data.frame(lsoda(x0, StochOut$t,SIRdet,parms)) #run ODE

#--------------
#plot everything
mf = par(mfrow=c(2,1)) #sets up 2x1 subplot
plot(StochOut$t,StochOut$N[1,],type='l',col='blue',xlab="time",ylab="Susceptibles",main="Stochastic vs deterministic SIR model with loss of immunity")
lines(StochOut$t,detOut$S,col='red')

plot(StochOut$t,StochOut$N[2,],type='l',col='blue',xlab="time",ylab="Infected",main="Stochastic vs deterministic SIR model with loss of immunity")
lines(StochOut$t,detOut$I,col='red')
par(mfrow=mf)
