#class notes week 8
#demographic stochasticity in continuous time
#how much does discrete time vs continuous time influence the model
#so we will start with an ODE model (ordinary differential equation) and then put in some demographic stochasicity 

#we will have: (see notes)
#total pop: N = S + I + R
#susceptible becomes infected: dS/dt =  -B*(S*I)/N
#infected: dI/dt = BSI/N - YI (amount recovered)
#recovered: dR/dt = YI 
#(See notes for above)

#We will be using non-dimensionalization to take units out and simplify the model
#so we will non-dimensionalize time, and turn it into disease generations (t)
#so tao=t(gamma+mew)

#TWO STEPS:
#1 Calculate wait time to next event
#2 calculate which event occurs (birth, death, infection, recovery)

#SIR (simple infection rate) model with and without demographic stochasticity

rm(list=ls()) #clear all
library(deSolve) #loads the ODE solving package

#--------------
#functions
#deterministic model
SIRdet = function(t,x,parms){
  with(as.list(c(parms,x)),{
    dS = p*(N-S) - R0*S*I/N #susceptibles
    dI = R0*S*I/N - I #infecteds
    list(c(dS,dI)) #return dS and dI in a vector
  })
}

#stochastic model
SIRstoch = function(R0,p,N,I0,S0,tf=100){
  #error checking
  #making sure pop sizes are whole numbers
  N = as.integer(N)
  I0 = as.integer(I0)
  S0 = as.integer(S0)
  if((I0+S0)>N) stop("Total population sizes S0+I0 are greater than total population size N")
  if(p>1) stop("Cannot have a birth/death rate p > 1")
  #assign remaining parameters and intial conditions
  alpha = 1-p #recovery rate
  Rec0 = N-S0-I0 #initial population of recoverd individuals
  vt = 0 #vector of time each event occured
  SIR = as.matrix(c(S0,I0,Rec0)) #build initial population sizes matrix to add to over time
  tstep = 1 #counter of events that have occured
  #outcome matrix
  dN = cbind(c(1,0,0),c(-1,0,0),c(-1,1,0),c(0,-1,1),c(0,-1,0),c(0,0,-1)) #the six possible outcomes of an event: birth, death of S, infection, recovery, death of I, and death of R
  nOut = ncol(dN) #number of possible outcomes
  while(vt[tstep]<tf) {#run until we hit tf
    #find rates of each event and their sum
    Nt = sum(SIR[,tstep])
    wt = c(p*Nt, p*SIR[1,tstep],R0*SIR[1,tstep]*SIR[2,tstep]/Nt, alpha*SIR[2,tstep],p*SIR[2,tstep],p*SIR[3,tstep]) #rate for each event -  make sure in same order as outcome matrix above
    sum_wt = sum(wt) #sum of wait times
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
R0= 8 #infection rate
p= 0.04 #birth/date rate
N= 1000 #population size
S0= floor(0.4*N) #initial number of susceptibles as a whole number
I0= floor(0.01*N) #floor cuts off everything after the decimal, rounding down and guranteeing a whole number initial infected


#--------------
#run everything
#run the stochastic model
StochOut = SIRstoch(R0=R0,p=p,N=N, I0=I0, S0=S0)

#run the deterministic model
parms = c(R0=R0,p=p,N=N) #parameters vector
x0= c(S=S0, I=I0) #initial population sizes
detOut=as.data.frame(lsoda(x0, StochOut$t,SIRdet,parms)) #run ODE (need to add time)

#--------------
#plot everything
mf = par(mfrow=c(2,1)) #sets up 2x1 subplot
plot(StochOut$t,StochOut$N[1,],type='l',col='blue',xlab="time",ylab="Susceptibles",main="Stochastic vs deterministic SIR model with births and deaths")
lines(StochOut$t,detOut$S,col='red')

plot(StochOut$t,StochOut$N[2,],type='l',col='blue',xlab="time",ylab="Infected",main="Stochastic vs deterministic SIR model with births and deaths")
lines(StochOut$t,detOut$I,col='red')
par(mfrow=mf)
