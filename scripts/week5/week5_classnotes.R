#Class Notes 4/30/2019 Continues-time models
source("scripts/ecl233_functions.R")

#goals for the day:
#1) learn how to integrate a continuous time model
#2) fit models to our paramecium data


rm(list=ls()) #clear env
library(deSolve) #load deSolve package for ordinary differential equations (ODE) integration, main function: lsoda

#for lsoda function, R requires a vector or a list of labelled parameters
v=c(x=1,y=4,z=100) #example of parameters as a list

v["x"] #this provides you the value labeled "x" in vector v

with(as.list(v),print(x)) #prints value x in list v

#------------------------------------------------------------
#function we want to integrate
#note our logistic function dN/dt=r*N(1-(Nt/K)) doesn't have time on right side, but we need it for lsoda so we put in a dummy variable
#note:if using multiple populations, label them differently in a vector
logInt = function(t,n,parms){
  with(as.list(parms), {
    dndt=r*n*(1-n/K) #logistic differential equation
    return(list(dndt)) #returns dndt as a list
  })
}

#function to minimize to fit data to the logistic and find best fit parameters
#takes in parameter values that we are guessing and then return what we want to minimize: the mean square error between simulate and original data
logIntMin=function(parms){
  out=as.data.frame(lsoda(n0,times,logInt,parms)) #runs ode with current guess for parameter values in parms
  mse = mean((out$n-nTrue)^2) #find mean square error between the simulated and original data
  return(mse) #return mean square error
}

#-------------------------------------------------------------
#next we actually need our parameter values and our time
r=0.5 #growth rate
K=100 #carrying capacity
parms = c(r=r,K=K) #creates parameter vector
tf=20 #final time
times =1:tf #vector of times to run over
n0=c(n=2) #initial population size labeled as n

#-------------------------------------------------------------
#run the model: use lsoda, inputs are always y, times, func, parms in that order
out=as.data.frame((lsoda(n0,times,logInt,parms)))

#plot output
plot(out$time, out$n, type="l", xlab="Time", ylab="population size")

#so lets fit some paremicium data to this to find the r and K of that population data set
#we are going to try different r and K values and look for the thing that gives me the minimum difference
# use the prebuild in function optum

#fit model to data to discover r and K for population
nTrue=as.matrix(read.table("data/paramecium1.txt")) #read actual data converted to a matrix
colnames(nTrue)=NULL #take out the column names
n0=c(n=nTrue[1]) #provide initial n0 and with labels used in the logIntMin function
tf=length(nTrue) #take final time from data
times=1:tf #creates full time series

#need to start with some guesses for our parameters
rguess=1 #guess for growth rate
Kguess=500 #capacity guess
parms0 = c(r=rguess,K=Kguess) #create a labeled parameters vector with initial guesses

#find fit: send optim our initial guess and function we want to minimize
optimOut= optim(parms0, logIntMin)
parms = optimOut$par #takes out best guess values for the parameters

#run model and compare to the data
nSim=as.data.frame(lsoda(n0,times,logInt,parms)) #run the model with our best guess parameters

#plot results
plot(times,nTrue,type="p",xlab="Time",ylab="Abundance")
lines(nSim$time,nSim$n)
