#Assignment week 5

#---------------------------------
# For this assignment, you have two choices for simulating two interacting species: (a) competition (with practice in fitting to data) or (b) predator-prey (with practice in considering multiple functional forms); pick one.

# Regardless of your choice, you'll need set up your ode function for multiple species:
# Initialize your population sizes with labels, e.g., n0=c(x=x0, y=y0, ...), then use these labels in your ode function
# odeFn = function(t, n, parms) {
# with(as.list(c(parms, n)), { # extract parameters from parms vector
# dx = ... # dx/dt
# dy = ... # dy/dt
# ...
# return(list(c(dx, dy, ...))) # return dn/dt as a list
# })
# }
# then running lsoda gets you each species as out$x, out$y, etc.

# Choice (a) competition
# fit paramecium data (from Leslie 1957 reading) to Lotka-Volterra competition model
# dn1/dt = r1*n1*(1-n1/K1-a12*n2/K1)
# dn2/dt = r2*n2*(1-n2/K2-a21*n1/K2)
# plot actual and simulated time series for both populations
# you get r's and K's - this is in the "parameciumrk.txt" file - don't forget to use colnames(your input matrix) = NULL to erase the default V1, V2, etc. column names
# you'll have to find the alpha's to the data in the "paramecium2.txt" file
# note that this means you have to reconstruct the parameters vector within optim - start with parmsrk=r and K, then do parms = c(parmsrk, a=a) to add a's within minimizing function

# Choice (b) predator-prey
# Create phase plane diagrams (predator density vs. prey density, where lines are trajectories over time) for the predator-prey mode:
# dN/dt = r*F(N)-b*G(N)*P # prey 
# dPdt = c*G(N)*P-m*P # predator 
# for four scenarios:
# (i) type I functional response (G(N)=N) and no density dependence in prey (F(N)=N)
# (ii) type I functional response (G(N)=N) and density dependence in prey (F(N)=N*(1-N/K))
# (iii) type II functional response (G(N)=N/(1+d*N)) and density dependence in prey (F(N)=N*(1-N/K))
# (iv) type III functional response (G(N)=N^2/(1+d*N^2)) and density dependence in prey (F(N)=N*(1-N/K))
# in addition to plotting trajectories, mark the starting and ending points. Bonus if you also put a point at the internal equilibrium (equilibrium with nonzero values for both prey and predator) for each model.
# some parameter recommendations, you can alter these as you see fit:
# b = 0.01 # predator attack rate
# c = 0.1*b # predator conversion of preation into reproduction
# m = 0.2 # predator mortality
# a = 1 # rate of prey capture/unit prey and time
# d = 0.0015 # handling time
# r = 0.5 # prey growth rate
# K = 1000 # prey carrying capacity
# N0=300, P0=50 for initial population sizes, although you might want to try a few different initial population sizes per plot



#------------------------------------
rm(list=ls()) #clear env
library(deSolve) #load deSolve package for ordinary differential equations (ODE) integration, main function: lsoda
#set up ode function
odeFn= function(t,n,parms){
  with(as.list(c(parms,n)),{ #extract parameters from parms vector
    dn1 = r1*n1*(1-n1/K1-a12*n2/K1)
    dn2 = r2*n2*(1-n2/K2-a21*n1/K2)
    return(list(c(dn1,dn2))) # return dn/dt as a list
    })
}

#function to minimize to fit data to the logistic and find best fit parameters
#takes in parameter values that we are guessing and then return what we want to minimize: the mean square error between simulate and original data
odeFnMin=function(parmsa){
  parms=c(parmsrk,parmsa) #need to reconstruct parms in function based on parmsa a guesses
  out=as.data.frame(lsoda(n0,times,odeFn,parms)) #runs ode with current guess for parameter values in parms
  mse1 = mean((out[,1]-nTrue[,1])^2)#find mean square error between the simulated and original data for pop1
  mse2 = mean((out[,2]-nTrue[,2])^2)#find mean square error between the simulated and original data for pop2
  return((mse1+mse2)/2) #return the mean of mse for both populations (why?)
}

#-------------------------------------
#load data sets
nLV=as.matrix(read.table("data/parameciumrk.txt")) #read actual data converted to a matrix
colnames(nLV)=NULL #take out the column names
r1=.7816
K1=559.6860
r2=.6283
K2=202.4931
parmsrk = c(r1=r1,r2=r2,K1=K1,K2=K2) #creates parameter vector
#next we need the alpha values a12 and a21
nTrue=as.matrix(read.table("data/paramecium2.txt"))
colnames(nTrue)=NULL
n1=nTrue[1,1]
n2=nTrue[1,2]
aguess12=2
aguess21=1
parms=c(parmsrk,a12=aguess12,a21=aguess21)
parmsa=c(a12=aguess12,a21=aguess21)
n0=c(n1=n1,n2=n2)
tf=length(nTrue) #take final time from data
times=1:tf #creates full time series


#-----------------------------------------

#find fit: send optim our initial guess and function we want to minimize
# optim takes in a function and gives you the parameter values that minimize MSE that function
optimOut= optim(parmsa, odeFnMin)
parms = optimOut$par

#run model and compare to the data
nSim=as.data.frame(lsoda(n0,times,odeFn,parms)) #run the model with our best guess parameters

#plot output
plot(nSim$time, nSim$n2, type="l", xlab="Time", ylab="population size",col="blue",ylim=c(0, 400))
lines(nSim$n1,col="red")
points(nTrue[,1],col="blue")
points(nTrue[,2],col="red")

