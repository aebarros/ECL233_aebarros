# two plots using Ricker model (N_t+1=N_t*exp(r*(1-N_t/K))):
# 1. plot N vs. t for four different values of r: 1.5, 2.3, 2.6, 3.0 (show different behaviors - one-cycle, two-cycle, four-cycle, chaos)
# 2. bifurcation plot for r from 1.5 to 3.6

# you'll have to use a for loop to iterate through time
# BUT you can do this most efficiently if you can figure out how to "vectorize" over multiple values of r

# other parameter suggestions (can use different values, your choice):
# K=100
# N0=50
# tf=100 in 1st plot, 500 in 2nd

# Watch out for any repeated code that you could make into a function.

# EXTRA CREDIT: Create a plot of the Lyapunov exponents associated with your Ricker bifurcation diagram. 

rm(list = ls()) #this clears the environment

#goal: use a Ricker model to plot N vs t for different values of r

Ricker<-function(nt,K=100,r=1.5){
  nt1=nt*exp(r*(1-nt/K)) #Ricker model over one time step
  return(nt1) # return the population size in the next time step
}


tf = 100 #set time limit


#try Ricker through time with different initial conditions
r=c(1.5,2.3,2.6,3)
rn0 = length(r) #how many inital conditions we are using

n=matrix(NA, nrow=tf, ncol=rn0) #sets up an empty matrix to fill in with our pop sizes


#want to apply function through time and loop through all initial conditons at once
n[1,] = c(50,50,50,50)
r
n
for(t in 1:(tf-1)) { #loop through time
  n[t+1,1] = Ricker(n[t,1],r=1.5) #Could not figure out the proper way to vectorize r for the calculations, but found this way that works, by breaking up each column/row fill into a different function in the for loop
  n[t+1,2] = Ricker(n[t,2],r=2.3)
  n[t+1,3] = Ricker(n[t,3],r=2.6)
  n[t+1,4] = Ricker(n[t,4],r=3)
}
matplot(1:tf,n,type="l",xlab="Time",ylab="Pop Size") #plots our matrix for us, plots all the initial conditions run at the same time

#how to do subplots:
#need to open a series of subplots, populate them and them close them
mf = par(mfrow=c(2,2)) #sets up two sub plots - two plots one row
for(stp in 1:4) plot(1:tf,n[,stp],type="l",xlab="Time",ylab="N",main=paste("r=",r[stp])) #plotting two initial conditions in each subplot with title saying what each initial condition was

# 2. Goal: create a bifurcation plot for r from 1.5 to 3.6
#note: much of this code was found from : http://rstudio-pubs-static.s3.amazonaws.com/880_1fc8964becf942f7a754468009603fa9.html and then played with to work
rm(list = ls()) #this clears the environment

rmax <- 3.6
plot(-1, -1, xlim = c(0, rmax), ylim = c(0, 500), xlab = "r", ylab = "N")
a <- 0.01
r <- seq(0, rmax, by = 0.1)
n <- 100

for (z in 1:length(r)) {
  nt <- vector()
  nt[1] <- 10
  for (i in 1:n) {
    
    nt[i+1] <- nt[i] * exp(r[z] *(1-nt[i]/100))
    
  }
  uval <- unique(nt[40:n])
  points(rep(r[z], length(uval)), uval, cex = 0.1, pch = 19)
}

