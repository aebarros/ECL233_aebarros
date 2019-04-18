# class 3 code: linear matrix models with Teasel data
#-----------------------
# 1. upload teasel data (read.table) from Werner and Caswell (1977)
# 2. run for 50 time steps starting with 2 individuals in stage 1
# 3. plot actual growth factor over time: N_{t+1}/N_t
# 4. plot expectation expected growth factor over time (leading eigenvalue): how much deviation do you get from not being in the stable age structure?

rm(list=ls()) # clear everything

# functions
leval = function(L) Re(eigen(L)$values[1]) # calculate leading eigenvalue of input matrix L

# parameters
TS = as.matrix(read.table("teasel_stage.txt")) # upload teasel data - check it out at the prompt as TS -- make sure you're in the same directiory as the data file (or add the path to the filename)
#------------
# run Teasel model, find total population size over time, also find eigenvalue
nstg = nrow(TS) # number of stage classes
N0 = c(2, rep(0, nstg-1)) # vector of initial population sizes
tf = 50 # final time step
Nt = matrix(0,nstg,tf) # initialize vector to fill in times
Nt[,1] = N0 # initial total population size
# run model
for(t in 1:(tf-1)) {
	Nt[,t+1] = TS%*%Nt[,t] # remember to use matrix multiplication!
	}
Ntot = colSums(Nt) # sum to get total population size over time
gf_time = Ntot[-1]/Ntot[-tf] # find growth factor for each time step: total population size each time step divided by the previous time step
# determine expectation with leading eigenvalue
lambda = leval(TS) # leading eigenvalue
#------------
# plot results - together
vt = 1:(tf-1) # vector of times
#pdf(file="teaselTime.pdf") # create file for saving plot 
plot(vt,gf_time,type="l",col="black",xlab="Time",ylab="Growth factor")
lines(vt,rep(lambda,1,tf-1),col="red")
#dev.off()
