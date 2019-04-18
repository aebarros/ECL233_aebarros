#Leslie Matrix pracitce

#load teasel data from WErner and Cassel (1977)
#run for 50 time steps starting with two individuals in stage 1
#then plot actual growth factor over time: Ntot[t+1]/N[t]
#compare to the expected long term growth factor (leading eigenvalue)

rm(list=ls())
#parameters
TS<-as.matrix(read.table("data/teasel_stage.txt")) #upload teasel data and convert from table into a matrix
nstg = nrow(TS) #tells us number of stages
n0 = c(2,rep(0,nstg-1)) #sets starting pop
n0
tf =500
vt = 1:tf #vector of times

#Run the model through time
Nt=matrix(NA,nstg,tf) #starting empty matrix to fill with pop sizes through time
Nt[,1]=n0 #fill in initial population vector
Nt

for(t in 1:(tf-1)) {
  Nt[,t+1] = TS%*%Nt[,t]
  }#runs the teasel matrix through time, giving us pop size for each age class for each time step up to 50
Nt[,500]
Ntot = colSums(Nt) #gives us total population size for each time step
Ntot

leval=function(L) Re(eigen(L)$values[1]) #function find leading eigen value of a matrix

#calculate expected long-term growth factor
lambda = leval(TS) #calculate leading eigenvalue

#plot results
plot(vt[-tf],Ntot[-1]/Ntot[-tf],type="l",col="black",xlab="time", ylab="Growth factor")# for y, want everything but the first entry, because...?
lines(vt[-tf],rep(lambda,tf-1),col="red")


