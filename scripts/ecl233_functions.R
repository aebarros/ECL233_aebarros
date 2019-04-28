#Script to track all ECL 233 functions that might be useful

#Ricker Growth Model
Ricker<-function(nt,K,r){
  nt1=nt*exp(r*(1-nt/K)) #Ricker model over one time step
  return(nt1) # return the population size in the next time step
}

#dispersal function following a laplacian dispersal kernel where a relates to the spread, probability an offspring born at x, disperses to y
prob.disp<-function(y,x){
  exp(-abs(y-m-x)/a)/(2*a)
}

#Leading eigen value function to provide us with lambda
leval=function(L) Re(eigen(L)$values[1]) #function find leading eigen value of a matrix
#calculate expected long-term growth factor

