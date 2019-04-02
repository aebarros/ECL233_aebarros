#Ricker model:
#N(t+1)=N(t)e^r(1-Nt/k)
#Beverton-Holt:
#N(t+1)=((e^r)Nt)/(1+((e^r-1)Nt)/K)
  #K= carrying capacity
  #t will not change for these plots/models
  #both of these models have high density dependence mortality, one is based off of canibalism the other on space

#note: I decided to do both models and the plot generation all in one function

rm(list = ls()) #this clears the environment
#goal: plot N_t+1 vs. N_t for Ricker and Beverton-Holt on the same plot for the same K and r values

#define functions
  #Ricker Model
growth=function(r,nt,K){
  n.Ricker=nt*exp(r*(1-nt/K)) #Ricker model
  n.BevHolt=(exp(r)*nt)/(1+(exp(r)-1)*nt/K)
  pdf(file="figures/exp1.pdf") #start a file to save
  plot(nt.base,n.BevHolt,type="l",col="blue",xlab="N(t)",ylab="N(t)+1") #include the plot in the function
  lines(nt.base,n.Ricker,col="black")
  legend("topleft", legend=c("Ricker model","Beverton-Holt Model"),col=c("blue","black"),lty = c(1,1))
  dev.off() #save the current plot to the pdf
}


#define inputs/parameters
r.base=.5
nt.base=seq(0,100, by=.1)
K.base=10

#run the functions for these different values
growth(r.base,nt.base,K.base) #calls the growth function with the above defined inputs

source("scripts/assignment_week_1.R") #runs everything above
