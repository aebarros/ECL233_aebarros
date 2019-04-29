#looking at a model of a population that lives on a 1 dimensional space
#along that transect, only a fixed proportion is habitable
#Each year t, Nt individuals have offspring in a density dependent matter
#then parents die and offspring disperse according to a dispersal kernal
#offspring that land withing [-L,L] survive
# we want to track how density of individuals within that space changes from year to year

#parameters
rm(list=ls())
source("scripts/ecl233_functions.R") #source function script
r=1
K=100
a=2
m=0
tf=100
L=20
bins=41

#Plot1
#lets think about my final plot:
#My X values will be -20=<x<=20
#my y values will be density values for each point x for time t, or N[t,x]
#so my final matrix would look something like a 41x100, with time on the rows, location at the columns, N[t,x] as the values
#To do this we need a few things:
  #Growth at x based on N[t,x]
  #the sum of N[t+1] from other locations * probability of dispersal to x
  #surivorship, which is 1 if individual disperses within -L=<x<=L (I think this might be possible to ignore?)

#Ricker model for growth function (note: I changed this from my base "Ricker" function, so that r wasn't a preset)
growth<-function(nt,K,r){
  nt1=exp(r*(1-nt/K)) #Ricker model over one time step
  return(nt1) # return the population size in the next time step
} #use Ricker function from ecl233_functions script
  
#dispersal function following a laplacian dispersal kernel where a relates to the spread, or probability an offspring born at x, disperses to y (this was provided in the assignment)
prob.disp<-function(y,x){
  exp(-abs(y-m-x)/a)/(2*a)
}

#Create matrix A which is the probability of dispersal between each site
#function to create A based on different L
dx=(2*L)/bins #set the delta of x

A.create=function(L){ #this funcion creates a dispersal probability matrix based on a value of L
  x=seq(-L+dx/2,L-dx/2,length=bins) #create vector x from last midpoint to first midpoint
  At=outer(x,x,prob.disp) #use the outer function (I'm still not positive on how it works) to create matrix At that is x by x in length, filled by probability of dispersal from any x to anyother x
  return(At)
}
A.create(20)

#create emptry n matrix
n=matrix(NA,tf,bins)
n[1,1:19]=0 #fill the first vector/timestep (is there a simpler way of doing this?)
n[1,20:22]=10
n[1,23:41]=0
  
A=A.create(20) #create matrix A of dispersal proability using the A.create function from earlier
t=1 
for(t in 1:(tf-1)){ #for loop that fills in matrix n
  n[t+1,]=A%*%(growth(n[t,],K,r)*n[t,]*dx) #need to fill in a matrix "n" here with # of rows=tf, # of columns = bins
}

#now we have our matrix (n) of density over time!!!

#plotting time
#first set a row 1 for matrix column numbers
x.points<-c(seq(from=-20,to=20,by=1)) #sets up points to use for x-axis
b=t(n) #take tranverse of the n matrix to plot it (not sure why I have to do this, but it seems to work)
matplot(b,type="l",xlab="location",ylab="density",col = "black",axes=F) #axes=F drops axis labels to replace
axis(2) #sets back y-axis
axis(side=1,at=1:nrow(b),labels = x.points) #sets x-axis


#part 2 eigen value question
n[t+1,]=A%*%(dx*n[t,]*(exp(r))) #removed the density dependence
#get eigen value from A*dx*e^r redifining A many different times for different Ls (why is n[t,] dropped?) because without density dependence growth doesn't change?
#so I need a for loop probably for different values of L, that creates A matrix of dispersal probability, and then gets an eigen value from it?

#first need a matrix to fill, dominant eigen value in one column, L in the other

Le <- seq(0.01,5,by=0.01) 
Lmat <- cbind(Le, NA)#create new matrix to populate

for(i in 1:length(Le)){
  Li=Le[i]
  dx=(2*Li)/bins #have to recreate dx based on the new L values
  Lmat[i,2]=leval(A.create(Li)*dx*exp(r)) #fills the out matrix
}

#make the plot
plot(Lmat, xlab="Transect Length (L)", ylab="Lambda")
abline(a=1,b=0, col="red")
