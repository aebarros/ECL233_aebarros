#load teasel data from WErner and Cassel (1977)
# A(N) is given by the teasel matrix TS with the 1,7 entry of 431 replaced by f/(1+a*sum(N)). 
##so we need to create a matrix of lambda for different fecudnity values based on f=c(1:75)


rm(list=ls())
#parameters
TS<-as.matrix(read.table("data/teasel_stage.txt")) #upload teasel data and convert from table into a matrix
nstg = nrow(TS) #tells us number of stages
n0 = c(2,rep(0,nstg-1)) #sets starting pop
n0
tf =50
vt = 1:tf #vector of times
a=.01
f=1
fmax=75

#Run the model through time
Nt=matrix(NA,nstg,tf) #starting empty matrix to fill with pop sizes through time
Nt[,1]=n0 #fill in initial population vector
Nt

leval=function(L) Re(eigen(L)$values[1]) #function find leading eigen value of a matrix
#calculate expected long-term growth factor

#maybe make a new TS_fecundity vector based on different values of f?
f=1
TS_fecun=matrix(NA,fmax,2)
TS_fecun[1,1]=f
Ntot=2
for(f in 1:(fmax-1)){
  TS_fecun[f+1,1]=f+1
  TS_fecun[f+1,2]=(f/(1+a*Ntot))
}

TS_fecun

#maybe try to make a for loop to fill in f_lambda matrix?
#make a matrix with column1=f, column2=lambda

f=1
f_lambda=matrix(NA,fmax,2)
f_lambda[1,1]=f
f_lambda

for(f in 1:(fmax-1)){
  f_lambda[f+1,1]=f+1 #this populates the f column with all the f values 1:75
  TS[1,7]=(TS_fecun[f+1,2])
  f_lambda[f+1,2]= leval(TS)#here we need to populate the lambda values, probably need to rebuild our leval function? but we need some way to update the TS matrix with a new fecundity based on f?
}
f_lambda


#############Part2################
#parameters
TS<-as.matrix(read.table("data/teasel_stage.txt")) #upload teasel data and convert from table into a matrix
nstg = nrow(TS) #tells us number of stages
n0 = c(2,rep(0,nstg-1)) #sets starting pop
n0
tf =50
vt = 1:tf #vector of times
a=.01
fmax=75
tf=500
a=.01
n0=c(2,0,0,0,0,0,0)
#Run the model through time
Nt=matrix(NA,nstg,tf) #starting empty matrix to fill with pop sizes through time
Nt[,1]=n0 #fill in initial population vector
Nt
#build a second matrix
Ntot_f=matrix(NA,fmax,2)
Ntot_f[1,1]=1
Ntot_f

###Double for loops### Is this working?
for(f in 1:fmax-1){
   Ntot_f[f+1,1]=f+1 #filling in f values
      for(t in 1:(tf-1)){ 
     Ntot=sum(Nt[,t]) #calculates Ntot for current time t
     TS[1,7]=(f/(1+a*Ntot)) #fills TS[1,7] with function value for set f
     Nt[,t+1] = TS%*%Nt[,t] #fills in rest of Nt matrix
     }
   Ntot_f[f+1,2]=sum(Nt[,500])
   } #runs the teasel matrix through time, giving us pop size for each age class for each time step up to 50

###PLOTS###
#1st plot of lambda vs f
mf = par(mfrow=c(1,2))
line_of_growth=1
plot(f_lambda[,1],f_lambda[,2],type="l",col="black",xlab="f-value", ylab="Growth factor")
lines(f_lambda[,1],rep(line_of_growth,fmax),col="red")
plot(Ntot_f[,1],Ntot_f[,2],type="l",col="black",xlab="f-value", ylab="Total population at stable equilibrium (t=500)")


