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

#Build first matrix for Nt (age strucutre numbers through time tf=50)
Nt=matrix(NA,nstg,tf) #starting empty matrix to fill with pop sizes through time
Nt[,1]=n0 #fill in initial population vector
Nt

leval=function(L) Re(eigen(L)$values[1]) #function find leading eigen value of a matrix
#calculate expected long-term growth factor

#maybe make a new TS_fecundity vector based on different values of f?
f=1
TS_fecun=matrix(NA,fmax,2)
TS_fecun[1,1]=f
TS_fecun
Ntot=0 #for density independence
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
#here we need to make a matrix of Ntot at t=500 for different values of f 1:75
#initial parameters
TS<-as.matrix(read.table("data/teasel_stage.txt")) #upload teasel data and convert from table into a matrix
nstg = nrow(TS) #tells us number of stages
n0 = c(2,rep(0,nstg-1)) #sets starting pop
n0
vt = 1:tf #vector of times
a=.01
fmax=75
tf=500
a=.01
n0=c(2,0,0,0,0,0,0)
#build a first matrix to collect population amounts through time
Nt=matrix(NA,nstg,tf) #starting empty matrix to fill with pop sizes through time
Nt[,1]=n0 #fill in initial population vector
Nt
#build a second matrix to collect total population of age group at t=500 across f
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
mtext("The affect of change of f-value on lambda and stable equilibrium pop size", line =-2, outer = TRUE)



#Questions:
#Are my plots right?
#Am I understanding this right?: f=fecundity, dividing it by (1+a*Ntot) at time t is adding in density dependence
#Why does the TS Leslie matrix have more values in 1st column?
#The initial starting Nt of (2,0,0,0,0) must be unrealistic, as there is high likelihood of both individuals not reaching maturity. Can we model the system so that fractions of individuals are rounded down, so if 0 reach maturity the population crashes to 0, and not just near zero?
  #Answer: view 2 as a density for unit of area


#notes/script from open lab:
#objective: function that creates matrix A as a function of f and ntot
tf=500
fmax=75
fun=function(f=75,ntot=0){#ntot=0 to be density independent
  A=TS
  A[1,7]=f/(1+a*ntot)
  return(A)
}
#objective:function that produces lambda for a given f
lam.f=function(f=75){
  A=fun(f=f,ntot = 0)
  return(leval(A))
}
#objective:input f,n0,tf output: ntot at tf (second plot)
ntot.tf=function(f=75,n0=c(2,0,0,0,0,0,0),tf=500){
  n=n0
  for(t in 1:(tf-1)){
    ntot=sum(n)
    n=fun(f,ntot)%*%n
  }
  return(sum(n))
}

#matrix for plot1
#parameters
fmax=75
f.lam=matrix(NA,fmax,2)
f.lam[1,1]=0
f.lam

for(f in (1:fmax-1)){
  f.lam[f+1,1]=f+1
  f.lam[f+1,2]=lam.f(f)
}

#matrix for plot2
tf.ntot=matrix(NA,fmax,2)
tf.ntot[1,1]=2
tf.ntot

for(f in(1:fmax-1)){
  tf.ntot[f+1,1]=f+1
  tf.ntot[f+1,2]=ntot.tf(f)
}
##Looks right!!!!
#Make plots!

#1st plot of lambda vs f
mf = par(mfrow=c(1,2))
line_of_growth=1
plot(f.lam[,1],f.lam[,2],type="l",col="black",xlab="f-value", ylab="Growth factor")
lines(f.lam[,1],rep(line_of_growth,fmax),col="red")
plot(tf.ntot[,1],tf.ntot[,2],type="l",col="black",xlab="f-value", ylab="Total population at stable equilibrium (t=500)")
mtext("The affect of change of f-value on lambda and stable equilibrium pop size", line =-2, outer = TRUE)

#lets try an extra plot
#objective: calculate Ntot over time for four different values of f
#objective:input f,n0,tf output: ntot at tf (second plot)
#build empty matrix
t.ntot=matrix(NA,tf,5)
t.ntot[1,]=c(1,2,2,2,2)
n=n0
t=1
n=matrix(NA,tf,7)
n[1,]=n0

ntot.fun=function(f=75,n0=c(2,0,0,0,0,0,0),tf=500){
  for(t in 1:(tf-1)){
    ntot=sum(n[t,])
    n[t+1,]=fun(f,ntot)%*%n[t,]
    #return(sum(n[t+1,]))
  }
    return(sum(n[t+1,]))
}
for(t in (1:tf-1)){
t.ntot[t+1,1]=t+1
t.ntot[t+1,2]=ntot.fun(f=1,t)
t.ntot[t+1,3]=ntot.fun(f=25)
t.ntot[t+1,4]=ntot.fun(f=50)
t.ntot[t+1,5]=ntot.fun(f=75)
}
