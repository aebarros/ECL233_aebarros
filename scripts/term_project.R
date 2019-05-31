#Term Project
#model Laplacian dispersal of two hypothetical sessile aquatic organisms species along a salinity gradient fit over a 1-dimensional linear transect. Both species are then fit into a discrete time Lotka-Volterra competition model [1]. Species X will have a lower salinity tolerance, and growth rate will decrease as salinity increases. The growth rate of species Y will be unaffected by higher salinity values but will be competitively excluded by Species X at lower salinity values.
#Species x is more competitive, with a higher intitial growth rate and higher competitive value

#things to add if time:
  #stochastic salinity shifts DONE!
rm(list=ls())
library(deSolve)
source("scripts/ecl233_functions.R") #source function script

#--------------
#functions
#dispersal kernal k(y,x) = exp(-abs(y-m-x)/a)/(2*a)
prob.disp<-function(y,x){
  exp(-abs(y-m-x)/a)/(2*a)
}
#where m is the mean displacement (0 in this case)
#where 2a^2 is the variance of displacement

#dispersal for each species:
Xprob.disp<-function(y,x){
exp(-abs(y-m-x)/Xa)/(2*Xa)
}
Yprob.disp<-function(y,x){
  exp(-abs(y-m-x)/Ya)/(2*Ya)
}

#competition equations
#Xn[t+1] = Xr*Xn[t]*(1-Xn[t]/Kx-a*Yn[t]/Kx)
Xgrowth=function(Xr,Xn,Kx,a,Yn,t){
  Xn1=Xn+Xr*Xn*(1-(Xn/Kx)-(a*Yn/Kx))
  return(Xn1)
}

#Yn[t+1] = Yr*Yn[t]*(1-Yn[t]/Ky-b*Xn[t]/Ky)
Ygrowth=function(Yr,Yn,Ky,b,Xn,t){
  Yn1=Yn+Yr*Yn*(1-(Yn/Ky)-(b*Xn/Ky))
  return(Yn1)
}

#Growth rate of species X (Xr) is decreased in higher salinity so that 
Xr.f=function(Xr0,ms,S)Xr0*(1-(S/ms)) #where ms is maximum salinity tolerance and S is local salinity
#Yr is not impacted by salinity

#---------------
#parameters
m=0 #default value used to create dispersal matrix
Xa=2 #value used to determine dispersal distance of X
Ya=2 #value used to determine dispersal distance of Y
Xr0=.6 #initial growth rate of X at 0 salinity
Yr=.3 #inisial growth rate of Y, not affected by salinity
Kx=100 #carrying cap of X
Ky=100 #carrying cap of Y
a=.3 #competition on X from Y
b=1 #competition on Y from X, X is more competitive at lower salinity values
ms=15 #maximum salinity tolerance of X, before growth rate stops
tf=1000 #length of run
L=20 #sets transect length
bins=41 #may need to change this
S=0 #default salinity value to test

#-----------------
#build dispersal matrixes
#function to create A based on different L
dx=(2*L)/bins #set the delta of x
A.create=function(L,prob.a){ #this funcion creates a dispersal probability matrix based on a value of L and using a probability dispersal function
  x=seq(-L+dx/2,L-dx/2,length=bins) #create vector x from last midpoint to first midpoint
  At=outer(x,x,prob.a) #use the outer function  to create matrix At that is x by x in length, filled by probability of dispersal from any x to anyother x
  return(At)
}
disp_X=A.create(L,Xprob.disp) #for now they are the same because Xa=Ya
disp_Y=A.create(L,Yprob.disp)

#create salinity values spread across transect for normal, dry, and wet years across different months
norm.salinity=seq(from=0, to=32,length.out = 2*L+1)
dry.salinity=rep(32,times=2*L+1)
dry.salinity[1:21]=seq(from=0, to=32, length.out = 21)
wet.salinity=rep(0,times=2*L+1)
wet.salinity[21:41]=seq(from=0, to=32, length.out = 21)

norm.y.salinity=matrix(NA,12,bins)
norm.y.salinity[0:2,]=matrix(rep(wet.salinity,times=2),ncol=bins,byrow = T)
norm.y.salinity[3:8,]=matrix(rep(norm.salinity,times=6),ncol=bins,byrow = T)
norm.y.salinity[9:11,]=matrix(rep(dry.salinity,times=3),ncol=bins,byrow = T)
norm.y.salinity[12,]=matrix(rep(wet.salinity,times=1),ncol=bins,byrow = T)

dry.y.salinity=matrix(NA,12,bins)
dry.y.salinity[0:2,]=matrix(rep(wet.salinity,times=2),ncol=bins,byrow = T)
dry.y.salinity[3:6,]=matrix(rep(norm.salinity,times=4),ncol=bins,byrow = T)
dry.y.salinity[7:12,]=matrix(rep(dry.salinity,times=6),ncol=bins,byrow = T)

wet.y.salinity=matrix(NA,12,bins)
wet.y.salinity[0:6,]=matrix(rep(wet.salinity,times=6),ncol=bins,byrow = T)
wet.y.salinity[7:10,]=matrix(rep(norm.salinity,times=4),ncol=bins,byrow = T)
wet.y.salinity[11:12,]=matrix(rep(wet.salinity,times=2),ncol=bins,byrow = T)


S.sal=rbind(norm.salinity,dry.salinity,wet.salinity) #create matrix of seasonal salinity options

#build Xr based on salinity spread for fixed salinity
Xr=Xr.f(Xr0,ms,S=norm.salinity)

#lets first just try and do the two populations without salinity
#create emptry Xn and Yn matrices
Xn=matrix(NA,tf,bins)
Xn[1,1:3]=10 #fill the first vector/timestep (is there a simpler way of doing this?)
Xn[1,4:41]=0

Yn=matrix(NA,tf,bins)
Yn[1,1:38]=0 #fill the first vector/timestep (is there a simpler way of doing this?)
Yn[1,39:41]=10

t=1
for(t in 1:(tf-1)){ #for loop that fills in matrices Xn and Yn
  Xn[t+1,]=disp_X%*%(Xgrowth(Xr,Xn[t,],Kx,a,Yn[t,],t)*dx)
  Yn[t+1,]=disp_Y%*%(Ygrowth(Yr,Yn[t,],Ky,b,Xn[t,],t)*dx)
}

#------------------
#plot matrices
#first set a row 1 for matrix column numbers
x.points<-norm.salinity #sets up points to use for x-axis
Yb=t(Yn)
Xb=t(Xn)#take tranverse of the n matrix to plot it (not sure why I have to do this, but it seems to work)
matplot(cbind(Yb,Xb),type="l",xlab="salinity",ylab="density",col="black",axes=F) #axes=F drops axis labels to replace
matlines(Xb[,100], col="green", lwd=2)
matlines(Yb[,100], col="red", lwd=2)
axis(2) #sets back y-axis
axis(side=1,at=1:nrow(Yb),labels = x.points) #sets x-axis

#plot nonstochastic model for total population of each species
Xb.sum=rowSums(Xn)
Yb.sum=rowSums(Yn)
plot(1, type="n", xlab="Months", ylab="N",main="Fixed Salinity", xlim=c(0, tf), ylim=c(0, max(Yb.sum,Xb.sum)))
matlines(Xb.sum,col="green")
matlines(Yb.sum,col='red')

#---------------
#So we have our normal, dry and wet year salinitiy spreads set above, now would like to add some yearly stochasticity over the transect. 
#Goal: recreate the above model so that every 12 generations (months in this case) the salinity vector has a chance to change between one of the three options, or stay the same.
#note: the final plot will have to have transect distance on the x-axis as salinity will be changing throughout the run
#note: I would also like to have a way of filling a matrix that tells me if its a normal, dry, or wet season every month

#create emptry Xn and Yn matrices
Xn.s=matrix(NA,tf,bins)
Xn.s[1,1:3]=10 #fill the first vector/timestep (is there a simpler way of doing this?)
Xn.s[1,4:41]=0

Yn.s=matrix(NA,tf,bins)
Yn.s[1,1:38]=0 #fill the first vector/timestep (is there a simpler way of doing this?)
Yn.s[1,39:41]=10

#probability matrix of moving between season types
#1=norm, 2=dry, 3=wet
P=rbind(c(0.7, 0.15, 0.15),
        c(0.8, 0.15, 0.05),
        c(0.8, 0.05, 0.15))
Z=numeric(tf) #create list to track seasonal salinity shifts
Z[1]=1

Year=numeric(tf) #Create Year matrix to populate and use to calculate months later
Year[1]=0

months=seq(from=1,to=12,length.out = 12)

k=dim(P)[1]#create k as the number of dimmesnions in seasonal shift matrix
for(t in 2:tf){ #for loop that fills in matrices Xn and Yn
  Year[t+1]=ifelse(t%%12==0,Year[t]+1,Year[t]) #this populates the Year matrix with what year we are in for the run
  month=(t-(12*Year[t])) #tells the loop what month we are on (1-12)
  Z[t]=ifelse(t%%12==0,sample(1:k,size=1,prob=P[Z[t-1],]),Z[t-1])#pick salinity year index (1,2, or 3)
  if(Z[t]==1){ #nested if else statements to select the proper water year
    S.sal=norm.y.salinity
  }else{
    if(Z[t]==2){
      S.sal=dry.y.salinity
    }else{
      S.sal=wet.y.salinity
    }
  }
  Xr=Xr.f(Xr0,ms,S=S.sal[month,]) #builds the X growth rate off of new Salinity spread
  Xn.s[t,]=disp_X%*%(Xgrowth(Xr,Xn.s[t-1,],Kx,a,Yn.s[t-1,],t)*dx)
  Yn.s[t,]=disp_Y%*%(Ygrowth(Yr,Yn.s[t-1,],Ky,b,Xn.s[t-1,],t)*dx)
}

#IT WORKED! PLOT IT UP!
x.points<-seq(from=1, to=41,length.out = 2*L+1) #sets up points to use for x-axis
Yb.s=t(Yn.s)
Xb.s=t(Xn.s)#take tranverse of the n matrix to plot it (not sure why I have to do this, but it seems to work)
matplot(cbind(Yb.s,Xb.s),type="l",xlab="River km",ylab="density",col="black",axes=F) #axes=F drops axis labels to replace
matlines(Xb.s[,tf], col="green", lwd=2)
matlines(Yb.s[,tf], col="red", lwd=2)

axis(2) #sets back y-axis
axis(side=1,at=nrow(Yb):1,labels = x.points) #sets x-axis

#next goal, plot total populations at each time step for each species and show plot through time
Xn.t=rowSums(Xn.s)
Yn.t=rowSums(Yn.s)
y.max=ifelse(max(Xn.t)>max(Yn.t),max(Xn.t),max(Yn.t))
plot(Xn.t,type='l',col="green",ylim=c(0,y.max),ylab="N",xlab="Month",main="Water Year Variation")
lines(Yn.t,col="red")

#determine what percent of years are norm, dry, wet
pct.n=100*sum(Z==1)/tf
pct.d=100*sum(Z==2)/tf
pct.w=100*sum(Z==3)/tf

#can I run the stochastic for loop for many many reps?
reps=100
XY.t=matrix(NA,tf,2*reps)

for(r in 1:reps){
  for(t in 2:tf){ #for loop that fills in matrices Xn and Yn
    Year[t+1]=ifelse(t%%12==0,Year[t]+1,Year[t]) #this populates the Year matrix with what year we are in for the run
    month=(t-(12*Year[t])) #tells the loop what month we are on (1-12)
    Z[t]=ifelse(t%%12==0,sample(1:k,size=1,prob=P[Z[t-1],]),Z[t-1])#pick salinity year index (1,2, or 3)
    if(Z[t]==1){ #nested if else statements to select the proper water year
      S.sal=norm.y.salinity
    }else{
      if(Z[t]==2){
        S.sal=dry.y.salinity
      }else{
        S.sal=wet.y.salinity
      }
    }
    Xr=Xr.f(Xr0,ms,S=S.sal[month,])
    Xn.s[t,]=disp_X%*%(Xgrowth(Xr,Xn.s[t-1,],Kx,a,Yn.s[t-1,],t)*dx)
    Yn.s[t,]=disp_Y%*%(Ygrowth(Yr,Yn.s[t-1,],Ky,b,Xn.s[t-1,],t)*dx)
    Xn.t=rowSums(Xn.s) #catches total pop of each time step
    Yn.t=rowSums(Yn.s)
    }
  XY.t[,r]=Xn.t #stores total pop of each time step in first half of matrix for Xn.s
  XY.t[,r+reps]=Yn.t #stores total pop of each time step in second half of matrix for Yn.s
  }
max(XY.t)
plot(1, type="n", xlab="Months", ylab="N",main="Water Year Variation 100 reps", xlim=c(0, tf), ylim=c(0, max(XY.t)))
matlines(XY.t[,1:reps],col="green")
matlines(XY.t[,(reps+1):(2*reps)],col='red')


#One final thing: Try and figure out what percentage of the reps each species is on top
species.avg=matrix(NA,1,(2*reps)) #create 1 row matrix to fill with column averages from XY.t
for(i in 1:(2*reps)){
  species.avg[,i]=mean(XY.t[,i])
}
species.winner=rep(0,times=reps)
for(i in 1:reps){
  species.winner[i]=ifelse(species.avg[,i]>species.avg[,i+reps],"X","Y")
}
100*sum(species.winner=="X")/reps
