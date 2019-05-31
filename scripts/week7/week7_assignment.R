#week 7 assignment
#Recreate intrinsic Mean Extinction Time (MTE) curves similar to Figure 3 of the Melbourne and Hastings (2008) reading. Specifically, let the mean field model be given by N[t+1]=N[t]*R*exp(-alpha*N[t]).Create a plot where R varies between 1.1 and 30 on the horizontal axes and the log of MET is plotted on the vertical axes. Plot 6 curves corresponding to the Poisson distribution and the negative binomial with size parameter values 0.5,1,10,50,100. As in Figure 3 of the paper, choose the alpha values so that the positive equilibrium of the model is always 30. alpha has to change with R to keep equilibrium population at 30

#for this we will use dpois and dnbinom

#formula to change alpha with R to keep equilibrium around 30
#alpha= -(1/30)*log(1/R)
aR=function(R) {
  a=log(R)/30
}

#field mean model set up with a static R
R=2 #if R is too high it gets to stochastic and crashes
a=aR(R)
#N[2]=N[1]*R*exp(-a*N[1])
fit=function(R,a,N){ #ricker model
  N*R*exp(-a*N)
}

N1=50 #let's try this with a random starting population
tf=1000
N=numeric(tf)
N[1]=N1
for(t in 2:tf){
  N[t]=rpois(n=1,lambda= fit(R,a,N[t-1]))
}
N[1:20]
plot(N,type="l") #with these numbers, this population doesn't often go extinct

#calculating MTE for fixed R value
chop=200
Q=matrix(NA,chop,chop)
for(i in 1:chop){
  Q[i,]=dpois(x=1:chop,lambda=fit(R,a,i))
}
stuff=eigen(t(Q))
v=Re(stuff$vectors[,1])
v=v/sum(v)
v
sum(v) # =1 so it's right

hist(N,freq=FALSE,col=rgb(1,0,0,0.5))
lines(v,lwd=2)

#we already found 1, where 1/(1-lambda) where lambda is the leading eigenvalue of Q

1/(1-Re(stuff$values[1])) #"intrinsic mean time to extinction" gives us expected time to extinction 700 million years
#here lambda is the probability of not going extinct which is VERY close to 1

#now want to find mean time to extinction given i initial individuals
#this is given by (Id-Q)^{-1}%*%ones where ones is a column vector of ones

extinct=solve(diag(chop)-Q)%*%matrix(1,chop,1) #gives inverse of ID matrix?
plot(1:100,extinct[1:100],type="b") #gives us a plot of mean time to extinction based on initial starting population

#So what we need is a function to calculate MTE from different R values
MTE=function(R){
  a=aR(R)
  chop=200
  Q=matrix(NA,chop,chop)
  for(i in 1:chop){
    Q[i,]=dpois(x=1:chop,lambda=fit(R,a,i))
  }
  stuff=eigen(t(Q))
  ext=1/(1-Re(stuff$values[1]))
  return(ext)
}
MTE(2)



#----------------
#now lets try for dbinom, we are doing a seperate function for each size variation. There is probably a more efficient way of doing this, but this is the classic coding dilemna for me: do I take a little extra time to repeat myself and do many different functions? Or do I take a really long time to figure out a more efficient method? I chose the former.
MTEdnbinom.5=function(R){
  a=aR(R)
  chop=200
  Q=matrix(NA,chop,chop)
  for(i in 1:chop){
    Q[i,]=dnbinom(x=1:chop,size=.5,mu=fit(R,a,i))
  }
  stuff=eigen(t(Q))
  ext=1/(1-Re(stuff$values[1]))
  return(ext)
}
MTEdnbinom1=function(R){
  a=aR(R)
  chop=200
  Q=matrix(NA,chop,chop)
  for(i in 1:chop){
    Q[i,]=dnbinom(x=1:chop,size=1,mu=fit(R,a,i))
  }
  stuff=eigen(t(Q))
  ext=1/(1-Re(stuff$values[1]))
  return(ext)
}
MTEdnbinom10=function(R){
  a=aR(R)
  chop=200
  Q=matrix(NA,chop,chop)
  for(i in 1:chop){
    Q[i,]=dnbinom(x=1:chop,size=10,mu=fit(R,a,i))
  }
  stuff=eigen(t(Q))
  ext=1/(1-Re(stuff$values[1]))
  return(ext)
}
MTEdnbinom50=function(R){
  a=aR(R)
  chop=200
  Q=matrix(NA,chop,chop)
  for(i in 1:chop){
    Q[i,]=dnbinom(x=1:chop,size=50,mu=fit(R,a,i))
  }
  stuff=eigen(t(Q))
  ext=1/(1-Re(stuff$values[1]))
  return(ext)
}
MTEdnbinom100=function(R){
  a=aR(R)
  chop=200
  Q=matrix(NA,chop,chop)
  for(i in 1:chop){
    Q[i,]=dnbinom(x=1:chop,size=100,mu=fit(R,a,i))
  }
  stuff=eigen(t(Q))
  ext=1/(1-Re(stuff$values[1]))
  return(ext)
}
MTEdnbinom(2)

#Now a for loop to populate a matrix for just the poison
Rs=seq(from=1.1,to=30, by=0.1)
RMTE=matrix(NA,length(Rs),7)
for(j in 1:length(Rs)){
  R=Rs[j]
  RMTE[j,1]=R
  RMTE[j,2]=(MTE(R))
  RMTE[j,3]=(MTEdnbinom.5(R))
  RMTE[j,4]=(MTEdnbinom1(R))
  RMTE[j,5]=(MTEdnbinom10(R))
  RMTE[j,6]=(MTEdnbinom50(R))
  RMTE[j,7]=(MTEdnbinom100(R))
}
matplot(RMTE[,1],log(RMTE[,-1]),type="l",xlab="R value",ylab="MTE log years")
legend('topright', legend=c("Poisson", "NegBin 0.5","NegBin 1","NegBin 10","NegBin 50","NegBin 100"), col=c("black", "hotpink", "cyan","blue","green","red"), lty=c(1,1)) #build the legend
