
#clear variables etc
rm(list=ls())
library(fitdistrplus)
load("scripts/week6/SARS.Rdata")
mf = par(mfrow=c(2,2))
#week 6 assignment
#p[k]=probability of having k offspring

#for including quarintine effects: 
#p.PAC[k]=(1-q)*p[k] for k>0 and p.PAC[0]=p[0]+q*(1-p[0])
#for HPC effects 
#g.HPC(s)=g(1-q+q*s) if g(s) is the probability of generating function without control

#---------------
#Part 1
#set parameters

q=0.4
tf=10
p.pac=numeric(34)
for(i in 1:34)p.pac[i]=(1-q)*sum(SARS==(i-1))/34 #looking at this setup:
#for p[i], it sums all the values of SARS that equal i-1, then divides by 34. this gives the probability for of seeing i offspring?
p.pac #probability of seeing 1,2,3 etc offspring...

outbreak.pac=function(tf=10){
  N=matrix(1,tf,1)
  for(t in 2:tf){
    if(N[t-1]>0){
      N[t] = sum(sample(0:33,size=N[t-1],replace=T,prob=p.pac))
    }else{N[t]=0}
  }
  return(N)
}
#plot(outbreak.pac(),type="l")

reps=10000
all.p.pac=matrix(NA,tf,reps)

for(i in 1:reps)all.p.pac[,i]=outbreak.pac(tf=tf)
#matplot(all.p.pac+1,type="l",lty=1,log="y")

#need to set hist up for all.p.pac where 10th generation is above 0
N=which(all.p.pac[10,]>0)
P.PAC=all.p.pac[,N]
P.PAC=colSums(P.PAC)
hist(P.PAC, xlab="Population at generation 10")

#for HPC at each step we know how many people are infected, for each of those we sample randomly from the SARS vector (which tells us how many they contacted). Add those together and take a rbinomial across that total (of 1-q, or probability of each met individual getting infected).
q=0.4
tf=10

p=numeric(34) #what is this doing? Making a vector of zeros length 34
for(i in 1:34)p[i]=sum(SARS==(i-1))/34
p #probability of seeing 1,2,3 etc offspring...

outbreak.hpc=function(tf=10){
  N=matrix(1,tf,1)
  for(t in 2:tf){
    if(N[t-1]>0){
      temp = N[t-1]+sum(sample(0:33,size=N[t-1],replace=T,prob=p.pac))
      N[t]=rbinom(n=1,size=temp,prob=.6)
        }
    else{N[t]=0}}
  return(N) #I think this is working?
}
#plot(outbreak.hpc(),type="l")

reps=10000
all.p.hpc=matrix(NA,tf,reps)

for(i in 1:reps)all.p.hpc[,i]=outbreak.hpc(tf=tf)
#matplot(all.p.hpc+1,type="l",lty=1,log="y") #something isn't right here, none of the trials go extinct

#need to set hist up for all.p.pac where 10th generation is above 0
N=which(all.p.hpc[10,]>0)
P.HPC=all.p.hpc[,N]
P.HPC=colSums(P.HPC)
hist(P.HPC,xlab="Population at generation 10")

#----------------
#Part 2
##Note for part B: notes on site are wrong: g.HPC(s)=g(q+(1-q)*s)
#g.PAC(s)=sum of 0-33 for p.pac(k)s^k


g.PAC=function(s)sum(p.pac*s^c(0:33))

q.ppac=matrix(NA,tf,2)
q.ppac[1,2]=0
qe <- seq(0.1,1,by=0.1)
#for(i in 1:34)p.pac[i]=(1-q)*sum(SARS==(i-1))/34
#i think i need a double for loop here?
for(j in (1:length(qe))){
  q=qe[j]
  for(i in 1:34){
    p.pac[i]=(1-q)*sum(SARS==(i-1))/34
    }#what horror have I wrought?
  temp=q.ppac[j,2]
  q.ppac[j+1,2]=g.PAC(temp)# I'm having issues for some reason populating this
  q.ppac[j+1,1]=q}
q.ppac #Why is my chance for extinction decreasing as my q increases?

plot(x=q.ppac[,1],y=q.ppac[,2],xlab="q value",ylab="probability of extinction by 10th generation") 


#HPC

g=function(s)q+(1-q)*s

q.phpc=matrix(NA,tf,2)
q.phpc[1,2]=0
qe <- seq(0.1,1,by=0.1)
#for(i in 1:34)p.pac[i]=(1-q)*sum(SARS==(i-1))/34
#i think i need a double for loop here?
for(j in (1:length(qe))){
  q=qe[j]
  temp=q.phpc[j,2]
  q.phpc[j+1,2]=g(temp)# I'm having issues for some reason populating this
  q.phpc[j+1,1]=q}
q.phpc #Why is my chance for extinction decreasing as my q increases?

plot(x=q.phpc[,1],y=q.phpc[,2],xlab="q value",ylab="probability of extinction by 10th generation") #this one seems right
