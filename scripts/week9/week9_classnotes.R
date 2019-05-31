#Class Notes Week 9
#Stochastic population dynamics of Bay Checkerspot Butterflies
rm(list=ls())
load("checkerspot.Rdata")
coef=c(0.932171,-4.348866e-05,-3262.42)

#lecture ntoes: the last two weeks we've been thinkign about demographic stochasticity, finit number of individuals with uncorrelated fates
#now we are going to focus on env stochasticity
#so in class we will be looking at a simpler single species model, for hw we will be looking at multiple species using a lottery model
#increase in precipitation variation is thought to have been why butterfly populations went extinct locally
#the larvae of the butterfly eat a key host plant, so in dry years the host plant does poorly and can't sustain larvae, in really wet years plants grow to fast for larvae to take advantage
#So we are going to test the precipiation variation hypothesis

#look at data
plot(pop$year,pop$N,type='b',pch=21,bg='red',log='y')
plot(rainfall$year,rainfall$precip,type='b',pch=21,bg='red',log='y')

#build a model that will make use of the rain data
f=function(N,W)N*exp(coef[1]+coef[2]*N+coef[3]*W^(-2)) #update rule saying what the effect of rainfall W is on the next years population

#low density (0) per-capita growth rate
r=function(W)(coef[1]+coef[2]+coef[3]*W^(-2))

#based on chalkboard notes which washed over me we should be able to fugure out when the population is doing well or poorly at low densities

index.pre=which(rainfall$year<1971)
index.post=which(rainfall$year>1970)

#compute low density per-capita growth rate for both scenarios

r.pre=mean(r(rainfall$precip[index.pre])) #this gives us a positive number, meaning at low density growth rates the population grows
r.post=mean(r(rainfall$precip[index.post])) #negative growth rate

#lets go to the full model

#goal is a function that takes initial N, length of run tf, some seqment of weather data W, and number of reps to run the model
  #output a matrix of simulations of length tf

N0=pop$N[1]
tf=100
reps=2
W=rainfall$precip[index.pre]

N=matrix(N0,tf,reps) #matrix to fill
for(t in 1:(tf-1)){
  W.sample=sample(x=W,size=reps,replace=TRUE)
  N[t+1,]=f(N[t],W.sample)
}

matplot(N,type='l',log='y')

#turn above into function
mclaughlin=function(N0=pop$N[1],tf=100,reps=2,W=rainfall$precip[index.pre]){
  N=matrix(N0,tf,reps) #matrix to fill
  for(t in 1:(tf-1)){
    W.sample=sample(x=W,size=reps,replace=TRUE)
    N[t+1,]=f(N[t,],W.sample)
  }
  return(N)
}

matplot(mclaughlin(reps=100,W=rainfall$precip[index.pre]),type='l',log="y") #can switch between index.pre and index.post to compare precip values before and after 1/1/1971

#compare single run vs ensemble views

ensemble=mclaughlin(reps=10000,tf=200)
single=mclaughlin(reps=1,tf=10000)

#plot historgrams of log densities and comparte
par(mfrow=c(1,2))
hist(log10(ensemble[tf,]))
hist(log10(single),xlim=c(-15,5))
