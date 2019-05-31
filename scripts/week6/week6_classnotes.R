###Week 6 class notes
#Goals:
#1) basic random umber commands
#2) model selection with maximum likelihood
#3) simulating branching processes

#clear variables etc
rm(list=ls())

#binomial distribution
#repeating KPs experiment of flipping a coin 24,000 times
rbinom(n=1,size=24000,prob=.5) #11884 heads on first run
#now we want to know what is the probability of getting 12,012 heads?
dbinom(x=12012,size=24000,prob=.5) #dbinom tells us probability of a certain outcome, .005 in this case
#what is the chance of getting within 12 of 12000?
pbinom(q=12012,size = 24000,prob=0.5)-pbinom(q=11988,size=24000,prob=.5) #tells us probability of getting between 11988-12012

#Central Limit Theorem illustration
me=c(6,42,11,0,1/5)

#sample X times from this distribution with replacement
X=100
many.me=sample(x=me,size=X,replace=T)
mean(many.me) #sample mean
mean(me) #expected mean or population mean

#do the above over many repititions
Y=10000 #number of reps
my.means=numeric(Y) #makes a vector of length Y
for(i in 1:Y)my.means[i]=mean(sample(x=me,size=X,replace=T)) #fills the vector with means of the samples

hist(my.means,50,freq=F) #fre1=F gives us the density of the distribution
#this looks like a normal distribution
#the Central Limit Theorem says that if you take the means of a large enough sample size, the distribution of those means will be a normal distribution

#plot an appropriate normal density over this
M=mean(me)
V=mean((me-M)^2) #true variance
SD=sqrt(V) #standard deviation of population
SD.sample.mean = SD/sqrt(X) #SD of a sample of the population (our samples are size X=100 in this case)

xs=seq(M-3*SD.sample.mean,M+3*SD.sample.mean,length=200)
ys=dnorm(xs,mean=M,sd=SD.sample.mean)

lines(xs,ys,lwd=3,col="red")

#poisson distributions
#a poisson is a distribution of a binomial variable where there are many 0s?

#fit poisson and negbinomial distributions to some data
load("scripts/week6/SARS.Rdata") #these numbers represend an outbreak of SARS in Beijing with a misdiagnosed patient. You have a patient 0 that infected 33 individuals, each which infected a certain number. of individuals

library(fitdistrplus)
model1=fitdist(data=SARS,distr="pois",method="mle") #this asks what lambda maximizes something
model2=fitdist(data=SARS,distr="nbinom",method="mle")

#plot the fits against the data
cdfcomp(list(model1,model2),lwd=3,legendtext=c("pois","nbinom")) #compares cumulative distributions against fitted distributions?
#looking ag above plot, the negative binomial model fits the data better

#check AICs? Akaki Information Criteria
model1$aic
model2$aic #much lower than model 1, so thats better apparently
#so select model2

#simulate disease outbreaks (or a population growing)
#branching process
#idea is N[t] = # infected at generation t
# offspring distribution p[k] = prob. of k infected
# N[t+1] = sum of N[t] random draws from the infected distribution

outbreak = function(tf = 10){
  #assuming N[1]=1
  N=matrix(1,tf,1)
  for(t in 2:tf){
    if(N[t-1]>0){
    N[t] = sum(sample(SARS,size=N[t-1],replace=T))
    }else{N[t]=0}
  }
  return(N)
}
plot(outbreak(),type="l") #if you plot many times you get different values

tf=10
reps=100
all.outbreak=matrix(NA,tf,reps)

for(i in 1:reps)all.outbreak[,i]=outbreak(tf=tf)
matplot(all.outbreak+1,type="l",lty=1,log="y")

#so we can ask, what is the probability no more infected by time 10?
sum(all.outbreak[tf,]==0)/reps #92% of runs disease went extinct

#the analytic approach to extinction probability
#probability generating functions (PGF)
#g(s)=sum[p[k]*s^k] where k is the number of offspring
#theory tells us if N[1]=1 then probability of extinction by time tf+1 is g(g(...g(0))) with tf compositions?

#do the above in R here:
#we need our pgf

p=numeric(34) #what is this doing? Making a vector of zeros length 34
q=.1 #why, as q increases, is my chance of extinction decreasing?
for(i in 1:34)p[i]=(1-q)*sum(SARS==(i-1))/34
p #probability of seeing 1,2,3 etc offspring...
g=function(s)sum(p*s^c(0:33))
#lets compute the probability of extinction by time tf
temp=0
for(i in 1:9)temp=g(temp)
temp #0.8977957
