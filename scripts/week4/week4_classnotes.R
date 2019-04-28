#class notes 4/23/19
rm(list=ls())

#Goal: build an integral projection model (IPM) using the provided data

#read in data
D<-read.csv("data/IPM-data.csv")

#look at the data
plot(D$size,D$surv) #looks like survival increases with size
plot(D$size,D$sizeNext) #looks like the bigger you are, the bigger you are next time step, although variance increases with size
plot(D$size,D$fec.flower) #as plants get bigger, more likely to have flowers

parents=which(D$fec.flower==1)
plot(D$size[parents],D$fec.seed[parents]) #plots size of those that made seeds with the number of seeds produced (bigger plants had more seeds)

#build model components
surv.model=glm(surv~size,data=D,family=binomial) #making a glm for survival vs size, using binomial model because survival is binary
summary(surv.model) #look at main thing, the coefficients. Intercept=-2.6, slope is positive

#make a function that uses x value (size) and surival model we made to predict survival,type="response" keeps data in original units
survival=function(x)predict(surv.model,data.frame('size'=x),type="response")

#plot function against the data using choosen x values
temp=range(cbind(D$size,D$sizeNext),na.rm=T)
temp #min=.5, max=8.85
a=temp[1] #min size for use
b=temp[2] #max size for use
rm(temp) #clear temporary variable for future use maybe

xs=seq(a,b,length=100)
xs
plot(D$size,jitter(D$surv,factor=0.05),pch=21,bg=rgb(0,0,1,.5),col=rgb(0,0,1,0.5))
lines(xs,survival(xs),lwd=2)

#now to growth
growth.model=glm(sizeNext~size,data=D) #predict sizeNext from size
summary(growth.model) #even though it is significant, we know it doesn't satisfy assumptions under linear regression (normal distribution)
hist(growth.model$residuals) #roughly normal residual distribution

growth.mean=function(x)predict(growth.model,data.frame('size'=x),type="response")
plot(D$size,D$sizeNext,pch=21)
lines(xs,growth.mean(xs),lwd=2,col="red") #looks like it does a reasonable job of predicting

#now we want to put together the components of the model
#create full growth kernel
#assuming (even though its not true) that Standard Deviation doesn't change with size
sigma=sd(growth.model$residuals)
sigma #0.83

#make growth kernel?
growth=function(y,x)dnorm(y,mean=growth.mean(x),sd=sigma)
growth(1,1) #testing above, probability of growing from size 1 to size 1 is .01459

data=D
source("scripts/week4/IPM-Fecundity.R") #this is giving us funtions made earlier in this script

#put it all together
K=function(y,x)survival(x)*growth(x,y)+F.all(y,x) #this is a function combining all the other functions

#discritize the model into a matrix
bins=100
dx=(b-a)/bins #gives us midpoints

xs=seq(a+dx/2,b-dx/2,length=bins) #makes a matrix of our midpoints
xs

A=outer(xs,xs,K)*dx #outer takes a function of x and y and valutes for many xs and ys, and valuates for all pairs
filled.contour(xs,xs,A)

max(Re(eigen(A)$values))
