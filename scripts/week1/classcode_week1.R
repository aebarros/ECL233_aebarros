#script lesson1 4/2/2019
#introduction: philosophy of mathematical modeling, basics of R

##########Basics of Modeling###########
#what is a model?: some simplyfied view of reality. mathmatical models formalize conceptual models
#why do we model?: to help us think clearly, also used to help make predictions
  #we can control everything in our model
  #a model should be as simple as possible but no simpler
  #explain it in algebra

#Models say what CAN happen, observational data says what DOES happen

#Make sure to comment while you command
rm(list = ls()) #this clears your environment

#sequence function
x=seq(0,1,by=0.1)
x
x[5]
x[5:7]
length(x)
max(x)
y=exp(x)
y+2
x*y #mutliplies element by element
plot(x,y,type='l',xlab='x',ylab='exp(y)')

#functions and functional breakdown make your code cleaner
#preplanning of the function:
  #goal: plot exponential growth for different growth rates and starting population sizes
  #n(t)=n(0)*e^(rt) is the solution to dn/dt=rn. this is for unregulated pop growth in continues time
  #r = instantaneous growth rate
   #so we will need a function of the above equation
   #inputs: growth rates, starting pop size, and time series
  #make plots

#plotExp for different growth rates
