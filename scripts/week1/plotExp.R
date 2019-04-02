#Make sure to comment while you command
rm(list = ls()) #this clears your environment

#functions and functional breakdown make your code cleaner
#preplanning of the function:
#goal: plot exponential growth for different growth rates and starting population sizes
#n(t)=n(0)*e^(rt) is the solution to dn/dt=rn. this is for unregulated pop growth in continues time
#r = instantaneous growth rate
#so we will need a function of the above equation
#inputs: growth rates, starting pop size, and time series
#make plots

#plotExp for different growth rates and initial pop sizes

#define funtion
#calculates exponential growth
  #inputs: pop growth rate (r), time (t), and initial pop size (n0)
expGrowth<-function(r,t,n0){
  n=n0*exp(r*t) #calculate exponential growth, r defaults to throwing back the last thing calculated (in this case n)
  return(n) #used to tell r what to throw back at us
}

expGrowth(.2,20,5) #caculates n after 20 cycles for growth rate of .2 with starting pop of 5

#define parameters/inputs, this is useful example of setting input values for use with the function
t=seq(0,10,by=0.1) #time series
r.base =1.2 #baseline growth rate
n0.base = 2 #baseline initial population size
r.double = 2*r.base #doubles the baseline growth rate
n0.double=2*n0.base

#run the function for different values
n.base=expGrowth(r.base,t,n0.base) #base population trajectory
n.r2=expGrowth(r.double,t,n0.base) #pop trajectory with r doubled
n.n02=expGrowth(r.base,t,n0.double) #pop trajectory with douple starting pop

#plot
#pdf(file="exp.pdf") #start a file to save
plot(t,n.r2,type="l",col="blue",xlab="time", ylab="pop size",log="y")#plot line that takes the most space first, note the log scale on the y-axis
lines(t,n.base,col="black") #adds a line for the baseline case to the first plot
lines(t, n.n02, col="red")
legend("topleft",legend=c("Base","r doubled","n0 doubled"),col=c("black","blue","red"),lty=c(1,1,1))
#dev.off() #save the current plot to the pdf

#time to give the assignment a shot!