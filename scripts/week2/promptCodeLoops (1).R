# class 2: learning loops
# goals for the day: learn different types of loops, how to simulate discrete-time models

# there are three main types of loops

# if loops: if an action is true, preform it
x = -1
x>=0
if(x>=0) {y = x} # there can always be multiple actions in the curly brackets
y
x = 1
x>=0
if(x>=0) {y = x}
y
# can also put in an else
x = -3
if(x>=0) {y = x} else {y = -x} # this is the code for absolute value
y
# there are also switch statements, which are basically a formalized version of doing repetitive if/else statements
# you can put multiple criteria in your if statement, for example:
if(x>=0 & y>=0) {z=1} # &=and: only do the action if all statements are true
z
if(x>=0 | y>=0) {z=1} # |=or: do the action if at least one of the statments is true
z
# also, you can apply booleans (TRUE/FALSE) to vectors, for example:
x = -5:5
x>0
# And TRUE=1 and FALSE=0 in expressions, for example:
(x>0)*4
# this can come in handy with vectorization (we'll get to than with for loops)
# The any/all commands are like or/and applied to an entire vector
any(x>0) # TRUE if the test is true for at least one entry, FALSE otherwise
all(x>0) # TRUE if the test is true for all entries, FALSE otherwise

# while loops: preform an action while something is true
# say, run the Beverton-Holt model N_t+1 = alpha*N_t/(1+beta*N_t)
# we're going to program it as a function because we're going to keep using it
# here's how to put in default values (do in a script)
bh = function(nt, alpha = 1.2, beta = 0.001) {
  nt1 = alpha*nt/(1+beta*nt)
  return(nt1)
}
# enter function and at prompt:
bh(2) # can use function with just nt
bh(2,alpha=2) # can change any or all of default values
bh(2,alpha=2,beta=0.01)
bh(2,beta=0.01,alpha=2) # with changed default values, can put in different order
# make sure it's working: try it out for values where you know the answer -- here, the equilibrium
alp = 2
bet = 0.01
Neq = (alp-1)/bet # expected equilibrium value
bh(Neq,beta=bet, alpha=alp)
Neq # should be the same
# plot n_t+1 vs. n_t to get sense of function
beta=0.01
vn = 1:10000 
plot(vn,bh(vn),type="l",xlab="n_t",ylab="n_t+1") 

# now let's use this in a while loop: run bh through time to equilibrium
# in a script:
minDiff = 5#0.01# # minimum difference that signals reached equilibrium (start w/5 to demonstrat with print, then change to 0.01 and run through plot)
#tmax = 10#1000 # maximum time allowed to run -- first demonstrate without -- then first do with 10 then 1000 to show first one reached breaks the loop
n0 = 100 # initial population size
n = rep(NA,2) # initialize n vector over time
n[1] = n0 # start with two individuals
n[2] = bh(n[1]) # do next one -- you'll see why in a sec
t = 2 # start at time 2 (done two steps)
while((n[t]-n[t-1])>minDiff){# & t<tmax) { # if greater than minDiff, keep running (at first -- then add & t<tmax part)
  n[t+1]=bh(n[t]) # find bh-value at next time step
  t=t+1 # update time step
  print(n) # show what's happening (comment out after first run)
  print(t)
}
# change minDiff to 0.01, comment out print statements, then re-run and plot to show:
plot(1:t,n,xlab="Time",ylab="Population size")
# BEWARE OF INFINITE LOOPS!
# here: in case don't reach equilibrium, should put in maximum time (change in script to include tmax and add to loop, re-run)

# for loops: perform an action for each value of x from 1 to xf
rm(n) # clear n to start over
tf = 5#100# # final time (instead of having a max n)
n = rep(NA,tf) # it's always a good idea to initiate vectors at full length for most efficient memory usage - note that we can do this with a for loop, as we know how many times it'll run, but not for a while loop, as we're uncertain of the number of iterations
n[1] = n0 # we know the first entry
for(t in 1:(tf-1)) { # for each time step (start t=1, does what's in the loop, then t=1, all the way up to tf-1)
  n[t+1]= bh(n[t]) # run function
  print(t) # showing you t and n for demonstration
  print(n)
}
# run with tf=5, then comment out print statements, run to tf=100
plot(1:tf,n,xlab="Time",ylab="Population size")

# for loops will probably be the ones you end up using the most
# they'll also be the ones you can avoid with "vectorization"
# remember that r can preform functions like exp to vectors:
exp(1)
exp(1:10)
# you could have also gotten to the exp(1:10) with a for loop:
y=rep(NA,10)
for(x in 1:10){y[x]=exp(x)}
y
# but it's much faster just to apply exp to a vector
# you can vectorize if the next one doesn't rely on the previous one -- e.g., exp(10) doesn't depend on exp(9)
# but when the next one depends on the previous one -- like in your Beverton-Holt model -- you can't vectorize
# say, for the Beverton-Holt model we wanted to look at trajectories with a variety of different initial population sizes
# -- back to the script --
vn0 = seq(10,300,by=10) # vector of initial population sizes
Nn0 = length(vn0) # how many population sizes looking over
# we're going to need to initialize our n trajectories as a matrix, with a row for each initial population size
n = matrix(NA,nrow=tf,ncol=Nn0)
# this gives us a matrix
# --- at the prompt - a little about matrices
M = matrix(1:10,nrow=2,byrow=TRUE)
M
M[1,] # first row
M[,2] # second column
M[1,1:5] # first row, columns 1-5
M[,-1] # everything but the first column
M[,1] = 10 # re-assign the first column
M
# we could do a double for loop: first loop through n0's, then loop through time
# i.e.:
ntest = n
for(stepn0 in 1:Nn0) {
  ntest[1,stepn0] = vn0[stepn0]
  for(t in 1:(tf-1)) {
    ntest[t+1,stepn0] = bh(ntest[t,stepn0])
  }
}
# but instead let's take a "vectorized" shortcut: 
# -- back to the script:
# just loop through time, but apply model to vector, with entry for each n0
n[1,] = vn0 # initialize population size for all columns, first row
for(t in 1:(tf-1)) { # loop through time
  n[t+1,] = bh(n[t,]) # fill in next time step for next row, all columns
}
matplot(1:tf,n,type="l",xlab="Time",ylab="Population size") # matplot is a way to plot multiple lines at once: can have to matrices plotted against each other (each x,y entry) or one x vector and a matrix of columns for each line corresponding to that x vector
# or, for cases with fewer lines, we could do subplots
mf = par(mfrow=(c(2,1))) # set up subplots
for(stp in 1:2) {plot(1:tf, n[,stp], type="l", xlab="Time", ylab="N",main=paste("n0=",vn0[stp]))} # note paste to put together value with string in title
par(mfrow=mf) # stop doing subplots

# your task:
# two plots using Ricker model (N_t+1=N_t*exp(r*(1-N_t/K))):
# 1. plot N vs. t for four different values of r: 1.5, 2.3, 2.6, 3.0 (show different behaviors - one-cycle, two-cycle, four-cycle, chaos)
# 2. bifurcation plot for r from 1.5 to 3.6
# class discussion: what is a bifurcation plot? what is chaos? why get more chaotic with faster growth?
# biological insight: vast array of possible behaviors from simple model

# you'll have to use a for loop to iterate through time
# BUT you can do this most efficiently if you can figure out how to "vectorize" over multiple values of r

# ideas on the general approach? (pseudocode together)

# other parameter suggestions:
# K=100
# N0=50
# tf=100 in 1st plot, 500 in 2nd
# tp=how many time points at end to extract for 2nd


# ------------------------------
# some add-ons
# different types of matrix multiplication are useful for vectorization in different ways:
v = 1:2 
w = 3:4
#  matrix product (or inner product)
v%*%w
# outer product -- creates array of multiplied combinations:
v%o%w
# Kronecker product -- creates vector of multiplied combinations: 
v%x%w
# can see how would be useful to apply calculation across different conditions or construct arrays to send into vectorized calculations
# can use with vector of ones to repeat a vector -- useful for when doing combinations of multiple parameters
v%x%rep(1,3) # as vector
v%x%matrix(1,1,3) # as matrix
matrix(rep(v,3),ncol=3,byrow=FALSE) # same thing
# you might even want to formalize the formation of these types of matrices (a vector repeated a certain number of times) with the function (name taken from Matlab):
repmat = function(a,m) {matrix(rep(a,m), ncol=m, byrow=FALSE)}

