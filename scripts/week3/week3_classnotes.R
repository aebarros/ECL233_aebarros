#week3 Classnotes Discrete Time Models II
#Talking about Leslie matrixes
#start with a (Leslie matrix) multiplied by a (vector of population at time i) that will end up giving us a (vector of number of individuals at each age group for time i+1)

v=1:2
w=3:4

v*w #this just multiplies one to one, not the matrix algebra we will want to use

v%*%w #

v%o%w #creates the outer multiplaction into a matrix array?

v%x%w #same as above, but as a vector

v%x%rep(1,3) #take v and multiplies by (1,1,1) and creates an array

#making a Leslie Matrix
L = rbind(c(1/2,1,3/4),c(2/3,0,0),c(0,1/3,0)) #top row is fecundity, second and third rows are surivoship to next age class
L
n0=c(2,0,0) #initial population size

L%*%n0 #this runs a Leslie Matrix giving us a number of individuals at next time step

tf=10 #set how many times steps we want to run this
Nt= matrix(NA,nrow(L),tf) #set up our future matrix to fill in, currently blank
Nt

Nt[,1] = n0 #fill in first column of matrix with initial conditions

for(t in 1:(tf-1)) Nt[,t+1] = L%*%Nt[,t] #this runs a for loop up to tf
Nt

colSums(Nt) #sums each of the age group columns and tells us total populations size through time

#Leslie matrix times a vector = eigen value times the vector?
#what is the eigen value?
#if a L matrix is 100x100, you have 100 eigen value eigen vector pairs?

#The biggest eigen value will swamp out other values, which will eventurally determine long term expected growth rate. This will also settle us into the eigen vector that corresponds to leading Eigen vector, which corresponds to our stable age structure.

#So if largest Eigenvalue is <1, expect population to decline in long term, if largest Eigen value is >1, expect to increase in long term

#fortunately we don't have to create a function to calculate eigen values and eigen vectors
ev=eigen(L) #eigen is function in base r
ev #$values gives us all eigen values (with largest first), $vectors gives us all eigen vectors (in this case three of each)

ev$values [1] #pulls out largest eigen value
Re(ev$values[1]) #pulls out largest eigen value without the annoying imaginary numbers
Re(ev$vectors[,1]) #pulls out leading eigen vector without imaginary numbers
