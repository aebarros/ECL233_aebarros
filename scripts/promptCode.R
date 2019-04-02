# Class 1 code at the prompt:
# introduction to R
# goals for the day:
# - basic arithmetic
# - plotting
# - scripts & functions

# at it's most basic level, R works as a high-functioning calculator
# for example:
2+2
exp(1)
# to get help...
help(exp)
# ...or use the menu
# note: another thing R is very handy for is statistical analysis, but we won't be covering that very deeply in this class

# you can also assign values
x = 1
x
exp(x)
# also: sqrt, abs, sin, cos, tan, log, log10

# and create plots.  To do this, well need lists of values: vectors
# there are three main ways to make vectors: c, rep, seq
c(1, 2, 3) # the most basic way to make a vector: put the values straight in
rep(1, 3) # make a vector of a value (here 1) repeated n times (here 3)
1:10 # sequence: for spacings of 1
x = seq(0, 1, by=0.1) # sequence" for spacings other than 1 - zero to one in steps of 0.1
x
x[5] # you can find a particular value of a vector
length(x) # you can also determine properties of a vector
max(x) # and find it's maximum and minimum with max and min
y = exp(x) # this evaluates each entry of the vector
y
x*y # this gives you entry-by-entry multiplication, you can do the same with +-/^
x+1 # this adds 1 to every entry
plot(x, y, type="l", xlab="x", ylab="y") # type=line
