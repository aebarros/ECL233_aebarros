# The fecundity kernel
# first, did you have kids or not
model.flower=glm(fec.flower~size,data=data,family=binomial)
# create the flowering function
F.flower=function(x)predict(model.flower,data.frame('size'=x),type='response')
# second, how many kids given that I flowered?
# create a new version of the data that eliminates the non-flowering
data2=data.frame(size=data$size[parents],fec.seed=data$fec.seed[parents])
model.kids=glm(fec.seed~size,data=data2,family=poisson)
F.kids=function(x)predict(model.kids,data.frame('size'=x),type='response')
# plot these
par(mfrow=c(1,2))
plot(data$size,data$fec.flower,pch=21,bg=rgb(1,0,0,0.5))
lines(xs,F.flower(xs))
plot(data$size[parents],data$fec.seed[parents],pch=21,bg=rgb(1,0,0,0.5))
lines(xs,F.kids(xs))
#estimate the establishment probability
est.prob=sum(is.na(data$size))/sum(data$fec.seed,na.rm=TRUE)
mean.size=mean(data$sizeNext[is.na(data$size)])
# estimate the distribution of sizes in first year
sd.size=sd(data$sizeNext[is.na(data$size)])
# Putting together the fecundity function
F.all=function(y,x)F.flower(x)*F.kids(x)*est.prob*dnorm(y,mean=mean.size,sd=sd.size)