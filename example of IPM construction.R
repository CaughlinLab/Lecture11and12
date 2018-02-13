#simulate fish datasets

sample.size=10000

#Drawing values from a uniform distribution for
#fish size, a predictor variable
fish.size=runif(sample.size,min=0.1,max=150)

#Setting the true values for the parameters in this model
true.slope=1
true.intercept=1
true.sd.val=6

#Use rnorm to draw random values of fish growth from a normal distribution
fish.sizeT2=rnorm(n=sample.size,mean=true.intercept+fish.size*true.slope,sd=true.sd.val)

#plot a histogram of the data
#note that breaks determines how many bins (histogram "boxes") there are in the plot
hist(fish.sizeT2,col="light blue",breaks=100)

#plot the data
plot(fish.sizeT2~fish.size,pch=21,bg="light blue",xlab="Fish weight (kg)", ylab="Fish growth (kg)",cex.axis=2,cex.lab=2)
#
#fit a statistical model using least squares regression:
m1=lm(fish.sizeT2~fish.size)

#extract the values from the regression
#extract the values from the regression
intercept=coef(m1)[1]
slope=coef(m1)[2]
sd.val=summary(m1)$sigma

#add the statistical results to the plot and observe the difference:
curve(coef(m1)[1]+coef(m1)[2]*x,add=T,col="pink",lwd=8,lty=2)


#simulating binomial data
#true values
true.binom.intercept=2
true.binom.slope=0.1

fish.survival=rbinom(n=sample.size,size=1,p=plogis(true.binom.intercept+true.binom.slope*fish.size))

#jitter is a useful trick to make the data more visible:
plot(jitter(fish.survival)~jitter(fish.size),xlab="Fish weight (kg)",ylab="Fish survival probability",cex.axis=2,cex.lab=2)

m1=glm(fish.survival~fish.size,family="binomial")

#extract coefficients
#extract coefficients
binom.intercept=coef(m1)[1]
binom.slope=coef(m1)[2]

curve(plogis(binom.intercept+x*binom.slope),add=T,col="red",lwd=5)

###########count data
true.poisson.intercept=-5
true.poisson.slope=0.04

fish.babies=rpois(n=sample.size,lambda=exp(true.poisson.intercept+true.poisson.slope*fish.size))

hist(fish.babies,breaks=100,xlab="Number of offspring",
     cex.lab=1.2,cex.axis=1.2,col="forestgreen")

plot(fish.babies~fish.size,xlab="Fish size",ylab="Number of offspring",cex.axis=2,cex.lab=2,
     pch=19,bg="forestgreen",cex=2)
m1=glm(fish.babies~fish.size,family="poisson")

##extract coefficients
##extract coefficients
poisson.intercept=coef(m1)[1]
poisson.slope=coef(m1)[2]

curve(exp(poisson.intercept+poisson.slope*x),add=T,col="blue",lwd=2)

#####ONE LAST PIECE
### WHAT IS THE SIZE OF NEW INDIVIDUALS?

#assume that they are 5 kg 
intercept.recruits=5
sd.recruits=3

############survival functions!
s_x<-function(size.x,binom.intercept,binom.slope){
  pred.survival=plogis(binom.intercept+binom.slope*size.x)
  return(pred.survival)
}

#try with some values!
s_x(100,binom.intercept,binom.slope)
s_x(10,binom.intercept,binom.slope)
s_x(50,binom.intercept,binom.slope)



#growth function
g_yx<-function(size.y,size.x,intercept,slope,sd.val){
  
  growth.prob=dnorm(size.y,mean=intercept+slope*size.x,sd=sd.val)
  return(growth.prob)
}
##try to make a curve!

curve(g_yx(x,100,intercept,slope,sd.val),from=1,to=200,
      xlab="New size",ylab="Probability density",
      main="New size of a 100 kg fish")

##########Pyx

P_yx=function(size.y,size.x,intercept,slope,sd.val,binom.intercept,binom.slope) {
  
  P=s_x(size.x,binom.intercept,binom.slope)*g_yx(size.y,size.x,intercept,slope,sd.val)
  return(P)
}

############### function for reproduction (poisson)
F_yx=function(size.x,size.y,intercept.recruits,sd.recruits,poisson.intercept,poisson.slope) {
  size.y=dnorm(size.y,mean=intercept.recruits,sd=sd.recruits)
  recruit.num=exp(poisson.intercept+poisson.slope*size.x)
  return(size.y*recruit.num)
}

F_yx(size.y=2,size.x=NULL,intercept.recruits,sd.recruits,poisson.intercept,poisson.slope)
F_yx(size.y=1,size.x=NULL,intercept.recruits,sd.recruits,poisson.intercept,poisson.slope)
F_yx(size.y=10,size.x=NULL,intercept.recruits,sd.recruits,poisson.intercept,poisson.slope)

Kyx<-function(size.y,size.x,intercept,slope,sd.val,binom.intercept,binom.slope,
              intercept.recruits,sd.recruits,poisson.intercept,poisson.slope) { #the integral projection kernel
  return(P_yx(size.y,size.x,intercept,slope,sd.val,binom.intercept,binom.slope)+
           F_yx(size.y,size.x,intercept.recruits,sd.recruits,poisson.intercept,poisson.slope)) #NOTE THAT FXY IS THE RECRUITMENT FUNCTION HERE. I have commented it out since you don't have this.
}

#how to construct the kernel?
#how to construct the kernel?
#how to construct the kernel?

#######kernel 1
BINS=100

#U is upper bound on size
#L is lower bound on size


bigmatrix<-function(intercept,slope,sd.val,binom.intercept,binom.slope,
                    intercept.recruits,sd.recruits,poisson.intercept,poisson.slope,L,U) {
  
  
  # boundary points b and mesh points y
  h <- (U - L)/BINS #boundary points
  
  midpts <- L + ((1:BINS) - 1/2) * h
  #set a matrix to loop over
  P1=matrix(NA,nrow=BINS,ncol=BINS)
  
  #set a matrix to loop over
  F1=matrix(NA,nrow=BINS,ncol=BINS)
  
  #loop to construct the matrix
  for(i in 1:BINS){
    for(j in 1:BINS){
      
      P1[i,j]<-h*P_yx(midpts[i],midpts[j],intercept,slope,sd.val,binom.intercept,binom.slope)
      
      F1[i,j]<-h*F_yx(midpts[i],midpts[j],intercept.recruits,sd.recruits,poisson.intercept,poisson.slope)
      
    }
  }
  
  K<-P1+F1
  
  colnames(K)<-midpts
  rownames(K)<-midpts
  
  return(K)
}

#source reads in source code. For this to work,
#the .R file must be in your working directory!
source("IPM_calcs.R")

ipm_kernel=bigmatrix(intercept,slope,sd.val,binom.intercept,binom.slope,
                intercept.recruits,sd.recruits,poisson.intercept,poisson.slope,L=0.1,U=200)

myimage(ipm_kernel,sp="fake fish:")

