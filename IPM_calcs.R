
############survival functions!
s_x<-function(size.x,binom.intercept,binom.slope){
  pred.survival=plogis(binom.intercept+binom.slope*size.x)
  return(pred.survival)
}


#growth function
g_yx<-function(size.y,size.x,intercept,slope,sd.val){
  
  growth.prob=dnorm(size.y,mean=intercept+slope*size.x,sd=sd.val)
  return(growth.prob)
}
#
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




get.lambda<-function(mat) {
  
  if(any(is.nan(mat)==T)) {return(NA)}
  else{
    
    if(any(is.infinite(mat)==T)) {return(NA)}
    
    else{
      if(any(is.na(mat)==T)) {return(NA)}
      
      else{
        
        return(Re(eigen(mat)$values[1]))}}}}
# 
# 

#code for class:
myimage<-function(x,sp,colpal=rainbow, ...){
  MIN <- min(x,na.rm=T)
  MAX <- max(x,na.rm=T)
  yLabels <- round(as.numeric(rownames(x)),2)
  xLabels <- round(as.numeric(colnames(x)),2)
  title <-c()
  
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rev(colpal(256)[1:200])
  ColorLevels <- seq(MIN, MAX, length=length(ColorRamp))
  
  # Reverse Y axis
  #reverse <- nrow(x) : 1
  yLabels <- yLabels #[reverse]
  x <- x #[reverse,]
  
  # Data Map
  #t(x)
  par(mar = c(5.3,5.3,2.5,2))
  image(1:length(xLabels), 1:length(yLabels),x, col=ColorRamp,
        axes=FALSE, zlim=c(MIN,MAX),main=paste(sp," pop growth rate","=",Re(round(eigen(x)$values[1],2))),
        xlab="Size at time t",ylab="Size at time t+1")
  
  abline(0,1,lwd=3)
  
  box()
  
  
  if( !is.null(title) ){
    title(main=title)
  }
  xR<-xLabels[seq(from=1,to=length(xLabels),by=30)]
  yR<-yLabels[seq(from=1,to=length(xLabels),by=30)]
  
  #for some reason, this is going from 1:192. I guess its by pixels, not by 
  #rownames
  
  axis(BELOW<-1, at=seq(from=1,to=dim(x)[1],length=length(xR)), labels=xR, cex.axis=1.09,las=2)
  axis(LEFT <-2, at=seq(from=1,to=dim(x)[1],length=length(xR)), labels=yR, las= HORIZONTAL<-1,
       cex.axis=1.09)
  
  
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
  par(mfrow=c(1,1))
}