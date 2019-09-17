library(truncnorm)
library(cubature)
library(metRology)


##################################################################
#####################Example 2.2##################################
##################################################################

#In the absence of the selection event, with a constant prior
#the posterior follows N(ybar,1/n)
mu <- (-1) #true value of mean
n <- 10 #sample size
nn <- 100000 #number of simulations

yn<-rnorm(nn,mean=mu,sd=sqrt(1/n))
qtile<-qnorm(0.9,mean=yn,sd=sqrt(1/n))

#proportion of the simulated '90% posterior quantile'>= -1 

mean(qtile >= -1)

#An alternative way to compute this proportion
#table 2.1

nn <- 100000  #number of simulations

Q<- rep(0, times = nn)  #posterior quantile

ybar<-rnorm(nn,mean=mu,sd=sqrt(1/n))
Q<-pnorm(-1,mean=ybar,sd=sqrt(1/n))

mean(Q <= 0.9)
mean(Q <= 0.95)
mean(Q <= 0.99)
hist(Q,main="Non-selective, n=10")

##################################################################
#####################Example 2.3##################################
##################################################################
#In the absence of the selection event, with a prior 1/sigma
#the marginal posterior of mu follows t(n-1, bary, S^2)
mu <- (-1) #true value of mean
n <- 10 #sample size

s <- 100000 #number of simulations
ybar<-rnorm(s,mean=mu,sd=sqrt(1/n))
yvar<-rgamma(s,shape=(n-1)/2, rate = (n-1)/2)
qtile<-qt.scaled(0.9,df=n-1, mean=ybar,sd=sqrt(yvar)/sqrt(n))

mean(qtile >= -1)
hist(qtile)

#An alternative way to compute this proportion
#table 2.2

Quan <- rep(0, times = s)
ybar<-rnorm(s,mean=mu,sd=sqrt(1/n))
yvar<-rgamma(s,shape=(n-1)/2, rate = (n-1)/2)
Quan<-pt.scaled(-1,df=n-1, mean=ybar,sd=sqrt(yvar)/sqrt(n))

mean(Quan <= 0.9)
hist(Quan,main="Non-selective, n=10")

##################################################################
#####################Example 2.4##################################
##################################################################
#A standard normal prior leads to the posterior N(n*ybar/(n+1),1/(n+1))
#Linear regression on computing the convergence rate 

Sdnormalprior<-function(n,q){
  s=1000000;mu=-1
  yn<-rnorm(s, mean=mu,sd=sqrt(1/n))
  ProbSet<-pnorm(-1,mean=(n*yn/(1+n)),sd=sqrt(1/(1+n)))
  return(mean(ProbSet<= q)-q)
}
#Figure 2.1
a<-log(seq(5,500,5))
for(t in seq(5,500,5)){b[t/5]<-Sdnormalprior(t,0.9)}
plot(a,b,xlab="log(n)",ylab="log(error)")
convergence <- lm(b ~ a)
summary(convergence)
abline(lm(b ~ a),col="blue")

##################################################################
#####################Example 2.5##################################
##################################################################

#known sigma^2=1;without random responses #constant prior;selection ybar>0 
Integral.Posterior<-function(mu,n,s,quantile){ 
  ybar<-rtruncnorm(s, a=0, b=Inf, mean = mu, sd = sqrt(1/n))
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for( i in 1:s){
    trunclili2<-function(mu)  {
      exp(dnorm(sqrt(n)*(ybar[i]-mu),0,1, log = T) - (pnorm(sqrt(n)*(mu), 0, 1, log.p = T)))
    }
    u[i]<-integrate(trunclili2,lower=-100,upper=mu)$value
    d[i]<-integrate(trunclili2,lower=-100,upper=mu)$value+integrate(trunclili2,lower=mu,upper=Inf)$value 
    p[i]<-u[i]/d[i]
  }
  print(summary(p))
  print(mean(p <=quantile))
  hist(p)
}

Integral.Posterior(-0.5,100,100000,0.99)


##################################################################
#####################Example 2.6##################################
##################################################################
#unknown sigma=1; without random responses #prior 1/sigma
par(mfrow=c(1,3))

mu=-1;n=10 #selection ybar1>0
ybar1<-rtruncnorm(1, a=0, b=Inf, mean = mu, sd = sqrt(1/n))

trunclili8<-function(m, sig){
  (1/sig)*exp(dnorm(ybar1, m, sig/sqrt(n), log = T)-(pnorm(sqrt(n)*m/sig, 0, 1, log.p = T)))}

x <- seq(-300, 50, length = 70)#mu
y <- seq(0,50, length = 50)#sigma
z <- outer(x, y, function(a, b) trunclili8(a, b))
nrz <- nrow(z);ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("blue", "pink") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
persp(x, y, z, col = color[facetcol], phi = 30, theta = 0, ticktype = "detailed",xlab = expression(mu), ylab = expression(sigma),main="mu=-1")
persp(x, y, z, col = color[facetcol], phi = 30, theta = 90, ticktype = "detailed",xlab = "mu", ylab = "sigma",main="mu=-1")
persp(x, y, z, col = color[facetcol], phi = 15, theta = 30, ticktype = "detailed",xlab = "mu", ylab = "sigma",main="mu=-1")


Integral.Posterior_U3<-function(mu,n,s,quantile){
  ybar1<-rtruncnorm(s, a=0, b=Inf, mean = mu, sd = sqrt(1/n))
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for(i in 1:s){
    trunclili8<-function(m,sig){
      (1/(sig^2))*exp(dnorm(sqrt(n)*(ybar1[i]-m)/sig, 0, 1, log = T)-(pnorm(sqrt(n)*(m)/sig, 0, 1, log.p = T)))
    }
    uu <- function(mu){ #calculate the double integral
      sapply(mu, function(a){
        integrate(function(sigma){
          trunclili8(a, sigma)
        }, lower = 0.1, upper = 1000)$value})
    }
    u[i]<-integrate(uu,lower=-100,upper=mu)$value
    d[i]<-u[i]+integrate(uu,lower=mu,upper=100)$value 
    p[i]<-u[i]/d[i]
    print(i)
  }
  print(summary(p))
  print(mean(p <= quantile))
  hist(p)
}

#unknown sigma=1; without random responses #prior exp(-sigma)
par(mfrow=c(1,3))

mu=-1;n=10
ybar1<-rtruncnorm(1, a=0, b=Inf, mean = mu, sd = sqrt(1/n))

trunclili8<-function(m,sig){
  exp(-sig)*exp(dnorm(ybar1, m, sig/sqrt(n), log = T)-(pnorm(sqrt(n)*(m)/sig, 0, 1, log.p = T)))}

x <- seq(-80, 20, length = 100)#mu
y <- seq(0,30, length = 100)#sigma
z <- outer(x, y, function(a, b) trunclili8(a, b))
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette( c("blue", "pink") )
nbcol <- 100
color <- jet.colors(nbcol)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
persp(x, y, z, col = color[facetcol], phi = 30, theta = 0, ticktype = "detailed",xlab = expression(mu), ylab = expression(sigma),main="mu=-1")
persp(x, y, z, col = color[facetcol], phi = 30, theta = 90, ticktype = "detailed",xlab = "mu", ylab = "sigma",main="mu=-1")
persp(x, y, z, col = color[facetcol], phi = 15, theta = 30, ticktype = "detailed",xlab = "mu", ylab = "sigma",main="mu=-1")


Integral.Posterior_U1<-function(mu,n,s,quantile){
  ybar1<-rtruncnorm(s, a=0, b=Inf, mean =mu, sd = sqrt(1/n))
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for(i in 1:s){
    trunclili9<-function(m,sig){
      (exp(-sig))*exp(dnorm(ybar1[i],m,sig/sqrt(n),log = T)-(pnorm(sqrt(n)*m/sig, 0, 1, log.p = T)))
    }
    uu <- function(mu){ #calculate the double integral
      sapply(mu, function(a){
        integrate(function(sigma){
          trunclili9(a, sigma)
        }, lower = 0, upper = Inf)$value})
    }
    u[i]<-integrate(uu,lower=-Inf,upper=mu)$value
    d[i]<-u[i]+integrate(uu,lower=mu,upper=Inf)$value 
    p[i]<-u[i]/d[i]
  }
  print(summary(p))
  print(mean(p <= quantile))
  hist(p)
}

##################################################################
#####################Example 2.7##################################
##################################################################
#standard normal prior, selective case

Integral.Posterior_sn<-function(mu,n,s,quantile){ #An integral of posterior from -Inf to true mu
  yn_1<-rtruncnorm(s, a=0, b=Inf, mean = mu, sd = sqrt(1/n))
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for( i in 1:s){
    trunclili4<-function(m)  {
      exp(dnorm(sqrt(1+n)*(n*yn_1[i]/(1+n)-m),0,1, log = T) - (pnorm(sqrt(n)*(m), 0, 1, log.p = T)))
    }
    u[i]<-integrate(trunclili4,lower=-Inf,upper=mu)$value
    d[i]<-integrate(trunclili4,lower=-Inf,upper=mu)$value+integrate(trunclili4,lower=mu,upper=Inf)$value 
    p[i]<-u[i]/d[i]
  }
  print(summary(p))
  print(mean(p <=quantile))
  hist(p)
}
Integral.Posterior_sn(-1,100000000000,10000,0.9)

par(mfrow=c(1,2))
s=15;mu=-1;n=10
yn_1<-rtruncnorm(s, a=0, b=Inf, mean = mu, sd = sqrt(1/n))
plot(seq(-10,5,0.01),trunclili4(seq(-10,5,0.01)),type="l",main="",xlab=expression(mu), ylab="Posterior density")
summary(yn_1)
for( i in 1:s){
  trunclili4<-function(m)  {
    exp(dnorm(sqrt(1+n)*(n*yn_1[i]/(1+n)-m),0,1, log = T) - (pnorm(sqrt(n)*(m), 0, 1, log.p = T)))
  }
  lines(seq(-10,5,0.01),trunclili4(seq(-10,5,0.01)),type="l",col=i+11)
}

abline(v=-1)

y11<-rnorm(s, mean = mu, sd = sqrt(1/n))
plot(seq(-10,5,0.01),trunc(seq(-10,5,0.01)),type="l",main="",xlab=expression(mu), ylab="Posterior density")

for( i in 1:s){
  trunc<-function(m)  {
    exp(dnorm(sqrt(1+n)*(n*y11[i]/(1+n)-m),0,1, log = T))
  }
  lines(seq(-10,5,0.01),trunc(seq(-10,5,0.01)),type="l",col=i+11)
}
##################################################################
#####################Example 3.1##################################
##################################################################
Integral.Posterior_trun<-function(mu,n,s,quantile){ 
  #An integral of posterior from -Inf to true mu
  yn_1<-rtruncnorm(s, a=0, b=Inf, mean = mu, sd = sqrt(1/n))
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for( i in 1:s){
    posterior<-function(mu)  {
      A<-(n^(3/2)*mu*exp(-(n*mu^2)/2)/(sqrt(2*pi)))/(exp(pnorm(sqrt(n)*(mu), 0, 1, log.p = T)))
      B<-(n/(2*pi))*exp(-n*mu^2)/(exp(pnorm(sqrt(n)*(mu), 0, 1, log.p = T)))^2
      prior<-sqrt(n-(A+B))
      return(prior*exp(dnorm(sqrt(n)*(yn_1[i]-mu),0,1, log = T) - (pnorm(sqrt(n)*(mu), 0, 1, log.p = T))))
    }
    u[i]<-integrate(posterior,lower=-8.5,upper=mu)$value
    d[i]<-u[i]+integrate(posterior,lower=mu,upper=Inf)$value 
    p[i]<-u[i]/d[i]
  }
  return(mean(p <=quantile))
  hist(p)
  #print(summary(p))
}


##################################################################
#####################Example 4.1##################################
##################################################################
#constant prior; selection on randomised data
#beta denotes the variance of each y_i
Integral.Posterior_random<-function(mu,alpha,beta,n,s,quantile){ 
  U<-rtruncnorm(s, a=0, b=Inf, mean = mu, sd = sqrt((beta+alpha)/n))
  V<-rnorm(s,mean = mu, sd = sqrt((beta+alpha)/n))
  Y<-(1/2)*(U+V)
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for( i in 1:s){
    trunclili2<-function(x)  {
      exp(dnorm(sqrt(n/beta)*(Y[i]-x),0,1, log = T) - (pnorm(sqrt(n/(beta+alpha))*(x), 0, 1, log.p = T)))
    }
    u[i]<-integrate(trunclili2,lower=-Inf,upper=mu)$value
    d[i]<-integrate(trunclili2,lower=-Inf,upper=mu)$value+integrate(trunclili2,lower=mu,upper=Inf)$value 
    p[i]<-u[i]/d[i]
  } 
  return(mean(p <= quantile))
}
Integral.Posterior_random(-1,1,1,5,10000,0.9)
##################################################################
#####################Example 4.2##################################
##################################################################

#In non-selective context, with prior exp(-sigma)
Integral.Posterior_U00<-function(mu,n,s,quantile){
  y1<-rnorm(s, mean = mu, sd = sqrt(1/n))
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for(i in 1:s){
    trunclili8<-function(m,sig){ 
      (exp(-sig))*exp(dnorm(y1[i],m,sqrt(n)/sig, log = T))
    }
    uu <- function(mu){ #calculate the double integral
      sapply(mu, function(a){
        integrate(function(sigma){
          trunclili8(a, sigma)
        }, lower = 0, upper = Inf)$value})
    }
    u[i]<-integrate(uu,lower=-Inf,upper=mu)$value
    d[i]<-u[i]+integrate(uu,lower=mu,upper=Inf)$value 
    p[i]<-u[i]/d[i]
    if(i==5000){print(p[i])}
  }
  print(summary(p))
  print(mean(p <= 0.9))
  print(mean(p <= 0.95))
  print(mean(p <= quantile))
  hist(p)
}

Integral.Posterior_U00(-1,50,2000,0.99)

# Another way of computing the above integral
Integral.Posterior_U11<-function(mu,n,s){
  y1<-rnorm(s, mean = mu, sd = sqrt(10/n))
  u<-rep(0,s)
  for(i in 1:s){
    trunclili88<-function(sig){
      (exp(-sig))*exp(pnorm(mu, y1[i],sig/sqrt(n), log = T))
    }
    u[i]<-integrate(trunclili88,lower=0,upper=Inf)$value
  }
  print(summary(u))
  print(mean(u <= 0.9))
  print(mean(u <= 0.95))
  print(mean(u <= 0.5))
  hist(u)
}
Integral.Posterior_U11(-1,1000,100000)

#In selective context, with prior exp(-sigma) and noises with variance alpha/n
Integral.Posterior_U5<-function(mu,alpha,n,s){
  U<-rtruncnorm(s, a=0, b=Inf, mean = mu, sd = sqrt((1+alpha)/n))
  V<-rnorm(s,mean = mu, sd = sqrt((1+alpha)/n))
  Y<-(1/2)*(U+V)
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for(i in 1:s){
    trunclili9<-function(m,sig){
      (exp(-sig))*exp(dnorm(Y[i],m,sig/sqrt(n),log = T)-(pnorm(sqrt(n)*(m)/sqrt(sig^2+alpha), 0, 1, log.p = T)))
    }
    uu <- function(mu){ #calculate the double integral
      sapply(mu, function(a){
        integrate(function(sigma){
          trunclili9(a, sigma)
        }, lower = 0, upper = Inf)$value})
    }
    u[i]<-integrate(uu,lower=-Inf,upper=mu)$value
    d[i]<-u[i]+integrate(uu,lower=mu,upper=Inf)$value 
    p[i]<-u[i]/d[i]
  }
  print(summary(p))
  print(mean(p <= 0.9))
  print(mean(p <= 0.95))
  print(mean(p <= 0.99))
  hist(p)
}
Integral.Posterior_U5(1,1,100,5000)
Integral.Posterior_U5(0,1,10,5000)
Integral.Posterior_U5(-1,1,10,5000)
Integral.Posterior_U5(-1,18,10,1000)

re1<-Integral.Posterior_U5(-1,0.9,10,1000)


#In non-selective context, with prior exp(-sigma)
# and the likelihood by the product of i.i.d. N(\mu,sigma^2)
Integral.Posterior_U12<-function(mu,n,s){
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for(i in 1:s){
    y <- rnorm(n, mean = mu, sd = 1)
    trunclili88<-function(sig){
      (exp(-sig)) * sig ^ (-n + 1) * exp(-n / (2 * sig ^ 2) * (sum(y ^ 2) / n - mean(y) ^ 2)) * pnorm(mu, mean = mean(y), sd = sig / sqrt(n))
    }
    
    trunclili89<-function(sig){
      (exp(-sig)) * sig ^ (-n + 1) * exp(-n / (2 * sig ^ 2) * (sum(y ^ 2) / n - mean(y) ^ 2))
    }
    
    u[i]<-integrate(trunclili88,lower=0,upper=Inf)$value
    d[i]<-integrate(trunclili89,lower=0,upper=Inf)$value
    p[i]<-u[i]/d[i]
  }
  print(summary(p))
  print(mean(p <= 0.9))
  print(mean(p <= 0.95))
  print(mean(p <= 0.99))
  hist(p)
  #return(mean(p <= 0.9))
}
Integral.Posterior_U12(-1,50,100000)# First order probability matching


#In selective context, with prior exp(-sigma)
#and the likelihood by the product of i.i.d. N(\mu,sigma^2)
#This integral diverges !
Integral.Posterior_U13<-function(mu,n,s,quantile){
  u<-rep(0,s)
  d<-rep(0,s)
  p<-rep(0,s)
  for(i in 1:s){
    y <- rnorm(n, mean = mu, sd = 1)
    count <- 0
    while(mean(y) <= 0){
      y <- rnorm(n, mean = mu, sd = 1)
      count <- count + 1
    }
    ybar1<-mean(y) #compute the sample mean bar{y}
      trunclili13<-function(m,sig){
      (exp(-sig)) * sig ^ (-n) *exp(-n / (2 * sig ^ 2) *(m-ybar1)^2) *exp(-n / (2 * sig ^ 2) * (sum(y ^ 2) / n - mean(y) ^ 2))/ pnorm(sqrt(n)*m/sig, 0, 1)
    }
    uu <- function(mu){ #calculate the double integral
      sapply(mu, function(a){
        integrate(function(sigma){
          trunclili13(a, sigma)
        }, lower = 1, upper = 100)$value})
    }
    u[i]<-integrate(uu,lower=-100,upper=mu)$value
    d[i]<-u[i]+integrate(uu,lower=mu,upper=100)$value 
    p[i]<-u[i]/d[i]
    print(i)
  }
  print(summary(p))
  print(mean(p <= quantile))
  hist(p)
}
Integral.Posterior_U13(-1,10,100,0.9)

#plot its joint posterior density
mu=-1;n=10

  y <- rnorm(n, mean = mu, sd = 1)
  count <- 0
  while(mean(y) <= 0){
    y <- rnorm(n, mean = mu, sd = 1)
    count <- count + 1
}
  ybar1<-mean(y)
  trunclili13<-function(m,sig){
  (exp(-sig)) * sig ^ (-n) *exp(-n / (2 * sig ^ 2) *(m-ybar1)^2) *exp(-n / (2 * sig ^ 2) * (sum(y ^ 2) / n - mean(y) ^ 2))/ pnorm(sqrt(n)*m/sig, 0, 1)
}

x <- seq(-100, 100, length = 100)#mu
yy <- seq(0,100, length = 100)#sigma
z <- outer(x, yy, function(a, b) trunclili13(a, b))
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette( c("blue", "pink") )
nbcol <- 100
color <- jet.colors(nbcol)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
persp(x, yy, z, col = color[facetcol], phi = 30, theta = 0, ticktype = "detailed",xlab = expression(mu), ylab = expression(sigma),main="mu=-1")
persp(x, yy, z, col = color[facetcol], phi = 30, theta = 90, ticktype = "detailed",xlab = "mu", ylab = "sigma",main="mu=-1")
persp(x, yy, z, col = color[facetcol], phi = 15, theta = 30, ticktype = "detailed",xlab = "mu", ylab = "sigma",main="mu=-1")


