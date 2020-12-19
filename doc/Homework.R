## ----echo=FALSE---------------------------------------------------------------
n <- 1000
u <- runif(n)
# F(x) = 1-(b/x)^a, x\geqb>0,a>0. 代入Pareto(2,2)得,x = 2*(1-u)^(-1/2)
x <- 2*(1-u)^(-1/2)
hist(x, prob = TRUE, main = expression(f(x) == 8 * x^-3))
y <- seq(0,100,.00001)
lines(y,8*y^(-3) )


## ----echo=FALSE---------------------------------------------------------------
#建立f_e函数，输入n，返回n个符合f_e分布的样本
f_e<-function(n){
  u1<-runif(n,-1,1)
  u2<-runif(n,-1,1)
  u3<-runif(n,-1,1)
  for (i in 1:n) {
    if ( abs(u3[i]) >= abs(u1[i]) && abs(u3[i]) >= abs(u2[i]))
      u[i]=u2[i] else
      u[i]=u3[i]
  }
  return(u)
}
hist(f_e(10000),prob = TRUE)
y<-seq(-1,1,0.01)
lines(y,3/4*(1-y^2))

## ----echo=FALSE---------------------------------------------------------------
r<-4
beta<-2
n<-1000
u<-runif(n)

#由F(y)=1-(beta\beta+y)^r, y>=0.求得以下y的表达式
y<-beta*(1/(1-u)^(1/r)-1)  

#由F(y)=1-(beta\beta+y)^r求得f(x)=r*beta^r/(beta+y)^(r+1)，代入r=4,beta=2可作图
hist(y, prob = TRUE, main = expression(f(x)==64/(2+y)^5)) 

z <- seq(0, 100, .0001)
lines(z, 64/(2+z)^5)

## -----------------------------------------------------------------------------
set.seed(12345)
m <- 1e8; 
x <- runif(m, min=0, max={pi/3})
theta.hat <- mean(sin(x)*pi/3)
print(c(theta.hat, -(cos({pi/3})-cos(0))))


## -----------------------------------------------------------------------------
MC.Phi <- function(x, R = 10000, antithetic = TRUE) {
  u <- runif(R/2,0,x)
  if (!antithetic) 
    v <- runif(R/2) 
  else v <- 1 - u
  u <- c(u, v)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
    g <- x[i] * exp(-(u * x[i])^2 / 2)
    cdf[i] <- mean(g) / sqrt(2 * pi) + 0.5
  }
  cdf
}
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1.95
for (i in 1:m) {
  MC1[i] <- MC.Phi(x, R = 1000, anti = FALSE)
  MC2[i] <- MC.Phi(x, R = 1000)
}


antithetic<-function(x,m=1e4,anti=TRUE){
        theta<-numeric(m)
         if(!anti) {####simple MC method
              t<-runif(m,min=0,max=x)
              theta<-exp(t)*x
         }
        else{  ###antithetic variate approach
                u<-runif(m/2,min=0,max=x)
                v<-x-u
                t<-c(u,v)
                theta<-exp(t)*x
        }
        theta
}
a<-seq(0.1,1,length=10)
theta1<-theta2<-PR<-numeric(length(a))
for(i in 1:length(a)){
        t1<-antithetic(a[i])##antithetic estimator
        t2<-antithetic(a[i],anti=F)##simple MC estimator
        theta1[i]<-mean(t1)
        theta2[i]<-mean(t2)
        PR[i]<-100*(var(t2)-var(t1))/var(t2)###variance reduction
}
print(round(rbind(x, theta1, theta2, PR), 10))

## ----echo=FALSE---------------------------------------------------------------
g<-function(x)x^2/sqrt(2*pi)*exp(-x^2/2)
f1<-function(x)2/(sqrt(2*pi))*exp(-(x-1)^2/2)
f2<-function(x)exp(1-x)

x<-seq(1,10,0.001)
plot(x,g(x),type="l",lwd="2",col="black",main="curves of three functions")  #绘制g(x),f1(x),f2(x)曲线
lines(x,f1(x),lwd="2",col="red")   
lines(x,f2(x),lwd="2",col="blue")   
legend("topright",
       legend =c('g(x)','f1(x)',"f2(x)") ,
       lty=1,
       col=c("black","red","blue"))   



## ----echo=FALSE---------------------------------------------------------------
set.seed(123)
n<-10000
u<-rnorm(n)
v<-rexp(n)
a<-abs(u)+1  
b<-v+1 
theta_hat1<-mean(g(a)/f1(a))
theta_hat2<-mean(g(b)/f2(b))
est<-c(theta_hat1,theta_hat2)
sd<-c(sd(g(a)/f1(a)),sd(g(b)/f2(b)))
func<-c("f1","f2")
mydata<-data.frame(func,est,sd)
mydata

## -----------------------------------------------------------------------------
set.seed(1234)
n<-10000;k<-5
g<-function(x){
  exp(-x)/(1+x^2)*(x>0)*(x<1)
}
f<-function(x){
  exp(-x)/(1-exp(-1))
}

#Method1: Stratified Importance Sampling
s<-n/k
est<-numeric(k)
var<-numeric(k)
for (i in 1:k) {
  u<-runif(s,(i-1)/k,i/k)
  x<--log(1-u*(1-exp(-1))) #inverse transform
  est[i]<-mean(g(x)/f(x)) #estimate in each subinterval
  var[i]<-sd(g(x)/f(x))^2*(s-1)/s #variance in each subinterval
}
lambda_hat1.est<-mean(est)  #estimate of method1
lambda_hat1.sd<-sqrt(sum(var)/s) #sd of method1

#Method2: Importance Sampling
v<-runif(n)
y<--log(1-v*(1-exp(-1))) #inverse transform
lambda_hat2.est<-mean(g(y)/f(y))  #estimate of method2
lambda_hat2.sd<-sd(g(y)/f(y)) #sd of method2

#Compare
c(lambda_hat1.est,lambda_hat2.est)
c(lambda_hat1.sd,lambda_hat2.sd)


## -----------------------------------------------------------------------------
set.seed(1227)
n <- 20
mu <- 1;sigma <- 2 
alpha <- 0.05
CL <- replicate(1000,expr = {
  x <- rlnorm(n,mu,sigma)
  abs(sqrt(n/var(log(x)))*(mean(log(x))-mu)) 
})
sum(CL<qt(0.975,n-1))
mean(CL<qt(0.975,n-1))

## -----------------------------------------------------------------------------
set.seed(12345)
n<-20
lambda<-2;alpha<-0.05
CL<-replicate(1000,expr = { 
  x <- rchisq(n, df = 2)
  abs(sqrt(n/var(x))*(mean(x)-lambda))
})
sum(CL<qt(0.975,n-1))
mean(CL<qt(0.975,n-1))

## -----------------------------------------------------------------------------
alpha <- .05
n<- 30 
m <- 2500

a<- c(seq(.1,3.6,.7),seq(6,49,1))#parameters of beta(a,a) distribution
t<-seq(1,50,1)#degrees of freedom of t distribution

sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

N <- length(a) 
pwr1 <- pwr2 <- numeric(N) 
#critical value for the skewness test 
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
sktests1 <- sktests2 <- numeric(m)
for (j in 1:N) {
  for (i in 1:m) {
    x <- rbeta(n,a[j],a[j])
    y <- rt(n,t[j])
    sktests1[i] <- as.integer(abs(sk(x)) >= cv) 
    sktests2[i] <- as.integer(abs(sk(y)) >= cv) 
    }
  pwr1[j] <- mean(sktests1) 
  pwr2[j] <- mean(sktests2) 
}

#plot power vs parameters of beta distribution
plot(a, pwr1, type = "b", xlab = "parameters of beta distribution", ylim = c(0,0.08))
se1 <- sqrt(pwr1* (1-pwr1) / m) 
#add standard errors 
lines(a, pwr1+se1, lty = 3) 
lines(a, pwr1-se1, lty = 3)

#plot power vs parameters of t distribution
plot(t, pwr2, type = "b", xlab = "parameters of t distribution", ylim = c(0,1),col="blue")
se2 <- sqrt(pwr2* (1-pwr2) / m) 
#add standard errors 
lines(t, pwr2+se2, lty = 3,col="blue") 
lines(t, pwr2-se2, lty = 3,col="blue")

#when put them together
plot(a, pwr1, type = "b", xlab = "parameters of distribution", ylab = "pwr", ylim = c(0,1))
lines(t, pwr2, type = "b", xlab = "t", ylim = c(0,1),col="blue")
abline(h = .05, lty = 3,col="blue")
legend("topright",
       legend =c("beta","t") ,
       lty=3,
       col=c("black",'blue'))

## -----------------------------------------------------------------------------
# generate samples under H1 to estimate power
n <- c(20,100,1000)#sample size
m <- 10000
sigma1 <- 1
sigma2 <- 1.5
power <- power2 <-  numeric(length(n))
count5test <- function(x,y){
  X <- x-mean(x)#centered by sample mean
  Y <- y-mean(y)
  outx <- sum(X>max(Y))+sum(X<min(Y))
  outy <- sum(Y>max(X))+sum(Y<min(X))
  return(as.integer(max(c(outx,outy))>5))
}
for(i in 1:length(n)){
  power[i] <- mean(replicate(m,expr = {
  x <- rnorm(n[i],0,sigma1)
  y <- rnorm(n[i],0,sigma2)
  count5test(x,y)
  }))
  pvalues <- replicate(m,expr={
    x <- rnorm(n[i],0,sigma1)
    y <- rnorm(n[i],0,sigma2)
    Ftest <- var.test(x, y, ratio = 1,
                      alternative = c("two.sided", "less", "greater"),
                      conf.level = 0.945, ...)
    Ftest$p.value})
  power2[i] <- mean(pvalues<=0.055)
}
power

power2



## -----------------------------------------------------------------------------
library(MASS)
set.seed(123)
alpha<-0.05
d<-2
n <-c(10,20,30,50,100,500)  #sample sizes 
cv <- qchisq(1-alpha,d*(d+1)*(d+2)/6) #crit. values for each n

sk <-function(x) { 
  n<-nrow(x)
  for (i in 1:d) {
    x[,i]<-x[,i]-mean(x[,i])
  }
  s<-solve(cov(x))
  b<-mean((x%*%s%*%t(x))*(x%*%s%*%t(x))*(x%*%s%*%t(x)))
  return(b*n/6)
}
#n is a vector of sample sizes 
#we are doing length(n) different simulations
p.reject <- numeric(length(n)) #to store sim. results 
m <- 1000
#num. repl. each sim.
mu<-rep(0,d)
sigma<-diag(rep(1,d))
for (i in 1:length(n)) { 
  sktests <- numeric(m)
  for (j in 1:m) { 
    x <- mvrnorm(n[i],mu,sigma) 
    sktests[j] <- as.integer(sk(x) >= cv) 
    }
  p.reject[i] <- mean(sktests)
}
p.reject


## -----------------------------------------------------------------------------
library(MASS)
set.seed(0)
alpha<-0.05
d<-2
n <-20 #sample sizes 
cv <- qchisq(1-alpha,d*(d+1)*(d+2)/6) #crit. values for each n
sk <-function(x) { 
  n<-nrow(x)
  for (i in 1:d) {
    x[,i]<-x[,i]-mean(x[,i])
  }
  s<-solve(cov(x))
  b<-mean((x%*%s%*%t(x))*(x%*%s%*%t(x))*(x%*%s%*%t(x)))
  return(b*n/6)
}

m <- 1000
epsilon <- c(seq(0, .15, .05), seq(.15, 0.9, .15)) 
N <- length(epsilon) 
pwr <- numeric(N) #critical value for the skewness test 

for (j in 1:N) {
  e <- epsilon[j] 
  sktests <- numeric(m) 
  for (i in 1:m) {
    sig <- sample(c(1,10), replace = TRUE, size = n, prob = c(1-e, e))
    x <- mvrnorm(n,rep(0,d),diag(rep(sig[1],d)))
    for (k in 2:n) {
      sigma<-diag(rep(sig[k],d))
      x <- rbind(x,mvrnorm(n,rep(0,d),sigma))
    }
    sktests[i] <- as.integer(sk(x) >= cv) 
  }
  pwr[j] <- mean(sktests) 
}
#plot power vs epsilon 
plot(epsilon, pwr, type = "b", xlab = bquote(epsilon), ylim = c(0,1))
se <- sqrt(pwr * (1-pwr) / m) #add standard errors 
lines(epsilon, pwr+se, lty = 3) 
lines(epsilon, pwr-se, lty = 3)

## ----warning=TRUE-------------------------------------------------------------
data(law, package = "bootstrap")
n <- nrow(law)
y <- law$LSAT
z <- law$GPA
theta.hat <- cor(y, z)
#compute the jackknife replicates, leave-one-out estimates
theta.jack <- numeric(n)
for (i in 1:n)
theta.jack[i] <- cor(y[-i],z[-i])
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
print(bias) #jackknife estimate of bias
se <- sqrt((n-1) *
mean((theta.jack - mean(theta.jack))^2))
print(se)


## ----warning=TRUE-------------------------------------------------------------
library(boot)
data('aircondit')
boot.obj <- boot(aircondit,R=2000,
                 statistic = function(x,i){mean(x[i,1])})
print(boot.ci(boot.obj,type = c("basic","norm","perc")))

boot.BCa <- function(x,th0,th,stat,conf=.95){
  x <- as.matrix(x)
  n <- nrow(x)
  N <- 1:n
  alpha <- (1+c(-conf,conf))/2
  zalpha <- qnorm(alpha)
  z0 <- qnorm(sum(th<th0)/length(th))
  th.jack <- numeric(n)
  for (i in 1:n) {
    J <- N[1:(n-1)]
    th.jack[i] <- stat(x[-i, ],J)
  }
  L <- mean(th.jack)-th.jack
  a <- sum(L^3)/(6*sum(L^2)^1.5)
  adj.alpha <- pnorm(z0+(z0+zalpha)/1-a*(z0+zalpha))
  limits <- quantile(th,adj.alpha,type = 6)
  return(list("est"=th0,"BCa"=limits))
}
n <- nrow(aircondit)
B <- 2000
x <- aircondit$hours
theta.b <- numeric(B)
theta.hat <- mean(x)
for (b in 1:B) {
  i <- sample(1:n,size = n,replace = TRUE)
  x <- aircondit$hours[i]
  theta.b[b] <- mean(x)
}
stat <- function(dat,index){
  mean(dat[index])
}
boot.BCa(x,th0 = theta.hat,th=theta.b,stat=stat)

library(boot)
data('aircondit')
boot.obj <- boot(aircondit,R=2000,
                 statistic = function(x,i){mean(x[i,1])})
print(boot.ci(boot.obj,type = c("basic","norm","perc","bca")))

## ----warning=TRUE-------------------------------------------------------------
set.seed(0)
library(boot)
library(bootstrap)
data(scor)
lambda_hat <- eigen(cov(scor))$values
theta_hat <- lambda_hat[1] / sum(lambda_hat)

n <- nrow(scor) # number of rows (data size)
# Jackknife
theta_j <- rep(0, n)
for (i in 1:n) {
x <- scor [-i,]
lambda <- eigen(cov(x))$values
theta_j[i] <- lambda[1] / sum(lambda)
# the i-th entry of theta_j is the i-th "leave-one-out" estimation of theta
}
bias_jack <- (n - 1) * (mean(theta_j) - theta_hat)
# the estimated bias of theta_hat, using jackknife
se_jack <- (n - 1) * sqrt(var(theta_j) / n)
# the estimated se of theta_hat, using jackknife
# print the answers

bias_jack
se_jack

## ----message=FALSE, warning=TRUE----------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic) 
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n)

for (k in 1:(n-1)) {#leave the first pair of data
  mag <- magnetic[-k] 
  che <- chemical[-k]
  for (i in k:(n-1)) {#leave the second pair of data
    y <- mag[-i] 
    x <- che[-i]
    z <- c(chemical[k],chemical[i+1])
    J1 <- lm(y ~ x) 
    yhat1 <- J1$coef[1] + J1$coef[2] * z
    e1[k,i+1] <- (magnetic[k] - yhat1[1])^2 + (magnetic[i+1] - yhat1[2])^2
    
    J2 <- lm(y~x+I(x^2)) 
    yhat2 <- J2$coef[1] + J2$coef[2] * z + J2$coef[3] * z * z
    e2[k,i+1] <- (magnetic[k] - yhat2[1])^2+(magnetic[i+1] - yhat2[2])^2
    
    J3 <- lm(log(y) ~ x) 
    logyhat3 <- J3$coef[1] + J3$coef[2] * z
    yhat3 <- exp(logyhat3) 
    e3[k,i+1] <- (magnetic[k] - yhat3[1])^2+(magnetic[i+1] - yhat3[2])^2
    
    J4 <- lm(log(y) ~ log(x)) 
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(z) 
    yhat4 <- exp(logyhat4) 
    e4[k,i+1] <- (magnetic[k] - yhat4[1])^2+(magnetic[i+1] - yhat4[2])^2     
  }
}
c(mean(e1),mean(e2),mean(e3),mean(e4))


## -----------------------------------------------------------------------------
set.seed(12345)
maxout <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}

alpha <- 0.05    #significance level
n1 <- 30;n2 <-50 #two different sample size
mu1 <- mu2 <- 0;sigma1 <- sigma2 <- 1
m <- 1000   # the times of Monte Carlo experiments

p_value <- replicate(m,expr={
  x1 <- rnorm(n1,mu1,sigma1)
  x2 <- rnorm(n2,mu2,sigma2)
  ts <- numeric(199+1)
  ts[1] <- maxout(x1,x2)
  for(i in 1:199){
    ind <- sample(1:(n1+n2),size = n1,replace = FALSE)
    x.perm <- c(x1,x2)[ind]; y.perm <- c(x1,x2)[-ind]
    ts[i+1] <- maxout(x.perm,y.perm)
  }
  mean(ts >= ts[1])
})
#estimate the type1 error
print(mean(p_value < alpha))


## -----------------------------------------------------------------------------
library(RANN) # implementing a fast algorithm
# for locating nearest neighbors
# (alternative R package: "yaImpute")
library(boot)
library(energy)
library(Ball)

x <- rnorm(n1,mu1,sigma1)
y <- rnorm(n2,mu2,sigma2)
z <- c(x, y) # pooled samplec8

#R code implementing NN test
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}


# Power comparison 
m <- 100; k<-3; p<-2; mu <- 0.8; set.seed(1245)
n1 <- n2 <- 10; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,3),ncol=p);
y <- matrix(rnorm(n2*p,0,2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*123)$p.value
}
alpha <- 0.3;
pow1 <- colMeans(p.values<alpha)
pow1


## -----------------------------------------------------------------------------
m <- 100; k<-3; p<-2; mu <- 0.5; set.seed(12345)
n1 <- n2 <- 10; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,1),ncol=p);
y <- matrix(rnorm(n2*p,1,2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.05;
pow2 <- colMeans(p.values<alpha)
pow2


## -----------------------------------------------------------------------------
m <- 1e2; k<-3; p<-2; mu <- 0.3
n1 <- n2 <- 10; R<-999; n <- n1+n2; N = c(n1,n2)
p.values <- matrix(NA,m,3)

for(i in 1:m){
 x <- matrix(rt(n1*p,1),ncol = p);
  y <- cbind(rnorm(n2),rnorm(n2,mean = 1))
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
  p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow3 <- colMeans(p.values<alpha)
pow3


## -----------------------------------------------------------------------------
m <- 100; k<-3; p<-2; mu <- 0.3; set.seed(345)
n1 <- 30; n2 <- 20; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,3.5),ncol=p);
y <- matrix(rnorm(n2*p,0,2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*1)$p.value
}
alpha <- 0.1;
pow4 <- colMeans(p.values<alpha)
pow4


## ----warning=FALSE, eval=FALSE------------------------------------------------
#  # density of the standard Laplace distribution
#  r<-function(x)0.5*exp(-abs(x))
#  rw.Metropolis <- function(sigma, x0, N) {
#    x <- numeric(N)
#    x[1] <- x0
#    u <- runif(N)
#    k<-0
#    for (i in 2:N) {
#      y <- rnorm(1, x[i-1], sigma)
#      if (u[i] <= (r(y) / r(x[i-1]))) x[i] <- y else {
#        x[i] <- x[i-1]
#        k<-k+1
#      }
#    }
#  return(list(x=x, k=k))
#  }
#  N <- 5000
#  sigma <- c(.05, .5, 2, 16)
#  x0 <- 25
#  rw1 <- rw.Metropolis(sigma[1], x0, N)
#  rw2 <- rw.Metropolis(sigma[2], x0, N)
#  rw3 <- rw.Metropolis(sigma[3], x0, N)
#  rw4 <- rw.Metropolis(sigma[4], x0, N)
#  
#  #the acceptance rates of each chain
#  print(1-c(rw1$k, rw2$k, rw3$k, rw4$k)/2000)
#  
#  plot(1:N,rw1$x,type = "l")
#  plot(1:N,rw2$x,type = "l")
#  plot(1:N,rw3$x,type = "l")
#  plot(1:N,rw4$x,type = "l")

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  r<-function(x)0.5*exp(-abs(x))
#  Gelman.Rubin <- function(psi) {
#    # psi[i,j] is the statistic psi(X[i,1:j])
#    # for chain in i-th row of X
#    psi <- as.matrix(psi)
#    n <- ncol(psi)
#    k <- nrow(psi)
#  
#    psi.means <- rowMeans(psi)
#    B <- n * var(psi.means)
#    psi.w <- apply(psi, 1, "var")
#    W <- mean(psi.w)
#    v.hat <- W*(n-1)/n + (B/n)
#    r.hat <- v.hat / W
#    return(r.hat)
#  }
#  
#  normal.chain <- function(sigma, N, X1) {
#    #generates a Metropolis chain for standard Laplace distribution
#    #with Normal(X[t], sigma) proposal distribution
#    #and starting value X1
#    x <- rep(0, N)
#    x[1] <- X1
#    u <- runif(N)
#    for (i in 2:N) {
#      xt <- x[i-1]
#      y <- rnorm(1, xt, sigma)
#      if (u[i] <= (r(y) / r(x[i-1]))) x[i] <- y else
#        x[i] <- x[i-1]
#    }
#    return(x)
#    }
#  sigma <- 2 #parameter of proposal distribution
#  k <- 4
#  n <- 15000
#  b <- 1000
#  #number of chains to generate
#  #length of chains
#  #burn-in length
#  #choose overdispersed initial values
#  x0 <- c(-10, -5, 5, 10)
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k) X[i, ] <- normal.chain(sigma, n, x0[i])
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi)) psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  #plot psi for the four chains
#  for (i in 1:k) plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
#  #plot the sequence of R-hat statistics
#  rhat <- rep(0, n)
#  for (j in (b+1):n) rhat[j] <- Gelman.Rubin(psi[,1:j])
#  plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#  abline(h=1.1, lty=2)
#  

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  c_k <- function(k,a){
#    sqrt(a^2*k/(k+1-a^2))
#  }
#  res <- function(k,a){
#    pt(c_k(k-1,a),df = k-1) - pt(c_k(k,a),df = k)
#  }
#  root.curve <- sapply(c(4:25,100,500,1000),function(k){uniroot(res,interval = c(0.0001,sqrt(k)-0.0001),k=k)$root})#因为是开区间，所以下界比0大一点，上界比根号k小一点
#  
#  result <- cbind(c(4:25,100,500,1000),root.curve)
#  knitr::kable(result)
#  

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  set.seed(1227)
#  n_A. <- 444
#  n_B. <- 132
#  n_AB <- 63
#  n_OO <- 361
#  p0 <- runif(1,0,1)
#  q0 <- runif(1,0,1-p0)
#  #
#  likelihood_e <- function(prob,p0,q0){
#    r0 <- 1-p0-q0
#    p <- prob[1]; q <- prob[2] ; r <- 1-p-q
#    - n_A. * (2*log(p)*(p0^2/(p0^2+2*p0*r0)) + log(2*p*r)*(2*p0*r0/(p0^2+2*p0*r0))) -
#      n_B. * (2*log(q)*(q0^2/(q0^2+2*q0*r0)) + log(2*q*r)*(2*q0*r0/(q0^2+2*q0*r0))) -
#      n_AB * log(2*p*q) - 2*n_OO * log(r)
#  }
#  #
#  iter <- 0;E1 <- 0;E2 <- 1
#  while(iter < 200 && abs(E1-E2)> 1e-6){
#    output <- optim(par = c(0.1,0.1),likelihood_e,p0 = p0,q0 = q0)
#    E1 <- E2;E2 <- output$value
#    p0 <- output$par[1]
#    q0 <- output$par[2]
#    iter <- iter + 1
#  }
#  estimate <- data.frame(p0,q0,iter)
#  colnames(estimate) <- c("p","q","iteration times")
#  knitr::kable(estimate)

## ----warning=FALSE------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  formulas <- list(
#    mpg ~ disp,
#    mpg ~ I(1 / disp),
#    mpg ~ disp + wt,
#    mpg ~ I(1 / disp) + wt
#  )
#  
#  #Use for loops
#  models.loop <- list()
#  for(i in 1:length(formulas)){
#    models.loop <- c(models.loop,list(lm(formulas[[i]],data = mtcars)))
#  }
#  
#  #use lapply
#  models <- lapply(formulas,lm,data = mtcars)
#  models

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  trials <- replicate(
#    100,
#    t.test(rpois(10, 10), rpois(7, 10)),
#    simplify = FALSE
#  )

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  trials <- replicate(
#    100,
#    t.test(rpois(10, 10), rpois(7, 10)),
#    simplify = FALSE
#  )
#  round(sapply(1:100,function(i){trials[[i]]$p.value}),3)
#  
#  #Use the[[ to get rid of anonymous function.
#  round(sapply(trials,"[[","p.value"),3)

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  testlist <- list(iris, mtcars, cars)
#  lmapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
#  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
#  if(simplify == TRUE){return(simplify2array(out))}
#  out
#  }
#  lmapply(testlist, mean, numeric(1))

## ----echo=TRUE, warning=FALSE, eval=FALSE-------------------------------------
#  library(Rcpp)
#  set.seed(1227)
#  
#  lap_f = function(x) exp(-abs(x))
#  
#  rw.Metropolis = function(sigma, x0, N){
#    x = numeric(N)
#    x[1] = x0
#    u = runif(N)
#    k = 0
#    for (i in 2:N) {
#      y = rnorm(1, x[i-1], sigma)
#      if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y
#      else {
#        x[i] = x[i-1]
#        k = k+1
#      }
#    }
#    return(list(x = x, k = k))
#  }
#  
#  cppFunction('List rw_Metropolis_cpp(double sigma, double x0, int N) {
#    List out(2);
#    NumericVector x(N);
#    x[0] = x0;
#    DoubleVector u = runif(N);
#    int k=0;
#    for(int i=1;i < N; i++) {
#      double y = as<double>(rnorm(1, x[i-1], sigma));
#      if (u[i] <= exp(abs(x[i-1])-abs(y))) {
#        x[i] = y;
#      }
#      else{
#        x[i] = x[i-1];
#        k = k + 1;
#      }
#    }
#    out[0] = x;
#    out[1] = k;
#    return (out);
#  }')
#  
#  N = 2000
#  sigma = c(.05, .5, 2, 16)
#  x0 = 25
#  par(mfrow = c(2,2))
#  
#  rej = numeric()
#  rej_c = numeric()
#  #chains
#  for (i in 1:length(sigma)) {
#    par(mfrow=c(1,2))
#    rw = rw.Metropolis(sigma[i],x0,N)$x
#    rw_c = rw_Metropolis_cpp(sigma[i],x0,N)[[1]]
#    plot(rw, type="l",
#         xlab=bquote(sigma == .(round(sigma[i],3))),
#         ylab="from R", ylim=range(rw))
#    plot(rw_c, type="l",
#         xlab=bquote(sigma == .(round(sigma[i],3))),
#         ylab="from Rcpp", ylim=range(rw_c))
#    rej[i] = rw.Metropolis(sigma[i],x0,N)$k
#    rej_c[i] = rw_Metropolis_cpp(sigma[i],x0,N)[[2]]
#  }
#  
#  #accept rate
#  acc = round((N-rej)/N,4)
#  acc_c = round((N-rej_c)/N,4)
#  res = rbind(acc, acc_c)
#  rownames(res) = c("Accept rates from R","Accept rates from Rcpp")
#  colnames(res) = paste("sigma",sigma)
#  knitr::kable(res)
#  
#  
#  #qqplot
#  for (i in 1:length(sigma)) {
#    rw = rw.Metropolis(sigma[i],x0,N)$x
#    rw_c = rw_Metropolis_cpp(sigma[i],x0,N)[[1]]
#    qqplot(rw, rw_c, xlab = "from R", ylab = "from Rcpp", main = bquote(sigma == .(sigma[i])))
#    f <- function(x) x
#    curve(f, col = 'red',add = TRUE)
#  }
#  
#  #Campare the computation time of the two functions with the function “microbenchmark”
#  library(microbenchmark)
#  for (i in 1:length(sigma)) {
#    ts <- microbenchmark(rw <- rw.Metropolis(sigma[i],x0,N),
#                         rw_c <- rw_Metropolis_cpp(sigma[i],x0,N))
#    time <- data.frame(summary(ts)[,c(3,5,6)])
#    rownames(time) <- c(paste("sigma",sigma[i],"from R"),paste("sigma",sigma[i],"from Rcpp"))
#    print(knitr::kable(time))
#  }

