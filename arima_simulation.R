## generate model
arma21 <- arima.sim(model = list(ar=c(0.8, 0.1), ma=0.8), n = 200, sd = 0.2)
acf(arma21)
pacf(arma21)

### Invertible or stationary ?
(fit.arma21 <- arima(arma21, order =c(2,0,1), include.mean = FALSE))
names(fit.arma21)
fit.arma21$coef
polyroot(rev(c(1, 0.8)))
polyroot(rev(c(1, -0.8, -0.1))) ## constraint

### simulation
nb.sim <- 10
nb.obs <- 200
sim.mat <- mat.or.vec(nb.obs, nb.sim)
for(i in 1:nb.sim){
  sim.mat[,i] <- arima.sim(model = list(ar=c(0.8, 0.1), ma=0.8), n = nb.obs, sd = 0.2)
}

### plot the realization of the process
plot.sim <- function(simulations){
  n1 <- dim(simulations)[1]
  n2 <- dim(simulations)[2]
  plot(simulations[,1], type = "l", main = "Process simulation")
  for(i in 2:n2){
    matlines(1:n1, simulations[,i], lty=1, col = 1)
  }
}
plot.sim(sim.mat)


### estimation and plots of the acf
nb.lag <- 50
par(mfrow=c(1,1))
acf(sim.mat[,1], lag.max = nb.lag)
for (i in 2:nb.sim){
  lines(1:(nb.lag+1), acf(sim.mat[,i], lag.max = nb.lag, plot = FALSE)$acf, lty = 1, col = i, type = "h")
}

### estimation and plots of the pacf
par(mfrow=c(1,1))
pacf(sim.mat[,1], lag.max = nb.lag)
for (i in 2:nb.sim){
  lines(pacf(sim.mat[,i], lag.max = nb.lag, plot = FALSE)$acf, lty = 1, col = i, type = "h")
}


### variance of each realization
data.frame("realization rank"= 1:nb.sim, "variance"= apply(sim.mat, 2, var))



## seasonality
log.c02 <- c(6.153, 6.308, 6.327, 6.344, 6.266, 6.121, 6.096, 6.081,6.178,
             6.411, 6.442, 6.505, 6.195, 6.085, 6.086, 6.079)
sigma <- 0.065 ## with six days data
mu <- 6.24
plot(log.c02, type = "l")

### seasonal model because the observations are made every 2nd of each month
plot(log.c02[9:16], type = "l")
lines(1:8, log.c02[1:8],lty = 1, col = "red", type = "l")

### predictions 

### trying to simulate 
sim.season <- arima((log.c02-mu),order=c(1,0,0), seasonal = list(order = c(1, 0, 0), period = 8), method = "CSS"
                    ,include.mean = FALSE)
sim.season

### prdictions at t+1 and t+2 where t=16
# let's write x = log.cO2 - mu
season.coef <- c(0.547, 0.86)
x <- log.c02 - mu
x17 <- season.coef[1]*x[9] + season.coef[2]*season.coef[1]*x[8] - season.coef[1]*x[16]
x18 <- season.coef[1]*x[10] + season.coef[2]*season.coef[1]*x[9] - season.coef[1]*x17
log.c02.17 <- x17 + mu
log.c02.18 <- x18 + mu


## predictions interval
log.c02.17.lower <- (x17 - 1.96*sigma) + mu
log.c02.17.upper <- (x17 + 1.96*sigma) + mu
log.c02.18.lower <- x18 - 1.96*sqrt(season.coef[1])*sigma*sqrt(1+(season.coef[1])^2) + mu
log.c02.18.upper <- x.18 + 1.96*sqrt(season.coef[1])*sigma*sqrt(1+(season.coef[1])^2) + mu
### for log.co2.17
par(mfrow=c(1,1))
plot(log.c02, type = "l", xlim = c(1,18))
segments(x0 = 16, y0 = log.c02[16], x1 = 17, y1 = log.c02.17, col = "red")
points(17, log.c02.17, lty = 1, col = "blue", lwd = 3, type = "b")
lines(x = c(17, 17), y = c(log.c02.17.upper,log.c02.17.lower), col='red',lwd=3, lty = "99")
### for log.co2.18
par(mfrow=c(1,1))
plot(log.c02, type = "l", xlim = c(1,18))
segments(x0 = 16, y0 = log.c02[16], x1 = 17, y1 = log.c02.17, col = "red")
points(17, log.c02.17, lty = 1, col = "blue", lwd = 3, type = "b")
points(18, log.c02.18, lty = 1, col = "blue", lwd = 3, type = "b")
segments(x0 = 17, y0 = log.c02.17, x1 = 18, y1 = log.c02.18, col = "red")
lines(x = c(18, 18), y = c(log.c02.18.upper,log.c02.18.lower), col='red',lwd=3, lty = "99")

## simulating seasonal processes

## simuation 1
sim.1 <- arima.sim(model = list(ar=-0.8, order=c(1,1,0)), n = 500)
arima(sim.1)
plot(sim.1)
acf(sim.1)
pacf(sim.1)

## simulation 2
# sim.2 <- arima.sim(list(order = c(12,0,0), ar = c(rep(0,11),0.85)), n = 500)
sim.2 <- arima.sim(model=list(ar=c(rep(0,11),0.85)), n=500)
plot(sim.2)
acf(sim.2, lag.max = 50)
pacf(sim.2, lag.max = 50)

## simulation 3
sim.3 <- arima.sim(model=list(ar=0.7, ma = c(rep(0,11),0.8)), n=500)
plot(sim.3)
acf(sim.3, lag.max = 50)
pacf(sim.3, lag.max = 50)


## simulation 4
library(sarima)
sim.4 <- arima.sim(n = 500, list( ar=c(0.8,rep(0,10), -0.7, 0.56)))
plot(sim.4)
acf(sim.4, lag.max = 50)
pacf(sim.4, lag.max = 50)

## simulation 5

## 100 simulations of the four arma processes
nb.sim <- 100
nb.obs <- 300


### functions of simulations with theta, phi and sigma as parameters
sim.arma <- function(phi, theta, sigma, nb.sim, nb.obs){
  sim.matrix <- mat.or.vec(nb.obs, nb.sim)
  for(i in 1:nb.sim){
    sim.matrix[,i] <- arima.sim(list(ar = phi, ma = theta), n = nb.obs, sd = sigma)
  }
  return(sim.matrix)
}

coef.estimation <- function(sim.matrix, method = "ML"){
  p <- dim(sim.matrix)[2]
  coef.list <- numeric(p)
  for(i in 1:p){
    coef.list[i] <- (arima(sim.matrix[,i], order=c(1,0,1), method = method, include.mean = FALSE)$coef)[1]
  }
  return(coef.list)
}

sigma.phi.estimation <- function(sim.matrix, method = "ML"){
  p <- dim(sim.matrix)[2]
  var.list <- numeric(p)
  for(i in 1:p){
    cov.var.mat <- arima(sim.matrix[,i], order=c(1,0,1), method = method, include.mean = FALSE)$var.coef
    var.list[i] <- cov.var.mat[1][1]
  }
  return(var.list)
}

process.description <- function(sim.matrix, method = "ML"){
  oldpar <- par(mfrow=c(3,1), mgp=c(2,0.7,0), mar=c(3,3,1.5,1))
  on.exit(par(oldpar))
  plot.sim(sim.matrix)
  coef.list <- coef.estimation(sim.matrix, method)
  var.list <- sigma.phi.estimation(sim.matrix, method)
  hist(coef.list, main = "histogram plot of the estimates phi", col = 2)
  plot(coef.list, var.list, main = "evolution of the variance per estimate of phi", type = "h")
}

### Process 1
sim.mat.1 <- sim.arma(phi = 0.9, theta = -0.4, sigma = 0.1, nb.sim = nb.sim, nb.obs = nb.obs)
plot.sim(sim.mat.1)
process.description(sim.matrix = sim.mat.1)
plot(sigma.phi.estimation(sim.matrix = sim.mat.1), type = "l", main = "estimated variance of phi hat")
abline(b=0, a=mean(sigma.phi.estimation(sim.matrix = sim.mat.1)), col ="red", lwd = "2")
#### get quantiles

### Process 2
sim.mat.2 <- sim.arma(phi = 0.9, theta = -0.4, sigma = 5, nb.sim = nb.sim, nb.obs = nb.obs)
process.description(sim.matrix = sim.mat.2)
plot(sigma.phi.estimation(sim.matrix = sim.mat.2), type = "l", main = "estimated variance of phi hat")
abline(b=0, a=mean(sigma.phi.estimation(sim.matrix = sim.mat.2)), col ="red", lwd = "2")

### Process 3
sim.mat.3 <- sim.arma(phi = 0.995, theta = -0.4, sigma = 0.1, nb.sim = nb.sim, nb.obs = nb.obs)
plot.sim(simulations = sim.mat.3)
process.description(sim.matrix = sim.mat.3, method = "CSS")
###  minimizing the sum of squared residuals (CSS) for estimation of the parameters cuz non stationarity
plot(sigma.phi.estimation(sim.matrix = sim.mat.3, method = "CSS"), type = "l", main = "variance of phi hat")
abline(b=0, a=mean(sigma.phi.estimation(sim.matrix = sim.mat.3, method = "CSS")), col ="red", lwd = "2")

### Process 4
sim.mat.4 <- sim.arma(phi = 0.995, theta = -0.4, sigma = 5, nb.sim = nb.sim, nb.obs = nb.obs)
process.description(sim.matrix = sim.mat.4, method = "CSS")
###  minimizing the sum of squared residuals (CSS) for estimation of the parameters cuz non stationarity
plot(sigma.phi.estimation(sim.matrix = sim.mat.4, method = "CSS"), type = "l", main = "variance of phi hat")
abline(b=0, a=mean(sigma.phi.estimation(sim.matrix = sim.mat.4, method = "CSS")), col ="red", lwd = "2")

###Last process
sim.mat.5 <- sim.arma(phi = 0.9, theta = -0.4, sigma = 1, nb.sim = 100, nb.obs = 100)
process.description(sim.matrix = sim.mat.5)
plot(sigma.phi.estimation(sim.matrix = sim.mat.5), type = "l", main = "estimated variance of phi hat")
abline(b=0, a=mean(sigma.phi.estimation(sim.matrix = sim.mat.5)), col ="red", lwd = "2")
