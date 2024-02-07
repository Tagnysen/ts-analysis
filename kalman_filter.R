### Read the data
rm(list=ls())
data.bulding = read.csv("data_building.csv")
View(data.bulding)

### Testing data and training data
nb.obs <- dim.data.frame(data.bulding)[1]
train.data <- data.bulding[1:(nb.obs-3), ]
test.data <- data.bulding[(nb.obs-2):nb.obs,]
n1 <- dim.data.frame(train.data)[1]
n2 <- dim.data.frame(test.data)[1]

### Plotting
plot(1:nb.obs, data.bulding$yTi, type = "l", col = 2, main = "Indoor Temperature", xlab = "time")
plot(1:nb.obs, data.bulding$Ta, type = "l", col = 3, main = "Ambient Temperature", xlab = "time")
plot(1:nb.obs, data.bulding$Ph, type = "l", col = 4, main = "Heating Input", xlab = "time")
plot(1:nb.obs, data.bulding$Gv, type = "l", col = 6, main = "Solar radiation", xlab = "time")

library(MASS)
A <- matrix(data = c(0.755, 0.1, 0.24, 0.9), ncol = 2, nrow = 2)
B <- matrix(data = c(0.005, 0, 0.127, 0, 0.335, 0), ncol = 3, nrow = 2)
C <- matrix(c(1,0),nrow=1)
Sigma1 <- matrix(data = c(0.5, 0, 0, 0.5), ncol = 2, nrow = 2)
Sigma2 <- matrix(0.5)

### Estimate using the Kalman filter
### model
X <- matrix(nrow=2,ncol=n1)
X[1, ] <- train.data$yTi
X[2, ] <- rep(mean(train.data$yTi), n1)   ### I assume that the first value of Tm = mean(Ti) at t=0
Y <- numeric(n1)
Y[1] <- C%*%X[,1]+sqrt(Sigma2) %*% rnorm(1)
U.init <- rbind(train.data[1, -c(1,2)])
## Simulation
for (I in 2:n1){
  X[,I] <- A%*%X[,I-1,drop=FALSE] + B%*%t(U.init) + mvrnorm(n = 1, mu = c(0,0), Sigma = Sigma1)
  Y[I] <- C%*%X[,I] + sqrt(Sigma2)%*% rnorm(1)
}

## Loading library
library("FKF")
## Running the Kalman filter with the parameters and initial values used for the simulation
kf <- fkf(a0=c(X[1,1], mean(train.data$yTi)), P0 = diag(0,2), dt = matrix(c(0,0),ncol=1), 
          Tt = A, ct=0, Zt = C, HHt = Sigma1, GGt = Sigma2, yt = matrix(Y,nrow=1))
str(kf)
plot(kf$at[1,], type = "l", col = 2, xlab="Time", ylab="Predicted values", main = "One step prediction of the states")
lines(kf$at[2,], col = 1)
legend(1, 45, legend=c("Ti", "Tm"), col=c("red", "black"), lty=1, cex=0.8)


### Predictions of the last third observations
pred.test <- matrix(nrow = 2, ncol = 3)
X.pred <- rbind(test.data$yTi, rep(0,3))
pred.test[,1] <- A%*%kf$att[,537] + B%*%t(test.data[1, -c(1,2)])
Pt.537 <- 
for (i in 2:3) {
  pred.test[,i] <- A%*%pred.test[,i-1] + B%*%t(test.data[i, -c(1,2)])
}
plot(test.data$yTi, type = "p")
lines(pred.test[1,] + cbind(0,-1.96*sqrt(Pt[1,1,]),1.96*sqrt(Pt[1,1,])),type="l", lty=c(1,2,2), col=3, lwd=2)






