library('survival')
library('msm')

###################### *** Data generation *** ######################
set.seed(2017)

# Number of patients
n=60

# Kendall's Tau
tau=0.3

# True ratio of the means
et = 1.33

# Shape of the Weibull baseline hazards
Weib_shape = 2

g = rgamma(n, shape=.5*(1/tau-1), rate=.5*(1/tau-1))
data <- data.frame(
  id = 1:n,
  TTP1 = 12 * rweibull(n, Weib_shape, 1/(g*Weib_shape)),
  TTP2 = 12 * rweibull(n, Weib_shape, et/(g*Weib_shape)))

C = runif(n, min=21, max=99)
data$status <- 1 * (data$TTP2 <= C)
data$TTP2[data$status == 0] <- C[data$status == 0]

data$ratio <- data$TTP2/data$TTP1
########################################################################


#################### *** Parametric estimation *** #####################
reg <- survreg(Surv(ratio, status)~1, data = data, dist="loglogistic")
shape <- exp(-reg$coef)
scale <- 1/(reg$scale)
estmean <- c(reg$coef,reg$scale)
estim <- 1/(1+(shape)^(scale))
sd <- deltamethod(~1/(1+(exp(-x1))^(1/x2)), estmean, reg$var) 
mean <- (pi/(shape*scale*sin(pi/scale)))

RES <- rbind('TRUE VALUE' = c(ESTIM=1/(1+(et)^(-Weib_shape)), SD=NA),
             'PARAMETRIC ESTIMATION'=c(estim, sd))
########################################################################
 


################## *** Non parametric estimation *** ###################
rankmax <- Vectorize(function(a, vec)  rev(which(vec <= a))[1], "a")
rankmin <- Vectorize(function(a, vec)  which(vec >= a)[1], "a")

data$TTP1m <- 1*data$TTP1
dataMOD <- data[, c('TTP1m', 'TTP2')]
dataMOD$TTP2[data$status == 0] <- Inf

lsort = sort(unlist(data[, c('TTP1m', 'TTP2')]))
rsort = sort(unlist(dataMOD))

mid1 = (rankmin(dataMOD$TTP1m, rsort) + rankmax(dataMOD$TTP1m, lsort)) /2
mid2 = (rankmin(data$TTP2, rsort) + rankmax(dataMOD$TTP2, lsort)) /2

prop = as.numeric((mid2 - mid1) >= 0)

estim <- mean(prop)
sd <- sqrt((estim*(1-estim))/n)

RES <- rbind(RES, 'NON PARAMETRIC ESTIMATION'=c(estim, sd))
########################################################################


# RESULTS
RES <- cbind(
  RES, LCI = RES[,1] + qnorm(.025) * RES[,2],
  UCI = RES[,1] + qnorm(.975) * RES[,2])
print(RES, na.print='', digits=2)