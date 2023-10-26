# This is a quick simulation demonstrating that 
# unmeasured confounding between the instrument 
# and treatment can bias IV (for the ATE) when we relax the
# constant-effects assumption. Here we make the
# treatment effect and first-stage effect correlated
# but leave all other effects constant. Every variable
# is continuous and all effects are linear. 

rm(list = ls())

library(AER)

n.obs <- 5000
n.sims <- 1000

sims <- data.frame(tau_est = rep(NA, n.sims),
                   f_test = rep(NA, n.sims))

zb1 <-.2
zb2 <- .2

db1 <- .4
db2 <- .4

for (i in 1:n.sims) {

  uzd <- rnorm(n.obs)
  udy <- rnorm(n.obs)
  
  tau <- rnorm(n.obs, .5, sd = .5)

  pi <- tau*.4 + rnorm(n.obs, sd = .5)

  z <- rnorm(n.obs, uzd*zb1)
  d <- rnorm(n.obs, z*pi + uzd*zb2 + udy*db1)
  y <- rnorm(n.obs, d*tau + udy*db2)
  dat <- data.frame(z,d,y)

  for_f <- summary(lm(d ~ z, dat))
  
  sims$f_test[i] <- for_f$fstatistic[1]
  
  for_est <- ivreg(y ~ d | z, data = dat)

  sims$tau_est[i] <- for_est$coefficients[2]
  
}

mean(sims$tau_est)

mean(sims$f_test)
