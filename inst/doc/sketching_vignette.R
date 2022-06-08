## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(sketching)
seed <- 220526  
set.seed(seed)  

## -----------------------------------------------------------------------------
Y <- AK$LWKLYWGE
intercept <- AK$CNST
X_end <- AK$EDUC
X_exg <- AK[,3:11]
X <- cbind(X_exg, X_end)
Z_inst <- AK[,12:(ncol(AK)-1)]
Z <- cbind(X_exg, Z_inst)
fullsample <- cbind(Y,intercept,X)
n <- nrow(fullsample)
d <- ncol(X)

## -----------------------------------------------------------------------------
# choice of m (data-oblivious sketch size)
target_size <- 0.05
target_power <- 0.8
S_constant <- (stats::qnorm(1-target_size) + stats::qnorm(target_power))^2
tau_limit <- 10
m_ols <- floor(n*S_constant/tau_limit^2) 
print(m_ols)

## -----------------------------------------------------------------------------
ys <- fullsample[,1]
reg <- as.matrix(fullsample[,-1])
fullmodel <- lm(ys ~ reg - 1)
# use homoskedasticity-only asymptotic variance
ztest <- lmtest::coeftest(fullmodel, df = Inf)
est <- ztest[(d+1),1] 
se <- ztest[(d+1),2]
print(c(est,se))
# use heteroskedasticity-robust asymptotic variance
ztest_hc <- lmtest::coeftest(fullmodel, df = Inf, 
            vcov = sandwich::vcovHC, type = "HC0")
est_hc <- ztest_hc[(d+1),1] 
se_hc <- ztest_hc[(d+1),2]
print(c(est_hc,se_hc))

## -----------------------------------------------------------------------------
subsample <- sketch(fullsample, m_ols, method = "bernoulli")
ys <- subsample[,1]
reg <- subsample[,-1]
submodel <- lm(ys ~ reg - 1) 
# use homoskedasticity-only asymptotic variance
ztest <- lmtest::coeftest(submodel, df = Inf)
est <- ztest[(d+1),1] 
se <- ztest[(d+1),2]
print(c(est,se))
# use heteroskedasticity-robust asymptotic variance
ztest_hc <- lmtest::coeftest(submodel, df = Inf, 
            vcov = sandwich::vcovHC, type = "HC0")
est_hc <- ztest_hc[(d+1),1] 
se_hc <- ztest_hc[(d+1),2]
print(c(est_hc,se_hc))

## -----------------------------------------------------------------------------
subsample <- sketch(fullsample, m_ols, method = "unif")
ys <- subsample[,1]
reg <- subsample[,-1]
submodel <- lm(ys ~ reg - 1) 
# use homoskedasticity-only asymptotic variance
ztest <- lmtest::coeftest(submodel, df = Inf)
est <- ztest[(d+1),1] 
se <- ztest[(d+1),2]
print(c(est,se))
# use heteroskedasticity-robust asymptotic variance
ztest_hc <- lmtest::coeftest(submodel, df = Inf, 
            vcov = sandwich::vcovHC, type = "HC0")
est_hc <- ztest_hc[(d+1),1] 
se_hc <- ztest_hc[(d+1),2]
print(c(est_hc,se_hc))

## -----------------------------------------------------------------------------
subsample <- sketch(fullsample, m_ols, method = "countsketch")
ys <- subsample[,1]
reg <- subsample[,-1]
submodel <- lm(ys ~ reg - 1) 
# use homoskedasticity-only asymptotic variance
ztest <- lmtest::coeftest(submodel, df = Inf)
est <- ztest[(d+1),1] 
se <- ztest[(d+1),2]
print(c(est,se))
# use heteroskedasticity-robust asymptotic variance
ztest_hc <- lmtest::coeftest(submodel, df = Inf, 
            vcov = sandwich::vcovHC, type = "HC0")
est_hc <- ztest_hc[(d+1),1] 
se_hc <- ztest_hc[(d+1),2]
print(c(est_hc,se_hc))

## -----------------------------------------------------------------------------
subsample <- sketch(fullsample, m_ols, method = "srht")
ys <- subsample[,1]
reg <- subsample[,-1]
submodel <- lm(ys ~ reg - 1) 
# use homoskedasticity-only asymptotic variance
ztest <- lmtest::coeftest(submodel, df = Inf)
est <- ztest[(d+1),1] 
se <- ztest[(d+1),2]
print(c(est,se))
# use heteroskedasticity-robust asymptotic variance
ztest_hc <- lmtest::coeftest(submodel, df = Inf, 
            vcov = sandwich::vcovHC, type = "HC0")
est_hc <- ztest_hc[(d+1),1] 
se_hc <- ztest_hc[(d+1),2]
print(c(est_hc,se_hc))

## -----------------------------------------------------------------------------
fullsample <- cbind(Y,intercept,X,intercept,Z)
n <- nrow(fullsample)
p <- ncol(X)
q <- ncol(Z)
# choice of m (data-oblivious sketch size)
target_size <- 0.05
target_power <- 0.8
S_constant <- (qnorm(1-target_size) + qnorm(target_power))^2
tau_limit <- 5
m_2sls <- floor(n*S_constant/tau_limit^2) 
print(m_2sls)

## -----------------------------------------------------------------------------
ys <- fullsample[,1]
reg <- as.matrix(fullsample[,2:(p+2)])
inst <- as.matrix(fullsample[,(p+3):ncol(fullsample)]) 
fullmodel <- ivreg::ivreg(ys ~ reg - 1 | inst - 1) 
# use homoskedasticity-only asymptotic variance
ztest <- lmtest::coeftest(fullmodel, df = Inf)
est <- ztest[(d+1),1] 
se <- ztest[(d+1),2]
print(c(est,se))
# use heteroskedasticity-robust asymptotic variance
ztest_hc <- lmtest::coeftest(fullmodel, df = Inf, 
            vcov = sandwich::vcovHC, type = "HC0")
est_hc <- ztest_hc[(d+1),1] 
se_hc <- ztest_hc[(d+1),2]
print(c(est_hc,se_hc))

## -----------------------------------------------------------------------------
# sketching methods for 2SLS
methods <- c("bernoulli","unif","countsketch","srht")
results_2sls <- array(NA, dim = c(length(methods),3))
for (met in 1:length(methods)){
  method <- methods[met]
    # generate a sketch
    subsample <- sketch(fullsample, m_2sls, method = method)
    ys <- subsample[,1]
    reg <- as.matrix(subsample[,2:(p+2)])
    inst <- as.matrix(subsample[,(p+3):ncol(subsample)]) 
    submodel <- ivreg::ivreg(ys ~ reg - 1 | inst - 1) 
    # use homoskedasticity-only asymptotic variance
    ztest <- lmtest::coeftest(submodel, df = Inf)
    est <- ztest[(d+1),1] 
    se <- ztest[(d+1),2]
    # use heteroskedasticity-robust asymptotic variance
    ztest_hc <- lmtest::coeftest(submodel, df = Inf, 
            vcov = sandwich::vcovHC, type = "HC0")
    est_hc <- ztest_hc[(d+1),1] 
    se_hc <- ztest_hc[(d+1),2]
  results_2sls[met,] <- c(est, se, se_hc)
}
rownames(results_2sls) <- methods
colnames(results_2sls) <- c("est", "non-robust se","robust se")
print(results_2sls)

