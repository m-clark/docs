library(mlogit)
data('Train', package='mlogit')
Tr <- mlogit.data(Train, choice = "choice", shape = "wide",
                  varying = 4:11, alt.levels = c(1,2), sep = "", id='id')
head(Tr)

library(dplyr)
Tr$price= Tr$price/100*2.20371
Tr$time = Tr$time/60

ml.train = mlogit(choice ~ price + time + change + comfort | 0, data=Tr)
summary(ml.train)



# standard 'multinomial logistic' model

ml <- foreign::read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
summary(ml)

# The data set contains variables on 200 students. The outcome variable is prog,
# program type. The predictor variables are social economic status, ses, a
# three-level categorical variable and writing score, write, a continuous
# variable.

ml$prog2 <- relevel(ml$prog, ref = "academic")
m1 = nnet::multinom(prog2 ~ ses + write, data = ml)
summary(m1)
coef(m1)

# compare to multiple binary logistic
ml$vocation = ml$prog2=='vocation'
ml$general = ml$prog2=='general'
glm(vocation ~ ses + write, data = ml[ml$prog2!='general',], family=binomial)
glm(general ~ ses + write, data = ml[ml$prog2!='vocation',], family=binomial)

ml2 = mlogit.data(ml, choice='prog2', id.var='id', shape='wide')
head(ml2, 12)

# from the help on

m2 = mlogit(prog2 ~ 0|ses + write, data = ml2)
summary(m2)
coef(m1)
coef(m2)


multinomregML <- function(par, X, y) {
  levs = levels(y)
  ref = levs[1]

  y0 = y==ref
  y1 = y==levs[2]
  y2 = y==levs[3]


  par = matrix(par, ncol=2)

  # V1 = X %*% par[1:4]
  # V2 = X %*% par[5:8]
  # ll = sum(-log(1 + exp(V1)+exp(V2))) + sum(V1[y1],V2[y2])  # more like mnlogit depiction

  V = X %*% par
  baseProbVec <- 1/(1 + rowSums(exp(V)))

  loglik = sum(log(baseProbVec))  + crossprod(c(V), c(y1, y2))
  loglik
}

# debugonce(multinomregML)
multinomregML(runif(8, -1,1),
              X=model.matrix(prog2 ~ ses + write, data = ml),
              y=ml$prog2)

test = optim(runif(8, -.1, .1), multinomregML, X=model.matrix(prog2 ~ ses + write, data = ml),
             y=ml$prog2, control=list(maxit=1000, reltol=1e-12, ndeps=rep(1e-8, 8),
                                      trace=T, fnscale=-1, type=3),
             method='BFGS')
test
coef(m2)

# debugonce(mlogit)
m2 = mlogit(prog2 ~ 0|ses + write, data = ml2)


library(mnlogit)
m3 = mnlogit(prog2 ~ 0|ses + write, data = ml2)
m3
# debugonce(mnlogit)
m3 = mnlogit(prog2 ~ 1|ses + write, data = ml2)
coef(m3)


cbind(coef(m2)[c(1,3,5,7,2,4,6,8)], coef(m3), test$par) %>% round(5)
cbind(logLik(m2), logLik(m3), test$value)



# Alternative specific and constant variables -----------------------------

# in this example, price is alternative invariant (Z) income is
# individual/alternative specific (X), and catch is alternative specific (Y)

data(Fish)
head(Fish)
fm <- formula(mode ~ price | income | catch)
# debugonce(mnlogit)
fit <- mnlogit(fm, Fish)
summary(fit)

# X dim nrow(Fish)/K x p + 1 (intercept)
# Z, Y nrow(N); Y has alt specific coefs; then for Z ref group dropped so nrow = nrow*(K-1)/K
# for ll everything through previous X the same
# then calc probmat for Y and Z, add to X probmat, and add to base

multinomregML2 <- function(par, X, Y, Z, respVec, choice) {

  N = sum(choice)
  K = length(unique(respVec))
  levs = levels(respVec)

  xpar = matrix(par[1:6], ncol=3)
  ypar = matrix(par[7:10], ncol=4)
  zpar = matrix(par[length(par)], ncol=1)

  # Calc X
  Vx  = X %*% xpar

  # Calc Y (mnlogit finds N x 1 results by going through 1:N, N+1:N*2 etc; then
  # makes 1 vector, then subtracts the first 1:N from whole vector, then makes
  # Nxk-1 matrix with N+1:end values (as 1:N are just zero)); creating the
  # vector and rebuilding the matrix is unnecessary though
  Vy = sapply(1:K, function(alt) Y[respVec == levs[alt],, drop=F] %*% ypar[alt])
  Vy = Vy[,-1] - Vy[,1]

  # Calc Z
  Vz = Z %*% zpar
  Vz = matrix(Vz, ncol=3)

  # all Vs must fit into N x K -1 matrix where N is nobs (i.e. individuals)
  V = Vx + Vy + Vz

  ll0 = crossprod(c(V), choice[-(1:N)])
  baseProbVec <- 1/(1 + rowSums(exp(V)))
  loglik = sum(log(baseProbVec)) + ll0
  loglik

  # note fitted values via
  # fitnonref = exp(V) * matrix(rep(baseProbVec, K-1), ncol = K-1)
  # fitref = 1-rowSums(fitnonref)
  # fits = cbind(fitref, fitnonref)
}


inits = runif(11, -.1, .1)
mdat = mnlogit(fm, Fish)$model  # this data already ordered!

# X has a constant value across alternatives; the coefficients regard selection of alternative relative to reference
X = cbind(1, mdat[mdat$`_Alt_Indx_`=='beach', 'income']); dim(X); head(X)

# Y will use the complete data to start; coefficients will be differences from reference alternative coefficient
Y = as.matrix(mdat[,'catch', drop=F]); dim(Y)

# Z are difference scores from reference group
Z = as.matrix(mdat[mdat$`_Alt_Indx_`!='beach','price', drop=F])
Z = Z-mdat[mdat$`_Alt_Indx_`=='beach','price']; dim(Z);

respVec = mdat$`_Alt_Indx_` # first 10 should be 0 0 1 0 1 0 0 0 1 1 after beach dropped

# debugonce(multinomregML2)
multinomregML2(inits, X, Y, Z, respVec, choice=mdat$mode)
out = optim(par=rep(0,11), multinomregML2, X=X, Y=Y, Z=Z, respVec=respVec, choice=mdat$mode,
            control=list(maxit=1000, reltol=1e-12, ndeps=rep(1e-8, 11),
                         trace=T, fnscale=-1, type=3),
            method='BFGS')
# out
# round(out$par, 3)
round(cbind(out$par, coef(fit)), 3)
cbind(logLik(fit), out$value)
