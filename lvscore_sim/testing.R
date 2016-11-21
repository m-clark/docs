nitems = 6
loadings = .8
effect =.5
sampsize = 1000
factorvar = 1

# lambda = rnorm(nitems-1, mean=loadings, sd=.1)  # unstandardized
# lambda = matrix(c(1, lambda), nrow=1)

# lambda = rnorm(nitems, mean=loadings, sd=.1)    # standardized (used in simulations)
lambda = rep(loadings, nitems)                    # just for this demo

# factors and some noise
f = matrix(rnorm(sampsize, mean=rep(0), sd=factorvar), ncol=1)
e = mvtnorm::rmvnorm(sampsize, sigma=diag(1-c(lambda)^2, nitems))

# observed responses
x = 0 + f%*%lambda + e

y = 2 + effect*scale(f) + rnorm(sampsize)
d = data.frame(x, y); colnames(d) = tolower(colnames(d))


library(lavaan)
population.model <- '
y ~ 2*1 + .5*f1

f1 =~ .8*x1 + 0.8*x2 + .8*x3

x1 ~~ (1-.8^2)*x1
x2 ~~ (1-.8^2)*x2
x3 ~~ (1-.8^2)*x3
'
population.model <- '
y ~ 2*1 + .5*f1

f1 =~ .8*x1 + 0.8*x2 + .8*x3 + .8*x4 + 0.8*x5 + .8*x6

x1 ~~ (1-.8^2)*x1
x2 ~~ (1-.8^2)*x2
x3 ~~ (1-.8^2)*x3
'

# generate data; note, standardized lv is default
myData <- simulateData(population.model, sample.nobs=sampsize)

psych::describe(myData)
psych::describe(d)


test.model <- '
y ~ f1
f1 =~ x1 + x2 + x3
'
test.model <- '
y ~ f1
f1 =~ x1 + x2 + x3 + x4 + x5 + x6
'

lavmod = sem(test.model, data=myData, meanstructure=T, std.lv=T)
summary(lavmod, standardized=T)

mymod = sem(test.model, data=d, meanstructure=T, std.lv=T)
summary(mymod, standardized=T)

cor(lavPredict(mymod), f)
plot(lavPredict(mymod), f)


library(rstan)

datalist = list(N=nrow(x), K=ncol(x), X=x, y=c(y))
bf0 <- stan(file = 'bayesfscore.stan',  data = datalist, iter = 10, verbose = FALSE)

bf <- stan(file = 'bayesfscore.stan', fit=bf0, data = datalist, iter=7000, warmup=2000, thin=20, cores=4)

print(bf, pars=c('b0', 'beta', 'sigma', 'lambda', 'u'))
# shinystan::launch_shinystan(bf)

cor(cbind(get_posterior_mean(bf, 'f')[,5],  lavPredict(mymod), f))
head(cbind(get_posterior_mean(bf, 'f')[,5], lavPredict(mymod), f))
head(cbind(get_posterior_mean(bf, 'lambda')[,5], parTable(mymod)[2:7, 'est']))