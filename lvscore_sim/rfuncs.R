
# Functions for data generation and model building ------------------------

create_data_1factor <- function(loadings=.5, nitems=3,
                                factorvar=1, sampsize=500,
                                effect=.5) {

  # loading matrix
  # lambda = rnorm(nitems-1, mean=loadings, sd=.1)
  # lambda = matrix(c(1, lambda), nrow=1)

  lambda = rnorm(nitems, mean=loadings, sd=.1)
  lambda = ifelse(lambda >= 1, .99, lambda)

  # factors and some noise
  f = matrix(rnorm(sampsize, mean=rep(0), sd=factorvar), ncol=1)
  e = mvtnorm::rmvnorm(sampsize, sigma=diag(1-lambda^2, nitems))

  # observed responses
  x = 0 + scale(f)%*%lambda + e

  y = 2 + effect*scale(f) + rnorm(sampsize)

  res = list(Indicators=x, Target=y)
  res
}


runLM <- function(dataList, type='fscore') {
  if(!type %in% c('fscore', 'sumscore', 'randomitem')) stop("type must be one of 'fscore', 'sumscore', or 'randomitem'")
  if(type == 'sumscore'){
    score = scale(rowSums(dataList$Indicators))[,1]
  }
  if(type == 'randomitem'){
    score = scale(dataList$Indicators[, sample(1:ncol(dataList$Indicators), 1)])[,1]
  }
  if(type == 'fscore'){
    fres = psych::fa(dataList$Indicators, covar=F, fm='ML')
    fres = psych::factor.scores(dataList$Indicators, fres, method='Thurstone')
    score = fres$scores
  }
  lmMod = lm(dataList$Target~score)
  sumLmMod = summary(lmMod)
  out = list(coefficients = coef(lmMod),
             se = sumLmMod$coefficients[,'Std. Error'],
             sigma = sumLmMod$sigma)
  return(out)
}

runLavaan <- function(dataList) {
  nitems = ncol(dataList$Indicators)
  inds = paste('X', 1:nitems, sep='+')
  pastefunc = function(nitems){
    if(nitems==1) return(paste0('X', nitems))
    paste0(paste0('X', nitems, ' + '), pastefunc(nitems-1))
  }

  mod = paste0('f =~ ', pastefunc(nitems))
  mod = paste0(mod, '\n y ~ f')
  modres = lavaan::cfa(mod, data=data.frame(y=dataList$Target, dataList$Indicators), std.lv=T, meanstructure=T)
  dplyr::filter(lavaan::parTable(modres), (lhs=='y' & rhs=='f') | (lhs=='y' & rhs=='y'))
}


summarizeResults <- function(lmres_fscore, lmres_sumscore, lmres_ranitem, lvres, raw=F) {
  coefResult = data.frame(lmfscore = sapply(lmres_fscore, function(x) x$coefficients[2]),
                          lmsumcore = sapply(lmres_sumscore, function(x) x$coefficients[2]),
                          lmitem = sapply(lmres_ranitem, function(x) x$coefficients[2]),
                          lav = sapply(lvres, function(x) x$est[1]))
  seResult = data.frame(lmfscore = sapply(lmres_fscore, function(x) x$se[2]),
                        lmsumcore = sapply(lmres_sumscore, function(x) x$se[2]),
                        lmitem = sapply(lmres_ranitem, function(x) x$se[2]),
                        lav = sapply(lvres, function(x) x$se[1]))
  sigmaResult = data.frame(lmfscore = sapply(lmres_fscore, function(x) x$sigma),
                           lmsumcore = sapply(lmres_sumscore, function(x) x$sigma),
                           lmitem = sapply(lmres_ranitem, function(x) x$sigma),
                           lav = sapply(lvres, function(x) x$est[2]^.5))
  if(raw) return(list(coef=coefResult, se= seResult, sigma=sigmaResult))
  data.frame(Parameter = c('coef','se','sigma'),
             rbind(colMeans(coefResult, na.rm=T),
                   colMeans(seResult, na.rm=T),
                   colMeans(sigmaResult, na.rm=T)))
}