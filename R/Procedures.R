### Main procedures

#' Global SSR minimizer for structural change model
#'
#'`doglob()` identify if the structural change model is i) pure or ii)
#' partial change model. The procedure then calls appropriate functions \code{\link[mbreaks]{dating}} to estimate
#' the pure change model and \code{\link[mbreaks]{nldat}} to estimate the partial change model.
#'
#'@param y matrix of dependent variable
#'@param z matrix of independent variables with coefficients allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m number of breaks in the structural change model
#'@param h Minimum segment length of regime considered in estimation. If users want to specify a particular value, please set `eps1=0`
#'@param eps convergence criterion for iterative recursive computation. (For partial change model ONLY)
#'@param eps1 trimming level
#'@param maxi maximum number of iterations. (For partial change model ONLY)
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation (Must be a `p x 1` matrix, where p is number of x variables)
#'@param printd  Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#'@return A list containing the following components:
#'\describe{
#'\item{glb}{Minimum global SSR.}
#'\item{datevec}{Vector of dates (optimal minimizers).}
#'\item{bigvec}{Associated SSRs with possible break dates combination.}
#'}
#'@export

doglob = function (y,z,x,m,eps,h,maxi,fixb,betaini,printd,eps1){

#check if model is pure or partial change
if (is.null(x)) {p = 0}
else {p = dim(x)[2]}
q = dim(z)[2]
bigT = dim(y)[1]

#check user specifications
conditions = check_trimming(bigT,eps1,m,h,p,q)
h=conditions$h
m=conditions$m


if(p == 0) {
  if (printd == 1){
  cat('This is a pure structural change model with the following specifications:\n')
  cat(paste(q,'regressors z with allowed to change coefficients\n'))
  cat(paste('maximum number of breaks:',m),'\n') }
  out = dating(y,z,h,m,q,bigT)
  }
else{
  if (printd == 1){
  cat('This is a partial structural change model with the following specifications:\n')
  cat(paste(q,'regressors z with allowed to change coefficients\n'))
  cat(paste(p,'regressors x with fixed coefficients\n'))
  cat(paste('maximum number of breaks:',m),'\n')
  if(fixb == 1) {
    check_beta0(betaini,p)
    cat('initial regime-wise coefficients: \n')
    cat(betaini)}
  cat(paste('convergence criterion:',eps),'\n')
  cat(paste('print iteration option (1)=TRUE/(0)=FALSE',printd),'\n')}
  out = nldat(y,z,x,h,m,p,q,bigT,fixb,eps,maxi,betaini,printd)
}
glb = out$glb
datevec = out$datevec
bigvec = out$bigvec
#printing results
if (printd == 1) {
for (i in 1:m){
  cat(paste('Model with',i,'breaks has SSR:',glb[i,1]),'\n')
  cat('The dates of breaks are:\n')
  cat(datevec[1:i,i])
  }
}
  return (out)
}

#' SupF, UDMax & WDMax testing procedure
#'
#' `dotest()` compute the test statistics and report the critical values of
#' the 2 main supF tests below:
#' \itemize{
#' \item SupF test of 0 vs m breaks
#' \item Double Max test proposed by Perron and Bai, 1998
#'}
#'
#'@param y_name matrix of dependent variable
#'@param z_name matrix of regressors which coefficients are allowed to change
#'across regimes.
#'@param x_name matrix of regressors which coefficients are constant across
#' regimes.
#'@param data the data set for estimation
#'@param m maximum number of breaks
#'@param const indicates whether the regression model include an
#'intercept changing across regimes. Default value is 1
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation (Must be a `p x 1` matrix, where p is number of x variables)
#'@param printd  Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#'@param prewhit option to use AR(1) for prewhitening
#'@param robust set to \code{1} to allow for heterogeneity
#' and autocorrelation in the residuals, \code{0} otherwise.
#' The method used is \emph{Andrews(1991)} automatic bandwidth with AR(1) approximation with quadratic
#' kernel. Note: Do not set to \code{1} if lagged dependent variables are
#' included as regressors.
#' @param hetdat option for the construction of the F tests. Set to 1 if want to
#' allow different moment matrices of the regressors across segments.
#' If \code{hetdat} = \code{0}, the same moment matrices are assumed for each segment
#' and estimated from the ful sample. It is recommended to set
#' \code{hetdat}=\code{1} if number of regressors \code{x} > \code{0}.
#' @param hetvar option for the construction of the F tests.Set to \code{1}
#' if users want to allow for the variance of the residuals to be different across segments.
#' If \code{hetvar}=\code{0}, the variance of the residuals is assumed constant
#' across segments and constructed from the full sample. \code{hetvar}=\code{1} when \code{robust} =\code{1})
#' @param hetomega used in the construction of the confidence intervals for the break
#' dates. If \code{hetomega}=\code{0}, the long run covariance matrix of zu is
#' assumed identical across segments
#' (the variance of the errors u if \code{robust}=\code{0})
#' @param hetq used in the construction of the confidence intervals for the break
#' dates. If \code{hetq}=\code{0}, the moment matrix of the data is assumed identical
#' across segments
#' @param printd Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#'@return A list that contains following:
#'\describe{
#'\item{ftest}{SupF test of 0 vs m (1 to maximum) breaks statistics}
#'\item{cv_supF}{Critical values for Sup F test }
#'\item{cv_Dmax}{Critical values for Double Max test}
#'\item{supF1}{table summarizing the SupF test (for viewing purposes)}
#'\item{UDMax}{table summarizing the Double Max test (including UDMax statistics and CVs)}
#'}
#'@export

dotest = function(y_name,z_name=NULL,x_name=NULL,data,
                  m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,prewhit=1,robust=1,
                  hetdat=1,hetvar=1,hetq=1,hetomega=1,const=1){
  siglev=matrix(c(10,5,2.5,1),4,1)
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x
  p = df$p
  q = df$q
  y_name = df$y_name
  z_name = df$z_name
  x_name = df$x_name


  if(eps1 == 0){ stop('The tests are undefined for 0% trimming level')}
  bigT = dim(y)[1]
  h = floor(eps1*bigT)
  #check user specifications
  conditions = check_trimming(bigT,eps1,m,h,p,q)

  h=conditions$h
  eps1=conditions$eps1
  m=conditions$m

  #check for collinearity/duplication of regressors
  check_duplicate(z,x,z_name,x_name)
  v_eps1 = c(0.05,0.10,0.15,0.20,0.25)
  if(!eps1 %in% v_eps1){
    stop('No matching trimming level available for testing')
  }

  out = doglob(y,z,x,m,eps,h,maxi,fixb,betaini,printd,eps1)
  datevec = out$datevec

  #procedure for F test
  #print('a) supF tests against a fixed number of breaks')

  ftest = matrix(0L, nrow = m, ncol = 1)
  wftest = matrix(0L, nrow = m, ncol = 1)

  for (i in 1:m){
    ftest[i,1] = pftest(y,z,i,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar)
    if(printd==1){
    print(paste('supF test for 0 versus',i,'breaks (scaled by q):',ftest[i,1]))}
  }

  cv_supF = matrix(0L,4,m)
  for (sl in 1:4){
    #critical values for supF test
    cv = getcv1(sl,eps1)
    #pad missing values with NA
    if(dim(cv)[2]<m){
      tmp_cv = dim(cv)[2]
      cv_supF[sl,1:tmp_cv] = cv[q,1:tmp_cv,drop=FALSE]
      cv_supF[sl,tmp_cv+1:m] = matrix(NA,1,m - tmp_cv)
    }else{
    cv_supF[sl,] = cv[q,1:m,drop=FALSE]}
    if (printd==1){
    print(paste('The critical values at the',siglev[sl,1],'% level are (for k = 1 to',m,'):'))
    print(cv_supF[sl,1:m,drop=FALSE])}
  }

  #procedure for Dmax and UDmax test

  #print('b) Dmax test against an unknown number of breaks')
  #print(paste('The UDmax test is:',max(ftest)))
  cv_Dmax = matrix(0L,4,1)
  for (sl in 1:4) {
    #critical values for Dmax test
    cvm = getdmax(sl,eps1)
    cv_Dmax[sl,1] = cvm[q,1]
    if(printd==1){
    print(paste('The critical values at the',siglev[sl,1],'% level is:',
                cvm[q,1]))}
  }


  for (sl in 1:4){
    #computation of WDmax test
    cv = getcv1(sl,eps1)
    for( i in 1:m){
      wftest[i,1] = cv[q,1] * ftest[i,1] / cv[q,i]
    }
    if (printd==1){
    print(paste('WDmax test at the',siglev[sl,1],'% level is:',max(wftest)))}
  }
  rownames(cv_supF) = siglev
  rownames(cv_Dmax) = siglev
  if (m > 5){
  UDmax = max(ftest[1:5,])}
  else{
    UDmax = max(ftest)
  }

  out = list('ftest' = ftest, 'cv_supF' = cv_supF,
             'cv_Dmax' = cv_Dmax, 'UDmax' = UDmax)
  out$mbreak = m
  class(out) = 'sbtests'

  out = compile_sbtests(out)
  return(out)
}


#' Sequential Sup F tests
#'
#' `doseqtests()` computes the sequential sup F tests of l versus l+1 for l from
#' 1 to m with each corresponding null hypothesis of maximum number of break is l
#' and alternative hypothesis is l+1. The l breaks under the null hypothesis
#' are taken from the global minimization estimation
#'
#' @param y_name name of dependent variable in the data set
#' @param z_name name of independent variables in the data set which coefficients are allowed to change
#' across regimes. \code{default} is a vector of 1 (Mean-shift model)
#' @param x_name name of independent variables in the data set which coefficients are constant across
#' regimes. \code{default} is `NULL`
#' @param data name of data set used
#' @param const indicates whether the regression model include an
#' intercept changing across regimes. Default value is 1
#' @param m maximum number of breaks
#' @param eps1 value of trimming (in percentage) for the construction
#' and critical values. Minimal segment length `h` will be set
#' at \code{default} = int(\code{eps1}*T) (T is total sample size).
#'  There are five options:
#' \itemize{
#' \item `eps1=0.05` Maximal value of `m` = 10.
#' \item `eps1=0.10` Maximal value of `m` = 8.
#' \item `eps1=.15` Maximal value of `m` = 5.
#' \item `eps1=.20` Maximal value of `m` = 3.
#' \item `eps1=.25` Maximal value of `m` = 2.
#' \item `eps1=0` is not allowed. The test is undefined for no trimming level.
#' }
#' @param prewhit set to \code{1} to apply AR(1) prewhitening prior to estimating
#' the long run covariance matrix.
#' @param robust set to \code{1} to allow for heterogeneity
#' and autocorrelation in the residuals, \code{0} otherwise.
#' The method used is \emph{Andrews(1991)} automatic bandwidth with AR(1) approximation with quadratic
#' kernel. Note: Do not set to \code{1} if lagged dependent variables are
#' included as regressors.
#' @param hetdat option for the construction of the F tests. Set to 1 if want to
#' allow different moment matrices of the regressors across segments.
#' If \code{hetdat} = \code{0}, the same moment matrices are assumed for each segment
#' and estimated from the ful sample. It is recommended to set
#' \code{hetdat}=\code{1} if number of regressors \code{x} > \code{0}.
#' @param hetvar option for the construction of the F tests.Set to \code{1}
#' if users want to allow for the variance of the residuals to be different across segments.
#' If \code{hetvar}=\code{0}, the variance of the residuals is assumed constant
#' across segments and constructed from the full sample. \code{hetvar}=\code{1} when \code{robust} =\code{1})
#' @param hetomega used in the construction of the confidence intervals for the break
#' dates. If \code{hetomega}=\code{0}, the long run covariance matrix of zu is
#' assumed identical across segments
#' (the variance of the errors u if \code{robust}=\code{0})
#' @param hetq used in the construction of the confidence intervals for the break
#' dates. If \code{hetq}=\code{0}, the moment matrix of the data is assumed identical
#' across segments
#' @param maxi number of maximum iterations for recursive calculations of finding
#' global minimizers.\code{default} = 10 (For partial change model ONLY)
#' @param eps convergence criterion for recursive calculations (For partial change model ONLY)
#' @param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#' the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#' @param betaini Initial \eqn{beta_0} to use in estimation (Must be a `p x 1` matrix, where p is number of x variables)
#' @param printd Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#' @return A list that contains following:
#'\describe{
#'\item{supfl}{SupF(l+1|l) test statistics.}
#'\item{cv}{Critical values for SupF(l+1|l) test.}}
#'
#'@examples
#'doseqtests('inf',c('inflag','lbs','inffut'),data=nkpc,prewhit=0)
#'
#'@export


doseqtests = function(y_name,z_name=NULL,x_name=NULL,data,
                    m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,
                    robust=1,hetdat=1,hetvar=1,hetq=1,hetomega=1,const=1) {

  siglev=matrix(c(10,5,2.5,1),4,1)
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x
  p = df$p
  q = df$q
  if(eps1 == 0){ stop('The tests are undefined for 0% trimming level')}
  bigT = dim(y)[1]
  h = floor(eps1*bigT)

  #check user specifications
  conditions = check_trimming(bigT,eps1,m,h,p,q)

  h=conditions$h
  eps1=conditions$eps1
  m=conditions$m

  #check for collinearity/duplication of regressors
  check_duplicate(z,x,z_name,x_name)

  v_eps1 = c(0.05,0.10,0.15,0.20,0.25)
  if(!eps1 %in% v_eps1){
    stop('No matching trimming level available for testing')
  }

  out = doglob(y,z,x,m,eps,h,maxi,fixb,betaini,printd,eps1)
  datevec = out$datevec
  bigvec = out$bigvec
  supfl = matrix (0L,m,1)
  ndat = matrix (0L,m,1)
  for (i in seq(1,m-1,1)){
    out1 = spflp1(bigvec,datevec[1:i,i,drop=FALSE],i+1,y,z,h,q,prewhit,robust,x,p,hetdat,hetvar)
    supfl[i+1,1] = out1$maxf
    #print(paste('The supF(',i+1,'|',i,') test is',supfl[i,1]))
    #print(paste('It corresponds to a new break at:',ndat[i,1]))
  }
  supfl[1,1] = pftest(y,z,1,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar)
  cv_supFl = matrix(0L,4,m)

  for (c in 1:4){
    cv = getcv2(c,eps1)
    cv_supFl[c,] = cv[q,1:m,drop=FALSE]
  }

  rownames(cv_supFl) = siglev


  out = list('supfl' = supfl, 'cv' = cv_supFl)
  out$mbreak = m
  class(out) = 'seqtests'
  out = compile_seqtests(out)
  return(out)

}

#' Estimating number of breaks via information criterion
#'
#' `doorder()` estimates the number of breaks
#'  using one of the following information criteria:
#' \itemize{
#'  \item modified Bayesian information criterion by Kurozumi and Tuvaandorj, 2011,
#'  \item modified Schwarz information criterion by Liu, Wu and Zidek, 1997,
#'  \item Bayesian information criterion by Yao, 1988} 
#'  and the structural break model corresponding to estimated number of breaks.
#'
#' @param y_name name of dependent variable in the data set
#' @param z_name name of independent variables in the data set which coefficients are allowed to change
#' across regimes. \code{default} is vector of 1 (Mean-shift model)
#' @param x_name name of independent variables in the data set which coefficients are constant across
#' regimes. \code{default} is `NULL`
#'@param data name of data set used
#'@param const indicates whether the regression model include an
#'intercept changing across regimes. Default value is 1
#'@param m maximum number of breaks
#'@param ic indicator which information criterion is used in selecting number of breaks:
#'\itemize{
#'\item \code{KT}
#'\item \code{BIC}
#'\item \code{LWZ} 
#'}
#'The default value is \code{KT}
#' @param eps1 value of trimming (in percentage) for the construction
#' and critical values. Minimal segment length `h` will be set
#' at \code{default} = int(\code{eps1}*T) (T is total sample size).
#  There are six options:
#' There are five options:
#' \itemize{
#' \item `eps1=0.05` Maximal value of `m` = 10.
#' \item `eps1=0.10` Maximal value of `m` = 8.
#' \item `eps1=.15` Maximal value of `m` = 5.
#' \item `eps1=.20` Maximal value of `m` = 3.
#' \item `eps1=.25` Maximal value of `m` = 2.
#' \item `eps1=0` This option allows users to explicitly specify
#' minimum segment length `h` parameters}
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd  Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#'@param prewhit set to \code{1} to apply AR(1) prewhitening prior to estimating
#' the long run covariance matrix.
#' @param robust set to \code{1} to allow for heterogeneity
#' and autocorrelation in the residuals, \code{0} otherwise.
#' The method used is \emph{Andrews(1991)} automatic bandwidth with AR(1) approximation with quadratic
#' kernel. Note: Do not set to \code{1} if lagged dependent variables are
#' included as regressors.
#' @param hetdat option for the construction of the F tests. Set to 1 if want to
#' allow different moment matrices of the regressors across segments.
#' If \code{hetdat} = \code{0}, the same moment matrices are assumed for each segment
#' and estimated from the ful sample. It is recommended to set
#' \code{hetdat}=\code{1} if number of regressors \code{x} > \code{0}.
#' @param hetvar option for the construction of the F tests.Set to \code{1}
#' if users want to allow for the variance of the residuals to be different across segments.
#' If \code{hetvar}=\code{0}, the variance of the residuals is assumed constant
#' across segments and constructed from the full sample. \code{hetvar}=\code{1} when \code{robust} =\code{1})
#' @param hetomega used in the construction of the confidence intervals for the break
#' dates. If \code{hetomega}=\code{0}, the long run covariance matrix of zu is
#' assumed identical across segments
#' (the variance of the errors u if \code{robust}=\code{0})
#' @param hetq used in the construction of the confidence intervals for the break
#' dates. If \code{hetq}=\code{0}, the moment matrix of the data is assumed identical
#' across segments
#' @param h Minimum segment length of regime considered in estimation. If users want to specify a particular value, please set `eps1=0`
#'@return A list of class `model` that contains one of the following:
#'\describe{
#'\item{mBIC}{change model with number of breaks selected by BIC}
#'\item{mLWZ}{change model with number of breaks selected by LWZ}
#'\item{mKT}{change model with number of breaks selected by KT}
#'}
#'
#'@examples
#'doorder('rate',data=real,ic=c('BIC'))
#'
#'@references Liu J, Wu S, Zidek JV (1997). \emph{"On Segmented Multivariate Regressions"},
#'Statistica Sinica, 7, 497-525.
#'Yao YC (1988). \emph{"Estimating the Number of Change-points via Schwartz Criterion"},
#' Statistics and Probability Letters, 6, 181-189.
#'Kurozumi E, Tuvaandorj P (2011). \emph{"Model Selection Criteria in Multivariate Models with
#'Multiple Structural Changes"}, Journal of Econometrics 164, 218-238.
#'@export

doorder = function(y_name,z_name = NULL,x_name = NULL,data,
                   m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,
                   betaini=0,printd=0,ic='KT',const=1,h=NULL,
                   prewhit=1,hetdat=1,hetq=1,hetomega=1,
                   hetvar=1,robust=1) {

  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x
  p = df$p
  q = df$q
  bigT = dim(y)[1]

  #check user specifications
  conditions = check_trimming(bigT,eps1,m,h,p,q)

  h=conditions$h
  eps1=conditions$eps1
  m=conditions$m

  #check for collinearity/duplication of regressors
  check_duplicate(z,x,z_name,x_name)

  if (p == 0){zz = z}
  else{zz = cbind(z,x)}
  #if(eps1 == 0){
    #recompute eps1_star with matching h when eps1 = 0 by user input
   # eps1_s = bigT/h/100
    #temp = doglob(y,z,x,m,eps,eps1_s,maxi,fixb,betaini,printd)
  #}else{
  temp = doglob(y,z,x,m,eps,h,maxi,fixb,betaini,printd,eps1)
  #}
  glb = temp$glb
  ssr0 = nssr(y,zz)
  delta0 = 0.1 #optimal parameters in LWZ paper
  c0 = 0.299
  glob= matrix(0L, nrow = m+1, ncol=1)
  glob[1,1] = ssr0
  glob[seq(2,m+1),1] = glb
  datevec = temp$datevec

  bic = matrix(0L,nrow = m+1, ncol = 1)
  lwz = matrix(0L,nrow = m+1, ncol = 1)
  kt  = matrix(0L,nrow = m+1, ncol = 1)
  for (i in seq(1,m+1)){
    #BIC criterion
    bic [i,1] = log(glob[i,1]/bigT) + log(bigT)*(i-1)*(q+1)/bigT
    #LWZ criterion
    lwz[i,1] = log(glob[i,1]/(bigT-i*q-i+1)) +
      ((i-1)*(q+1)*c0*(log(bigT))^(2+delta0))/bigT
    #Kurozumi and Tuvaandori (2011)
    if (i==1){
      bd=c(0,bigT)}
    else{
      bd=c(0,datevec[1:i-1,i-1],bigT)}
   for (l in seq(1,i)){
    segy   = y[seq(bd[l]+1,bd[l+1],1),,drop=FALSE];
    segz   = z[seq(bd[l]+1,bd[l+1],1),,drop=FALSE];
    segres = segy-segz%*%solve(t(segz)%*%segz)%*%t(segz)%*%segy;
    dt     = bd[l+1]-bd[l];
    kt[i,1]= kt[i,1]+(dt*log(t(segres)%*%segres/dt)+q*log(dt));
    }

    kt[i,1]    = kt[i,1]+2*i*log(bigT);
  }

  mBIC = which.min(bic) - 1
  mLWZ = which.min(lwz) - 1
  mKT = which.min(kt) - 1


  if (ic == 'BIC'){
    mSEL=mBIC
    p_name = 'BIC'
  }else if(ic == 'LWZ') {
    mSEL=mLWZ
    p_name = 'LWZ'
  }else if(ic == 'KT') {
    mSEL=mKT
    p_name = 'KT'
  }else{
    stop('No such criterion found. Please select either BIC, LWZ or KT')
  }

  if (mSEL == 0){
    message('\nThere are no breaks selected by ',p_name,' and estimation is skipped\n')
    out = list()
    out$p_name = p_name
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
  }
  else{
    date = temp$datevec[seq(1,mSEL,1),mSEL,drop=FALSE]
    out = estim(mSEL,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = p_name
    out$nbreak = mSEL
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out$const = const
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out = compile_model(out)
    return(out)
  }

}




#' Estimating number of breaks using sequential tests
#'
#' `dosequa()` sequentially increases the number of breaks from `1` to `m`
#' until the sequential tests reject and estimate the structural change model
#' with corresponding estimated breaks. The procedure is proposed by
#' Bai and Perron, 1998.
#'
#' @param y_name name of dependent variable in the data set
#' @param z_name name of independent variables in the data set, which coefficients are allowed to change
#' across regimes. Default value is vector of 1 (Mean-shift model).
#' @param x_name name of independent variables in the data set, which coefficients are constant across
#' regimes. Default value is \code{NULL}.
#' @param data name of the data set used
#' @param m maximum number of breaks
#' @param eps1 value of trimming (in percentage) for the construction
#' and critical values. Minimal segment length `h` will be set
#' at default value = \code{int(eps1 * T)} (T is total sample size).
#' There are five options:
#' \itemize{
#'   \item \code{eps1 = 0.05} Maximal value of \code{m} = 10.
#'   \item \code{eps1 = 0.10} Maximal value of \code{m} = 8.
#'   \item \code{eps1 = 0.15} Maximal value of \code{m} = 5.
#'   \item \code{eps1 = 0.20} Maximal value of \code{m} = 3.
#'   \item \code{eps1 = 0.25} Maximal value of \code{m} = 2.
#'   \item \code{eps1 = 0} This option is not allowed.
#' }
#' @param eps convergence criterion for iterative recursive computation
#' @param maxi maximum number of iterations
#' @param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#' the model will use values given in \code{betaini}. If \code{0}, \code{betaini} is skipped
#' @param betaini Initial \eqn{\beta_0} to use in estimation
#' @param printd Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#' @param prewhit set to \code{1} to apply AR(1) prewhitening prior to estimating
#' the long run covariance matrix.
#' @param robust set to \code{1} to allow for heterogeneity
#' and autocorrelation in the residuals, \code{0} otherwise. 
#' The method used is \emph{Andrews(1991)} automatic bandwidth with AR(1) approximation with quadratic
#' kernel. Note: Do not set to \code{1} if lagged dependent variables are
#' included as regressors.
#' @param hetdat option for the construction of the F tests. Set to 1 if you want to
#' allow different moment matrices of the regressors across segments. 
#' If \code{hetdat} = \code{0}, the same moment matrices are assumed for each segment
#' and estimated from the full sample. It is recommended to set
#' \code{hetdat = 1} if number of regressors \code{x} > \code{0}.
#' @param hetvar option for the construction of the F tests. Set to \code{1}
#' if users want to allow for the variance of the residuals to be different across segments.
#' If \code{hetvar} = \code{0}, the variance of the residuals is assumed constant
#' across segments and constructed from the full sample. \code{hetvar} = \code{1} when \code{robust} = \code{1}
#' @param hetomega used in the construction of the confidence intervals for the break
#' dates. If \code{hetomega} = \code{0}, the long run covariance matrix of zu is
#' assumed identical across segments (the variance of the errors u if \code{robust} = 0).
#' @param hetq used in the construction of the confidence intervals for the break
#' dates. If \code{hetq} = \code{0}, the moment matrix of the data is assumed identical
#' across segments.
#' @param signif significance level used in the sequential test to select number of breaks.
#' \itemize{
#'   \item 4: 1\% level
#'   \item 3: 2.5\% level
#'   \item 2: 5\% level
#'   \item 1: 10\% level
#' }
#' @param const indicates whether the regression model includes an
#' intercept changing across regimes. Default value is 1
#'
#' @return A list of `model` class with the number of breaks selected by sequential tests.
#' 
#' @examples
#' dosequa('rate', data = real, signif = 1)
#' 
#' @export

dosequa = function(y_name,z_name=NULL,x_name=NULL,data,
                   m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                   prewhit=1,robust=1,hetdat=1,hetvar=1,hetq=1,hetomega=1,const=1,signif=2) {


  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x
  p = df$p
  q = df$q
  bigT = dim(y)[1]

  #check user specifications
  conditions = check_trimming(bigT,eps1,m,h,p,q)
  h=conditions$h
  eps1=conditions$eps1
  m=conditions$m

  #check for collinearity/duplication of regressors
  check_duplicate(z,x,z_name,x_name)

  nbreak = 0
  siglev=matrix(c(10,5,2.5,1),4,1)

if(eps1 == 0){
  stop('The estimation procedure is invalid since it requires testing, but eps1 is set to 0')
}

   # print(paste('Output from the sequential procedure at significance level',
   #              siglev[j,1],'%'))
    out_seq = sequa(m,signif,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
    nbr = out_seq$nbreak


    #print(paste('The sequential procedure estimated the number of breaks at:',nbr))
    if (nbr > 0) {datese = as.matrix(out_seq$dv0)}
    else{cat("\nThere are no breaks selected by sequential procedure and estimation is skipped\n")
      out = list()
      out$p_name = 'dosequa'
      out$nbreak = nbr
      class(out) = 'model'
      return(out)
    }

    nbreak = nbr
    if (nbr!=0){
      dateseq = t(datese)
    }

  mSEL = nbreak

  if (mSEL == 0){
    message('\nThere are no breaks selected by sequential and estimation is skipped\n')
    out = list()
    out$p_name = 'dosequa'
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
    }
  else{
    date = dateseq
    date = t(date)
    out = estim(mSEL,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = 'dosequa'
    out$nbreak = mSEL
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out$const = const
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out$signif = signif
    out = compile_model(out)
    return(out)
  }


}


#' Estimating number of breaks using repartition procedure
#'
#' `dorepart()` computes the repartition estimates of the breaks obtained
#' by the sequential method by Bai, 1995. It allows estimates that have the same
#' asymptotic distribution as those obtained by global minimization. Otherwise,
#' the output from the procedure "estim" below does not deliver asymptotically
#' correct confidence intervals for the break dates.
#'
#' @param y_name name of dependent variable in the data set
#' @param z_name name of independent variables in the data set, whose coefficients are allowed to change
#' across regimes. \code{default} is a vector of 1 (Mean-shift model).
#' @param x_name name of independent variables in the data set whose coefficients are constant across
#' regimes. \code{default} is \code{NULL}.
#' @param data name of the data set used
#' @param m maximum number of breaks
#' @param eps1 value of trimming (in percentage) for the construction
#' and critical values. Minimal segment length \code{h} will be set
#' at \code{default} = \code{int(eps1 * T)} (T is total sample size).
#' There are five options:
#' \itemize{
#'   \item \code{eps1 = 0.05} Maximal value of \code{m} = 10.
#'   \item \code{eps1 = 0.10} Maximal value of \code{m} = 8.
#'   \item \code{eps1 = 0.15} Maximal value of \code{m} = 5.
#'   \item \code{eps1 = 0.20} Maximal value of \code{m} = 3.
#'   \item \code{eps1 = 0.25} Maximal value of \code{m} = 2.
#'   \item \code{eps1 = 0} This option is not allowed.
#' }
#' @param eps convergence criterion for iterative recursive computation
#' @param maxi maximum number of iterations
#' @param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#' the model will use values given in \code{betaini}. If \code{0}, \code{betaini} is skipped
#' @param betaini Initial \eqn{\beta_0} to use in estimation
#' @param printd Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#' @param m Maximum number of structural changes allowed. If not specified,
#' \code{m} will be set to \code{default} value matching the \code{eps1} input
#' @param const indicates whether the regression model includes an
#' intercept changing across regimes. Default value is 1
#' @param prewhit set to \code{1} to apply AR(1) prewhitening prior to estimating
#' the long run covariance matrix.
#' @param robust set to \code{1} to allow for heterogeneity
#' and autocorrelation in the residuals, \code{0} otherwise. 
#' The method used is \emph{Andrews(1991)} automatic bandwidth with AR(1) approximation with quadratic
#' kernel. Note: Do not set to \code{1} if lagged dependent variables are
#' included as regressors.
#' @param hetdat option for the construction of the F tests. Set to 1 if you want to
#' allow different moment matrices of the regressors across segments. 
#' If \code{hetdat} = \code{0}, the same moment matrices are assumed for each segment
#' and estimated from the full sample. It is recommended to set
#' \code{hetdat = 1} if number of regressors \code{x} > \code{0}.
#' @param hetvar option for the construction of the F tests. Set to \code{1}
#' if users want to allow for the variance of the residuals to be different across segments.
#' If \code{hetvar} = \code{0}, the variance of the residuals is assumed constant
#' across segments and constructed from the full sample. \code{hetvar} = \code{1} when \code{robust} = \code{1}
#' @param signif significance level used to sequential test to select number of breaks.
#' \itemize{
#'   \item 4: 1\% level
#'   \item 3: 2.5\% level
#'   \item 2: 5\% level
#'   \item 1: 10\% level
#' }
#'
#' @return A list of class \code{model} for the structural break model estimated by
#' the repartition procedure.
#' 
#' @examples
#' dorepart('inf', 'inflag', 'inffut', data = nkpc)
#' 
#' @references
#' Bai, J. 1995, \emph{"Estimating Breaks One at a Time"}, Econometric Theory, 13,
#' 315-352
#' 
#' @export


dorepart = function(y_name,z_name = NULL,x_name = NULL,data,
                    m=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,robust=1,hetdat=1,hetvar=1,const=1,signif=2){

  #check dataset conformity
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x
  p = df$p
  q = df$q
  bigT = dim(y)[1]

  #check user specifications
  conditions = check_trimming(bigT,eps1,m,h,p,q)

  h=conditions$h
  eps1=conditions$eps1
  m=conditions$m

  #check for collinearity/duplication of regressors
  check_duplicate(z,x,z_name,x_name)

  if(eps1 == 0){
    stop('The estimation procedure is invalid since it requires testing, but eps1 is set to 0')
  }

  reparv = matrix (0L,4,m)
  siglev=matrix(c(10,5,2.5,1),4,1)


  temp = sequa(m,signif,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
  nbreak = temp$nbreak
    if (temp$nbreak == 0){
      message('\nThere are no breaks selected by sequential procedure and the repartition procedure is skipped.\n')
      out=list()
      out$p_name = 'dorepart'
      out$nbreak = 0
      class(out) = 'model'
      return(out)
    }
    else {
      repartda = preparti(y,z,temp$nbreak,
                          temp$dv0,
                          h,x,p)
      reparv = repartda
    }


  #estimate the date at 5% significant level
  mSEL = nbreak


  if (mSEL == 0){
    message('\nThere are no breaks selected by sequential procedure and estimation is skipped\n')
    out = list()
    out$p_name = 'dorepart'
    out$nbreak = mSEL
    class(out) = 'model'
    return(out)
  }
  else{
    date = reparv
    hetq=1
    hetomega=1
    out = estim(mSEL,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = 'dorepart'
    out$nbreak = mSEL
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out$const = const
    out$y_name = y_name
    out$z_name = z_name
    out$x_name = x_name
    out$y = y
    out$x = x
    out$z = z
    out$signif = signif
    out = compile_model(out)
    return(out)
  }


  return (reparv)
}

#'Estimate a model with pre-specified number of breaks
#'
#'`dofix()` compute a structural change model with pre-specified number
#'of breaks.
#'
#' @param y_name name of dependent variable in the data set
#' @param z_name name of independent variables in the data set which coefficients are allowed to change
#' across regimes. \code{default} is vector of 1 (Mean-shift model)
#' @param x_name name of independent variables in the data set which coefficients are constant across
#' regimes. \code{default} is `NULL`
#'@param data name of data set used
#'@param const indicates whether the regression model include an
#'intercept changing across regimes. Default value is 1
#'@param fixn number of breaks specified
#' @param eps1 value of trimming (in percentage) for the construction
#' and critical values. Minimal segment length `h` will be set
#' at \code{default} = int(\code{eps1}*T) (T is total sample size).
#  There are five options:
#' \itemize{
#'   \item \code{eps1 = 0.05} Maximal value of \code{m} = 10.
#'   \item \code{eps1 = 0.10} Maximal value of \code{m} = 8.
#'   \item \code{eps1 = 0.15} Maximal value of \code{m} = 5.
#'   \item \code{eps1 = 0.20} Maximal value of \code{m} = 3.
#'   \item \code{eps1 = 0.25} Maximal value of \code{m} = 2.
#'   \item \code{eps1=0} This option allows users to explicitly specify
#' minimum segment length `h` parameters.
#' }
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd  Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#'@param prewhit set to \code{1} to apply AR(1) prewhitening prior to estimating
#' the long run covariance matrix.
#' @param robust set to \code{1} to allow for heterogeneity
#' and autocorrelation in the residuals, \code{0} otherwise.
#' The method used is \emph{Andrews(1991)} automatic bandwidth with AR(1) approximation with quadratic
#' kernel. Note: Do not set to \code{1} if lagged dependent variables are
#' included as regressors.
#' @param hetdat option for the construction of the F tests. Set to 1 if want to
#' allow different moment matrices of the regressors across segments.
#' If \code{hetdat} = \code{0}, the same moment matrices are assumed for each segment
#' and estimated from the ful sample. It is recommended to set
#' \code{hetdat}=\code{1} if number of regressors \code{x} > \code{0}.
#' @param hetvar option for the construction of the F tests.Set to \code{1}
#' if users want to allow for the variance of the residuals to be different across segments.
#' If \code{hetvar}=\code{0}, the variance of the residuals is assumed constant
#' across segments and constructed from the full sample. \code{hetvar}=\code{1} when \code{robust} =\code{1})
#' @param hetomega used in the construction of the confidence intervals for the break
#' dates. If \code{hetomega}=\code{0}, the long run covariance matrix of zu is
#' assumed identical across segments
#' (the variance of the errors u if \code{robust}=\code{0}).
#' @param hetq used in the construction of the confidence intervals for the break
#' dates. If \code{hetq}=\code{0}, the moment matrix of the data is assumed identical
#' across segments
#' @param h Minimum segment length of regime considered in estimation. If users want to specify a particular value, please set `eps1=0`
#'
#'@examples
#'dofix('rate',data=real,fixn=3)
#'
#'@return out A list of class `model` contains all information about the
#'estimated structural change model with `fixn` breaks
#'
#'@export

dofix = function(y_name,z_name = NULL,x_name=NULL,data,
                    fixn=5,eps=0.00001,eps1=0.15,maxi=10,fixb=0,betaini=0,printd=0,
                    prewhit=1,robust=1,hetdat=1,hetvar=1,hetq=1,hetomega=1,const=1,h=NULL){

  if(fixn<0){
    warning('\nThe maximum number of breaks cannot be negative, set prespecified breaks = 2\n')
    fixn = 2
  }

  #check dataset conformity
  df = process_data(y_name = y_name,z_name = z_name,x_name = x_name,data=data,const)
  y = df$y
  z = df$z
  x = df$x
  p = df$p
  q = df$q
  y_name = df$y_name
  z_name = df$z_name
  x_name = df$x_name
  bigT = dim(y)[1]

  #check user specifications
  conditions = check_trimming(bigT,eps1,fixn,h,p,q)

  h=conditions$h
  eps1=conditions$eps1
  fixn=conditions$m

  #check for collinearity/duplication of regressors
  check_duplicate(z,x,z_name,x_name)

  #if(eps1 == 0){
    #recompute eps1_star with matching h when eps1 = 0 by user input
    #eps1_s = bigT/h
   # t_out = doglob(y,z,x,fixn,eps,eps1_s,maxi,fixb,betaini,printd)
  #}else{
    t_out = doglob(y,z,x,fixn,eps,h,maxi,fixb,betaini,printd,eps1)
  #}
  datevec = t_out$datevec
if(length(datevec) == 0){
  message('\nThere are no breaks selected by the procedure\n')
  out = list()
  out$p_name = 'dofix'
  out$nbreak = length(datevec)
  class(out) = 'model'
  return(out)
}else{
  date = datevec[,fixn,drop=FALSE]
  hetq=1
  hetomega=1
  fix_mdl = estim(fixn,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
  out = fix_mdl
  out$p_name = 'fix'
  out$nbreak = fixn
  class(out) = 'model'
  out$numz = q
  out$numx = p
  out$const = const
  out$y_name =y_name
  out$x_name =x_name
  out$z_name =z_name
  out$y = y
  out$z = z
  out$x = x
  out = compile_model(out)
return(out)}
}


#'Format output of n break model
#'
#'\code{compile_model()} compiles the information of \code{model} class object \code{x} into three
#'main tables:
#'\describe{
#'\item{date_tab}{Table for estimated break date in the model
#'with 90\% and 95\% confidence intervals based on \code{robust},\code{hetomega}, \code{hetq} options
#'for errors and \code{prewhit} option.}
#'\item{RS_tab}{Table for estimated coefficients for \code{z} regressors
#'with corrected standard errors based on \code{robust},\code{hetdat},\code{hetvar}
#'options for errors and \code{prewhit} option.}
#'\item{FS_tab}{Table for estimated coefficients for \code{x} regressors
#'with corrected standard errors based on \code{robust}, \code{hetdat}, and \code{hetvar} options for errors and \code{prewhit} option.}
#'}
#'@param x the \code{model} class to format
#'@param digits number of digits displayed in console. Default value is 3
#'@return Formatted \code{x} with the following appended tables:
#'\describe{
#'\item{`date_tab`}{A data frame storing the break date estimated by the model,
#'and their corresponding confidence intervals.}
#'\item{`RS_tab`}{A data frame storing the estimated coefficients which allowed to
#'change across regimes with corrected standard errors.}
#'\item{`FS_tab`}{A data frame storing the estimated coefficients which is constant
#'across regimes with corrected standard errors.}
#'}
#' @note
#' \itemize{
#' \item If \code{x} returns 0 number of estimated break,
#' the function will return \code{NULL} value instead of the list in \code{Value}.
#' \item If \code{x} is a pure structural break, the `FS_tab` will return \code{NULL} in
#' \code{Value}.
#' }
#'@export
compile_model = function (x,digits=3){
  if (x$nbreak == 0){return(NULL)}
  else {
    #format date estimation
    CI_95 = c()
    CI_90 = c()
    coln = c()
    for (i in 1:x$nbreak){
      coln = c(coln, paste('Break',i,sep=''))
      CI_95 = c(CI_95,paste('(',x$CI[i,1],',',x$CI[i,2],')',sep=''))
      CI_90 = c(CI_90,paste('(',x$CI[i,3],',',x$CI[i,4],')',sep=''))
    }
    date_tab = data.frame(date = t(x$date),stringsAsFactors = FALSE)
    date_tab = rbind(date_tab,CI_95)
    date_tab = rbind(date_tab,CI_90)
    colnames(date_tab) = coln
    rownames(date_tab) = c('Date','95% CI','90% CI')


    #format regime-varying coefficients
    rnameRS = c()
    for (i in 1:x$numz){
        rnameRS = cbind(rnameRS,x$z_name[i])
      }

    cnameRS = c()
    coefRS = c()
    for (i in 1:(x$nbreak+1)){
    cnameRS = cbind(cnameRS,paste('Regime',i))
    }
    for (j in 1:x$numz){
      coefRSj = c()
    for (i in 1:(x$nbreak+1)){
      coefRSj = cbind(coefRSj,paste(format(round(x$beta[(j-1)*(x$nbreak+1)+i,1],digits),nsmall=digits),
                              paste('(',format(round(x$SE[(j-1)*(x$nbreak+1)+i,1],digits),nsmall=digits),')',sep='')))
      }
      coefRS = rbind(coefRS,coefRSj)
    }

  }
    rnameRSf = paste(rnameRS,'(SE)',sep=' ')
    coef_tabRS=data.frame(co = coefRS)


    rownames(coef_tabRS) = rnameRSf
    colnames(coef_tabRS) = cnameRS


    #format full sample coefficients
    if(x$numx == 0){coef_tabRW = NULL}
    else{
      rnameRW = c()
      for (i in 1:x$numx){
        rnameRW = cbind(rnameRW,x$x_name[i])
      }

    cnameRW = 'Full sample'
    coefRW = c()
    for (j in 1:x$numx){
      coefRW = cbind(coefRW,paste(format(round(x$beta[x$numz*(x$nbreak+1)+j,1],digits),nsmall=digits),
                                  paste('(',format(round(x$SE[x$numz*(x$nbreak+1)+j,1],digits),nsmall=digits),')',sep='')))}

    rnameRWf = paste(rnameRW,'(SE)',sep=' ')

    coef_tabRW=data.frame(co = coefRW)
    colnames(coef_tabRW) = rnameRWf
    rownames(coef_tabRW) = cnameRW}

    table = list('date_tab' = date_tab,'RS_tab' = coef_tabRS, 'FS_tab' = coef_tabRW)
    x$tab = table
    return(x)
  }



#'Summary output of a structural breaks `model`
#'
#'`print` the output of the S3 class `model` with all relevant information:
#'\itemize{
#'\item name of procedure used to obtain number of breaks in the model
#'\item print a table summarizing the break date estimation
#'(including confidence interval for the estimated date)
#'\item print a table summarizing the estimated coefficients for `z` regressors
#'\item print a table summarizing the estimated coefficients for `x` regressors (if any)
#'}
#'
#'@param x object of S3 class `model`
#'@param ... further arguments passed to or from other methods.
#'@return
#'No return value, called for printing to console the following information in `x`:
#'\itemize{
#'\item Basic details of the model: name of prodecures invoked,
#'number of estimated breaks, pure/partial structural change model,global min SSR
#'\item `date_tab` summarizes estimated break dates, see \code{\link{compile_model}}
#'\item `RS_tab` summarizes estimated coefficients allowed to change
#'across regimes, see \code{\link{compile_model}}
#'\item `FS_tab` summarizes estimated coefficients constant across regimes,
#' see \code{\link{compile_model}}
#'}
#'
#'@export

print.model <- function(x,...)
{
  #print procedure used to select number of breaks
  proc = switch(x$p_name,'dosequa'='sequential procedure', 'BIC' = 'BIC', 'LWZ' = 'LWZ',
                'KT'='KT',
                'dorepart' = 'repartition procedure', 'fix' = 'specified number of breaks')
  digits = max(3L, getOption("digits") - 3L)

  if (x$nbreak == 0){
    cat(paste('\nNo breaks were found using',proc),'\n')
  }else{
    if(x$p_name == 'dosequa' || x$p_name == 'dorepart'){
      lev = switch(x$signif,'10%','5%','2.5%','1%')
      cat(paste('\nThe number of breaks is estimated by',proc,'at',lev,'significance level\n'))
    }
      else{
    cat(paste('\nThe number of breaks is estimated by',proc,'\n'))}

  if(x$numx == 0){
    cat(paste('Pure change model with',x$nbreak,'estimated breaks.',sep=' '))
  }else if(x$numx > 0){
    cat(paste('Partial change model with', x$nbreak,'estimated breaks.',sep=' '))
  }
  cat('\nMinimum SSR =',
              format(round(x$SSR,3),nsmall=3),'\n')

  cat('\nEstimated date:\n')
  print(x$tab$date_tab,quote=FALSE)

  cat('\nEstimated regime-specific coefficients:\n')
  print(x$tab$RS_tab,quote=FALSE)


  if(x$numx == 0) {cat('\nNo full sample regressors\n')}
  else{
    cat('\nEstimated full-sample coefficients:\n\n')
    print(x$tab$FS_tab,quote='FALSE')}}


  invisible(x)
}


#' Compile the Output of Sup Wald Test
#'
#' `compile_sbtests` formats the output of `sbtests` into two tables.
#'
#' @param x An `sbtests` class object.
#' @param digits The number of decimal places displayed.
#'
#' @return A modified `sbtests` object, `x`, with two appended data frames:
#' \describe{
#'   \item{supF1}{A data frame containing SupF test statistics for testing 0 versus m breaks,
#'     where m is the maximum number of breaks considered in `x`. It includes critical values
#'     at the \emph{10\%, 5\%, 2.5\%, and 1\%} levels.}
#'   \item{UDMax}{A data frame containing Double Max test statistics with critical values
#'     at the \emph{10\%, 5\%, 2.5\%, and 1\%} levels.}
#' }
#'
#' @export

compile_sbtests <- function(x,digits = 3)
{
  if(x$mbreak == 0){
    return(x)
  }
  else{
  cnames1 = c()
  for (i in 1:x$mbreak){
    if (i == 1){
    cnames1 = cbind(cnames1, paste(i,'break'))}
    else
    {
      cnames1 = cbind(cnames1,paste(i,'breaks'))
    }
  }
  ftest = t(x$ftest)


  supF1 = data.frame(ftest = format(round(ftest,3),nsmall=digits))
  cv_supF = format(round(x$cv_supF,3),nsmall = digits)
  colnames(cv_supF) = colnames(supF1)
  supF1 = rbind(supF1,cv_supF)

  rownames(supF1) = c('Sup F','10% CV','5% CV','2.5% CV','1% CV')
  colnames(supF1) = cnames1

  UDmax = data.frame(UDmax = format(round(x$UDmax,3),nsmall = digits),
                     cv = format(round(t(x$cv_Dmax),3),nsmall = digits))

  colnames(UDmax) = c('UDMax','10% CV','5% CV','2.5% CV','1% CV')


  x$supF1 = supF1
  x$UDmax = UDmax
  return(x)}
}


#'Print Sup F and UDMax tests
#'
#'`print` prints the following information from a `sbtests` class object:
#'\describe{
#'\item{`supF1`}{A table reports sup F tests of 0 versus `1` upto `m` breaks with critical values for
#' \emph{1\%, 2.5\%, 5\%, and 10\%} significance levels.}
#'\item{`UDmax`}{A table reporting Double Max tests with critical values for
#' \emph{1\%, 2.5\%, 5\%, and 10\%} significance levels.}
#'}
#'
#'@param x class `sbtests` object
#'@param ... further arguments passed to or from other methods
#'
#'@return No return value, only for printing formatted `sbtests` class object to console
#'@examples
#'supF = dotest('inf','inflag',data=nkpc)
#'print(supF)
#'
#'@export

print.sbtests <- function(x,...)
{ if(x$mbreak == 0){
  warning('\nThe test is undefined for no break model\n')
}else{
  cat('\na) SupF tests against a fixed number of breaks\n\n')
  print(x$supF1,quote=FALSE)
  cat('\nb) UDmax tests against an unknown number of breaks\n\n')
  print(x$UDmax,quote=FALSE)}
  invisible(x)
}


#'Compile the output of sequential Sup Wald test
#'
#'`compile_seqtests` formats the output of the `seqtests` class object to 1 table
#'\describe{
#'\item{`sfl`}{A table containing sequential sup F tests statistics
#' of `l` versus `l+1` for `l` in `1` to `m` breaks
#'with critical values of the corresponding tests at
#' \emph{1\%, 2.5\%, 5\%, and 10\%} significance levels.}
#'}
#'
#'@param x `seqtests` class object
#'
#'@return class `seqtests` list `x` with appended data frame `sfl` containing the
#'sequential SupF test statistics with critical values at
#'\emph{10\%, 5\%, 2.5\%, and 1\%} level.
#'
#'
#'@export
compile_seqtests = function(x){
  if(x$mbreak==1){
    #message('\nThe test is exactly 0 versus 1 break, hence the sequential test is not repeated\n')
    #x$sfl = NULL
    return(x)
  }else if (x$mbreak==0){
    warning('\nThe test is undefined for maximum break = 0\n')
    x$sfl = NULL
    return(x)
  }
  else{
  nbreak = x$mbreak
  cnames = c()
  for (i in seq(0,nbreak-1)){
    cnames = cbind(cnames, paste('supF(',i+1,'|',i,')',sep=''))
  }
  temp_supfl = format(round(x$supfl,3),nsmall = 3)
  #temp_ndat = format(x$ndat,nsmall=0)
  temp_cv = format(round(x$cv,3),nsmall = 3)

  sfl =data.frame(supfl = t(temp_supfl[seq(1,nbreak,1),1,drop=FALSE]))
  #ndat = t(temp_ndat[seq(1,nbreak,1),1,drop=FALSE])
  cv = temp_cv[,seq(1,nbreak,1),drop=FALSE]
  #colnames(ndat) = colnames(sfl)
  colnames(cv) = colnames(sfl)
  sfl = rbind(sfl,cv)
  colnames(sfl) = cnames
  rownames(sfl) = c('Seq supF','10% CV','5% CV', '2.5% CV', '1% CV')
  x$sfl = sfl
  return (x)}
}

#'Print sequential SupF tests
#'
#'`print` prints the object of class `seqtests` with the following information
#'\itemize{
#'\item Maximum number of breaks `m` in the tests
#'\item `sfl` table with sequential sup F tests statistics of l versus
#'l+1 breaks up to `m` breaks}
#'
#'@param x `seqtests` class object.
#'@param ... further arguments passed to or from other methods.
#'
#'@return No return value, only for printing formatted `seqtests` class object to console
#'@examples
#'seq_supF = doseqtests('inf','inflag',data=nkpc)
#'print(seq_supF)
#'@export
print.seqtests = function(x,...){
  if(x$mbreak==1){
      cat('\nThe test is exactly 0 versus 1 break, hence the sequential test is not repeated\n')
    }else if (x$mbreak==0){
      cat('\nThe test is undefined for maximum break = 0\n')
    }
  else{
  cat('\nsupF(l+1|l) tests using global optimizers under the null\n\n')
  print(x$sfl,quote=FALSE)}
  invisible(x)
}



#' Plot structural change model
#'
#' `plot_model()` visualizes any object of class `model` with comparison between real,
#' fitted values between  model of `m` breaks and null model of `0` breaks with options
#' for confidence interval of break date.
#' @param title title of the graph
#' @param CI confidence intervals for break date and coefficient estimates visualize in terms of
#' fitted values
#' @param model object of class `model` in `mbreaks` package
#' @importFrom ggplot2 ggplot aes annotate geom_segment geom_line ggtitle .data coord_cartesian scale_color_manual geom_ribbon
#' @examples
#' rate = dofix('rate',data=real,fixn=2)
#' plot_model(rate,title='Ex-post US exchange rate')
#'
#' @return No return value, called for plotting class `model` object. For more details
#' on `model` class, see [compile_model]
#'
#' @export

plot_model = function(model,CI=0.95,title=NULL){
  m = model$nbreak
  if(m==0){
    warning('The model has no break. Visualization for comparison between structural breaks
        versus no breaks is skipped')
    return(NULL)
  }


  zreg = model$z
  xreg = model$x
  p = model$numx
  q = model$numz
  date = model$date
  beta = model$beta

  #comparison between structural break vs no break
  y = model$y
  ypred_break = model$fitted
  fixreg = cbind(xreg,zreg)
  fixbeta = OLS(y,fixreg)
  ypred_fix = fixreg%*%fixbeta
  x_t = seq(1,dim(y)[1],1)
  tbl = data.frame(x_t,y,ypred_break,ypred_fix,stringsAsFactors = TRUE)
  colnames(tbl) = c('time','y','ypred_break','ypred_fix')


  #labels and annotations for date's CIs
  date_lab = c()
  for (i in 1:m){
    date_lab = c(date_lab,paste('Date',i))
  }
  vline_seg = data.frame(date,date,rep(Inf,m),rep(-Inf,m))
  colnames(vline_seg) = c('x','xend','y','yend')
  y_pos=c()
  for (i in 1:m){
    y_pos = rbind(y_pos,(10.2+i/2)/10*min(y))
  }
  model$CI[model$CI<0] = 0
  model$CI[model$CI>dim(y)[1]] = dim(y)[1]
  CI_seg95 = data.frame(model$CI[,1],model$CI[,2],y_pos)
  CI_seg90 = data.frame(model$CI[,3],model$CI[,4],y_pos)

  #compute CIs for estimation y
  zbar = diag_par(zreg,m,date)
  if (p == 0){
    reg = zbar
  }
  else{
    reg = cbind(xreg,zbar)
  }

  if(!CI==0.95&&!CI==0.90){
    warning('Not available CI level, set to 95%')
    CI = 0.95
  }
  if(CI == 0.95){
   CI_seg = CI_seg95
   sd = 1.960*model$SE
   beta_lb = beta-sd
   beta_ub = beta+sd
  }else if (CI==0.90){
    CI_seg = CI_seg90
    sd = 1.645*model$SE
    beta_lb = beta-sd
    beta_ub = beta+sd
  }
  tbl$ypred_ub = reg%*%beta_ub
  tbl$ypred_lb = reg%*%beta_lb
  colnames(CI_seg) = c('lb','ub','y')

  if(is.null(title)){
  if(model$numx==0){
    if (m>1){
    title = paste('Pure change with',m,'breaks')}
    else{
      title = paste('Pure change with',m,'break')
    }
  }else{
    if (m>1){
      title = paste('Partial change with',m,'breaks')}
    else{
      title = paste('Partial change with',m,'break')
    }
  }}

  grph = ggplot2::ggplot(data = tbl, ggplot2::aes(x=.data$time))+
    ggplot2::coord_cartesian(xlim=c(0,max(x_t)+1))+
    ggplot2::geom_line(ggplot2::aes(y=.data$ypred_break,color='y_break'),size=0.3)+
    ggplot2::geom_line(ggplot2::aes(y=.data$ypred_fix,color='y_fix'),size=0.3)+
    ggplot2::geom_line(ggplot2::aes(y=.data$y,color='y'),size=0.2)+
    ggplot2::geom_ribbon(ggplot2::aes(ymax=.data$ypred_ub, ymin=.data$ypred_lb), fill="gray", alpha=.35)+
    #ggplot2::geom_ribbon(ggplot2::aes(ymax=.datax.upper, ymin=x.lower), fill="pink", alpha=.5)
    ggplot2::scale_color_manual(name=paste('Legends'),
                                breaks = c('y','y_break','y_fix'),
                                values = c('y'='black','y_break' = 'blue', 'y_fix' = 'red'),
      labels = c(expression(y),expression(hat(y)[m]), expression(hat(y)[0]))
      )+
    ggplot2::annotate('text',x = model$date[,1]-3*max(x_t)/100, y=10/11*max(y),label= date_lab,angle = 90)+
    ggplot2::geom_segment(data=vline_seg,
                          ggplot2::aes(x=.data$x,y=.data$y,
                                       xend=.data$xend,yend=.data$yend),
                          alpha =0.85,colour='purple',linetype='dashed')+
    ggplot2::geom_segment(data=CI_seg,
                          ggplot2::aes(x=.data$lb,xend=.data$ub,y=.data$y,yend=.data$y),
                          colour='red',alpha=0.6)+
  ggplot2::ggtitle(title)+ggplot2::ylab(model$y_name)
  grph
}


