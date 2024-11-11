######## Model estimation given date vector #########

#' Structural change model estimation
#'
#' `estim()` estimates the structural change model by OLS given specified vector of break dates
#' It also computes and reports confidence intervals for the
#' break dates based on asymptotic distributions of break date  and
#' corrected standard errors of coefficients estimates given the structure of covariance matrix
#' for model errors by specifying error options `robust`, `hetomega`, `hetq`, `hetdat` and `hetvar`
#'
#'@param y matrix of dependent variable
#'@param z matrix of regressors with coefficients are allowed to change across
#'regimes
#'@param x matrix of regressors with coefficients are constant across regimes
#'@param q number of `z` regressors z
#'@param p number of regressors x
#'@param m number of breaks
#'@param b vector of break dates
#'@param robust,hetomega,hetq,hetdat,hetvar options for assumptions on the error terms.
#'For more details, please refer to \code{\link{mdl}}.
#'@param prewhit option to use prewhitening process based on AR(1) approximation
#'@return A list containing the following components:
#'\describe{
#'\item{date}{List of estimated breaks.}
#'\item{CI}{List of Confidence Intervals for each corresponding break.}
#'\item{beta}{Estimated coefficients of the regression. The first
#'(\code{m}+1)*\code{q} are coefficients of \code{q} variables \code{z} that change across regimes.
#'The last \code{p} are coefficients of \code{p} variables \code{x}
#'that are constant across regimes.}
#'\item{SE}{Corrected standard errors for the coefficients' estimates}}
#'@export

estim = function(m,q,z,y,b,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar){
  if (m == 0){
    warning('There are no breaks in this model and estimation is skipped')
    return (NULL)}
  else{
    bigT = dim(z)[1]
    d = (m+1)*q + p
    vdel = matrix(0L,nrow = d, ncol = d)
    #construct zbar matrix. Diagonal partition of Z
    #at the estimated break date
    zbar = diag_par(z,m,b)

    #estimation and printing
    if (p == 0){
      reg = zbar
    }
    else{
      reg = cbind(x,zbar)
    }

    #estimation of beta and delta in pure/partial model
    beta = OLS(y,reg)
    vdel = pvdel(y,z,m,q,bigT,b,prewhit,robust,x,p,1,hetdat,hetvar)
    colnames(beta) = 'coefficients'

    SE = matrix(0L,d,1)
    colnames(SE) = 'corrected SE'
    for (i in 1:d){
      #print(paste('Corrected SE for coefficient',i,'is',sqrt(vdel[i,i])))
      SE[i,1] = sqrt(vdel[i,i])
    }

    if (robust == 0 && hetdat == 1 && hetvar == 0){
      warning('In this case robust=0, hetdat=1 and hetvar=0, the "corrected" are the same as that of the printout except for a different small sample correction.')
      SE[i,1] = sqrt(vdel[i,i])
    }

    #confidence interval for break date
    bound = interval(y,z,zbar,b,q,m,robust,prewhit,hetomega,hetq,x,p)
    CI_95 = bound[,c(1,2)]
    CI_90 = bound[,c(3,4)]

    for (i in 1:m){

      #print(paste('The 95% C.I for the',i,'th break is:',bound[i,1],' ',bound[i,2]))
      #print(paste('The 90% C.I for the',i,'th break is:',bound[i,3],' ',bound[i,4]))
    }
    CI = cbind(bound[,1],bound[,2],bound[,3],bound[,4])
    colnames(CI) = c('lower 95% CI','upper 95% CI','lower 90% CI','upper 90% CI')
    rownames(CI) = c(1:m)
    fitted = as.matrix(reg%*%beta)
    resid = as.matrix(y - fitted)

    colnames(resid) = 'residuals'
    colnames(fitted) = 'fitted.values'
    SSR = t(resid)%*%resid
    colnames(SSR) = 'global min SSR'
    rownames(SSR) = c()
    out = list('SE' = SE, 'CI' = CI, 'beta' = beta, 'date' = b,
               'SSR' = SSR, 'resid' = resid,'fitted.values' = fitted)
    return (out)
  }
}


##### Date estimation #######
# Functions and sub-functions to compute global
# minimizers of the model by dynamic programming approach given the SSR objective function

#' Calculate optimal point for segment partition
#'
#' `parti()` searchs for the optimal one break partition for a segment starting
#' from start to last using the vector storing the sum of squared residuals for any
#' segments between `b_start` and `b_last`
#'
#' @param start start date index of the segment
#' @param last end date index of the segment
#' @param b_start first possible break date
#' @param b_end last possible breakdate
#' @param bigT sample period T
#' @return A list containing the following components:
#' \item{ssrmin}{associated SSR of optimal break}
#' \item{dx}{optimal date (global minimizer)}
#'
#' @noRd
parti = function (start,b_start,b_end,last,bigvec,bigT) {
  #initialize array to store betas and value of index for procedure
  dvec = matrix(0L , nrow = bigT, ncol = 1)
  ini = (start-1)*bigT - (start-2)*(start-1)/2 + 1
  j = b_start

  #start the loop on the segment to find optimal break
  while (j <= b_end){
    jj = j - start
    k = j*bigT - (j-1)*j/2+last-j
    dvec[j,1] = bigvec[ini+jj,1] + bigvec[k,1]
    j = j+1
  }
  #get min SSR
  ssrmin = min( dvec[seq(b_start,b_end,1),1] )
  #get the date with lowest SSR
  dx = (b_start - 1) + which.min(dvec[seq(b_start,b_end,1)])
  out = list('ssrmin' = ssrmin, 'dx' = dx)
  return(out)
}

#' Optimal one break partitions with sequential method
#'
#' `partione()` calculates an optimal one break partitions for a segment that
#' starts at date start and ends at date last. It returns the optimal break
#' and the associated SSR. Procedure used with the sequential method when
#' the T*(T+1)/2 vector of SSR is not computed.
#'
#'@param b_start first possible break date
#'@param b_end last possible break date
#'@param last end of segment considered
#'@param vssr minimum SSRs of associated break date
#'@param vssrev recursive SSRs of the model
#'@return A list containing the following components:
#' \item{ssrmin}{associated SSR of optimal break}
#' \item{dx}{optimal date (global minimizer)}
#'
#' @noRd
partione = function(b1,b2,last,vssr,vssrev){
  dvec = matrix(0L,nrow = last, ncol = 1)
  j = b1;
  while (j <= b2) {
    dvec[j,1]=vssr[j,1]+vssrev[last-j,1];
    j=j+1;
  }
  ssrmin=min(dvec[seq(b1,b2),1]);
  dx=(b1-1)+which.min(dvec[seq(b1,b2),1]);
  out = list('ssrmin' = ssrmin, 'dx' = dx)
  return(out)
}


#' Optimal one break partition in partial structural change model
#'
#' `onebp()` computes the optimal one break partition in partial structural
#' change model by searching over all possible breaks given x regressors
#' have unchanged coefficients. Iteration to convergence is used to deal with
#' 2 sets of estimates needed to obtain: full-sample coefficients and regime-specific
#' coefficients
#'
#' @param y matrix of dependent variables
#' @param z matrix of regressors with coefficients
#' allowed to change across regimes
#' @param x matrix of regressors with constant coefficients across regimes
#' @param h minimal length of segment
#' @param start initial date to search
#' @param last last date to search
#'
#' @return A list containing the following components
#' \item{ssrmin}{associated SSR of optimal break}
#' \item{dx}{optimal date (global minimizer)}
#'
#' @noRd

onebp = function(y,z,x,h,start,last) {
  ssrind = 999999999999999
  i = matrix(h,ncol=1)

  while(i <= last - start + 1 - h){

    zb = diag_par(z[start:last,,drop=FALSE],1,i)
    y_reg = y[start:last,1]
    x_reg = cbind(x[start:last,],zb)
    bb = OLS(y_reg,x_reg)
    resid = y_reg - x_reg %*% bb
    ssrn = t(resid) %*% resid

    if (ssrn < ssrind){
      ssrind = ssrn
      bdat = i
    }

    i = i + 1
  }
  bd = bdat + start - 1
  out = list('ssrind' = ssrind,'bd' = bd) #convert code use ssrind to ssrmin, bd to dx
  return(out)
}



#'Computation of global minimizer for pure structural change model
#'
#'`dating()` computes break points that globally minimizes SSR via dynamic programming approach.
#'To avoid recursion depth increases as number of breaks in the model increases, a temporary
#'array is used to store optimal partition with corresponding SSR for all permissible
#'subsamples for all 1:m-1 breaks. For the m-th break, the problem becomes finding where to insert
#'the last feasible m+1-th segment into the sample partitioned by m-1 breaks
#'to obtain minimum SSR over the sample
#'
#'@param y matrix of dependent variable
#'@param z matrix of regressors with coefficients
#' allowed to change across regimes
#'@param h minimum length of segment
#'@param m maximum number of breaks
#'@param q number of `z` regressors
#'@param bigT sample period T
#'
#'@return A list containing the following components:
#'\item{glb}{minimum global SSR}
#'\item{datevec}{Vector of dates (optimal minimizers)}
#'\item{bigvec}{Associated SSRs}
#'
#'@export
dating = function(y,z,h,m,q,bigT){

  #initialize arrays to store results
  datevec = matrix(0L, nrow = m, ncol = m)
  optdat = matrix(0L, nrow = bigT, ncol = m)
  optssr = matrix(0L, nrow = bigT, ncol = m)
  dvec = matrix(0L, nrow = bigT, ncol = 1)
  glb = matrix(0L,nrow = m, ncol = 1)
  bigvec = matrix(0L,nrow = bigT*(bigT+1)/2,ncol = 1)

  #calculate all possible SSR and store
  for (i in 1:(bigT-h+1)) {
    vecssr = ssr(i,y,z,h,bigT)
    bigvec[seq((i-1)*bigT+i-(i-1)*i/2, i*bigT - (i-1)*i/2 ,1),1] = vecssr[seq(i,bigT,1),1]
  }

  #base case: 1 break
  if (m == 1) {
    out = parti(1,h,bigT-h,bigT,bigvec,bigT)
    datevec[1,1] = out$dx
    glb[1,1] = out$ssrmin
  }
  #
  #more than 1 break
  else {
    #change the end point from smallest to full sample T, with case m = 1
    for (j1 in seq(2*h,bigT,1)){
      out = parti(1,h,j1-h,j1,bigvec,bigT)
      optssr[j1,1] = out$ssrmin
      optdat[j1,1] = out$dx
    }
    glb[1,1] = optssr[bigT,1]
    datevec[1,1] = optdat[bigT,1]

    #with case m >= 2
    for (ib in 2:m){
      if (ib == m) {
        jlast = bigT
        for (jb in seq(ib*h,jlast-h,1)) {
          dvec[jb,1] = optssr[jb,ib-1] + bigvec[(jb+1)*bigT - jb*(jb+1)/2,1]
        }
        optssr[jlast,ib] = matrix(t(min(dvec[seq(ib*h,jlast-h,1)])))
        optdat[jlast,ib] = ib*h-1 + which.min(dvec[seq(ib*h,jlast-h,1)])
      }

      else {
        for (jlast in seq((ib+1)*h,bigT,1)){
          for (jb in seq(ib*h,jlast-h,1)){
            dvec[jb,1] = optssr[jb,ib-1] + bigvec[jb*bigT - jb*(jb-1)/2 + jlast -jb,1]
          }
          optssr[jlast,ib] = min(dvec[seq(ib*h,jlast-h,1),1])
          optdat[jlast,ib] = ib*h-1 + which.min(dvec[seq(ib*h,jlast-h,1),1])
        }
      }

      datevec[ib,ib] = optdat[bigT,ib]

      for (i in seq(1,ib-1,1)){
        xx = ib-i
        datevec[xx,ib] = optdat[datevec[xx+1,ib],xx]
      }
      glb[ib,1] = optssr[bigT,ib]
    }
  }

  out = list('glb' = glb, 'datevec' = datevec, 'bigvec' = bigvec)
  return (out)
}

#'Computation of global minimizer for partial structural change model
#'
#'`nldat()` computes the break dates of a partial structural change model
#'for a pre-specified number of breaks `m`. The procedure iterates between
#'estimating the invariant and changing coefficients of `x` and `z` regressors
#'until convergence, by noting that the residuals from linear regression model between
#'`y` and `x` regressors is a pure structural change model,
#' while the residuals from pure structural change model between `y` and `z` regressors
#' is a linear regression
#'
#'@param y dependent variable in matrix form
#'@param z matrix of regressors which coefficients are allowed to change across regimes
#'@param x matrix of regressors which coefficients are constant across regime
#'@param h minimum segment length
#'@param m number of breaks
#'@param p number of `z` regressors
#'@param q number of `x` regressors
#'@param bigT the sample size T
#'@param fixb option to use initial \eqn{\beta} If \code{1}, procedure requires \code{betaini}.
#'If \code{0}, procedure will not use initial beta values
#'@param betaini initial beta values. Required when use with option \code{fixb}
#'@param eps Convergence criterion (For partial change model ONLY)
#'@param maxi Maximum number of iterations (For partial change model ONLY)
#'@param printd option to print output from iterated estimations. If \code{1}, the results
#'for each iteration will be printed in console log. If \code{0}, no output will be printed
#'@return A list containing the following components:
#'\item{glb}{minimum global SSR}
#'\item{datevec}{Vector of dates (optimal minimizers)}
#'\item{bigvec}{Associated SSRs}
#'@references Bai J, Perron P (1998). \emph{"Estimating and Testing Linear Models with Multiple Structural
#'Changes"} Econometrica, 66, 47-78.
#'Bai J, Perron P (2003). \emph{"Computation and Analysis of Multiple Structural Change Models"}
#'Journal of Applied Econometrics 18, 1-22
#'@export

nldat = function(y,z,x,h,m,p,q,bigT,fixb,eps,maxi,betaini,printd){

  #initialize storage
  glb = matrix(0L , nrow = m, ncol = 1)
  globnl = matrix(0L , nrow = m, ncol = 1)
  datevec = matrix(0L, nrow = m, ncol = m)
  datenl = matrix(0L, nrow = m, ncol = m)

  #initialize current max break
  mi = 1
  while(mi <= m){
    if (printd == 1){
      print(paste('Breaks of model',mi))
    }
    if (fixb == 0){
      qq = p+q
      zz = cbind(x,z)

      #initial estimate of the model
      out = dating(y,zz,h,mi,qq,bigT)
      date = out$datevec[1:mi,mi,drop=FALSE]
      bigvec = out$bigvec

      #partition regressors with initial estimated date
      xbar = diag_par(x,mi,date)
      zbar = diag_par(z,mi,date)

      #calculate initial values of estimate
      teta = OLS(y,cbind(zbar,xbar))
      delta1 = teta[seq(1,q*(mi+1),1),1,drop=FALSE]
      beta1 = OLS(y - zbar %*% delta1, x)

      #calculate initial SSR of the model
      resid1 = y - x%*%beta1 - zbar%*%delta1
      ssr1 = t(resid1) %*% resid1

      if(printd==1) {
        print('The iterations are initialized with')
        prmatrix(delta1)
        prmatrix(beta1)
        print('With break date')
        print(date)}
    }
    else {
      check_beta0(betaini,p)
      beta1 = betaini
      ssr1 = -5
    }

    #start the iterations
    length = 999999999999
    i = 1

    while (length > eps) {

      out = dating(y-x%*%beta1,z,h,mi,q,bigT)
      #store the date vector for current max
      date = out$datevec[1:mi,mi,drop=FALSE]
      bigvec = out$bigvec
      zbar = diag_par(z,mi,date)

      #update estimates based on new partition
      teta1 = OLS(y,cbind(x,zbar))
      beta1 = teta1[seq(1,p,1),1,drop=FALSE]
      delta1 = teta1[seq(p+1,p+q*(mi+1),1),1,drop=FALSE]


      #check convergence condition
      resid_n = y - cbind(x,zbar) %*% teta1
      ssrn = t(resid_n) %*% resid_n
      length = abs(ssrn - ssr1)

      if(printd==1){
        print(paste('Iteration',i))
      }

      #check upper bound of iterations and update
      if (i >= maxi){
        print('Number of iterations has reached MAX ITER')
      }
      else {
        i = i+1
        ssr1 = ssrn
        glb[mi,1]=ssrn
        datevec[1:mi,mi] = date
      }

    }
    #finished current max breaks & update
    mi = mi + 1

  }


  out = list('glb' = glb, 'datevec' = datevec, 'bigvec' = bigvec)
  return(out)
}


#' Sequential procedure to obtain number of breaks and break dates
#'
#' `sequa()` compute the sequential Ftest(l+1|l) statistics and estimate the
#' corresponding break date based on the decision rule explained in Bai and Perron, 1998
#'
#' @param y dependent variable in matrix form
#' @param z matrix of regressors which coefficients are allowed to change across regimes
#' @param x matrix of regressors which coefficients are constant across regime
#' @param h minimum segment length
#' @param m number of breaks
#' @param p number of `z` regressors
#' @param q number of `x` regressors
#' @param bigT sample size T
#' @param signif significant level used in hypothesis test for decision
#' rule regarding continuation of estimating next break:
#' * 1: 10% significance level
#' * 2: 5% significance level
#' * 3: 2.5% significance level
#' * 4: 1% significance level
#' @param eps1 trimming level: `{0.05,0.10,0.15,0.20,0.25}`
#' @param q number of `z` regressors with changing coefficients across regimes
#' @param prewhit prewhitening procedure proposed by Andrews and Monahan (1992), that
#' is, a VAR(1) filter applied to wt*ut, where ut is regression residuals and the HAC
#' covariance matrix estimator is constructed based on the filtered series and the AR1
#' coefficient estimate are parametrically accounted for. The default value is 1.
#' @param robust allows heteroskedasticity and autocorrelation in estimated residual by using the HAC co-
#' variance matrix estimator in which the Quadratic Spectral kernel is used with the
#' bandwidth selected via the AR(1) approximation proposed by Andrews (1991). Default value is 1
#' @param hetvar allows for the variance of errors to be different across the segments
#' determined by the estimated breaks dates when constructing the F test statistics. If
#' set hetvar=0, the errors are assumed to have the same variance across the segments.
#' The default value is hetvar=1 (Note that hetvar=0 is not allowed when robust=1)
#' @param hetdat allows for the second moment matrices of wt to be di¤erent across the segments when constructing the F test. If you set hetdat=0,
#' wt is assumed to have the same second moment matrix across the segments. The default value is hetdat=1.
#' @references Bai J, Perron P (1998). \emph{"Estimating and Testing Linear Models with Multiple Structural
#' Changes"} Econometrica, 66, 47-78.
#' @noRd

sequa = function(m,signif,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1){

  dv = matrix(0L, nrow = m+2, ncol = 1)
  dv2 = matrix(0L, nrow = m+2, ncol = 1)
  ftestv = matrix(0L, nrow = m+1,ncol = 1)

  cv = getcv2(signif,eps1)
  dv[1,1] = 0

  if (p == 0){
    y_rev = rot90(rot90(y))
    z_rev = rot90(rot90(z))
    vssrev = ssr(1,y_rev,z_rev,h,bigT)
    vssr = ssr(1,y,z,h,bigT)
    out = partione(h,bigT-h,bigT,vssr,vssrev)
    datx = out$dx
    ssrmin = out$ssrmin
  }
  else{
    out = onebp(y,z,x,h,1,bigT)
    datx = out$bd
    ssrmin = out$ssrind
  }

  dv[2,1] = datx

  ftest=pftest(y,z,1,q,bigT,dv[2,1,drop=FALSE],prewhit,robust,x,p,hetdat,hetvar)

  if (ftest < cv[q,1]) {
    nbreak = 0
    dv[2,1] = 0
    #dv0 = dv[seq(2,nbreak+1,1),1]
    nseg = 1
  }
  else{
    #print(paste('First break found at:',datx))
    nbreak = 1
    nseg = 2
    dv[nseg+1,1] = bigT
  }

  while(nseg <= m){
    ds = matrix(0L,nseg+1,1)
    ftestv = matrix(0L,nseg+1,1)

    i_s = 1

    while(i_s <= nseg){
      length = dv[i_s+1,1] - dv[i_s,1]

      if(length >= 2*h){
        if(p==0){
          y_temp = y[seq(dv[i_s,1]+1,dv[i_s+1,1]),1,drop=FALSE]
          z_temp = z[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          vssr = ssr(1,y_temp,z_temp,h,length)
          y_temp_rev = rot90(rot90(y_temp))
          z_temp_rev = rot90(rot90(z_temp))
          vssrev = ssr(1,y_temp_rev,z_temp_rev,h,length)
          out = partione(h,length-h,length,vssr,vssrev)
          ds[i_s,1] = out$dx
          ftestv[i_s,1] = pftest(y_temp,z_temp,1,q,length,ds[i_s,1,drop=FALSE],prewhit,
                                 robust,0,p,hetdat,hetvar)
        }
        else{
          y_temp = y[seq(dv[i_s,1]+1,dv[i_s+1,1]),1,drop=FALSE]
          z_temp = z[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          x_temp = x[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          out = onebp(y,z,x,h,dv[i_s,1]+1,dv[i_s+1,1])
          ds[i_s,1] = out$bd - dv[i_s,1]
          ftestv[i_s,1] = pftest(y_temp,z_temp,1,q,length,ds[i_s,1],
                                 prewhit,robust,x_temp,p,hetdat,hetvar)
        }
      }
      else{
        ftestv[i_s,1] = 0.0
      }
      i_s = i_s + 1
    }

    maxf = max(ftestv[seq(1,nseg,1),1])

    if (maxf < cv[q,nseg]){
      #print(nbreak)
      #dv0 = dv[seq(2,nbreak+1,1),1]
    }
    else {
      newseg = which.max(ftestv[seq(1,nseg),1])
      dv[nseg+2,1] = ds[newseg,1] + dv[newseg,1]
      nbreak = nbreak + 1
      #check this sort
      dv2 = sort(dv[seq(2,nseg+2,1),1])
      dv2 = matrix(dv2, ncol = 1)
      dv[1,1] = 0
      dv[seq(2,nseg+2,1),1] = dv2

    }
    nseg = nseg + 1
  }

  #print('The sequential procedure has reached the upper limit')
  if (nbreak < 1) {dv0 = c()}
  else{
    dv0 = dv[seq(2,nbreak+1,1),1]}
  out = list('nbreak' = nbreak, 'dv0' = dv0)
}

#' Prepartion procedure
#'
#' `preparti()`
#'
#' @param y matrix of dependent variable
#' @param z matrix of `z` regressors with coefficients are allowed to change across regimes
#' @param nbreak number of break dates
#' @param dateseq vector of break date
#' @param h minimum segment length
#' @param x matrix of `x` regressors with coefficients do not change across regimes
#' @param p number of `x` regressors
#'
#' @return matrix of SSR corresponding to partitioned date
#'
#' @noRd
#'

preparti = function(y,z,nbreak,dateseq,h,x,p) {
  bigT = dim(z)[1]
  q = dim(z)[2]

  #care if nbreak is matrix or scalar
  dv = matrix(0L, nrow = nbreak+2, ncol = 1)
  dv[1,1] = 0
  dv[seq(2,nbreak+1,1),1] = dateseq

  dv[nbreak+2,1] = bigT
  ds = matrix(0L,nrow = nbreak, ncol = 1)
  dr = matrix(0L,nrow = nbreak, ncol = 1)

  for (is in 1:nbreak){
    length = dv[is+2,1] - dv[is,1]
    if (p == 0){
      index = seq(dv[is,1]+1,dv[is+2,1],1)
      y_temp = y[index,1,drop=FALSE]
      z_temp = z[index,,drop=FALSE]
      vssr = ssr(1,y_temp,z_temp,h,length)
      y_temp_rev = rot90(rot90(y_temp))
      z_temp_rev = rot90(rot90(z_temp))
      vssrev = ssr(1,y_temp_rev,z_temp_rev,h,length)
      out = partione(h,length-h,length,vssr,vssrev)
      ds[is,1] = out$dx
      dr[is,1] = ds[is,1] + dv[is,1]
    }
    else{
      out = onebp(y,z,x,h,dv[is,1]+1,dv[is+2,1])
      ds[is,1] = out$bd
      dr[is,1] = ds[is,1]
    }

  }
  return(dr)
}


