#Testing procedures

#'SupF test for 0 vs i breaks
#'
#'Function compute the supF test statistics of testing procedure with
#' null hypothesis: no break versus alternative hypothesis: \code{i} breaks.
#'
#'@param y dependent variables
#'@param z independent variables with coefficients are allowed to change across
#'regimes
#'@param q number of \code{z} regressors
#'@param x independent variables with constant coefficients across regimes
#'@param p number of \code{x} regressors
#'@param i number of breaks in the model
#'@param bigT sample period T
#'@param datevec \code{i} estimated dates from the model
#'@param robust,hetdat,hetvar options for assumptions on error terms
#'@param prewhit Options for prewhitening process
#'@return ftest SupF test results
#'@export

pftest = function(y,z,i,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar){
  #construct matrix R
  rsub = matrix(0L, nrow = i , ncol = i+1)
  datevec = as.matrix(datevec)
  j = 1
  while(j<=i){
    rsub[j,j] = -1
    rsub[j,j+1] = 1
    j=j+1
  }
  rmat = kron(rsub,diag(1,q))
  date = datevec[seq(1,i,1),i,drop=FALSE]
  zbar = diag_par(z,i,date)

  if (p==0){
    delta = OLS(y,zbar)
  }
  else {
    dbdel = OLS(y, cbind(zbar,x))
    delta = dbdel[seq(1,(i+1)*q) , 1,drop=FALSE]
  }

  vdel = pvdel(y,z,i,q,bigT,date,prewhit,robust,x,p,0,hetdat,hetvar)
  fstar = t(delta) %*% t(rmat) %*% solve(rmat %*% vdel %*% t(rmat)) %*%
    rmat %*% delta
  ftest = (bigT - (i+1)*q - p) %*% fstar / (bigT*i)

  return(ftest)
}

#' SupF(l+1|l) test
#'
#'Function computes the test statistics of supF(l+1|l) test with null hypothesis
#'is l=\code{nseg}-1 and alternative hypothesis is l+1.
#'The l breaks under the null hypothesis are taken from the global minimization.
#'
#'@param bigvec associated SSR of estimated break date under H0
#'@param dt estimated date under H0
#'@param nseg number of segment under H1
#'@param z independent variables with coefficients are allowed to change across
#'regimes
#'@param q number of \code{z} regressors
#'@param x independent variables with constant coefficients across regimes
#'@param p number of \code{x} regressors
#'@param prewhit,robust,hetdat,hetvar options on residuals/errors
#'@return A list that contains the following:
#'\itemize{
#'\item {maxf}{Maximum value of test}
#'\item{newd}{Additional date in alternative hypothesis }
#'}
#'@export
spflp1 = function(bigvec,dt,nseg,y,z,h,q,prewhit,robust,x,p,hetdat,hetvar){
  #
  ssr = matrix(0L,nrow = nseg, ncol = 1)
  ftestv = matrix(0L, nrow = nseg, ncol = 1)
  bigT = dim(z)[1]
  dv = matrix(0L,nrow = nseg+1, ncol = 1)
  dv[2:nseg,1] = dt
  dv[nseg+1,1] = bigT
  ds = matrix(0L, nrow = nseg, ncol = 1)

  i_n = 0
  for (is in 1:nseg){
    length = dv[is+1,1] - dv[is,1]

    if(length >= 2*h){
      if (p == 0){
        out = parti(dv[is,1]+1,dv[is,1]+h,dv[is+1,1]-h,dv[is+1,1],bigvec,bigT)
        ssr[is,1] = out$ssrmin
        ds[is,1] = out$dx
        y_test = y[seq(dv[is,1]+1,dv[is+1,1],1),1,drop=FALSE]
        z_test = z[seq(dv[is,1]+1,dv[is+1,1],1),,drop=FALSE]
        ftestv[is,1] = pftest(y_test,z_test,1,q,length,ds[is,1,drop=FALSE]-dv[is,1,drop=FALSE],
                              prewhit,robust,0,p,hetdat,hetvar)
      }
      else{
        out = onebp(y,z,x,h,dv[is,1]+1,dv[is+1,1])
        ssr[is,1] = out$ssrind
        ds[is,1] = out$bd
        y_test = y[seq(dv[is,1]+1,dv[is+1,1],1),1,drop=FALSE]
        z_test = z[seq(dv[is,1]+1,dv[is+1,1],1),,drop=FALSE]
        x_test = x[seq(dv[is,1]+1,dv[is+1,1],1),,drop=FALSE]
        ftestv[is,1] = pftest(y_test,z_test,1,q,length,ds[is,1,drop=FALSE]-dv[is,1,drop=FALSE],
                              prewhit,robust,x_test,p,hetdat,hetvar)
      }
    }
    else {
      i_n = i_n+1
      ftestv[is,1] = 0.0
    }
  }

  if (i_n == nseg) {
    #print(paste('Given the location of the breaks from the global optimization with',
    #            nseg,'breaks there was no more place to insert an additional breaks that satisfy the minimal length requirement.'))
  }

  maxf = max(ftestv[1:nseg,1])
  newd = ds[which.max(ftestv[1:nseg,1]),1]
  if (newd == 0){newd = NA}
  out = list('maxf' = maxf, 'newd' = newd)
  return(out)
}

