#Main program to handle users' input

#' Comprehensive structural change estimation and testing
#'
#' `mdl()` calls main functions of the `mbreaks` package to execute the following
#' estimation procedures:
#' \describe{
#' \item{`dotest()`}{ Function \code{\link{dotest}} conducts Sup F tests of `0` versus `m` breaks and Double Max tests.}
#' \item{`doseqtests()`}{Function \code{\link{doseqtests}} conducts the sequential Sup F tests of `l` versus `l+1` breaks.}
#' \item{`doorder()`}{Function \code{\link{doorder}} conducts the number of breaks selection from `1` to `m`
#'  breaks using the following information critera: \code{KT},\code{BIC}, and \code{LWZ}.}
#' \item{`dosequa()`}{Function \code{\link{dosequa}} conducts the number of breaks selection by sequential tests
#'  from `1` to `m` breaks using sequential Sup F tests.}
#' \item{`dofix()`}{Function \code{\link{dofix}} conducts structural break model estimation with `fixn` breaks.}
#' }
#' All the procedures automatically identify if the `model` is either i) pure structural
#' breaks model or ii) partial structural breaks model
#'
#' @param y_name name of dependent variable in the data set.
#' @param z_name name of independent variables in the data set which coefficients are allowed to change
#' across regimes. \code{default} is vector of 1 (Mean-shift model).
#' @param x_name name of independent variables in the data set which coefficients are constant across
#' regimes. \code{default} is \code{NULL}.
#' @param data the data set for estimation.
#' @param const indicates whether the regression model include an
#'intercept changing across regimes. Default value is 1.
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
#'   \item \code{eps1 = 0} This option allows users to explicitly specify
#' minimum segment length `h` parameters. However, this option will not
#' be allowed for testing and testing related functions.}
#' The default value is set at \code{eps1 = 0.15}.
#' @param m Maximum number of structural changes allowed. If not specify,
#' m will be set to \code{default} value matching `eps1` input.
#' @param prewhit set to \code{1} to apply AR(1) prewhitening prior to estimating
#' the long run covariance matrix.
#' @param robust set to \code{1} to allow for heterogeneity
#' and autocorrelation in the residuals, \code{0} otherwise.
#' The method used is Andrews(1991) automatic bandwidth with AR(1) approximation with quadratic
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
#' assumed identical across segments (the variance of the errors u if \code{robust}=\code{0})
#' @param hetq used in the construction of the confidence intervals for the break
#' dates. If \code{hetq}=\code{0}, the moment matrix of the data is assumed identical
#' across segments
#' @param maxi number of maximum iterations for recursive calculations of finding
#' global minimizers.\code{default} = 10 (For partial change model ONLY).
#' @param eps convergence criterion for recursive calculations (For partial change model ONLY)
#' @param fixn number of pre-specified breaks. \code{default} = -1. It will be replaced
#' automatically to 2 if no specification is given (For partial change model ONLY)
#' @param printd Print option for model estimation. \code{default} = 0, to
#' suppress intermediate outputs printing to console
#'@param signif significance level used to sequential test to select number of breaks.
#'\itemize{
#'   \item 4: 1\% level
#'   \item 3: 2.5\% level
#'   \item 2: 5\% level
#'   \item 1: 10\% level
#' }
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{\beta_0} to use in estimation (Must be a `p x 1` matrix, where `p` is number of x variables)
#'@param h Minimum segment length of regime considered in estimation. If users want to specify a particular value, please set `eps1=0`
#' @return  A list that contains the following:
#'\describe{
#'\item{sbtests}{A list of class `sbtests` representing Sup F tests of 0 versus m breaks and Double Max tests.}
#'\item{seqtests}{A list of class `seqtests` representing sequential Sup F test of l versus l+1 breaks.}
#'\item{BIC}{A list of class `model` with structural break model estimated by number of breaks by \code{BIC} criterion.}
#'\item{LWZ}{A list of class `model` with structural break model estimated by number of breaks by \code{LWZ} criterion.}
#'\item{KT}{A class `model` with structural break model estimated by number of breaks by \code{KT} criterion.}
#'\item{sequa}{A class `model` with structural break model estimated by number of breaks by sequential tests.}
#'\item{fix}{A class `model` with structural break model estimated by pre-specified `fixn` number of breaks.}
#'}
#'
#' Note: All \code{default} values of error assumptions (\code{robust},
#' \code{hetdat}, \code{hetvar}, \code{hetq}) are set to 1. The implications on
#' the structure of model\'s errors related to individual settings are explained within
#' the arguments section for each option.
#'
#' @seealso \code{\link{dotest}}, \code{\link{doseqtests}}, \code{\link{doorder}}, \code{\link{dosequa}}, and \code{\link{dofix}}
#' which are functions called by `mdl()`.
#' @examples
#'  US_rate = mdl('rate',data=real)
#'  nkpc_lbs = mdl('inf',c('inflag','lbs','inffut'),data=nkpc,prewhit = 0)
#'
#' @export

mdl <- function(y_name,z_name = NULL,x_name = NULL,data,eps1 = 0.15,m = 5,prewhit = 1,
                           robust = 1,hetdat = 1,hetvar = 1,hetomega = 1,hetq = 1,
    maxi = 10,eps = 0.00001,fixn=-1,fixb=0,betaini=0,printd = 0,const=1,signif=2,h=NULL){

  if(printd==1){
  cat('The options chosen are:\n')
  cat(paste(' i) hetdat = ',hetdat),'\n',paste('ii) hetvar = ',hetvar),'\n',
      paste('iii) hetomega = ',hetomega),'\n',paste('iv) hetq = ',hetq),'\n',
      paste('v) robust = ',robust),'\n',paste('vi) prewhite = ',prewhit),'\n')
  }

  mdl = list()

  #set significant level
  siglev=matrix(c(10,5,2.5,1),4,1)



  if (fixn == -1){
    fixn = 2
  }


  #Sup F test
  if (eps1!=0){
  mdl$sbtests = dotest(y_name=y_name,z_name=z_name,x_name=x_name,data=data,m=m,
                     eps=eps,eps1=eps1,maxi=maxi,fixb=fixb,
                     betaini=betaini,printd=printd,prewhit=prewhit,robust=robust,
                     hetdat=hetdat,hetvar=hetvar,const=const)
  #Sequential test
  mdl$seqtests = doseqtests(y_name=y_name,z_name=z_name,x_name=x_name,data=data,m=m,
                        eps=eps,eps1=eps1,maxi=maxi,fixb=fixb,
                        betaini=betaini,printd=printd,prewhit=prewhit,robust=robust,
                        hetdat=hetdat,hetvar=hetvar,const=const)
  #Sequential procedure to select num of breaks and estimate break dates
  mdl$SEQ = dosequa(y_name=y_name,z_name=z_name,x_name=x_name,data=data,
                    m=m,eps=eps,eps1=eps1,maxi=maxi,fixb=fixb,betaini=betaini,
                    printd=printd,prewhit=prewhit,robust=robust,
                    hetq=hetq,hetomega=hetomega,hetdat=hetdat,
                    hetvar=hetvar,const=const,signif=signif)}
  else{ #when eps1 set to 0, no test is included
    mdl$sbtests = NULL
    mdl$seqtests = NULL
    mdl$SEQ = NULL
  }


  #Information criteria to select num of breaks and estimate break dates
  mdl$BIC = doorder(y_name=y_name,z_name=z_name,x_name = x_name,data=data,
                    m=m,eps=eps,eps1=eps1,maxi=maxi,fixb=fixb,
                    betaini=betaini,printd=printd,ic = 'BIC',const=const,h=h,
                    prewhit=prewhit,robust=robust,
                    hetdat=hetdat,hetvar=hetvar,hetomega=hetomega,hetq=hetq)
  mdl$LWZ = doorder(y_name=y_name,z_name=z_name,x_name = x_name,data=data,
                    m=m,eps=eps,eps1=eps1,maxi=maxi,fixb=fixb,
                    betaini=betaini,printd=printd,ic = 'LWZ',const=const,h=h,
                    prewhit=prewhit,robust=robust,
                    hetdat=hetdat,hetvar=hetvar,hetomega=hetomega,hetq=hetq)
  mdl$KT = doorder(y_name=y_name,z_name=z_name,x_name = x_name,data=data,
                    m=m,eps=eps,eps1=eps1,maxi=maxi,fixb=fixb,
                    betaini=betaini,printd=printd,ic = 'KT',const=const,h=h,
                   prewhit=prewhit,robust=robust,
                   hetdat=hetdat,hetvar=hetvar,hetomega=hetomega,hetq=hetq)



  #Repartition procedure to select num of breaks and estimate break dates
  #mdl$repart = dorepart(y_name=y_name,z_name=z_name,x_name=x_name,data=data,
  #                      m=m,eps=eps,eps1=eps1,maxi=maxi,fixb=fixb,betaini=betaini,
  #                      printd=printd,prewhit=prewhit,robust=robust,hetdat=hetdat,hetvar=hetvar,const=const,signif=signif)

  #mdl$fix = dofix(y_name=y_name,z_name=z_name,x_name=x_name,data=data,
  #                fixn=fixn,eps=eps,eps1=eps1,maxi=maxi,fixb=fixb,betaini=betaini,
  #                printd=printd,prewhit=prewhit,robust=robust,hetdat=hetdat,hetvar=hetvar,const=const)
  #reorganize the results into the list
  class(mdl) <- 'mdl'
  return(mdl)
  }


#'Print information of `mbreaks` comprehensive procedure
#'
#'`print` prints the class `mdl` object with default showing only
#'certain procedures called
#'by `mdl()` function including: `seqtests` class object, `sbtests` class
#'object, and `model` class object using KT information criterion
#'
#'@param x class `mdl` object
#'@param ... further arguments passed to or from other methods
#'
#'@examples
#' rate = mdl('rate',data=real)
#' print(rate)
#'@return
#'No return value, only for printing `model`, `sbtests` and `seqtests` class objects
#'invoked during `mdl()`.
#'
#'@export
print.mdl <- function(x,...)
{
  digits = max(3L, getOption("digits") - 3L)
  cat('-----------------------------------------------------')
  if (is.null(x$KT)) {cat('\nNo breaks were founded by KT\n')}
  else{
  print(x$KT)}
  cat('-----------------------------------------------------')
  if(is.null(x$sbtests)){cat('\nTrimming level sets to 0, no supF tests included')}
  else{
  print(x$sbtests)}
  cat('-----------------------------------------------------')
  if(is.null(x$seqtests)){cat('\nTrimming level sets to 0, no sequential F tests included')}
  else{
  print(x$seqtests)}
  cat('-----------------------------------------------------')

  cat(paste('\nTo access additional information about specific procedures
  (not included above), type stored variable name + \'$\' + procedure name'))
  invisible(x)
}






