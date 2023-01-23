## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "")
library(mbreaks)

## ----estimate US rate---------------------------------------------------------
#the data set for the example is real.Rda
data(real)
#carry out all testing and estimating procedures via a single call
rate = mdl('rate',data=real,eps1=0.15)
#display the results from the procedures
rate

## ----graph_USr, fig.width=7---------------------------------------------------
#estimating the mean-shift model with BIC (the default option is ic=`'KT'`, which use modified BIC as criterion)
rate_mBIC = doorder('rate',data=real)
#NOTE: equivalent to rate$KT; type rate$KT to compare with new result


# visualization of estimated model with modified BIC (in the argument, we can replace rate$KT with rate_mBIC for exact same graph; recall that `$` is the operator to refer to field BIC in return list from mdl())
plot_model(rate$KT, title = 'US Exchange rate')

## ----reproduce_graph_USr, fig.width=7-----------------------------------------
  #collect model information
  m = rate_mBIC$nbreak           #number of breaks
  y = rate_mBIC$y                #vector of dependent var
  zreg = rate_mBIC$z             #matrix of regressors with changing coefs
  date = rate_mBIC$date          #estimated date
  fity = rate_mBIC$fitted.values #fitted values of model
  bigT = length(y)
  #compute the null model
  fixb = solve((t(zreg) %*% zreg)) %*% t(zreg) %*% y
  fity_fix = zreg%*%fixb    #fitted values of null model
  
  
  #plots the model
  tx = seq(1,bigT,1)
  range_y = max(y)-min(y);
  plot(tx,y,type='l',col="black", xlab='time',ylab="y", 
       ylim=c(min(y)-range_y/10,max(y)+range_y/10),lty=1)
  #plot fitted values series for break model
  lines(tx, fity,type='l', col="blue",lty=2)
  #plot fitted values series for null model
  lines(tx, fity_fix,type='l', col="dark red",lty=2)
  
  #plot estimated dates + CIs
  for (i in 1:m){
    abline(v=date[i,1],lty=2)
    if (rate_mBIC$CI[i,1] < 0){rate_mBIC$CI[i,1] = 0}
    if(rate_mBIC$CI[i,2]>bigT){ rate_mBIC$CI[i,2]=bigT}
    segments(rate_mBIC$CI[i,1],min(y)*(12+i/m)/10,rate_mBIC$CI[i,2],min(y)*(12+i/m)/10,lty=1,col='red')
  }
  
  legend(0,max(y)+range_y/10,legend=c("observed y",paste(m,'break y'),"0 break y"),
        lty=c(1,2,2), col=c("black","blue","red"), ncol=1)
  

## ----reproduce_table_PY-------------------------------------------------------
data(nkpc)
#x_t is GDP gap
  z_name = c('inflag','ygap','inffut')
  #we can invoke each test separately by using dotest() and doseqtests()
  supF_ygap = dotest('inf',z_name,data=nkpc,prewhit = 0, eps1 = 0.1,m=1)
  #z regressors' names are passed in the argument as an array, which equivalent to above argument call with z_name
  seqF_ygap = doseqtests('inf',c('inflag','ygap','inffut'),data=nkpc,prewhit = 0, eps1=0.1)
  #see test results
  supF_ygap
  seqF_ygap
  
  
#x_t is labor income share 
  #or invoke all tests using mdl() 
  nkpc_lbs = mdl('inf',c('inflag','lbs','inffut'),data=nkpc,prewhit = 0, eps1=0.1, m=5)
  nkpc_lbs$sbtests
  nkpc_lbs$seqtests
  

## ----re_estimate model given known breaks, echo = FALSE-----------------------
#only need to re-estimate model with output gap since we use mdl() for income share, we can obtain the estimated sequential model from SEQ (which is returned from mdl() as a list element)
# It is recommended to store desirable options as variables and set arguments = variables to avoid mistakes and save time
eps1 = 0.1
prewhit = 0
ygap_fixn = dofix('inf',z_name,data=nkpc,fixn=1,prewhit=prewhit,eps1=eps1)
#or use data-dependent sequential approach
ygap_SEQ = dosequa('inf',z_name,data=nkpc,prewhit=prewhit,eps1=eps1)


## ----index_date, include=FALSE------------------------------------------------
ygap_date = ygap_fixn$date

## ----reproduce sub-sample IV estimates, include=FALSE-------------------------
k=4
#list of instruments
instruments = c('inflag','lbslag','ygaplag','spreadlag','dwlag','dcplag')
#list of endogenous
regressors = c('inffut','inflag','lbs')

bigT = dim(nkpc)[1]
#independent variable
Y = as.matrix(nkpc[,'inf',drop=FALSE])

#form matrix of instruments
Z = as.matrix(nkpc[,instruments])
Z = cbind(rep(1,151),Z)
#endogenous variable
X_e = as.matrix(nkpc$inffut,drop=FALSE)
#first stage regression 
#X_res = (Z%*%solve(t(Z)%*%Z)%*%t(Z))%*%X_e
#2nd stage regressors
X = as.matrix(nkpc[,regressors])
X = cbind(rep(1,151),X)

#partition the regressors
T1 = seq(1,ygap_date)
T2 = seq(ygap_date+1,bigT)
Y1 = Y[T1,1,drop=FALSE]
Y2 = Y[T2,1,drop=FALSE]



#### R version #####
#multiplication difference
X_resR = (Z%*%solve(t(Z)%*%Z)%*%t(Z))%*%X_e

#2nd stage regressors
XR = as.matrix(nkpc[,regressors])
XR = cbind(rep(1,151),XR)
XhR = as.matrix(nkpc[,c('inflag','lbs')])
XhR = cbind(rep(1,151),XhR,X_resR)

Xh1R = XhR[T1,]
Xh2R = XhR[T2,]

#full sample estimate:
betaR = solve(t(XhR)%*%XhR)%*%t(XhR)%*%Y


#subsample estimates:
beta1R = solve(t(Xh1R)%*%Xh1R)%*%t(Xh1R)%*%Y1
beta2R = solve(t(Xh2R)%*%Xh2R)%*%t(Xh2R)%*%Y2

#compute variance
res1R = Y1 - Xh1R%*%beta1R
res2R = Y2 - Xh2R%*%beta2R
#no prewhitening to match paper
hac1R = mbreaks:::correct(Xh1R,res1R,0)
hac2R = mbreaks:::correct(Xh2R,res2R,0)
vhac1R = solve(t(Xh1R)%*%Xh1R)%*%hac1R%*%solve(t(Xh1R)%*%Xh1R) #4 regressors
vhac1R = vhac1R*(125-k)
vhac2R = solve(t(Xh2R)%*%Xh2R)%*%hac2R%*%solve(t(Xh2R)%*%Xh2R)
vhac2R = vhac2R*(bigT-125-k)
stdhac1R = sqrt(diag(vhac1R))
stdhac2R = sqrt(diag(vhac2R))


## ----display results with R results, echo=FALSE-------------------------------
colnames = c('$\\mu$','$\\gamma(\\pi_{t-1})$','$\\kappa(x_t)$','$\\beta(E_t\\pi_{t+1})$')
rownames = c('1960:Q1-1991:Q1','$SE_1$','1991:Q2-1997:Q4','$SE_2$')
IV_estimatesR = data.frame(round(t(beta1R),3))
IV_estimatesR = rbind(IV_estimatesR,round(stdhac1R,3))
IV_estimatesR = rbind(IV_estimatesR,round(t(beta2R),3))
IV_estimatesR = rbind(IV_estimatesR,round(stdhac2R,3))
colnames(IV_estimatesR) = colnames
rownames(IV_estimatesR) = rownames
knitr::kable(IV_estimatesR)

