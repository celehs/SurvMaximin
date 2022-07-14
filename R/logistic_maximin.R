##########################################################################
###### Generate beta_star
# Arguments
# - B: coefficients of L groups, dim (p, L)
# - Sigma: covariates matrix, dim (p, p)
# - delta: default 0; the ridge penalty parameter;
### if B^{T}\Sigma B is nearly singular, we may set delta to be non-zero, like 0.5.
##########################################################################
#' Generate transfer-learning survmaximin coefficients
#'
#' @param B_source A p (number of variables) * L (number of source sites) numeric matrix.
#' @param Sigma_target Data covariance matrix for the target site.
#' @param delta Default 0; the ridge penalty parameter.
#'
#' @return A list with the following components:
#'
#' \tabular{ll}{
#'  \code{beta.est} \tab Estimated beta coefficients for the target site. \cr
#'  \code{weight} \tab A vector containing trained weights for each source site. \cr
#' }
#' @export
#'
#' @examples
#' data(B_source); data(Sigma_target)
#' output <- survmaximin(B_source, Sigma_target, delta=0.5)
#' beta.maximin <- output$beta.est
survmaximin <- function(B_source, Sigma_target, delta=0){
  L <- dim(B_source)[2]
  Gamma <- t(B_source)%*%Sigma_target%*%B_source
  opt.weight <- rep(NA, L)
  # Problem Definition
  v <- CVXR::Variable(L)
  Diag.matrix <- diag(eigen(Gamma)$values)
  for (ind in 1:L){
    Diag.matrix[ind, ind] <- max(Diag.matrix[ind, ind], 0.001)
  }
  Gamma.positive <- eigen(Gamma)$vectors %*% Diag.matrix %*% t(eigen(Gamma)$vectors)
  objective <- CVXR::Minimize(CVXR::quad_form(v, Gamma.positive + diag(delta, L)))
  constraints <- list(v>=0, sum(v)==1)
  prob.weight <- CVXR::Problem(objective, constraints)
  if (CVXR::is_dcp(prob.weight)){
    # Problem Solution
    result <- CVXR::solve(prob.weight)
    opt.status <- result$status
    opt.sol <- result$getValue(v)
  }
  for (l in 1:L){
    opt.weight[l] <- opt.sol[l]*(abs(opt.sol[l])>10^{-8})
  }
  return(list(beta.est=B_source%*%opt.weight,weight=opt.weight))
}



#' Generate federate-learning coefficients
#'
#' @param B_all A p (number of variables) * L (number of sites) numeric matrix.
#' @param Sigma_all A list containing the data covariance matrix for each site.
#' @param delta Default 0; the ridge penalty parameter.
#'
#' @return A list with each element as below, being the estimated results for each site:
#'
#' \tabular{ll}{
#'  \code{beta.est} \tab Estimated beta coefficients for the target site. \cr
#'  \code{weight} \tab A vector containing trained weights for each source site. \cr
#' }
#' @export
#'
#' @examples
#' data(B_all); data(Sigma_all)
#' output <- survmaximin_fed(B_all, Sigma_all, delta=0.5)
#' beta.maximin <- output$beta.est
survmaximin_fed <- function(B_all, Sigma_all, delta=0){
  L <- dim(B_all)[2] - 1
  n.site = length(Sigma_all)
  out = c()
  for(i in 1:n.site){
    Sigma_target = Sigma_all[[i]]
    B_source = B_all[, -i]
    out_part = survmaximin(B_source, Sigma_target, delta=delta)
    out = c(out, list(out_part))
  }
  return(out)
}






#### Accuracy measures
#### yyi is risk score
#### Di event status
ROC.Est.FUN=function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F)
{
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0
  for(k in 1:pp)
  {
    yy = yy0;
    # if(!is.null(fpr0)){
    #   tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
    #   fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
    #   TPR = stats::approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y;
    #   TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR);
    #   yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))
    #   FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    # }else{
    #   TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
    #   FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    # }
    TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
    FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)


    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    #AUC <- sum((sum_I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum_I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
    #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
    AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
  out
}

S.FUN=function(yy,Yi,Di,yes.smooth=F)
{
  # if(yes.smooth){
  #   Y1i = Yi[Di==1]; n1 = sum(Di); bw = stats::bw.nrd(Y1i)/n1^0.6
  #   c(t(rep(1/n1,n1))%*%stats::pnorm((Y1i-VTM(yy,n1))/bw))
  # }else{
  #   return((sum_I(yy,"<",Yi,Vi=Di)+sum_I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  # }
  return((sum_I(yy,"<",Yi,Vi=Di)+sum_I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)



  ##sum_I(yy,"<=",Yi,Vi=Di)/sum(Di)
}

sum_I <- function(yy,FUN,Yi,Vi=NULL)
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}


######## Generate (X=min(T, C), delta=(T<=C)) from cox models ##########
# Z=z[[1]]; b=b[,1]; a=1
generate_survivaldata <- function(Z, b, a){
  n = dim(Z)[1]
  t=sqrt( (-log(1-stats::runif(n, 0, 1))/(exp(Z%*%b)))*8/a )
  c=stats::rexp(n, rate = 1)
  x=apply(cbind(t,c),1,min)
  delta=as.numeric(t<= c)

  return(list(x,delta))
}


ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  rho^exponent
}

generate_sim_data <- function(p, L, signal, rho, tau, nl, n.valid){
  #### coefficients b: same for the same setting
  set.seed(1234)
  b=matrix(NA,p,L)
  n = length(nl)
  if (signal=='moderate'){
    b[,L] = c(0.5,0.4,0.3,0.2,-0.2,-0.3,-0.4,-0.5,rep(0,p-8))
  } else if(signal=="low"){
    b[,L] = c(0.5,0.4,0.3,0.2,-0.2,-0.3,-0.4,-0.5,rep(0,p-8))/2
  }
  for (l in 1:(L-1)){
    # b[,l]=b[,L]+c(rnorm(8,mean=0,stats::sd=tau),rep(0,p-8))
    aa=tau*(3*(l>5)+(l<=5))
    b[,l]=b[,L]+c(-aa,aa,-aa,aa,-aa,aa,-aa,aa,rep(0,p-8))
  }
  # round(b[1:10,],2)

  #### Generate simulated dataset Z from multivariate normal: rho, sig
  set.seed(as.numeric(Sys.time()))


  if (rho=='compound'){
    Sigma=matrix(0.5,nrow=p,ncol=p)+diag(0.5,nrow=p,ncol=p);Sigma_Q=Sigma+0.1
  } else if (rho=='AR'){
    Sigma=ar1_cor(p, 0.5); # rho=0.3, 0.7 // correlation of Z
    Sigma_Q=Sigma+0.1 }

  z = z.valid.all = NULL
  for (l in 1:(L-1)){
    aa=MASS::mvrnorm(nl[l],rep(0,p), Sigma)
    aa[,3]=(aa[,3]>stats::quantile(aa[,3],.3)); aa[,3]=(aa[,3]-mean(aa[,3]))/stats::sd(aa[,3])
    aa[,2]=(aa[,2]>stats::quantile(aa[,2],.6)); aa[,2]=(aa[,2]-mean(aa[,2]))/stats::sd(aa[,2])
    z[[l]]=aa

    aa=MASS::mvrnorm(n.valid,rep(0,p), Sigma)
    aa[,3]=(aa[,3]>stats::quantile(aa[,3],.3)); aa[,3]=(aa[,3]-mean(aa[,3]))/stats::sd(aa[,3])
    aa[,2]=(aa[,2]>stats::quantile(aa[,2],.6)); aa[,2]=(aa[,2]-mean(aa[,2]))/stats::sd(aa[,2])

    z.valid.all[[l]]=aa
  }
  l=L
  aa=MASS::mvrnorm(nl[l],rep(0,p), Sigma_Q)
  aa[,3]=(aa[,3]>stats::quantile(aa[,3],.3)); aa[,3]=(aa[,3]-mean(aa[,3]))/stats::sd(aa[,3])
  aa[,2]=(aa[,2]>stats::quantile(aa[,2],.6)); aa[,2]=(aa[,2]-mean(aa[,2]))/stats::sd(aa[,2])
  z[[l]]=(aa-apply(aa,2,mean))/apply(aa,2,stats::sd)

  aa=MASS::mvrnorm(n.valid,rep(0,p), Sigma_Q)
  aa[,3]=(aa[,3]>stats::quantile(aa[,3],.3)); aa[,3]=(aa[,3]-mean(aa[,3]))/stats::sd(aa[,3])
  aa[,2]=(aa[,2]>stats::quantile(aa[,2],.6)); aa[,2]=(aa[,2]-mean(aa[,2]))/stats::sd(aa[,2])
  z.valid.all[[l]]=(aa-apply(aa,2,mean))/apply(aa,2,stats::sd)


  aa=MASS::mvrnorm(n.valid,rep(0,p), Sigma_Q)
  aa[,3]=(aa[,3]>stats::quantile(aa[,3],.3)); aa[,3]=(aa[,3]-mean(aa[,3]))/stats::sd(aa[,3])
  aa[,2]=(aa[,2]>stats::quantile(aa[,2],.6)); aa[,2]=(aa[,2]-mean(aa[,2]))/stats::sd(aa[,2])
  z.valid=(aa-apply(aa,2,mean))/apply(aa,2,stats::sd)

  x=NULL; delta=NULL; x.valid.all = delta.valid.all = NULL
  for (l in 1:L){
    temp=generate_survivaldata(z[[l]], b[,l], a=1+0.05*l)
    x[[l]]=temp[[1]]; delta[[l]]=temp[[2]]
    temp_valid = generate_survivaldata(z.valid.all[[l]], b[,l], a=1+0.05*l)
    x.valid.all[[l]]=temp_valid[[1]]; delta.valid.all[[l]]=temp_valid[[2]]
  }
  l=L
  temp=generate_survivaldata(z.valid, b[,l], a=1+0.05*l)
  x.valid=temp[[1]]; delta.valid=temp[[2]]
  return(list(`train` = list(`x` = x, `z` = z, `delta` = delta),
         `valid` = list(`x` = x.valid, `delta` = delta.valid,
                        `z` = z.valid),
         `valid.all` = list(`x` = x.valid.all, `z` = z.valid.all,
                            `delta` = delta.valid.all)))
}
