# #!/usr/bin/env Rscript
# ## setwd("/n/data1/hsph/biostat/celehs/lab/xix636/survmaximin/code")
# # source("logistic_maximin.R")
# library(MASS)
# library(survival)
# library(glmnet)
# library(CVXR)
# library(survC1)
# set.seed(1)
# p=25
# nQ=600
# n.valid=2000
# nl=sapply(1:15,function(l){300*ceiling(l/3)})
# nl=c(nl,nQ)
# n = sum(nl)
# L=length(nl)
# ind.site=NULL; for (l in 1:(L-1)){ind.site=c(ind.site,rep(l,nl[l])) }
# signal = 'moderate'; rho = 'AR'; tau = .1
#
# data = generate_sim_data(p, L, signal, rho, tau, nl, n.valid)
# x = data$train$x; z = data$train$z; delta = data$train$delta
# x.valid = data$valid$x; z.valid = data$valid$z; delta.valid = data$valid$delta
#
# ################ Analysis ################
# ### start below, we should wrap these codes into functions
#
# ### Local Cox model training
# ### START MAXIMIN TIMER
# timestart=Sys.time()
#
# ## b_hat from all sites
# b_hat=matrix(0,p,L)
# for (l in 1:(L)){
#   outLas <- glmnet::cv.glmnet(z[[l]], survival::Surv(x[[l]],delta[[l]]), family='cox', alpha=1)
#   b_hat[,l]=as.vector(coef(outLas, s=outLas$lambda.min))
# }
#
#
# ######## auc.maximin.site
# l=L
# B_source = b_hat[,-l]
# Sigma_target = t(z[[l]])%*%z[[l]] / nrow(z[[l]])
#
# usethis::use_data(B_source, Sigma_target, overwrite = TRUE)
# usethis::use_data(x.valid, z.valid, delta.valid, overwrite = TRUE)
# saveRDS(data$valid.all, '~/program/maximin/valid_all.rds')
# saveRDS(data$valid, '~/program/maximin/valid.rds')
#
#
# ######## auc.maximin.site
# B_all = b_hat
# Sigma_all = c()
# for(l in 1:length(z)){
#   Sigma_all = c(Sigma_all, list(t(z[[l]])%*%z[[l]] / nrow(z[[l]])))
# }
#
# usethis::use_data(B_all, Sigma_all, overwrite = TRUE)
#
#
# output= survmaximin(B_source, Sigma_target, delta=0.5)
# beta.maximin<-output$beta.est
#
# auc.maximin=ROC.Est.FUN(Di=delta.valid,yyi=z.valid%*%beta.maximin,yy0=0.5)[1]
# c.maximin=survC1::Est.Cval(data.frame(cbind(x.valid,delta.valid,z.valid%*%beta.maximin)),tau=5)$Dhat
#
#
# timeend=Sys.time()
# time.maximin<-as.numeric(as.POSIXct(timeend)-as.POSIXct(timestart), units="secs")
# print(paste0("time.maximin: ",time.maximin))
