# #!/usr/bin/env Rscript
# p=25
# nQ=600
# n.valid=2000
# nl=sapply(1:15,function(l){300*ceiling(l/3)})
# nl=c(nl,nQ)
# n = sum(nl)
#
# source("/n/data1/hsph/biostat/celehs/lab/haz364/maximin_simulations/code/logistic_maximin.R")
#
# nQ=600
# n.valid=2000
# nl=sapply(1:15,function(l){300*ceiling(l/3)})
# nl=c(nl,nQ)
# n = sum(nl)
# L=length(nl)
# ind.site=NULL; for (l in 1:(L-1)){ind.site=c(ind.site,rep(l,nl[l])) }
# signal = 'moderate'; rho = 'AR'; tau = .1
#
# ## replication
# simu=500
# for (i in start:end){
#
# for (rho in c('compound','AR')){
#   # rho='AR' # // correlation of Z #
#
# for (signal in c('moderate','low')){
#   # signal='moderate' # // correlation of Z #
#
# for (tau in c(0.05,.1,.2)){
#   # tau=0.2 # // perturbation of beta
#
# #### coefficients b: same for the same setting
# set.seed(1234)
# b=matrix(NA,p,L)
# if (signal=='moderate'){
#   b[,L] = c(0.5,0.4,0.3,0.2,-0.2,-0.3,-0.4,-0.5,rep(0,p-8))
#   } else if(signal=="low"){
#     b[,L] = c(0.5,0.4,0.3,0.2,-0.2,-0.3,-0.4,-0.5,rep(0,p-8))/2
#     }
# for (l in 1:(L-1)){
#   # b[,l]=b[,L]+c(rnorm(8,mean=0,sd=tau),rep(0,p-8))
#   aa=tau*(3*(l>5)+(l<=5))
#   b[,l]=b[,L]+c(-aa,aa,-aa,aa,-aa,aa,-aa,aa,rep(0,p-8))
#   }
# # round(b[1:10,],2)
#
#
# #### Generate simulated dataset Z from multivariate normal: rho, sig
# set.seed(as.numeric(Sys.time()))
#
# ar1_cor <- function(n, rho) {
#   exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
#                     (1:n - 1))
#   rho^exponent
# }
#
# if (rho=='compound'){
#   Sigma=matrix(0.5,nrow=p,ncol=p)+diag(0.5,nrow=p,ncol=p);Sigma_Q=Sigma+0.1
# } else if (rho=='AR'){
#   Sigma=ar1_cor(p, 0.5); # rho=0.3, 0.7 // correlation of Z
#   Sigma_Q=Sigma+0.1 }
#
# z=NULL
# for (l in 1:(L-1)){
#   aa=mvrnorm(nl[l],rep(0,p), Sigma)
#   aa[,3]=(aa[,3]>quantile(aa[,3],.3)); aa[,3]=(aa[,3]-mean(aa[,3]))/sd(aa[,3])
#   aa[,2]=(aa[,2]>quantile(aa[,2],.6)); aa[,2]=(aa[,2]-mean(aa[,2]))/sd(aa[,2])
#   z[[l]]=aa
#   # aa=mvlognormal(nlall[l], Mu = rep(1,p), Sigma = rep(scal[l], p), R = Sigma)
#   # aa[,1]=(aa[,1]>quantile(aa[,1],.3))
#   # aa[,2]=(aa[,2]>quantile(aa[,2],.6))
#   # Z[[l]]=(aa-apply(aa,2,mean))/apply(aa,2,sd)
# }
# l=L
# aa=mvrnorm(nl[l],rep(0,p), Sigma_Q)
# aa[,3]=(aa[,3]>quantile(aa[,3],.3)); aa[,3]=(aa[,3]-mean(aa[,3]))/sd(aa[,3])
# aa[,2]=(aa[,2]>quantile(aa[,2],.6)); aa[,2]=(aa[,2]-mean(aa[,2]))/sd(aa[,2])
# z[[l]]=(aa-apply(aa,2,mean))/apply(aa,2,sd)
# aa=mvrnorm(n.valid,rep(0,p), Sigma_Q)
# aa[,3]=(aa[,3]>quantile(aa[,3],.3)); aa[,3]=(aa[,3]-mean(aa[,3]))/sd(aa[,3])
# aa[,2]=(aa[,2]>quantile(aa[,2],.6)); aa[,2]=(aa[,2]-mean(aa[,2]))/sd(aa[,2])
# z.valid=(aa-apply(aa,2,mean))/apply(aa,2,sd)
#
# ######## Generate (X=min(T, C), delta=(T<=C)) from cox models ##########
# # Z=z[[1]]; b=b[,1]; a=1
# generate_survivaldata <- function(Z, b, a){
#   n = dim(Z)[1]
#   t=sqrt( (-log(1-runif(n, 0, 1))/(exp(Z%*%b)))*8/a )
#   c=rexp(n, rate = 1)
#   x=apply(cbind(t,c),1,min)
#   delta=as.numeric(t<= c)
#
#   return(list(x,delta))
# }
#
# x=NULL; delta=NULL
# for (l in 1:L){
#   temp=generate_survivaldata(z[[l]], b[,l], a=1+0.05*l)
#   x[[l]]=temp[[1]]; delta[[l]]=temp[[2]]
# }
# print( round(sapply(1:L,function(l){mean(delta[[l]])}),2) )
# print(sapply(1:L,function(l){sum(delta[[l]])}) )
#
# l=L
# temp=generate_survivaldata(z.valid, b[,l], a=1+0.05*l)
# x.valid=temp[[1]]; delta.valid=temp[[2]]
#
#
# ################ Analysis ################
# ######## pool of L-1 sites
# zpool=NULL;xpool=NULL; deltapool=NULL
# for (l in 1:(L-1)){
#   zpool=rbind(zpool,z[[l]])
#   xpool=c(xpool,x[[l]])
#   deltapool=c(deltapool,delta[[l]])
# }
#
# # missing
# z[[L]]=z[[L]][,2:p]
# z.valid=z.valid[,2:p]
#
# timestart=Sys.time()
#
# ## b from all sites
# b_hat=matrix(0,p,L)
# for (l in 1:(L-1)){
#   outLas <- cv.glmnet(z[[l]], Surv(x[[l]],delta[[l]]), family='cox', alpha=1)
#   b_hat[,l]=as.vector(coef(outLas, s=outLas$lambda.min))
#   # fitcox <-tryCatch(coxph(Surv(X[[l]],delta[[l]])~z[[l]]), error = function(e) list(coef=NA))
#   # b_hat[,l]=as.vector(fitcox$coefficients)
# }
# # round(b_hat[1:10,],2)
# l=L
# outLas <- cv.glmnet(z[[l]][1:200,], Surv(x[[l]][1:200],delta[[l]][1:200]), family='cox', alpha=1)
# b_hat[,l]=c(NA,as.vector(coef(outLas, s=outLas$lambda.min)))
#
#
# #### projection
# betatemp=b_hat[,L]
# ind.na=which(is.na(betatemp))
#
# beta.impute=matrix(0,nrow(b_hat)-length(ind.na),L-1)
# for (l in 1:(L-1)){
#   temp=ginv(t(z[[l]][,-ind.na]) %*% z[[l]][,-ind.na]) %*% (t(z[[l]][,-ind.na])%*%z[[l]][,ind.na])
#   beta.impute[,l]=b_hat[-ind.na,l]+temp*b_hat[ind.na,l]
# }
# # round(beta.impute[1:10,],2)
#
# B = beta.impute
# Sigma = t(z[[L]]) %*% z[[L]]/nrow(z[[L]])
# output=beta_star(B, Sigma, delta=0.5)
# beta.maximin<-output$beta.est
#
# c.maximin=Est.Cval(data.frame(cbind(x.valid,delta.valid,z.valid%*%beta.maximin)),tau=5)$Dhat
#
#
# ######## maximin p-1 beta, no projection
# ## b from all sites
# b_hat_1=matrix(0,p-1,L)
# for (l in 1:(L-1)){
#   outLas <- cv.glmnet(z[[l]][,2:p], Surv(x[[l]],delta[[l]]), family='cox', alpha=1)
#   b_hat_1[,l]=as.vector(coef(outLas, s=outLas$lambda.min))
#   # fitcox <-tryCatch(coxph(Surv(X[[l]],delta[[l]])~z[[l]]), error = function(e) list(coef=NA))
#   # b_hat[,l]=as.vector(fitcox$coefficients)
# }
# l=L
# outLas <- cv.glmnet(z[[l]][1:200,], Surv(x[[l]][1:200],delta[[l]][1:200]), family='cox', alpha=1)
# b_hat_1[,l]=as.vector(coef(outLas, s=outLas$lambda.min))
# # round(b_hat_1[1:19,],2)
#
# B = b_hat_1[,-L]
# Sigma = t(z[[L]]) %*% z[[L]]/nrow(z[[L]])
# output=beta_star(B, Sigma, delta=0.5)
# beta.maximin.less<-output$beta.est
# # round(cbind(beta.maximin,beta.maximin.less),2)
#
# c.maximin.less=Est.Cval(data.frame(cbind(x.valid,delta.valid,z.valid%*%beta.maximin.less)),tau=5)$Dhat
#
#
# cor(z[[l]]%*%beta.maximin,z[[l]]%*%beta.maximin.less)
# corr=cor.test(z[[l]]%*%beta.maximin, z[[l]]%*%beta.maximin.less, method = 'spearman')$estimate
#
# timeend=Sys.time()
# time.maximin<-as.numeric(as.POSIXct(timeend)-as.POSIXct(timestart), units="secs")
# print(paste0("time.maximin: ",time.maximin))
#
# ################ Distributed learning
# ######## meta analysis
# timestart=Sys.time()
#
#
# beta_meta=rep(0,p); wei=rep(0,p)
# for (l in 1:(L-1)){
#   fitcox <-tryCatch(coxph(Surv(x[[l]],delta[[l]])~z[[l]]), error = function(e) list(coef=NA))
#   if( !is.null(fitcox$coefficients) & !any(is.na(fitcox$coef)) ){
#     wt <- 1/( ( summary(fitcox)$coef[,3] )^2 )
#     beta_meta=beta_meta+fitcox$coef * wt
#     wei=wei+wt
#   }
# }
# beta_meta=beta_meta/wei # round(as.numeric(beta_meta),2)
# l=L
# # auc.meta=ROC.Est.FUN(Di=delta.valid,yyi=z.valid%*%beta_meta[-1],yy0=0.5)[1]
# c.meta=Est.Cval(data.frame(cbind(x.valid,delta.valid,z.valid%*%beta_meta[-1])),tau=5)$Dhat
#
# timeend=Sys.time()
# time.meta<-as.numeric(as.POSIXct(timeend)-as.POSIXct(timestart), units="secs")
# print(paste0("time.meta: ",time.meta))
#
# ######## Rui
# timestart=Sys.time()
#
# dat=data.frame('site'=ind.site,'time'=xpool,'status'=deltapool,zpool)
#
# fit.ODAC <-tryCatch(DistCox(mydata = dat, id.local = 1,init_est = beta_meta,strat = T,output.ODACO1 = T),
#                     error = function(e) list(coef=NA))
# if( !any(is.na(fit.ODAC$coef)) ){
#   l=L
#   # auc.dist=ROC.Est.FUN(Di=delta.valid,yyi=z.valid%*%fit.ODAC$beta_tilde,yy0=0.5)[1]
#   # auc.dist1=ROC.Est.FUN(Di=delta[[l]],yyi=z[[l]]%*%fit.ODAC$beta_tilde1,yy0=0.5)[1]
#   c.dist=Est.Cval(data.frame(cbind(x.valid,delta.valid,z.valid%*%fit.ODAC$beta_tilde[-1])),tau=5)$Dhat
# }else{
#   fit.ODAC <-tryCatch(DistCox(mydata = dat, id.local = 1,init_est = 'local',strat = T,output.ODACO1 = T),
#                       error = function(e) list(coef=NA))
#   l=L
#   # auc.dist=ROC.Est.FUN(Di=delta.valid,yyi=z.valid%*%fit.ODAC$beta_tilde,yy0=0.5)[1]
#   # auc.dist1=ROC.Est.FUN(Di=delta[[l]],yyi=z[[l]]%*%fit.ODAC$beta_tilde1,yy0=0.5)[1]
#   c.dist=Est.Cval(data.frame(cbind(x.valid,delta.valid,z.valid%*%fit.ODAC$beta_tilde[-1])),tau=5)$Dhat
# }
#
# timeend=Sys.time()
# time.dist<-as.numeric(as.POSIXct(timeend)-as.POSIXct(timestart), units="secs")
# print(paste0("time.dist: ",time.dist))
#
# result.c=data.frame( (cbind(p,rho,signal,tau,c.meta,c.dist,c.maximin,c.maximin.less,c.self,c.self.400,c.self.600,corr,time.maximin,
#                             time.dist,time.meta)))
# print(result.c)
#
# # res.beta=data.frame(cbind("true"=b[,L],'beta.meta'=beta_meta,'beta.dist'=fit.ODAC$beta_tilde,'beta.dist1'=fit.ODAC$beta_tilde1,
# #                           'beta.maximin.pool'=as.numeric(beta.maximin.pool),'beta.maximin'=as.numeric(beta.maxmin),
# #                           'beta.self'=as.numeric(beta.self) ))
# # round(res.beta,2)
#
# write.table(result.c, paste(dir.output,batch,"_temp_missing_p",p,".txt", sep =""),append=TRUE,col.names = FALSE, row.names= FALSE)
#
# } # tau
# } # signal
# } # rho
# print(i)
#
# }
#
#
#
#
# # ################ summarize ################
# # rm(list=ls())
# # library(tidyverse)
# # setwd("/n/data1/hsph/biostat/celehs/lab/haz364/maximin_simulations/Output_07012022/")
# # files=list.files()
# # p=50
# # result=NULL
# # for(ii in 1:5){
# #   junk=read.table(paste(getwd(),'/', ii,'_temp_missing_p',p,'.txt', sep =""),  header = FALSE, sep = "", dec = ".",fill = TRUE)
# #   result=rbind(result,junk)
# # }
# #
# # names(result)=c('p','rho','signal','tau','c.meta','c.dist','c.maximin','c.maximin.less',
# #                  'c.self','c.self.400','c.self.600','corr','time.maximin',
# #                 'time.dist','time.meta')
# #
# # res=result
# #
# # signal='moderate'
# # temp1=res[res$signal==signal,c("c.meta","rho",'signal',"tau")];
# # temp1$auc=colnames(temp1)[1]; colnames(temp1)[1]='effect_size'
# # temp1=temp1 %>% group_by(rho,signal,tau) %>%
# #   summarise(rho=rho[1],signal=signal[1],tau=tau[1],auc=auc[1],effect_size=mean(as.numeric(effect_size),na.rm=T) )
# #
# # temp2=res[res$signal==signal,c("c.dist","rho",'signal',"tau")];
# # temp2$auc=colnames(temp2)[1]; colnames(temp2)[1]='effect_size'
# # temp2=temp2 %>% group_by(rho,signal,tau) %>%
# #   summarise(rho=rho[1],signal=signal[1],tau=tau[1],auc=auc[1],effect_size=mean(as.numeric(effect_size),na.rm=T) )
# #
# # temp3=res[res$signal==signal,c("c.maximin","rho",'signal',"tau")];
# # temp3$auc=colnames(temp3)[1]; colnames(temp3)[1]='effect_size'
# # temp3=temp3 %>% group_by(rho,signal,tau) %>%
# #   summarise(rho=rho[1],signal=signal[1],tau=tau[1],auc=auc[1],effect_size=mean(as.numeric(effect_size),na.rm=T) )
# #
# # temp4=res[res$signal==signal,c("c.self","rho",'signal',"tau")];
# # temp4$auc=colnames(temp4)[1]; colnames(temp4)[1]='effect_size'
# # temp4=temp4 %>% group_by(rho,signal,tau) %>%
# #   summarise(rho=rho[1],signal=signal[1],tau=tau[1],auc=auc[1],effect_size=mean(as.numeric(effect_size),na.rm=T) )
# #
# # temp5=res[res$signal==signal,c("c.maximin.less","rho",'signal',"tau")];
# # temp5$auc=colnames(temp5)[1]; colnames(temp5)[1]='effect_size'
# # temp5=temp5 %>% group_by(rho,signal,tau) %>%
# #   summarise(rho=rho[1],signal=signal[1],tau=tau[1],auc=auc[1],effect_size=mean(as.numeric(effect_size),na.rm=T) )
# #
# # temp6=res[res$signal==signal,c("c.self.400","rho",'signal',"tau")];
# # temp6$auc=colnames(temp6)[1]; colnames(temp6)[1]='effect_size'
# # temp6=temp6 %>% group_by(rho,signal,tau) %>%
# #   summarise(rho=rho[1],signal=signal[1],tau=tau[1],auc=auc[1],effect_size=mean(as.numeric(effect_size),na.rm=T) )
# #
# # temp7=res[res$signal==signal,c("c.self.600","rho",'signal',"tau")];
# # temp7$auc=colnames(temp7)[1]; colnames(temp7)[1]='effect_size'
# # temp7=temp7 %>% group_by(rho,signal,tau) %>%
# #   summarise(rho=rho[1],signal=signal[1],tau=tau[1],auc=auc[1],effect_size=mean(as.numeric(effect_size),na.rm=T) )
# #
# # plot.res=rbind(temp1,temp2,temp3,temp4,temp5,temp6,temp7)
# #
# # plot.res.low=plot.res
# # plot.res.moderate=plot.res
# #
# # plot.res.20=rbind.data.frame(plot.res.low,
# #                              plot.res.moderate)
# # plot.res.20$p=p
# #
# # plot.res.final=plot.res.20
# # write.csv(plot.res.final,
# #           file=paste0("/n/data1/hsph/biostat/celehs/lab/haz364/maximin_simulations/results/plot_res_missing_scenario1.csv"),
# #           row.names = F)
# #
# # res.comp.t=result %>%
# #   group_by(p,rho,signal,tau) %>%
# #   summarize("time.maximin"=mean(time.maximin),
# #             "time.dist"=mean(time.dist),
# #             "time.meta"=mean(time.meta)) %>%
# #   ungroup()
# # write.csv(res.comp.t,
# #           file=paste0("/n/data1/hsph/biostat/celehs/lab/haz364/maximin_simulations/results/res_comp_time_missing_scenario1.csv"),
# #           row.names = F)
#
# # # save(plot.res,file='set3.rda')
# #
# # pdf(file=paste0('p_',p,'_signal_',unique(plot.res$signal),'.pdf'),height = 18, width =18)
# # par(mfrow=c(3,3))
# # plot.res %>% ggplot(aes(x=0, y=effect_size, fill=auc)) +
# #   scale_fill_manual(breaks = c('c..meta',"c.dist","c.maximin", "c.self","c.self.400","c.self.600", "c.truebeta"),
# #                     values=c( "grey",'chartreuse3','dodgerblue2','indianred1','indianred1','indianred1','red'))+
# #   geom_bar(stat="identity", color="black",width = 0.8, position = "dodge") +
# #   facet_grid(rows=vars(rho),cols=vars(tau),labeller = label_both) +
# #   xlab(paste('signal_',unique(plot.res$signal))) +ylab("c-stat") +ylim(0,1) +
# #   scale_y_continuous(breaks=seq(0,1,0.1))+scale_x_continuous(breaks=seq(1,20,1))+
# #   geom_hline(yintercept=1,linetype="dashed") +
# #   geom_hline(yintercept=.5,linetype="dashed") +
# #   theme(axis.text.x = element_text(angle=35, hjust = 1,size=11,face="bold"),
# #         strip.text.x = element_text(size = 13),
# #         strip.text.y = element_text(size = 13),
# #         axis.text=element_text(size=11,face="bold"),
# #         axis.title=element_text(size=13,face="bold"),
# #         plot.title = element_text(color="red", size=16, face="bold"),
# #         legend.text = element_text(size=15,face="bold"),
# #         # plot.margin=grid::unit(c(1,1,1,3), "cm"),
# #         legend.position="bottom")
# # dev.off()
# #
