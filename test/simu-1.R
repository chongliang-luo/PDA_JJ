


## J&J data simulation  
setwd('/Users/chongliang/Dropbox/R/DistCox/simu/JJ_pda')
# rm(list=ls()) 

require(devtools)
library(survival)
library(data.table)
# install_github('Penncil/pda')
library(pda)

# require(ODACO) 
# source('/Users/chongliang/Dropbox/R/DistCox/engine/DistCox_misc.R')


# source('/Users/chongliang/Dropbox/R/DistCox/engine/DistCox_misc.R')
# run_ODACH_with_pda: run ODACH or ODACH_CC with pda silently at local
# source('/Users/chongliang/Dropbox/R/DistCox/simu/ODAC_ODACH/simu/DistCox.R')  
source('/Users/chongliang/Dropbox/R/DistCox/simu/JJ_pda/misc.R')




## compare   
methods <- c('PooledSt', 'Pooled', 'Meta', 'ODACH', 'ODAC')
bt = log(c(0.77, 1.2, 1.5, 1.2, 1.8))

px <- 5 
K <- 10   
# n = c(rep(300,K/2), rep(200,K/2))
# N <- sum(n)
# rr  <- 0.34
nrep <- 100  
sites = paste0('site', 1:K)

BETA_hat <- array(NA, c(nrep, length(methods), px, 2))
simu_data = list()

tt <- Sys.time()
for(irep in 1:nrep){
  cat('\n----- irep = ', irep, '/', nrep, '-----', as.character(Sys.time())) 
      
  all_data = generate_Cox_data_all_sites( )
  # mean(all_data$all_data$status)  
  mydata <- data.table(all_data$all_data )
  myparam <- all_data$param 
  simu_data[[irep]] = all_data
  
  ## meta
  sum_K_wt <- rep(0, myparam$px)
  sum_K_beta_w <- rep(0, myparam$px) 
  for(sid in sites){ 
    # coxph sometimes error when # events too small
    fit_cox_i <- tryCatch(coxph(Surv(time, status) ~Trt+Age+Sex+RE+Mutation, data=mydata[site==sid,]) 
                          , error = function(e) list()) 
    if(!is.null(fit_cox_i$coef) & !any(is.na(fit_cox_i$coef))){
      wt <- 1/(summary(fit_cox_i)$coef[,3])^2  # inverse-var as weight
      sum_K_wt <- sum_K_wt + wt
      sum_K_beta_w <- sum_K_beta_w + fit_cox_i$coef * wt 
    } else {
      cat(sid)
    }
  }
  b.meta = sum_K_beta_w / sum_K_wt
  se.meta = sqrt(1 / sum_K_wt)
      
  #1# pooled, stratified by site
  fit_coxst_pkg <- coxph(Surv(time, status)~Trt+Age+Sex+RE+Mutation+ strata(site), data=mydata)
  BETA_hat[irep,1, ,1] <- fit_coxst_pkg$coef
  BETA_hat[irep,1, ,2] <- summary(fit_coxst_pkg)$coef[,3]
  
  #2# pooled, not stratified by site
  fit_cox_pkg <- coxph(Surv(time, status)~Trt+Age+Sex+RE+Mutation, data=mydata)
  BETA_hat[irep,2, ,1] <- fit_cox_pkg$coef
  BETA_hat[irep,2, ,2] <- summary(fit_cox_pkg)$coef[,3]
  
  #3# meta
  BETA_hat[irep,3, ,1] <- b.meta
  BETA_hat[irep,3,, 2] <- se.meta
  
  #4# ODACH
  mydata_split <- split(mydata[,-'t_event'], by='site', keep.by=F) 
  control <- list(project_name = 'J&J simulation',
                  step = 'initialize',
                  sites = sites,
                  heterogeneity = T,
                  model = 'ODAC',
                  family = 'cox',
                  outcome = "Surv(time, status)",
                  variables = c('Trt', 'Age', 'Sex', 'RE', 'Mutation'),
                  # xlev = list(sex=c('F', 'M')),  #levels of all categorical X's, with the first being the reference
                  optim_maxit = 300,
                  lead_site = 'site1',
                  init_method = "meta",
                  upload_date = as.character(Sys.time()) ) 
  fit_distcox <- run_ODACH_with_pda(control, mydir='pda/', mydata_split, upload_without_confirm=T, silent_message=T)
  BETA_hat[irep,4, ,1] <- fit_distcox$btilde
  BETA_hat[irep,4, ,2] <- fit_distcox$setilde
   
  # mydata0 = data.frame(mydata[,-'t_event'])
  # mydata0$site = c(rep(1,400), rep(2:7,each=300), rep(8:10,each=100))
  # fit_distcox <- DistCox(mydata0, id.local=1, init_est = b.meta, strat=T, verbose=F)
  # BETA_hat[irep, ik, ir,4, ] <- fit_distcox$beta_tilde
  
  #5# ODAC  
  control$heterogeneity = F
  fit_distcox <- run_ODACH_with_pda(control, mydir='pda/', mydata_split, upload_without_confirm=T, silent_message=T)
  BETA_hat[irep,5, ,1] <- fit_distcox$btilde
  BETA_hat[irep,5, ,2] <- fit_distcox$setilde
}
tt <- Sys.time() - tt

which(abs(BETA_hat[,3,1,1] -BETA_hat[,1,1,1]) > abs(BETA_hat[,4,1,1] -BETA_hat[,1,1,1]) )

# save(list = ls(), file='JJ_pda_simu_20251123.rda')

# prepare data for output, select the one with most X's coef to be better than meta (irep=100)
which(rowSums(abs(BETA_hat[,3,,1] -BETA_hat[,1,,1]) > abs(BETA_hat[,4,,1] -BETA_hat[,1,,1]) ) ==5)
# write.csv(all_data$all_data, file='JJ_pda_simu_data_20251123.csv')



boxplot(BETA_hat[,2:5,1,1] - BETA_hat[,1,1,1], ylim=c(-1, 1)/20)
abline(h=0)
# Trt, Age, Sex, Race/Ethnicity, Mutation

boxplot(BETA_hat[,1:5,1,1], ylim=c(-0.5, 0.2)); abline(h=bt[1])
boxplot(BETA_hat[,1:5,1,2], ylim=c(0.06, 0.09)) 
apply(BETA_hat[,1:5,1,1],2,mean)
apply(BETA_hat[,1:5,1,1]-1.96*BETA_hat[,1:5,1,2] < bt[1] & BETA_hat[,1:5,1,1]+1.96*BETA_hat[,1:5,1,2] > bt[1],2,mean,na.rm=T)

 



################################### ggplot  ####################################
require(ggplot2)
require(reshape2)
nm <- nn[1]
# BETA_hat[irep, unbalance, event.rate, method, mean/var, beta1/2]
beta.bias <- BETA_hat[, , ,3:6,2] 
for(ii in 1:4) beta.bias[,,,ii] <- beta.bias[,,,ii]-BETA_hat[, , ,1,2]  
dim(beta.bias) # 200   2   3   4

bias.df <- melt(beta.bias, varnames = c('rep', 'unbalance', 'event.rate', 'method'),
                value.name = "bias.to.pooled")
bias.df$unbalance <- as.factor(bias.df$unbalance)
# levels(bias.df$unbalance) <- paste0((2*nm-nn), '*', (K/2), '+', nn, '*', K/2)
levels(bias.df$unbalance) <- c('500*10', '750*5+250*5')

bias.df$event.rate <- as.factor(bias.df$event.rate)
levels(bias.df$event.rate) <- paste0('event rate = ', rr*100, '%')

bias.df$method <- as.factor(bias.df$method)
levels(bias.df$method) <- methods[3:6]   
head(bias.df)


## 4 methods shonw in supp  
# pdf('simu/ODACH/ODACH_simu_box_20220218_rev_m4.pdf', width = 12, height = 6)
# ggplot(bias.df, aes(x=unbalance, y=bias.to.pooled, fill=method)) +
ggplot(bias.df[bias.df$unbalance=='500*10',], aes(x=method, y=bias.to.pooled )) +
  geom_boxplot(notch = TRUE, outlier.size = .5) + 
  facet_grid(. ~ event.rate 
             # rows=vars(unbalance)
             # , cols=vars(event.rate)
             # , scales="free"
  ) +
  # geom_line() +
  # geom_hline(mapping = NULL, data = NULL, ..., yintercept,
  # na.rm = FALSE, show.legend = NA) + 
  # geom_point(size=2) +
  # xlim(-0.2, 4.2) + 
  ylim(-0.2, 0.1) +
  labs(title='', 
       # 'Boxplot: bias to pooled est, true effect size = 2 \n event rate = 20%, 5%, 2%, 1%, large site * 5 + small site * 5',
       x =  '',  
       y = "Relative bias to gold standard" ) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)
        , axis.ticks.x = element_blank()
        # , plot.title = element_text(family, face, colour, size)
        # , panel.background = element_blank()
        , panel.grid.minor = element_blank()
        # , panel.grid.major = element_blank()
        # , panel.border = element_rect()
        , axis.text.x = element_text(size = 10, face="bold")  
        , axis.text.y = element_text(size = 15) # angle=30,  hjust=1, size = 9
        , strip.text.x = element_text(size = 15)
        , legend.title = element_blank()
        , legend.position = c(0.1, 0.2)
        # , legend.key = element_rect(size=0.2)
        # , legend.key.height = unit(0.4, "cm")
  )  
# dev.off()


## 2 methods shonw in manuscript (meta and ODACH)
# pdf('simu/ODACH/ODACH_simu_box_20220218_rev_m2.pdf', width = 7, height = 6)
ggplot(bias.df[bias.df$unbalance=='500*10'&bias.df$method%in%methods[3:4],], aes(x=method, y=bias.to.pooled )) +
  geom_boxplot(notch = TRUE, outlier.size = .5) + 
  facet_grid(. ~ event.rate 
             # rows=vars(unbalance)
             # , cols=vars(event.rate)
             # , scales="free"
  ) + 
  ylim(-0.2, 0.1) +
  labs(title='', 
       # 'Boxplot: bias to pooled est, true effect size = 2 \n event rate = 20%, 5%, 2%, 1%, large site * 5 + small site * 5',
       x =  '',  
       y = "Relative bias to gold standard" ) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)
        , axis.ticks.x = element_blank()
        # , plot.title = element_text(family, face, colour, size)
        # , panel.background = element_blank()
        , panel.grid.minor = element_blank()
        # , panel.grid.major = element_blank()
        # , panel.border = element_rect()
        , axis.text.x = element_text(size = 10, face="bold")  
        , axis.text.y = element_text(size = 15) # angle=30,  hjust=1, size = 9
        , strip.text.x = element_text(size = 15)
        , legend.title = element_blank()
        , legend.position = c(0.1, 0.2)
        # , legend.key = element_rect(size=0.2)
        # , legend.key.height = unit(0.4, "cm")
  ) 
# dev.off()
