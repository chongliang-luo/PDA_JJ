
h0 <- function(t, p, lam){ ## weibull baseline hazard
  p * lam^p * t^(p-1)
}

f0 <- function(t, p, lam){ # Weibull pdf, p=shape, 1/lam=scale
  p * lam^p * t^(p-1) * exp(-(lam*t)^p)
}

S0inv_weibull <- function(u, p=weibull_shape, lam=1/weibull_scale){ # inv surv of weibull 
  lam^(-1) * (-log(1-u))^(1/p)
  # exp(-(lam*t)^p)
}


## generate data of survival time, using Cox PH with Weibull baseline hazard 
generate_Cox_data_all_sites <- function( 
    param = list(# N = c(rep(350,10/2), rep(150,10/2)), 
                 N = c(400, rep(300,6), rep(100,3)),
                 px = 5, 
                 beta = log(c(0.77, 1.2, 1.5, 1.2, 1.8)),
                 beta0 = -1, # intercept to adjust event rate
                 S0inv = S0inv_weibull,
                 weibull_scale_range = c(30,40),
                 weibull_shape_range = c(8,12),
                 censor = T,
                 event.rate = 0.34,
                 tie = F),
    X = NULL, verbose=F){ 
  # here are some averages and ranges from RCT and RW data for data simulation:
  # - Event rate: 0.34 (0.10–0.59)
  # - Baseline hazard: 0.023 (0.015–0.053)
  # - Median survival time: 34.4 months (13.2–47.5)
  # - Hazard ratio: 0.77 (0.75–0.80)
  # - Sample size: 2,500 total # of sites: 10
  # - Covariates: Age, Sex, Race/Ethnicity, Treatment regimen, Mutation subtype
  # Treatment: trt vs control, HR=0.77
  # Age, continuous, (age-60)/10, HR=1.2
  # Sex, M vs F,  HR=1.5
  # Race/Ethnicity, NHB vs NHW, HR=1.2
  # Mutation subtype, 1 vs 0, HR=1.8  
  N = param$N 
  K = length(N) # param$K
  px = param$px
  beta = param$beta
  beta0 = param$beta0
  S0inv = param$S0inv
  weibull_scale_range = param$weibull_scale_range
  weibull_shape_range = param$weibull_shape_range
  weibull_scale = runif(K, weibull_scale_range[1], weibull_scale_range[2])
  weibull_shape = exp(runif(K, log(weibull_shape_range[1]), log(weibull_shape_range[2]))) 
  param$weibull_scale = weibull_scale
  param$weibull_shape = weibull_shape
  censor = param$censor
  event.rate = param$event.rate
  tie = param$tie 
  
  # Covariates X, by sites
  # allow covar mean to be diff across sites
  site = rep(1:K, N)
  age.m = runif(K, -0.5, 0.5)
  sex.m = runif(K, 0.5, 0.6)  # % Male
  RE.m = runif(K, 0.05, 0.25) # % NHB
  mutation.m = runif(K, 0.1, 0.3) 
  param$age.m = age.m
  param$sex.m = sex.m
  param$RE.m = RE.m
  param$mutation.m = mutation.m
  if(is.null(X)){
    X = matrix(0, 0, 5)
    for(k in 1:K){ 
      X <- rbind(X,
                 cbind(Trt=rep(0:1,each=N[k]/2),
                       Age= round(rnorm(N[k],age.m[k],1),1),
                       Sex=rbinom(N[k], 1, sex.m[k]),
                       RE=rbinom(N[k], 1, RE.m[k]),
                       Mutation=rbinom(N[k], 1, mutation.m[k]) ))
    }
  } 
  Xb <- c(X %*% beta + beta0)
  expXb <- exp(Xb)
  # plot(density(Xb))
  
  # generate true event Time with Weibull baseline distribution
  U <- runif(sum(N), 0, 1)
  t_event = rep(NA, sum(N))
  for(k in 1:K){ 
    t_event[site==k] <- S0inv(U[site==k], p=weibull_shape[k], lam=1/weibull_scale[k]) / expXb[site==k]^(1/weibull_shape[k])  
  }
  t_event = round(t_event, 2)
  # summary(t_event) 
  # adjust weibull_scale_range and weibull_shape_range to match Median survival time: 34.4 months (13.2–47.5)
  
  ################################################################################################
  # NOTE: censoring has to be non-informative, 
  # and how to generating censoring time for given event.rate is tricky
  # better not to use: 
  # 1. one censoring time, e.g. quantile(t_event, event.rate)
  # 2. censoring time within a range, e.g. runif(N, 0.1, quantile(t_event, min(0.99,event.rate*2)))
  # suggest use:
  # Weibull r.v. with the same shape but scale weibull_scale_c = weibull_scale * ratio
  # notice ith sub has hazard w scale = weibull_scale / expXb[i]^(1/weibull_shape)
  # do some integral, know: 
  # censoring rate = mean(1/(1+(ratio)^weibull_shape*expXb))
  # thus do 1-d search for opt ratio
  ################################################################################################
  # if(event.rate < 1){
  #   ratio = ((1/(1-event.rate)-1)/mean(expXb))^(1/mean(weibull_shape))
  #   fn <- function(rr) abs(mean(1/(1+rr^mean(weibull_shape)*expXb)) - (1-event.rate) )
  #   ratio <- optimize(fn, interval = c(0.5, 1.5)*ratio)$minimum
  #   weibull_scale_c <- mean(weibull_scale) * ratio
  #   # weibull_scale_c <- weibull_scale * ((1/(1-event.rate)-1)/mean(expXb))^(1/weibull_shape)
  #   
  #   t_censor <- rweibull(N, mean(weibull_shape), weibull_scale_c)
  #   # t_censor <- runif(N, 0.1, quantile(t_event, min(0.99,event.rate*2)))
  #   # t_censor <- quantile(t_event, event.rate)
  #   ind_event <- ifelse(t_event < t_censor, 1, 0)
  #   # mean(ind_event) # to match Event rate: 0.34
  #   t_surv <- pmin(t_event, t_censor)
  #   t_surv = round(t_surv, 2)
  # } else {
  #   ind_event <- rep(1, length(t_event))
  #   t_surv <- t_event
  # }
  
  # administrative censoring (e.g. a fixed FU time)
  t_censor <- quantile(t_event, event.rate)
  t_censor
  ind_event <- ifelse(t_event < t_censor, 1, 0)
  mean(ind_event) # to match Event rate: 0.34
  t_surv <- pmin(t_event, t_censor)
  t_surv = round(t_surv, 2)
  
  if(verbose) cat('censoring rate = ', 1-mean(ind_event), '\n')
  if(verbose) plot(density(t_surv))
  if(tie) t_surv <- ceiling(t_surv)  # pmax(round(t_surv), 1)  # 
  
  all_data <- data.frame(site=paste0('site',site), t_event = t_event, time=t_surv, status=ind_event, X)
  return(list(all_data = all_data, param = param))
}




## run ODACH or ODACH_CC with pda silently at local 
run_ODACH_with_pda <- function(control, mydir, mydata, upload_without_confirm=T, silent_message=T){
  file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])
  N = sum(unlist(lapply(mydata, function(a) nrow(a))))
  K = length(control$sites)
  px = length(control$variables)
  # line #180 in ODAC.R 
  for(sid in 1:length(mydata)) mydata[[sid]][,1]=round(mydata[[sid]][,1],4)
  
  # ############################  STEP 1: initialize  ###############################
  ## assume lead site1: enter "1" to allow transferring the control file
  pda(site_id = control$lead_site, control = control, dir = mydir, 
      upload_without_confirm= upload_without_confirm, silent_message=silent_message)
  for(sid in K:1) pda(site_id = control$sites[sid], ipdata = mydata[[sid]], dir=mydir, 
                      upload_without_confirm =upload_without_confirm, silent_message=silent_message)
  
  if(control$heterogeneity ==F){ # ODAC
    # ############################  STEP 2: derivativeUWZ  ###############################
    for(sid in K:1) pda(site_id = control$sites[sid], ipdata = mydata[[sid]], dir=mydir, 
                        upload_without_confirm = upload_without_confirm, silent_message=silent_message)
  }
  
  # ############################  STEP 3: derivative  ###############################
  for(sid in K:1) pda(site_id = control$sites[sid], ipdata = mydata[[sid]], dir=mydir, 
                      upload_without_confirm = upload_without_confirm, silent_message=silent_message)
  
  # ############################  STEP 4: estimate  ###############################
  tt=tryCatch(pda(site_id = control$lead_site, ipdata = mydata[[1]], dir=mydir, 
                  upload_without_confirm = upload_without_confirm, silent_message=silent_message),error=function(e) NA)
  
  config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
  if(!is.na(tt[1])){
    fit.pda <- pdaGet(paste0(control$lead_site,'_estimate'), config = config)
    return(list(btilde = fit.pda$btilde, setilde=sqrt(diag(solve(fit.pda$Htilde))/N))) # 
  }else{
    return(list(btilde = NA, setilde=NA))
  }
}
