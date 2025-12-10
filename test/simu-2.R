
## demo(ODAC)
# demo(ODACT)
# demo(LATTE)
## run the example in local directory:

setwd('/Users/chongliang/Dropbox/R/DistCox/simu/JJ_pda')

require(data.table)
require(pda)
library(survival)

mydata = fread('JJ_pda_simu_data_20251123.csv', drop = 1)
mydata[,table(site)]
sites = unique(mydata$site) # site names
K = length(sites)  # number of sites
mydata_split <- split(mydata[,-'t_event'], by='site', keep.by=F) 

for(i in K:1){
  # write.csv(mydata_split[[i]], file=paste0('OTA-dry-run/data_site',i,'.csv'), row.names = F)
  if (!dir.exists(sites[i])) {
    dir.create(sites[i])
  } 
}
######################################################################
################################ ODACH ############################### 
######################################################################
         
## fit Cox model using pooled data, stratified by site
fit.pool <- coxph(Surv(time, status)~Trt+Age+Sex+RE+Mutation, data=mydata)
round(summary(fit.pool)$coef, 4)
#     Trt      Age      Sex       RE Mutation 
# -0.1792   0.1803   0.3664   0.1667   0.5959

## ODACH 
# specify your working directory, default is the tempdir
mydir <- 'pda/ODACH'   # tempdir()

# control 
control <- list(project_name = 'PDA J&J demo, ODACH',
                step = 'initialize', # ODACH start with initialize step
                sites = sites,
                heterogeneity = T,   # heterogeneous baseline hazards across sites
                model = 'ODAC',
                family = 'cox',
                outcome = "Surv(time, status)", 
                variables = c('Trt', 'Age', 'Sex', 'RE', 'Mutation'), # covariate var, Trt is the treatment/control assignment
                optim_maxit = 300,
                lead_site = 'site1',   # site1 is the "lead" site
                init_method = "meta",  # use the meta-estimates as the initial estimates for ODACH
                upload_date = as.character(Sys.time()) )

# lead site to create control file
pda(site_id = 'site1', control = control, dir = mydir, upload_without_confirm = T, silent_message = T)

# STEP 1: initialize  
for(i in K:1){
  pda(site_id = sites[i], ipdata = mydata_split[[i]], dir=mydir, upload_without_confirm = T, silent_message = T) 
}
 

# STEP 2: derive
for(i in K:1){
  pda(site_id = sites[i], ipdata = mydata_split[[i]], dir=mydir, upload_without_confirm = T, silent_message = T) 
}


# STEP 3: estimate 
pda(site_id = 'site1', ipdata = mydata_split[[1]], dir=mydir, upload_without_confirm = T, silent_message = T)

# the PDA ODACH is now completed! 

# compare the surrogate estimate with the pooled and meta estimates
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.odach <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet('control', config)
cbind(b.pool= round(fit.pool$coef,4),
      b.meta =control$beta_init,
      b.odach=fit.odach$btilde)

cbind(se.pool= round(summary(fit.pool)$coef[,3],4),
      # b.meta =control$beta_init,
      b.odach=fit.odach$setilde)
 



######################################################################
################################ ODACT ############################### 
######################################################################
mydata[,site:=as.numeric(gsub('site','',site))]
mydata = mydata[order(as.numeric(gsub('site','',site)), time), ]
# mydata = mydata[order(site,time), ]
mydata_split <- split(mydata[,-'t_event'], by='site', keep.by=F) 

summary(mydata$time)
evalt = seq(24,by=3,len=4) # time points to estimate beta(t)
h = 6   # window bandwidth  
px = 5    # 5 coef's 


## fit Cox with beta(t) using pooled data, stratified by site
fit.pool <- mycoxph_bt(mydata[,-c(1,2)], site=mydata$site, fn=llpl_st, times=evalt, h=h, betabar=rep(0,px), hessian=T)
fit.pool$b.pool

## ODACT 
# specify your working directory, default is the tempdir
mydir <- 'pda/ODACT'   # tempdir()

# control
control <- list(project_name = 'PDA J&J simulation, ODACT',
                step = 'initialize',
                sites = sites,
                heterogeneity = TRUE,
                model = 'ODACT',
                family = 'cox',
                outcome = "Surv(time, status)",
                variables = c('Trt', 'Age', 'Sex', 'RE', 'Mutation'),
                times = seq(24,by=3,len=4), # time points to estimate beta(t)  
                bandwidth = 6, 
                optim_maxit = 300,
                lead_site = 'site1',
                init_method = "meta",
                upload_date = as.character(Sys.time()) )
 
# lead site to create control file
pda(site_id = 'site1', control = control, dir = mydir, upload_without_confirm = T, silent_message = T)

# STEP 1: initialize  
for(i in K:1){
  pda(site_id = sites[i], ipdata = mydata_split[[i]], dir=mydir, upload_without_confirm = T, silent_message = T) 
} 


# STEP 2: derive
for(i in K:1){
  pda(site_id = sites[i], ipdata = mydata_split[[i]], dir=mydir, upload_without_confirm = T, silent_message = T) 
}


# STEP 3: estimate 
pda(site_id = 'site1', ipdata = mydata_split[[1]], dir=mydir, upload_without_confirm = T, silent_message = T)

# the PDA ODACT is now completed! 

# compare the surrogate estimate with the pooled and meta estimates
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.odact <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet('control', config)
cbind(b.pool= round(fit.pool$b.pool[1,],4),
      b.meta =control$beta_init[1,],
      b.odact=fit.odact$btilde[1,])

cbind(se.pool= round(summary(fit.pool)$coef[,3],4),
      # b.meta =control$beta_init,
      se.odact=fit.odact$setilde)





######################################################################
################################ LATTE ############################### 
######################################################################

# we will use the survival status (0=survived at EOS, 1=observed death before EOS) as binary outcome to conduct LATTE
 
# !!! specify your working directory, .json files will be written to this dir
mydir <- 'pda/LATTE'    
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))]) # clear any existing files
 
# simulate some nco outcomes 
# Negative Control Outcomes (NCOs) setup for LATTE calibration
# nco_outcomes = c("nco1", "nco2", "nco3")
# mydata$nco1 = rbinom(nrow(mydata), 1,0.3)
# mydata$nco2 = rbinom(nrow(mydata), 1,0.3)
# mydata$nco3 = rbinom(nrow(mydata), 1,0.3)
# 
# # simulate the time to event for nco outcomes
# mydata$nco1_time = sample(mydata$time, nrow(mydata), replace = TRUE)
# mydata$nco2_time = sample(mydata$time, nrow(mydata), replace = TRUE)
# mydata$nco3_time = sample(mydata$time, nrow(mydata), replace = TRUE)

# split data by sites
sites = unique(mydata$site)
K = length(sites)
mydata = data.frame(mydata) 
mydata_split <- split(mydata, mydata$site)[c(1,3:10,2)]
 
# control 
control <- list(project_name = 'PDA J&J demo, LATTE',
                # LATTE start with initialize step 
                step = "initialize", 
                sites = sites,
                model = "LATTE",
                family = "binomial",
                outcome = "status",
                # nco_outcomes = nco_outcomes,
                variables = c("Age", "Sex", "RE", "Mutation"),
                treatment = "treatment", # Trt
                lead_site = "site1",
                # balance covar distributions
                balancing_method = "stratification", 
                outcome_model = "logistic", 
                upload_date = as.character(Sys.time()) )

# lead site to create and write the control.json file
pda(site_id = 'site1', control = control, dir = mydir, upload_without_confirm = T, silent_message = T)


# STEP 1: initialize
for(i in K:1){ 
  pda(site_id = sites[i], ipdata = mydata_split[[i]], dir=mydir, upload_without_confirm = T, silent_message = T)
}
 

# STEP 2: estimate
pda(site_id = 'site1', ipdata = mydata_split[[1]], dir=mydir, upload_without_confirm = T, silent_message = T)

# the PDA LATTE is now completed!

# present the LATTE results
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.latte <- pdaGet(name = 'site1_estimate', config = config) 
# LATTE (with no NCO calibration) 
c(fit.latte$by_outcome[[outcome_id]]$coefficients, 
  fit.latte$by_outcome[[outcome_id]]$se )
 