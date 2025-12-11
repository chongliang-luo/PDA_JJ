# we will use the survival status (0=survived at EOS, 1=observed death before EOS) as binary outcome to conduct LATTE

# !!! specify your working directory, .json files will be written to this dir
setwd('/Users/chongliang/Dropbox/PDA_development/PDA_JJ/Local-dry-run/LATTE/')
# mydir <- "/Users/luli/Documents/developer/pda1210/pda/demo/pda_LATTE_demo"
# dir.create(mydir)
mydir = getwd()
# file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))]) # clear any existing files

# simulate some nco outcomes 
# Negative Control Outcomes (NCOs) setup for LATTE calibration
mydata = fread("../../Deliverable-dry-run/JJ_pda_simu_data_20251123.csv",drop = 1)
# mydata = data %>% select(-X)
nco_outcomes = c("nco1", "nco2", "nco3")
mydata$nco1 = rbinom(nrow(mydata), 1,0.3)
mydata$nco2 = rbinom(nrow(mydata), 1,0.3)
mydata$nco3 = rbinom(nrow(mydata), 1,0.3)

# simulate the time to event for nco outcomes
mydata$nco1_time = sample(mydata$time, nrow(mydata), replace = TRUE)
mydata$nco2_time = sample(mydata$time, nrow(mydata), replace = TRUE)
mydata$nco3_time = sample(mydata$time, nrow(mydata), replace = TRUE)

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
                treatment = "Trt", # Trt
                lead_site = "site1",
                # balance covar distributions
                balancing_method = "stratification", 
                outcome_model = "logistic", 
                upload_date = as.character(Sys.time()) )

pda(site_id = 'site1', control = control,  dir = mydir)
# lead site to create and write the control.json file
# pda(site_id = 'site1', control = control, dir = mydir, upload_without_confirm = T, silent_message = T)


# STEP 1: initialize
for(i in K:1){ 
  pda(site_id = sites[i], ipdata = mydata_split[[i]], dir=mydir, upload_without_confirm = T, silent_message = T)
}


# STEP 2: estimate
# pda(site_id = 'site1', ipdata = mydata_split[[1]], dir=mydir, upload_without_confirm = T, silent_message = T)
pda(site_id = 'site1', dir=mydir)

# the PDA LATTE is now completed!

# present the LATTE results
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.latte <- pdaGet(name = 'site1_estimate', config = config) 
# LATTE (with no NCO calibration) 
c(fit.latte$coefficients, 
  fit.latte$se )
