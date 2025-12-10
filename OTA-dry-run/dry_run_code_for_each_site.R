
# load required packages
require(data.table)
require(survival)
require(devtools)
# use the most recent GitHub version (1.3.1 as of 12/08/2025)
devtools::install_github("penncil/pda")
require(pda)

# !!! specify your working directory, and your site ID
setwd('/Users/chongliang/Dropbox/PDA_development/PDA_JJ/OTA-dry-run/')
mysite = 'site1' # 'site2'


# read in your site data
mydata = fread(paste0('data_',mysite,'.csv'))
sites = paste0('site', 1:10) # unique site names

 

######################################################################
################################ ODACH ############################### 
######################################################################

# !!! specify your working directory, .json files will be written to this dir
mydir <- 'site1_lead'   # 'site2'
if (!dir.exists(mydir)) {
  dir.create(mydir)
} 
file.remove(list.files(mydir,full.names=T)[grepl('.json', list.files(mydir))]) # clear any existing files



# STEP 0 (only lead site): create and write the control.json file  
control <- list(project_name = 'PDA J&J demo, ODACH',
                # ODACH start with initialize step
                step = 'initialize', 
                sites = sites,
                # assume heterogeneous baseline hazards across sites, i.e. stratified by site
                heterogeneity = T,   
                model = 'ODAC',
                family = 'cox',
                outcome = "Surv(time, status)", 
                # Trt is the treatment/control assignment 
                variables = c('Trt', 'Age', 'Sex', 'RE', 'Mutation'), 
                optim_maxit = 300,
                # site1 is assigned as the "lead" site
                lead_site = 'site1',   
                # use the meta-estimates as the initial estimates for ODACH
                init_method = "meta",  
                upload_date = as.character(Sys.time()) )
pda(site_id = 'site1', control = control, dir = mydir)
# upload control.json to OTA


# STEP 1: initialize
# download control.json from OTA to your working dir
pda(site_id = mysite, ipdata = mydata, dir=mydir)
# upload siteX_initialize.json to OTA
# (only lead site) also upload updated control.json to OTA


# STEP 2: derive
# download control.json from OTA to your working dir
pda(site_id = mysite, ipdata = mydata, dir=mydir)
# upload siteX_derive.json to OTA
# (only lead site) also upload updated control.json to OTA


# STEP 3 (only lead site): estimate 
pda(site_id = mysite, ipdata = mydata, dir=mydir)


# the PDA ODACH is now completed!


# (only lead site): present the ODACH results
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.odach <- pdaGet(name = 'site1_estimate', config = config)
data.frame(var=control$variables, 
           b.odach=fit.odach$btilde,   # beta coef's
           se.odach=fit.odach$setilde) # st errors of the beta's


# (only lead site):
# please collect all the .json files and upload to OTA!





######################################################################
################################ ODACT ############################### 
######################################################################

# !!! specify your working directory, .json files will be written to this dir
mydir <- 'pda/ODACT'    
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))]) # clear any existing files

 

# STEP 0 (only lead site): create and write the control.json file  
control <- list(project_name = 'PDA J&J demo, ODACT',
                # ODACT start with initialize step
                step = 'initialize', 
                sites = sites,
                # assume heterogeneous baseline hazards across sites, i.e. stratified by site
                heterogeneity = T,   
                model = 'ODACT',
                family = 'cox',
                outcome = "Surv(time, status)", 
                # Trt is the treatment/control assignment 
                variables = c('Trt', 'Age', 'Sex', 'RE', 'Mutation'), 
                # time points to estimate coefficients beta(t) 
                times = seq(24,by=3,len=4), 
                # bandwidth of data around the time points
                bandwidth = 6,
                optim_maxit = 300,
                # site1 is assigned as the "lead" site
                lead_site = 'site1',   
                # use the meta-estimates as the initial estimates for ODACT
                init_method = "meta",  
                upload_date = as.character(Sys.time()) )
pda(site_id = 'site1', control = control, dir = mydir)
# upload control.json to OTA


# STEP 1: initialize
pda(site_id = mysite, ipdata = mydata, dir=mydir)
# upload siteX_initialize.json to OTA
# (only lead site) also upload updated control.json to OTA


# STEP 2: derive
pda(site_id = mysite, ipdata = mydata, dir=mydir)
# upload siteX_derive.json to OTA
# (only lead site) also upload updated control.json to OTA


# STEP 3 (only lead site): estimate 
pda(site_id = mysite, ipdata = mydata, dir=mydir)


# the PDA ODACT is now completed!


# (only lead site): present the ODACT results
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.odact <- pdaGet(name = 'site1_estimate', config = config)
b.odact = fit.odact$btilde
se.odact = fit.odact$setilde
row.names(b.odact) = row.names(se.odact) = control$variables
colnames(b.odact) = colnames(se.odact) = paste0('time_', control$times)
b.odact   # beta coef's
se.odact  # s.e. of the beta's


# (only lead site):
# please collect all the .json files and upload to OTA!
 


