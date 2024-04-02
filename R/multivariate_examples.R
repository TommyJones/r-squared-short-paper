#### load libraries ####
library(mvrsquared)
library(tidyverse)

#### load data ####

##### amitriptyline #####
# source: M. V. Rudorfer (1982). Cardiovascular Changes and Plasma Drug Levels after
#         Amitriptyline Overdose. Journal of Toxicology-Clinical Toxicology, 19, 67-71.
# secondary source: R. A. Johnson and D. W. Wichern (2007). Applied Multivariate Statistical Analysis
#                   6th Ed. Pearson Prentice Hall, Upper Saddle River, NJ, USA.
# amit          <- read.table('./Data/T7-6.DAT', header = FALSE)

ami_data <- 
  read.table("data-raw/multivariate-example/T7-6.DAT", header = FALSE)

names(ami_data)  <- c('tot', 'ami', 'sex', 'amt', 'pr', 'diap', 'qrs')

ami_data <- as_tibble(ami_data)

# quick description: data containing 17 overdoses of the drug amitriptyline
# outcomes: total TCAD plasma level (tot), amount of amitriptyline present in TCAD plasma levels (ami)
# covariates: sex (1 for female, 0 for male), amount of antidepressants taken at time of overdose (amt)
#             PR wave measurement (pr), diastolic blood pressure (diap), QRS wave measurement (qrs)

##### pulp and paper properties #####
# primary source: R. A. Johnson and D. W. Wichern (2007). Applied Multivariate Statistical Analysis
#                 6th Ed. Pearson Prentice Hall, Upper Sadle River, NJ, USA.

pulp_data <- read.table("data-raw/multivariate-example/T7-7.dat")

names(pulp_data)   <- c('bl', 'em', 'sf', 'bs', 'afl', 'lff', 'fff', 'zst')



# quick description: Measurements of properties of pulp fibers and the paper made from them
# outcomes: breaking length (bl), elastic modulus (em), stress at failure (sf), burst strength (bs)
# covariates: arithmetic fiber length (afl), long fiber fraction (lff), fine fiber fraction (fff), zero span tensile (zst)



#### Amitriptyline Analysis ####
# ya            <- with(ami_data, cbind(tot, ami))

ya <- 
  ami_data |> 
  select(
    tot,
    ami
  ) |>
  as.matrix()

##### full model #####
ami_full <- 
  lm(cbind(tot, ami) ~ amt + sex + pr + diap + qrs, data = ami_data)

# ami_full    <- lm(cbind(tot, ami) ~ amt + sex + pr + diap + qrs, data = ami_data)

ypaf <- predict(ami_full)

calc_rsquared(ya, ypaf)

##### single covariate #####
ami_amt     <- lm(cbind(tot, ami) ~ amt, data = ami_data)
yp_amt      <- predict(ami_amt)

calc_rsquared(ya, yp_amt)

#### two covariates ####
ami_am_sex   <- lm(cbind(tot, ami) ~ amt + sex, data = ami_data)
yp_am_sex    <- predict(ami_am_sex)

ami_am_pr   <- lm(cbind(tot, ami) ~ amt + pr, data = ami_data)
yp_am_pr    <- predict(ami_am_pr)

ami_am_di   <- lm(cbind(tot, ami) ~ amt + diap, data = ami_data)
yp_am_di    <- predict(ami_am_di)

ami_am_qr   <- lm(cbind(tot, ami) ~ amt + qrs, data = ami_data)
yp_am_qr    <- predict(ami_am_qr)

calc_rsquared(ya, yp_am_sex)
calc_rsquared(ya, yp_am_pr)
calc_rsquared(ya, yp_am_di)
calc_rsquared(ya, yp_am_qr)

## percent change from one covariate #
(calc_rsquared(ya, yp_am_sex) - calc_rsquared(ya, yp_amt))/calc_rsquared(ya, yp_amt) # highest percent change
(calc_rsquared(ya, yp_am_pr) - calc_rsquared(ya, yp_amt))/calc_rsquared(ya, yp_amt)
(calc_rsquared(ya, yp_am_di) - calc_rsquared(ya, yp_amt))/calc_rsquared(ya, yp_amt)
(calc_rsquared(ya, yp_am_qr) - calc_rsquared(ya, yp_amt))/calc_rsquared(ya, yp_amt)

#### three covariates ####
ami_am_sex_pr   <- lm(cbind(tot, ami) ~ amt + sex + pr, data = ami_data)
yp_am_sex_pr    <- predict(ami_am_sex_pr)

ami_am_sex_di   <- lm(cbind(tot, ami) ~ amt + sex + diap, data = ami_data)
yp_am_sex_di    <- predict(ami_am_sex_di)

ami_am_sex_qr   <- lm(cbind(tot, ami) ~ amt + sex + qrs, data = ami_data)
yp_am_sex_qr    <- predict(ami_am_sex_qr)

calc_rsquared(ya, yp_am_sex_pr)
calc_rsquared(ya, yp_am_sex_di)
calc_rsquared(ya, yp_am_sex_qr)

## percent change from one covariate #
(calc_rsquared(ya, yp_am_sex_pr) - calc_rsquared(ya, yp_am_sex))/calc_rsquared(ya, yp_am_sex) # highest percent change
(calc_rsquared(ya, yp_am_sex_di) - calc_rsquared(ya, yp_am_sex))/calc_rsquared(ya, yp_am_sex)
(calc_rsquared(ya, yp_am_sex_qr) - calc_rsquared(ya, yp_am_sex))/calc_rsquared(ya, yp_am_sex)

#### four covariate model ####
ami_am_sex_pr_di  <- lm(cbind(tot, ami) ~ amt + sex + pr + diap, data = ami_data)
yp_am_sex_pr_di   <- predict(ami_am_sex_pr_di)

ami_am_sex_pr_qr  <- lm(cbind(tot, ami) ~ amt + sex + pr + qrs, data = ami_data)
yp_am_sex_pr_qr   <- predict(ami_am_sex_pr_qr)

calc_rsquared(ya, yp_am_sex_pr_di)
calc_rsquared(ya, yp_am_sex_pr_qr)

## percent change from one covariate #
(calc_rsquared(ya, yp_am_sex_pr_di) - calc_rsquared(ya, yp_am_sex_pr))/calc_rsquared(ya, yp_am_sex_pr) # highest percent change
(calc_rsquared(ya, yp_am_sex_pr_qr) - calc_rsquared(ya, yp_am_sex_pr))/calc_rsquared(ya, yp_am_sex_pr)


#### Pulp Analysis ####
Yf        <- with(pulp_data, cbind(bl, em, sf, bs))

pup_full  <- lm(cbind(bl, em, sf, bs) ~ afl + lff + fff + zst, data = pulp_data)
ypf       <- predict(pup_full)

calc_rsquared(Yf, ypf)


#### single covariate ####
pup_afl   <- lm(cbind(bl, em, sf, bs) ~ afl, data = pulp_data)
yp_afl    <- predict(pup_afl)

pup_lff   <- lm(cbind(bl, em, sf, bs) ~ lff, data = pulp_data)
yp_lff    <- predict(pup_lff)

pup_fff   <- lm(cbind(bl, em, sf, bs) ~ fff, data = pulp_data)
yp_fff    <- predict(pup_fff)

pup_zst   <- lm(cbind(bl, em, sf, bs) ~ zst, data = pulp_data)
yp_zst    <- predict(pup_zst)

calc_rsquared(Yf, yp_afl)
calc_rsquared(Yf, yp_lff)
calc_rsquared(Yf, yp_fff)
calc_rsquared(Yf, yp_zst) # best model

#### two covariate model ####
pup_zst_afl   <- lm(cbind(bl, em, sf, bs) ~ zst + afl, data = pulp_data)
yp_zst_afl    <- predict(pup_zst_afl)

pup_zst_lff   <- lm(cbind(bl, em, sf, bs) ~ zst + lff, data = pulp_data)
yp_zst_lff    <- predict(pup_zst_lff)

pup_zst_fff   <- lm(cbind(bl, em, sf, bs) ~ zst + fff, data = pulp_data)
yp_zst_fff    <- predict(pup_zst_fff)

calc_rsquared(Yf, yp_zst_afl) 
calc_rsquared(Yf, yp_zst_lff) 
calc_rsquared(Yf, yp_zst_fff)

(calc_rsquared(Yf, yp_zst_afl) - calc_rsquared(Yf, yp_zst))/calc_rsquared(Yf, yp_zst)
(calc_rsquared(Yf, yp_zst_lff) - calc_rsquared(Yf, yp_zst))/calc_rsquared(Yf, yp_zst)
(calc_rsquared(Yf, yp_zst_fff) - calc_rsquared(Yf, yp_zst))/calc_rsquared(Yf, yp_zst) # highest percent change

#### three covariate model ####
pup_zst_fff_afl   <- lm(cbind(bl, em, sf, bs) ~ zst + fff + afl, data = pulp_data)
yp_zst_fff_afl    <- predict(pup_zst_fff_afl)

pup_zst_fff_lff   <- lm(cbind(bl, em, sf, bs) ~ zst + fff + lff, data = pulp_data)
yp_zst_fff_lff    <- predict(pup_zst_fff_lff)

calc_rsquared(Yf, yp_zst_fff_afl)
calc_rsquared(Yf, yp_zst_fff_lff)

(calc_rsquared(Yf, yp_zst_fff_afl) - calc_rsquared(Yf, yp_zst_fff))/calc_rsquared(Yf, yp_zst_fff)
(calc_rsquared(Yf, yp_zst_fff_lff) - calc_rsquared(Yf, yp_zst_fff))/calc_rsquared(Yf, yp_zst_fff) # highest percent change


####  Construct tables of results for export ####

### AMI Table ----
output_ami <- tibble(
  model = c(
    "Full model",
    
    "Single covariate: amt",
    
    "Two covariates: amt, sex",
    "Two covariates: amt, pr",
    "Two covariates: amt, diap",
    "Two covariates: amt, qrs",
    
    "Three covariates: amt, sex, pr",
    "Three covariates: amt, sex, diap",
    "Three covariates: amt, sex, qrs",
    
    "Four covariates: amt, sex, qrs, diap",
    "Four covariates: amt, sex, qrs, qrs"
  ),
  
  r2 = c(
    calc_rsquared(ya, ypaf),
    
    calc_rsquared(ya, yp_amt),
    
    calc_rsquared(ya, yp_am_sex),
    calc_rsquared(ya, yp_am_pr),
    calc_rsquared(ya, yp_am_di),
    calc_rsquared(ya, yp_am_qr),
    
    calc_rsquared(ya, yp_am_sex_pr),
    calc_rsquared(ya, yp_am_sex_di),
    calc_rsquared(ya, yp_am_sex_qr),
    
    calc_rsquared(ya, yp_am_sex_pr_di),
    calc_rsquared(ya, yp_am_sex_pr_qr)
  ),
  
  pct_change_from_single = c(
    NA,
    
    NA,
    
    (calc_rsquared(ya, yp_am_sex) - calc_rsquared(ya, yp_amt))/calc_rsquared(ya, yp_amt), # highest percent change
    (calc_rsquared(ya, yp_am_pr) - calc_rsquared(ya, yp_amt))/calc_rsquared(ya, yp_amt),
    (calc_rsquared(ya, yp_am_di) - calc_rsquared(ya, yp_amt))/calc_rsquared(ya, yp_amt),
    (calc_rsquared(ya, yp_am_qr) - calc_rsquared(ya, yp_amt))/calc_rsquared(ya, yp_amt),
    
    (calc_rsquared(ya, yp_am_sex_pr) - calc_rsquared(ya, yp_am_sex))/calc_rsquared(ya, yp_am_sex), # highest percent change
    (calc_rsquared(ya, yp_am_sex_di) - calc_rsquared(ya, yp_am_sex))/calc_rsquared(ya, yp_am_sex),
    (calc_rsquared(ya, yp_am_sex_qr) - calc_rsquared(ya, yp_am_sex))/calc_rsquared(ya, yp_am_sex),
    
    (calc_rsquared(ya, yp_am_sex_pr_di) - calc_rsquared(ya, yp_am_sex_pr))/calc_rsquared(ya, yp_am_sex_pr), # highest percent change
    (calc_rsquared(ya, yp_am_sex_pr_qr) - calc_rsquared(ya, yp_am_sex_pr))/calc_rsquared(ya, yp_am_sex_pr)
    
  )
)


### Pulp Table ----
output_pulp <- tibble(
  model = c(
    "Full model",
    
    "Single covariate: afl",
    "Single covariate: lff",
    "Single covariate: fff",
    "Single covariate: zst",
    
    "Two covariates: zst, afl",
    "Two covariates: zst, lff",
    "Two covariates: zst, fff",

    "Three covariates: zst, fff, afl",
    "Three covariates: zst, fff, lff"
  ),
  
  r2 = c(
    calc_rsquared(Yf, ypf),
    
    calc_rsquared(Yf, yp_afl),
    calc_rsquared(Yf, yp_lff),
    calc_rsquared(Yf, yp_fff),
    calc_rsquared(Yf, yp_zst), # best model
    
    calc_rsquared(Yf, yp_zst_afl), 
    calc_rsquared(Yf, yp_zst_lff), 
    calc_rsquared(Yf, yp_zst_fff),
    
    calc_rsquared(Yf, yp_zst_fff_afl),
    calc_rsquared(Yf, yp_zst_fff_lff)
  ),
  
  pct_change_from_single = c(
    NA,
    
    NA,
    NA,
    NA,
    NA,
    
    (calc_rsquared(Yf, yp_zst_afl) - calc_rsquared(Yf, yp_zst))/calc_rsquared(Yf, yp_zst),
    (calc_rsquared(Yf, yp_zst_lff) - calc_rsquared(Yf, yp_zst))/calc_rsquared(Yf, yp_zst),
    (calc_rsquared(Yf, yp_zst_fff) - calc_rsquared(Yf, yp_zst))/calc_rsquared(Yf, yp_zst), # highest percent change
    
    (calc_rsquared(Yf, yp_zst_fff_afl) - calc_rsquared(Yf, yp_zst_fff))/calc_rsquared(Yf, yp_zst_fff),
    (calc_rsquared(Yf, yp_zst_fff_lff) - calc_rsquared(Yf, yp_zst_fff))/calc_rsquared(Yf, yp_zst_fff) # highest percent change
  )
)

#### Save results tables ####

write_rds(
  output_ami,
  file = "data-derived/multivariate-ami.rds"
)

write_rds(
  output_pulp,
  file = "data-derived/multivariate-pulp.rds"
)

