#### load libraries ####
library(mvrsquared)

#### load data ####

##### amitriptyline #####
# source: M. V. Rudorfer (1982). Cardiovascular Changes and Plasma Drug Levels after
#         Amitriptyline Overdose. Journal of Toxicology-Clinical Toxicology, 19, 67-71.
# secondary source: R. A. Johnson and D. W. Wichern (2007). Applied Multivariate Statistical Analysis
#                   6th Ed. Pearson Prentice Hall, Upper Saddle River, NJ, USA.
# amit          <- read.table('./Data/T7-6.DAT', header = FALSE)
amiData         <- read.table("http://static.lib.virginia.edu/statlab/materials/data/ami_data.DAT")
names(amiData)  <- c('tot', 'ami', 'sex', 'amt', 'pr', 'diap', 'qrs')
# quick description: data containing 17 overdoses of the drug amitriptyline
# outcomes: total TCAD plasma level (tot), amount of amitriptyline present in TCAD plasma levels (ami)
# covariates: sex (1 for female, 0 for male), amount of antidepressants taken at time of overdose (amt)
#             PR wave measurement (pr), diastolic blood pressure (diap), QRS wave measurement (qrs)

##### pulp and paper properties #####
# primary source: R. A. Johnson and D. W. Wichern (2007). Applied Multivariate Statistical Analysis
#                 6th Ed. Pearson Prentice Hall, Upper Sadle River, NJ, USA.
pulpData          <- read.table('./Data/T7-7.DAT', header = FALSE)
names(pulpData)   <- c('bl', 'em', 'sf', 'bs', 'afl', 'lff', 'fff', 'zst')
# quick description: Measurements of properties of pulp fibers and the paper made from them
# outcomes: breaking length (bl), elastic modulus (em), stress at failure (sf), burst strength (bs)
# covariates: arithmetic fiber length (afl), long fiber fraction (lff), fine fiber fraction (fff), zero span tensile (zst)



#### Amitriptyline Analysis ####
Ya            <- with(amiData, cbind(tot, ami))

##### full model #####
ami_full    <- lm(cbind(tot, ami) ~ amt + sex + pr + diap + qrs, data = amiData)
Ypaf        <- predict(ami_full)
calc_rsquared(Ya, Ypaf)

##### single covariate #####
ami_amt     <- lm(cbind(tot, ami) ~ amt, data = amiData)
Yp_amt      <- predict(ami_amt)

calc_rsquared(Ya, Yp_amt)

#### two covariates ####
ami_am_sex   <- lm(cbind(tot, ami) ~ amt + sex, data = amiData)
Yp_am_sex    <- predict(ami_am_sex)

ami_am_pr   <- lm(cbind(tot, ami) ~ amt + pr, data = amiData)
Yp_am_pr    <- predict(ami_am_pr)

ami_am_di   <- lm(cbind(tot, ami) ~ amt + diap, data = amiData)
Yp_am_di    <- predict(ami_am_di)

ami_am_qr   <- lm(cbind(tot, ami) ~ amt + qrs, data = amiData)
Yp_am_qr    <- predict(ami_am_qr)

calc_rsquared(Ya, Yp_am_sex)
calc_rsquared(Ya, Yp_am_pr)
calc_rsquared(Ya, Yp_am_di)
calc_rsquared(Ya, Yp_am_qr)

## percent change from one covariate #
(calc_rsquared(Ya, Yp_am_sex) - calc_rsquared(Ya, Yp_amt))/calc_rsquared(Ya, Yp_amt) # highest percent change
(calc_rsquared(Ya, Yp_am_pr) - calc_rsquared(Ya, Yp_amt))/calc_rsquared(Ya, Yp_amt)
(calc_rsquared(Ya, Yp_am_di) - calc_rsquared(Ya, Yp_amt))/calc_rsquared(Ya, Yp_amt)
(calc_rsquared(Ya, Yp_am_qr) - calc_rsquared(Ya, Yp_amt))/calc_rsquared(Ya, Yp_amt)

#### three covariates ####
ami_am_sex_pr   <- lm(cbind(tot, ami) ~ amt + sex + pr, data = amiData)
Yp_am_sex_pr    <- predict(ami_am_sex_pr)

ami_am_sex_di   <- lm(cbind(tot, ami) ~ amt + sex + diap, data = amiData)
Yp_am_sex_di    <- predict(ami_am_sex_di)

ami_am_sex_qr   <- lm(cbind(tot, ami) ~ amt + sex + qrs, data = amiData)
Yp_am_sex_qr    <- predict(ami_am_sex_qr)

calc_rsquared(Ya, Yp_am_sex_pr)
calc_rsquared(Ya, Yp_am_sex_di)
calc_rsquared(Ya, Yp_am_sex_qr)

## percent change from one covariate #
(calc_rsquared(Ya, Yp_am_sex_pr) - calc_rsquared(Ya, Yp_am_sex))/calc_rsquared(Ya, Yp_am_sex) # highest percent change
(calc_rsquared(Ya, Yp_am_sex_di) - calc_rsquared(Ya, Yp_am_sex))/calc_rsquared(Ya, Yp_am_sex)
(calc_rsquared(Ya, Yp_am_sex_qr) - calc_rsquared(Ya, Yp_am_sex))/calc_rsquared(Ya, Yp_am_sex)

#### four covariate model ####
ami_am_sex_pr_di  <- lm(cbind(tot, ami) ~ amt + sex + pr + diap, data = amiData)
Yp_am_sex_pr_di   <- predict(ami_am_sex_pr_di)

ami_am_sex_pr_qr  <- lm(cbind(tot, ami) ~ amt + sex + pr + qrs, data = amiData)
Yp_am_sex_pr_qr   <- predict(ami_am_sex_pr_qr)

calc_rsquared(Ya, Yp_am_sex_pr_di)
calc_rsquared(Ya, Yp_am_sex_pr_qr)

## percent change from one covariate #
(calc_rsquared(Ya, Yp_am_sex_pr_di) - calc_rsquared(Ya, Yp_am_sex_pr))/calc_rsquared(Ya, Yp_am_sex_pr) # highest percent change
(calc_rsquared(Ya, Yp_am_sex_pr_qr) - calc_rsquared(Ya, Yp_am_sex_pr))/calc_rsquared(Ya, Yp_am_sex_pr)


#### Pulp Analysis ####
Yf        <- with(pulpData, cbind(bl, em, sf, bs))

pup_full  <- lm(cbind(bl, em, sf, bs) ~ afl + lff + fff + zst, data = pulpData)
Ypf       <- predict(pup_full)

calc_rsquared(Yf, Ypf)


#### single covariate ####
pup_afl   <- lm(cbind(bl, em, sf, bs) ~ afl, data = pulpData)
Yp_afl    <- predict(pup_afl)

pup_lff   <- lm(cbind(bl, em, sf, bs) ~ lff, data = pulpData)
Yp_lff    <- predict(pup_lff)

pup_fff   <- lm(cbind(bl, em, sf, bs) ~ fff, data = pulpData)
Yp_fff    <- predict(pup_fff)

pup_zst   <- lm(cbind(bl, em, sf, bs) ~ zst, data = pulpData)
Yp_zst    <- predict(pup_zst)

calc_rsquared(Yf, Yp_afl)
calc_rsquared(Yf, Yp_lff)
calc_rsquared(Yf, Yp_fff)
calc_rsquared(Yf, Yp_zst) # best model

#### two covariate model ####
pup_zst_afl   <- lm(cbind(bl, em, sf, bs) ~ zst + afl, data = pulpData)
Yp_zst_afl    <- predict(pup_zst_afl)

pup_zst_lff   <- lm(cbind(bl, em, sf, bs) ~ zst + lff, data = pulpData)
Yp_zst_lff    <- predict(pup_zst_lff)

pup_zst_fff   <- lm(cbind(bl, em, sf, bs) ~ zst + fff, data = pulpData)
Yp_zst_fff    <- predict(pup_zst_fff)

calc_rsquared(Yf, Yp_zst_afl) 
calc_rsquared(Yf, Yp_zst_lff) 
calc_rsquared(Yf, Yp_zst_fff)

(calc_rsquared(Yf, Yp_zst_afl) - calc_rsquared(Yf, Yp_zst))/calc_rsquared(Yf, Yp_zst)
(calc_rsquared(Yf, Yp_zst_lff) - calc_rsquared(Yf, Yp_zst))/calc_rsquared(Yf, Yp_zst)
(calc_rsquared(Yf, Yp_zst_fff) - calc_rsquared(Yf, Yp_zst))/calc_rsquared(Yf, Yp_zst) # highest percent change

#### three covariate model ####
pup_zst_fff_afl   <- lm(cbind(bl, em, sf, bs) ~ zst + fff + afl, data = pulpData)
Yp_zst_fff_afl    <- predict(pup_zst_fff_afl)

pup_zst_fff_lff   <- lm(cbind(bl, em, sf, bs) ~ zst + fff + lff, data = pulpData)
Yp_zst_fff_lff    <- predict(pup_zst_fff_lff)

calc_rsquared(Yf, Yp_zst_fff_afl)
calc_rsquared(Yf, Yp_zst_fff_lff)

(calc_rsquared(Yf, Yp_zst_fff_afl) - calc_rsquared(Yf, Yp_zst_fff))/calc_rsquared(Yf, Yp_zst_fff)
(calc_rsquared(Yf, Yp_zst_fff_lff) - calc_rsquared(Yf, Yp_zst_fff))/calc_rsquared(Yf, Yp_zst_fff) # highest percent change

