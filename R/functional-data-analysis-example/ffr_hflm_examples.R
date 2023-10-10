#### load libraries ####
library(tidyverse)
library(mvrsquared)
library(refund)

#### load data ####
emg <- read_table("data-raw/functional-data-analysis-example/EMG.dat", col_names = FALSE) |>
  as.matrix() |>
  t()

acc <- read_table("data-raw/functional-data-analysis-example/LipAcc.dat", col_names = FALSE) |>
  as.matrix() |>
  t()

pos <- read_table("data-raw/functional-data-analysis-example/LipPos.dat", col_names = FALSE) |>
  as.matrix() |>
  t()

#### fit penalized FFR ####
modFFR  <- pffr(pos ~ ff(acc, xind = 1:ncol(acc), basistype = 'te',
                         splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                           k = c(10, 10))), yind = 1:ncol(pos))
Ypf     <- predict(modFFR)

Ypf_r2 <- calc_rsquared(pos, Ypf)


#### fit penalized HFLM ####
modHFP  <- pffr(pos ~ ff(emg, xind = 1:ncol(emg), basistype = 'te', limits = "s<=t",
                         splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                           k = c(10, 10))), yind = 1:ncol(pos))
Yp      <- predict(modHFP)

Yp_r2 <- calc_rsquared(pos, Yp)


#### fit penalized FFR+HFLM ####
modFH   <- pffr(pos ~ ff(acc, xind = 1:ncol(acc), basistype = 'te',
                         splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                           k = c(10, 10))) +
                      ff(emg, xind = 1:ncol(emg), basistype = 'te', limits = "s<=t",
                         splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                           k = c(10, 10))), 
                yind = 1:ncol(pos))

Ypfh    <- predict(modFH)

Ypfh_r2 <- calc_rsquared(pos, Ypfh)


#### Create table for use in paper ####
ffr_hflm_r2 <- 
  tibble(
    name = c("FFR", "HFP", "FH"),
    r2 = c(Ypf_r2, Yp_r2, Ypfh_r2)
  )

write_rds(
  ffr_hflm_r2,
  "data-derived/functional-data-analysis-r2.rds"
)

