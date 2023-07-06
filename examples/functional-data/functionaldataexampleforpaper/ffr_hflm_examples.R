#### load libraries ####
library(mvrsquared)
library(refund)

#### load data ####
emg   <- t(as.matrix(read.table('EMG.dat', header = FALSE)))
acc   <- t(as.matrix(read.table('LipAcc.dat', header = FALSE)))
pos   <- t(as.matrix(read.table('LipPos.dat', header = FALSE)))


#### fit penalized FFR ####
modFFR  <- pffr(pos ~ ff(acc, xind = 1:ncol(acc), basistype = 'te',
                         splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                           k = c(10, 10))), yind = 1:ncol(pos))
Ypf     <- predict(modFFR)
calc_rsquared(pos, Ypf)


#### fit penalized HFLM ####
modHFP  <- pffr(pos ~ ff(emg, xind = 1:ncol(emg), basistype = 'te', limits = "s<=t",
                         splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                           k = c(10, 10))), yind = 1:ncol(pos))
Yp      <- predict(modHFP)
calc_rsquared(pos, Yp)


#### fit penalized FFR+HFLM ####
modFH   <- pffr(pos ~ ff(acc, xind = 1:ncol(acc), basistype = 'te',
                         splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                           k = c(10, 10))) +
                      ff(emg, xind = 1:ncol(emg), basistype = 'te', limits = "s<=t",
                         splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                           k = c(10, 10))), 
                yind = 1:ncol(pos))
Ypfh    <- predict(modFH)
calc_rsquared(pos, Ypfh)
