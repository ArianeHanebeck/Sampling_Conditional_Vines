#########
# Setup #
#########

# Set working directory
setwd("")

# Libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rvinecopulib)
library(VineCopula)
library(kde1d)

# Loading files
source('Sample_Function.R') # Functions for sampling
source('Transformation.R') # Transforming a vine from VineCopula setup to rvinecopulib setup

# Loading Stan file (for univariate conditional distributions)
STAN = stan_model(file='STAN.stan')

###################################################
# 5 dim example - univariate conditional sampling #
###################################################

# Specify the vine of interest
bicop11 = bicop_dist("bb1", 0, c(1,2))
bicop12 = bicop_dist("gu", 0, 3)
bicop13 = bicop_dist("gaus",0, 0.9)
bicop14 = bicop_dist("clay",0, 3)

bicop21 = bicop_dist("clay",0, 1.4)
bicop22 = bicop_dist("frank",0, 2.8)
bicop23 = bicop_dist("bb7", 0, c(1,1.5))

bicop31 = bicop_dist("frank", 0, 1.8)
bicop32 = bicop_dist("gaus",0, 0.3)

bicop41 = bicop_dist("gaus", 0, 0.05)

vinecopdist = vinecop_dist(list(list(bicop11, bicop12, bicop13, bicop14), 
                              list(bicop21, bicop22, bicop23),
                              list(bicop31, bicop32),
                              list(bicop41)), 
                           matrix(c(2,3,4,5,1,
                                    3,4,5,2,0,
                                    4,5,3,0,0,
                                    5,4,0,0,0,
                                    5,0,0,0,0),5,5))

# If vine is defined in rvinecopulib setup, transform to VineCopula setup
RVM = transformation_to_RVM_all(vinecopdist)

# We are interested in the conditional density of 2|1345
ucon = c(0.2, FALSE, 0.45, 0.3, 0.78) # Define the conditioning values
indexcon = c(1, FALSE, 3, 4, 5) # Set the index to false in the conditioned values

# Sample from the conditional density 1000/5000/10000 times for each of the 4 chains
Sample1000 = sample_from_conditional(1000, STAN, STAN2, RVM, indexcon, ucon, thin=10, seed=555)
Sample5000 = sample_from_conditional(5000, STAN, STAN2, RVM, indexcon, ucon, thin=10, seed=12345)
Sample10000 = sample_from_conditional(10000, STAN, STAN2, RVM, indexcon, ucon, thin=10, seed=321)

# Extract the samples of the second variable
s1 = extract(Sample1000,permute=TRUE)$ucalculate
s5 = extract(Sample5000,permute=TRUE)$ucalculate
s10 = extract(Sample10000,permute=TRUE)$ucalculate

############
# Plotting #
############

x = seq(0,1,0.005)

# Transformation to z-scale
z1000 = qnorm(s1)
kde1000 = kde1d(z1000)

z5000 = qnorm(s5)
kde5000 = kde1d(z5000)

z10000 = qnorm(s10)
kde10000 = kde1d(z10000)

# Define the function DensIntegrate5 for comparison
# It obtains the conditional density of 2|1345 via integration
DensIntegrate5 = function(vinedistr, u1_int, u3_int,  u4_int,  u5_int){
  
  commondensity = function(u1_common, u2_common, u3_common, u4_common, u5_common){
    return(dvinecop(c(u1_common, u2_common, u3_common, u4_common, u5_common), vinedistr))
  }
  
  commondensityu2 = function(u2)
  {
    return(commondensity(u1_int, u2, u3_int, u4_int, u5_int))
  }
  
  commondensityu2 = Vectorize(commondensityu2)
  
  conditionaldensity = function(u2_cond){
    return(commondensityu2(u2_cond)/integrate(commondensityu2, 0, 1)$value)
  }
  
  return(conditionaldensity)
}

# Plot to compare the results
# Transformation based on Section 7.1 Transformation One-Dimensional Case in the paper
plot(x, DensIntegrate5(vinecopdist, ucon[1], ucon[3], ucon[4], ucon[5])(x), col = "black", 
     type = "l", lwd = 2, 
     xlab = bquote(u[2]),  ylab = "",
     sub = bquote(u[1]==.(round(ucon[1],2)) ~","~ u[3]==.(round(ucon[3],2))
                  ~","~ u[4]==.(round(ucon[4],2))
                  ~","~ u[5]==.(round(ucon[5],2))),
     cex.sub = 1.4, cex.lab = 1.4)
title(ylab = bquote(c(u[2] ~"|"~ u[1] ~","~ u[3] ~","~ u[4] ~","~ u[5])), line = 2.3, 
      cex.lab = 1.4)

lines(x, dkde1d(qnorm(x), kde1000)/dnorm(qnorm(x)), col = "red", lty = 2)
lines(x, dkde1d(qnorm(x), kde5000)/dnorm(qnorm(x)), col = "green", lty = 3)
lines(x, dkde1d(qnorm(x), kde10000)/dnorm(qnorm(x)), col = "blue", lty = 4)

legend("topright", legend = c("True", "n = 1000", "n = 5000", "n = 10000"), 
       col = c("black", "red", "green", "blue"), lty = 1:4, cex = 1.1)

