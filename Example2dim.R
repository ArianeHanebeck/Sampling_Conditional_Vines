#########
# Setup #
#########

# Set working directory
setwd("")

# Libraries
library(rstan)
rstan_options(auto_write = TRUE)
library(rvinecopulib)
library(VineCopula)
library(kde1d)
library(plotly)
library(MASS)

# Loading files
source('Sample_Function.R') # Functions for sampling
source('Transformation.R') # Transforming a vine from VineCopula setup to rvinecopulib setup

# Loading Stan file (for multivariate conditional distributions)
STAN2 = stan_model(file = 'STAN2.stan')

##################################################
# 5 dim example - bivariate conditional sampling #
##################################################

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

# We are interested in the conditional density of 24|135
ucon = c(0.2, FALSE, 0.45, FALSE, 0.78)
indexcon = c(1, FALSE, 3, FALSE, 5)

# Sample from the conditional density 5000 times for each of the 4 chains
Sample5000 = sample_from_conditional(5000, STAN, STAN2, RVM, indexcon, ucon,  thin = 10, seed = 12345)

# Extract the samples of the second and fourth variable
s = extract(Sample5000,permute = TRUE)$ucalculate
s1 = s[,1]
s2 = s[,2]

############
# Plotting #
############

x <- seq(0,1,length.out = 200)
y <- seq(0,1,length.out = 200)

# Transformation to z-scale
# Based on section 7.4 Transformation Two-Dimensional Case in the paper
z1 = qnorm(s1)
z2 = qnorm(s2)

kdz1 <- kde2d(z1, z2, n = 200)
dx = dnorm(kdz1$x)
dy = dnorm(kdz1$y)
d = dx%*%t(dy)

kdz1$z = kdz1$z / d
kdz1$x = pnorm(kdz1$x)
kdz1$y = pnorm(kdz1$y)

# Plot the 3d plot based on the MCMC samples
fig1 <- plot_ly(x = kdz1$x, y = kdz1$y, z = kdz1$z, showscale = FALSE, scene = 'scene1') %>% add_surface()%>% 
  layout(scene1 = list(xaxis = list(title = list(text = "u<sub>2</sub>", font = list(size = 20))),
                       yaxis = list(title = list(text = "u<sub>4</sub>", font = list(size = 20))),
                       zaxis = list(title = list(text = "c(u<sub>2</sub>, u<sub>4</sub> | u<sub>1</sub>, u<sub>3</sub>, u<sub>5</sub>)", font = list(size = 20))),
                       camera = list(eye = list(x = 1.5,y = 1.5,z = 1.25))),
         title = list(text = paste0("u<sub>1</sub> = ", round(ucon[1],2), 
                                  ", u<sub>3</sub> = ", round(ucon[3],2), 
                                  ", u<sub>5</sub> = ", round(ucon[5],2)), y = 0.9,
                      font = list(size = 20)))
fig1

# Define the function DensBivInteg5 for comparison
# It obtains the conditional density of 24|135 via integration
DensBivInteg5 = function(vinedistr, u1_int, u3_int, u5_int){
  
  commondensity = function(u1_common, u2_common, u3_common, u4_common, u5_common){
    return(dvinecop(c(u1_common, u2_common, u3_common, u4_common, u5_common), vinedistr))
  }
  
  commondensityu24 = function(u2,u4)
  {
    return(commondensity(u1_int, u2, u3_int, u4, u5_int))
  }
  
  commondensityu24 = Vectorize(commondensityu24)
  
  conditionaldensity = function(u2_cond, u4_cond){
    return(commondensityu24(u2_cond, u4_cond)/integral2(commondensityu24, 0, 1, 0, 1)$Q)
  }
  
  return(conditionaldensity)
}

# Plot the 3d plot based on the integrated function
zz <- outer(x,y, function(x,y) DensBivInteg5(vinecopdist, ucon[1], 
                                             ucon[3], ucon[5])(x,y))

figt <- plot_ly(x = x,y = y, z = zz, scene = 'scene1') %>% add_surface()%>% 
  layout(scene1 = list(xaxis = list(title = list(text = "u<sub>2</sub>", font = list(size = 20))),
                       yaxis = list(title = list(text = "u<sub>4</sub>", font = list(size = 20))),
                       zaxis = list(title = list(text = "c(u<sub>2</sub>, u<sub>4</sub> | u<sub>1</sub>, u<sub>3</sub>, u<sub>5</sub>)", 
                                               font = list(size = 15))),
                       camera = list(eye = list(x = 1.5,y = 1.5,z = 1.25))),
         title = list(text = paste0("u<sub>1</sub> = ", round(ucon[1],2), 
                                  ", u<sub>3</sub> = ", round(ucon[3],2), 
                                  ", u<sub>5</sub> = ", round(ucon[5],2)), y = 0.9,
                      font = list(size = 20)))
figt

# Compare the contour plots
par(pty = "s")
contour(zz, drawlabels = FALSE, col = "blue", nlevels = 5, xlim = c(0,1), ylim = c(0,1))
contour(kdz1, add = TRUE, col = "black", drawlabels = FALSE, nlevels = 5)
legend("bottomleft", legend = c("integ.", "estim."), col = c("blue", "black"), lty = 1)
