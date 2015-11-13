#
# script to run various simulations of character evolution and
# then estimate ancestral states using brownian motion
# plus a few options for visually summarizing the results
#

# load the libraries and custom functions needed
source("setup.R")

# setup simulation parameters
replicates <- 1
max_time <- 10
root_state <- 0.0

# a vector to hold the times at which each regime ends
regimes <- c(2.5, 5, 7.5, 10)

# for each regime, speciation and extinction are each calculated as the sum of two functions: 
# 1) a function of the character value x
# 2) a function of the number of lineages (diversity dependence)

# speciation as a function of the character value x
lambda <- list(function(x) noroptimal.x(x, y0=0, y1=2.0, xmid=0, s2=1),
               function(x) noroptimal.x(x, y0=0, y1=2.0, xmid=3, s2=1),
               function(x) noroptimal.x(x, y0=0, y1=2.0, xmid=6, s2=1),
               function(x) noroptimal.x(x, y0=0, y1=2.0, xmid=9, s2=1))
#lambda <- list(function(x) constant.x(x, 0.3))

# speciation as a function of the number of lineages d
lambda_d <- list(function(d) constant.x(d, 0.0), 
                 function(d) constant.x(d, 0.0),
                 function(d) constant.x(d, 0.0),
                 function(d) constant.x(d, 0.0))

# extinction as a function of the character value x
mu <- list(function(x) constant.x(x, 0.0), 
           function(x) constant.x(x, 0.0),
           function(x) constant.x(x, 0.0),
           function(x) constant.x(x, 0.0))

# extinction as a function of the number of lineages d
mu_d <- list(function(d) { rep(0.05 * (d-1), length(d)) }, # extinction increase linearly with slope 0.05
             function(d) { rep(0.05 * (d-1), length(d)) },
             function(d) { rep(0.05 * (d-1), length(d)) },
             function(d) { rep(0.05 * (d-1), length(d)) })

# character evolving with brownian motion, one function per regime
char <- list(make.brownian.with.drift(0, 1.0),
             make.brownian.with.drift(0, 1.0),
             make.brownian.with.drift(0, 1.0),
             make.brownian.with.drift(0, 1.0))

# run the simulations:
simulations <- simulate_trees(replicates, lambda, lambda_d, mu, mu_d, char, regimes, max_time, root_state, include_extinct=FALSE) 

# plot some results
plot_simulations(replicates, simulations)


