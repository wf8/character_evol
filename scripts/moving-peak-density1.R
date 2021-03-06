#
# script to run various simulations of character evolution and
# then estimate ancestral states using brownian motion
# plus a few options for visually summarizing the results
#

# setup simulation parameters
replicates <- 100
max_time <- 10
root_state <- 0.0
output_file <- "moving-peak-density1"

# parameters values to simulate under
birth <- c(0.0, 0.1, 0.2, 0.3, 0.2, 0.3, 0.4, 0.3, 0.4, 0.6, 0.0, 0.1, 0.2, 0.3, 0.2, 0.3, 0.4, 0.3, 0.4, 0.6, 0.0, 0.1, 0.2, 0.3, 0.2, 0.3, 0.4, 0.3, 0.4, 0.6)
death <- c(0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3)
beta <-  c(.01, .01, .01, .01, .01, .01, .01, .01, .01, .01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

peak = c(0.5, 1.0, 2.0)
sigma = c(0.1, 1.0, 10)
move = c(0.5, 1, 2) # final position of peak after 10 moves: 5, 10, 20

# a vector to hold the times at which each regime ends
regimes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

results <- list()
iter <- 1

for (i in 1:length(birth)) {

    for (j in 1:length(peak)) {
    
        for (k in 1:length(sigma)) {
        
            for (l in 1:length(move)) {

                print(paste("Simulating under parameter set:", iter, "out of", length(birth) * length(peak) * length(sigma) * length(move)))

                # for each regime, speciation and extinction are each calculated as the sum of two functions: 
                # 1) a function of the character value x
                # 2) a function of the number of lineages (diversity dependence)

                # speciation as a function of the character value x
                lambda <- list(function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=0, s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=1*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=2*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=3*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=4*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=5*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=6*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=7*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=8*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=9*move[l], s2=sigma[k]),
                               function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=10*move[l], s2=sigma[k]))

                # speciation as a function of the number of lineages d
                lambda_d <- list(function(d) constant.x(d, 0.0),
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0)) 

                # extinction as a function of the character value x
                mu <- list(function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]),
                           function(x) constant.x(x, death[i]))

                # extinction as a function of the number of lineages d
                mu_d <- list(function(d) { rep(0.1 * (d-1), length(d)) }, # extinction increases linearly with slope 0.1
                             function(d) { rep(0.1 * (d-1), length(d)) },
                             function(d) { rep(0.1 * (d-1), length(d)) },
                             function(d) { rep(0.1 * (d-1), length(d)) },
                             function(d) { rep(0.1 * (d-1), length(d)) },
                             function(d) { rep(0.1 * (d-1), length(d)) },
                             function(d) { rep(0.1 * (d-1), length(d)) },
                             function(d) { rep(0.1 * (d-1), length(d)) },
                             function(d) { rep(0.1 * (d-1), length(d)) },
                             function(d) { rep(0.1 * (d-1), length(d)) })

                # character evolving with brownian motion, one function per regime
                char <- list(make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]))

                # run the simulations:
                results[[iter]] <- simulate_trees(replicates, lambda, lambda_d, mu, mu_d, char, regimes, max_time, root_state, include_extinct=FALSE)
                iter <- iter + 1

                # plot some results
                #plot_simulations(replicates, simulations, est_type="lp")

                # save the results
                save(results, file=paste("output/", output_file, ".RData", sep=""))

            }

        }

    }
}
