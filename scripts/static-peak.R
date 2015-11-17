#
# script to run various simulations of character evolution and
# then estimate ancestral states using brownian motion
# plus a few options for visually summarizing the results
#

# setup simulation parameters
replicates <- 100
max_time <- 10
root_state <- 0.0
output_file <- "static-peak"

# parameters values to simulate under
birth <- c(0.0, 0.1, 0.2, 0.3, 0.2, 0.3, 0.4, 0.3, 0.4, 0.6, 0.0, 0.1, 0.2, 0.3, 0.2, 0.3, 0.4, 0.3, 0.4, 0.6, 0.0, 0.1, 0.2, 0.3, 0.2, 0.3, 0.4, 0.3, 0.4, 0.6)
death <- c(0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3)
beta <-  c(.01, .01, .01, .01, .01, .01, .01, .01, .01, .01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

peak = c(0.5, 1.0, 2.0)
sigma = c(0.1, 1.0, 10)

# a vector to hold the times at which each regime ends
regimes <- c(10)

results <- list()

for (i in 1:length(birth)) {

    for (j in 1:length(peak)) {
    
        for (k in 1:length(sigma)) {

            print(paste("Simulating under parameter set:", i * j * k, "out of", length(birth) * length(peak) * length(sigma)))

            # for each regime, speciation and extinction are each calculated as the sum of two functions: 
            # 1) a function of the character value x
            # 2) a function of the number of lineages (diversity dependence)

            # speciation as a function of the character value x
            lambda <- list(function(x) noroptimal.x(x, y0=birth[i], y1=peak[j], xmid=0, s2=sigma[k]))

            # speciation as a function of the number of lineages d
            lambda_d <- list(function(d) constant.x(d, 0.0)) 

            # extinction as a function of the character value x
            mu <- list(function(x) constant.x(x, death[i]))

            # extinction as a function of the number of lineages d
            mu_d <- list(function(d) constant.x(d, 0.0))

            # character evolving with brownian motion, one function per regime
            char <- list(make.brownian.with.drift(0, beta[i]))

            # run the simulations:
            results[[i]] <- simulate_trees(replicates, lambda, lambda_d, mu, mu_d, char, regimes, max_time, root_state, include_extinct=FALSE)

            # plot some results
            #plot_simulations(replicates, simulations, est_type="lp")

            # save the results
            save(results, file=paste("output/", output_file, ".RData", sep=""))

        }

    }
}
