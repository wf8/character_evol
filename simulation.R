#
# script to run various simulations of character evolution and
# then estimate ancestral states using brownian motion
# plus a few options for visually summarizing the results
#

# load the libraries and custom functions needed
source("setup.R")

# we will time the simulation
start_time <- Sys.time()

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
mu_d <- list(function(d) { rep(0.05 * d, length(d)) }, # extinction increase linearly with slope 0.05
             function(d) { rep(0.05 * d, length(d)) },
             function(d) { rep(0.05 * d, length(d)) },
             function(d) { rep(0.05 * d, length(d)) })

# character evolving with brownian motion, one function per regime
char <- list(make.brownian.with.drift(0, 1.0),
             make.brownian.with.drift(0, 1.0),
             make.brownian.with.drift(0, 1.0),
             make.brownian.with.drift(0, 1.0))

# a list to store each simulation's data in:
simulations <- list()

# run the simulations
for (i in 1:replicates) {

    # simulate tree and character
    sim_data <- tree.quasse.regimes(list(lambda, lambda_d, mu, mu_d, char), regimes=regimes, 
                                    max.t=max_time, x0=root_state, single.lineage=FALSE, include.extinct=FALSE)

    # get tip states and branching times from simulated data
    tip_states <- sim_data$tip.state
    branch_times <- as.vector(branching.times(sim_data))

    # reorder and extract the simulated ancestral states 
    sim_anc_states <- root_state
    for (j in (length(tip_states) + 1):max(sim_data$orig$idx2)) {

       sim_anc_states <- c(sim_anc_states, sim_data$orig$state[ sim_data$orig$idx2 == j ] )

    }

    # infer ancestral states
    est_data <- ace(tip_states, sim_data, type="continuous")
    est_anc_states <- as.vector(est_data$ace)

    # calculate contrasts between simulated and estimated ancestral states
    node_differences <- sim_anc_states - est_anc_states
    root_difference <- sim_anc_states[1] - est_anc_states[1]

    # save all the data for this simulation
    simulations[[i]] <- list(sim_data=sim_data, tip_states=tip_states, branch_times=branch_times, 
                             sim_anc_states=sim_anc_states, est_data=est_data, est_anc_states=est_anc_states, 
                             node_differences=node_differences, root_difference=root_difference)

}

# record the final time
processing_time <- Sys.time() - start_time

# various methods for summarizing the results:
par(mfrow=c(2,2))

# plot two traitgrams on top of one another
traitgram_sim(simulations[[1]]$tip_states, simulations[[1]]$sim_anc_states, simulations[[1]]$sim_data, xlab="trait")
traitgram_est(simulations[[1]]$tip_states, simulations[[1]]$est_anc_states, simulations[[1]]$sim_data)

# plot rescaled branch times against difference between simulated and estimated states
branch_times <- unlist( sapply(simulations, function(x){max(x$branch_times) - x$branch_times}) )
sim_anc_states <- unlist( sapply(simulations, function(x){x$sim_anc_states}) )
est_anc_states <- unlist( sapply(simulations, function(x){x$est_anc_states}) )
plot(branch_times, sim_anc_states[1:length(est_anc_states)] - est_anc_states, xlab="branching times", ylab="ancestral state differences")

# plot lineage through time curve
lineages <- attr(simulations[[1]]$sim_data$orig, "lineages_thru_time")
ages <- attr(simulations[[1]]$sim_data$orig, "ages")
plot(ages, lineages, type="l", xlab="time")

# view tree unbalance
plot(ladderize(simulations[[1]]$sim_data))
axisPhylo(backward=FALSE)

# plot the root differences for all simulations
#root_differences <- sapply(simulations, function(x){x$root_difference})
#hist(root_differences, breaks=20)


