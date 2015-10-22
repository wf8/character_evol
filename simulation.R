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
replicates <- 4
max_time <- 10
root_state <- 0.0

# a vector to hold the times at which each regime ends
regimes <- c(10)

# a list of speciation functions for the regime shifts, one function per regime
#lambda <- list(function(x) noroptimal.x(x, y0=0, y1=30, xmid=0, s2=10))
lambda <- list(function(x) constant.x(x, 0.3))

# a list of extinction functions for the regime shifts, one function per regime
mu <- list(function(x) constant.x(x, 0.0))

# character evolving with brownian motion, one function per regime
char <- list(make.brownian.with.drift(0, 1.0))

# a list to store each simulation's data in:
simulations <- list()

# run the simulations
for (i in 1:replicates) {

    # simulate tree and character
    sim_data <- tree.quasse.regimes(list(lambda, mu, char), regimes=regimes, 
                                    max.t=max_time, x0=root_state, single.lineage=FALSE)

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
                             rescaled_branch_times=rescaled_branch_times, sim_anc_states=sim_anc_states, 
                             est_data=est_data, est_anc_states=est_anc_states, 
                             node_differences=node_differences, root_difference=root_difference)

}

# record the final time
processing_time <- Sys.time() - start_time

# various methods for summarizing the results:

# plot two traitgrams on top of one another
traitgram_sim(simulations[[1]]$tip_states, simulations[[1]]$sim_anc_states, simulations[[1]]$sim_data)
traitgram_est(simulations[[1]]$tip_states, simulations[[1]]$est_anc_states, simulations[[1]]$sim_data)

# plot the root differences for all simulations
#root_differences <- sapply(simulations, function(x){x$root_difference})
#hist(root_differences, breaks=20)

# plot rescaled branch times against difference between simulated and estimated states
#branch_times <- unlist( sapply(simulations, function(x){x$branch_times}) )
#sim_anc_states <- unlist( sapply(simulations, function(x){x$sim_anc_states}) )
#est_anc_states <- unlist( sapply(simulations, function(x){x$est_anc_states}) )
#plot(branch_times, sim_anc_states - est_anc_states)

