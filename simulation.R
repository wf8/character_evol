#
# script to run various simulations of character evolution and
# then estimate ancestral states using brownian motion
#

# load the libraries and custom functions needed
source("setup.R")

# we will time the simulation
start_time <- Sys.time()

# setup simulation parameters starting with the tree height
max_time <- 10

# a vector to hold the times at which each regime ends
regimes <- c(10)

# a list of speciation functions for the regime shifts, one function per regime
#lambda <- list(function(x) noroptimal.x(x, y0=0, y1=30, xmid=0, s2=10))
lambda <- list(function(x) constant.x(x, 0.3))

# a list of extinction functions for the regime shifts, one function per regime
mu <- list(function(x) constant.x(x, 0.0))

# character evolving with brownian motion, one function per regime
char <- list(make.brownian.with.drift(0, 1.0))

# variables to store our simulated data in:
sim_data <- list()
tip_states <- list()
branch_times <- list()
sim_anc_states <- list()
rescaled_branch_time <- list()

# variables to store our estimated data in:
est_data <- list()
est_anc_states <- list()
node_differences <- list()
root_differences <- vector()

# run the simulation
replicates <- 1
for (i in 1:replicates) {

    # simulate tree and character
    sim_data[[i]] <- tree.quasse.regimes(list(lambda, mu, char), regimes=regimes, max.t=max_time, x0=0, single.lineage=FALSE)

    # get tip states and branching times from simulated data
    tip_states[[i]] <- sim_data[[i]]$tip.state
    branch_times[[i]] <- as.vector(branching.times(sim_data[[i]]))
    rescaled_branch_time[[i]] <- branch_times[[i]] / max(branch_times[[i]])

    # reorder and extract the simulated ancestral states 
    sim_anc_states[[i]] <- 0.0
    for (j in length(sim_data[[i]]$tip.state)+1:max(sim_data[[i]]$orig$idx2)) {

       sim_anc_states[[i]] <- c(sim_anc_states[[i]], sim_data[[i]]$orig$state[ sim_data[[i]]$orig$idx2==j ] )

    }

    # infer ancestral states
    est_data[[i]] <- ace(tip_states[[i]], sim_data[[i]], type="continuous")
    est_anc_states[[i]] <- as.vector(est_data[[i]]$ace)

    # calculate contrasts between simulated and estimated ancestral states
    node_differences[[i]] <- sim_anc_states[[i]] - est_anc_states[[i]]
    root_differences[i] <- sim_anc_states[[i]][1] - est_anc_states[[i]][1]
}

# record the final time
processing_time <- Sys.time() - start_time

#hist(tree_depths, breaks=20)
#hist(root_differences)

# plot rescaled branch times against difference between simulated and estimated states
#plot(rescaled_branch_time, sim_anc_states - est_anc_states)

# plot the two traitgrams on top of one another
traitgram_sim(tip_states[[1]], sim_anc_states[[1]], sim_data[[1]])
traitgram_est(tip_states[[1]], sim_data[[1]])


# plot the ancestral states as circles on the nodes
#plot(sim_data, label.offset=.05)
#nodelabels()
#tiplabels()
#tiplabels(pch=21, cex=tip_states, bg="yellow")
#nodelabels(pch=21, cex= anc_states$ace, bg="yellow")
