#
# This file loads the libraries and custom functions necessary to
# run simulations of character evolution and phylogeny under a moving
# selective regime.
#

library(diversitree)
library(phytools)
library(picante)


# Main function that runs simulations and saves results.
simulate_trees <- function(replicates, lambda, lambda_d, mu, mu_d, char, regimes, max_time, root_state, include_extinct=FALSE) {

    # progress bar
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0.0)

    # we will time the simulation
    start_time <- Sys.time()

    # a list to store each simulation's data in:
    simulations <- list()

    # run the simulations
    for (i in 1:replicates) {
    
        # simulate tree and character
        sim_tree <- tree.quasse.regimes(list(lambda, lambda_d, mu, mu_d, char), regimes=regimes,
                                        max.t=max_time, x0=root_state, single.lineage=TRUE, include.extinct=include_extinct)

        # check if all lineages went extinct
        if (is.null(sim_tree)) {

            simulations$data[[i]] <- list("extinct")
        
        } else {

            # remember the finishing time of this simulation
            finishing_time <- max(attr(sim_tree$orig, "ages"))

            # get tip states and branching times from simulated data
            tip_states <- sim_tree$tip.state
            branch_times <- c(finishing_time, as.vector(branching.times(sim_tree)))

            # reorder and extract the simulated ancestral states
            sim_anc_states <- root_state
            #sim_anc_states <- vector() # this will not include the root state
            for (j in (length(tip_states) + 1):max(sim_tree$orig$idx2)) {

               sim_anc_states <- c(sim_anc_states, sim_tree$orig$state[ sim_tree$orig$idx2 == j ] )

            }

            # infer ancestral states using
            # Brownian motion (BM) model fitted by residual maximum likelihood
            est_data_bm <- ace(tip_states, sim_tree, type="continuous")
            est_root_state_bm <- as.vector(est_data_bm$ace)[1]
            est_anc_states_bm <- c(est_root_state_bm, as.vector(est_data_bm$ace))

            # infer ancestral states using
            # linear parsimony (LP)
            est_data_lp <- ace_lp(tip_states, sim_tree)
            est_root_state_lp <- est_data_lp[1]
            est_anc_states_lp <- c(est_root_state_lp, est_data_lp)

            # calculate contrasts between simulated and BM estimated ancestral states
            node_differences_bm <- sim_anc_states[1:length(est_anc_states_bm)] - est_anc_states_bm
            root_difference_bm <- sim_anc_states[1] - est_anc_states_bm[1]
            
            # calculate contrasts between simulated and LP estimated ancestral states
            node_differences_lp <- sim_anc_states[1:length(est_anc_states_lp)] - est_anc_states_lp
            root_difference_lp <- sim_anc_states[1] - est_anc_states_lp[1]

            # save all the data for this simulation
            simulations$data[[i]] <- list(sim_tree=sim_tree, tip_states=tip_states, branch_times=branch_times,
                                          sim_anc_states=sim_anc_states, est_data_bm=est_data_bm, est_anc_states_bm=est_anc_states_bm,
                                          node_differences_bm=node_differences_bm, root_difference_bm=root_difference_bm, 
                                          est_data_lp=est_data_lp, est_anc_states_lp=est_anc_states_lp,
                                          node_differences_lp=node_differences_lp, root_difference_lp=root_difference_lp, 
                                          finishing_time=finishing_time)

        }

        # update progress bar
        setTxtProgressBar(pb, i/replicates)

    }

    # record the final time
    simulations$processing_time <- Sys.time() - start_time
    
    # remeber the simulations parameters
    simulations$replicates <- replicates 
    simulations$lambda <- lambda
    simulations$lambda_d <- lambda_d
    simulations$mu <- mu 
    simulations$mu_d <- mu_d 
    simulations$char <- char 
    simulations$regimes <- regimes
    simulations$max_time <- max_time
    simulations$root_state <- root_state

    writeLines("\n")

    simulations
}


# Plots various summaries of the results.
plot_simulations <- function(replicates, simulations, est_type="bm") {

    if (est_type == "bm") {
        est_anc_states = simulations[[1]]$est_anc_states_bm
        all_est_anc_states <- unlist( sapply(simulations[1:replicates], function(x){x$est_anc_states_bm}) )
    } else {
        est_anc_states = simulations[[1]]$est_anc_states_lp
        all_est_anc_states <- unlist( sapply(simulations[1:replicates], function(x){x$est_anc_states_lp}) )
    }

    par(mfrow=c(2,2))

    # plot two traitgrams on top of one another
    traitgram_given(simulations[[1]]$tip_states, simulations[[1]]$sim_anc_states, simulations[[1]]$sim_tree, simulations[[1]]$finishing_time, lab="trait", method="sim")
    traitgram_given(simulations[[1]]$tip_states, simulations[[1]]$est_anc_states_bm, simulations[[1]]$sim_tree, simulations[[1]]$finishing_time, color="red", method="est")
    traitgram_given(simulations[[1]]$tip_states, simulations[[1]]$est_anc_states_lp, simulations[[1]]$sim_tree, simulations[[1]]$finishing_time, color="orange", method="est")

    # plot rescaled branch times against difference between simulated and estimated states
    branch_times <- unlist( sapply(simulations[1:replicates], function(x){x$finishing_time - x$branch_times}) )
    sim_anc_states <- unlist( sapply(simulations[1:replicates], function(x){x$sim_anc_states}) )
    state_differences <- sim_anc_states[1:length(all_est_anc_states)] - all_est_anc_states
    plot(branch_times, state_differences, xlab="time", ylab="ancestral state differences")

    # plot lineage through time curve
    lineages <- attr(simulations[[1]]$sim_tree$orig, "lineages_thru_time")
    ages <- attr(simulations[[1]]$sim_tree$orig, "ages")
    plot(ages, lineages, type="l", xlab="time")

    # view tree unbalance
    plot(ladderize(simulations[[1]]$sim_tree), root.edge=TRUE)
    axisPhylo(backward=FALSE, root.time=round(simulations[[1]]$sim_tree$root.edge, 2))

    # plot the root differences for all simulations
    #root_differences <- sapply(simulations, function(x){x$root_difference})
    #hist(root_differences, breaks=20)

}


# Calculates ancestral states using linear parsimony as
# described in Swofford and Maddison 1987
ace_lp <- function(tip_states, tree) {

    # 1: pass down the tree and get the state sets (Farris intervals) for each internal node
    tree <- reorder(tree, "postorder")
    state_sets = list()
    for (i in 1:length(tree$edge[,1])) {
        node_id <- tree$edge[i,2]
        if (node_id <= length(tip_states))
            state_sets[[ node_id ]] <- c(tip_states[ node_id ])
        else {
            desc_id <- tree$edge[,2][ tree$edge[,1] == node_id ]
            mins <- c(min(state_sets[[ desc_id[1] ]]), min(state_sets[[ desc_id[2] ]]))
            maxs <- c(max(state_sets[[ desc_id[1] ]]), max(state_sets[[ desc_id[2] ]]))
            state_sets[[ node_id ]] <- c( min( maxs ), max( mins ) )
        }
    }
    root_id <- tree$edge[,1][length(tree$edge[,1])]
    desc_id <- tree$edge[,2][ tree$edge[,1] == root_id ]
    mins <- c(min(state_sets[[ desc_id[1] ]]), min(state_sets[[ desc_id[2] ]]))
    maxs <- c(max(state_sets[[ desc_id[1] ]]), max(state_sets[[ desc_id[2] ]]))
    state_sets[[ root_id ]] <- c( min( maxs ), max( mins ) )
    # 2: preorder traversal and for each node calculate median of ancestral node state and Farris intervals.
    anc_states <- c( median( state_sets[[ root_id ]] ) )
    for (i in length(tree$edge[,1]):1) {
        node_id <- tree$edge[i,2]
        if (node_id > root_id) {
            anc_id <- tree$edge[i,1]
            anc_state <- median( c(state_sets[[ node_id ]], anc_states[ anc_id - root_id + 1 ]) )
            anc_states[node_id - root_id + 1] <- anc_state
        }
    }
    anc_states
}


# Modified tree.quasse function to accept regimes. 
# lambda, mu, and char are each input as a list of functions where
# each function corresponds to a regime. The regime times are
# input as the vector 'regimes'.
tree.quasse.regimes <- function(pars, regimes=NA, max.taxa=Inf, max.t=Inf,
                        include.extinct=TRUE, x0=NA,
                        single.lineage=TRUE, verbose=FALSE) {
    if ( is.na(x0) )
        stop("x0 must be specified")
    else if ( length(x0) != 1 )
        stop("x0 must be of length 1")
    if ( is.na(regimes) || !is.vector(regimes) )
        stop("regimes must be a vector of times.")
    stopifnot(is.list(pars), all(sapply(pars, is.list)))
  
    info <- make.tree.quasse.regimes(pars, regimes, max.taxa, max.t, x0, single.lineage,
                             verbose)
    if ( single.lineage )
        info <- info[-1,]
    phy <- me.to.ape.quasse(info)
    if ( include.extinct || is.null(phy) )
        phy
    else {
        phy2 <- prune(phy)
        # after pruning extinct lineages, check if any survived
        if (class(phy2) != "phylo")
            NULL
        else {
            # add root stem
            finishing_time <- max(attr(phy2$orig, "ages"))
            phy2$root.edge <- finishing_time - max(branching.times(phy2))
            phy2
        }
    }
}


# Helper function to simulate quasse trees with regimes
# modified from package diversitree
make.tree.quasse.regimes <- function(pars, regimes, max.taxa=Inf, max.t=Inf, x0,
                             single.lineage=TRUE,
                             verbose=FALSE, k=500, ...) {
  lambda   <- pars[[1]]
  lambda_d <- pars[[2]]
  mu       <- pars[[3]]
  mu_d     <- pars[[4]]
  char     <- pars[[5]]
  
  if ( single.lineage ) {
    info <- data.frame(idx=1, len=1e-8, parent=0, state=x0,
                       extinct=FALSE, split=FALSE)
  } else {
    info <- data.frame(idx=1:2, len=1e-8, parent=0, state=x0,
                       extinct=FALSE, split=FALSE)
  }

  lineages <- which(!info$extinct & !info$split)
  n.taxa <- c(length(lineages))
  ages <- c(0)
  t <- 0
  t.left_regime <- regimes[1]
  t.left_total <- max.t
  i <- 1
  j <- 2
  while ( n.taxa[1] <= max.taxa && n.taxa[1] > 0 && t.left_total > 0 ) {
      verbose = FALSE
      while ( t.left_regime > 0 && t.left_total > 0 ) {
        x <- run.until.change(lineages, info, k, lambda[[i]], lambda_d[[i]], mu[[i]], mu_d[[i]], char[[i]], t.left_regime)
        lineages <- x[[1]]
        info <- x[[2]]
        n.taxa <- c(n.taxa, length(lineages))
        t <- t + x[[4]]
        ages <- c(ages, t)
        t.left_total <- t.left_total - x[[4]]
        t.left_regime <- t.left_regime - x[[4]]
        if ( verbose )
          cat(sprintf("%s: %d [%2.3f]\n",
                      c("-", " ", "+")[sign(x[[3]])+2], n.taxa[j], t))
        j <- j + 1
      }
      i <- i + 1
      t.left_regime <- regimes[i]
  }

  if ( n.taxa[i] > max.taxa ) {
    ## Drop final speciation event.
    drop <- info[nrow(info)-1,]
    info$split[drop$parent] <- FALSE
    info$state[drop$parent] <- drop$state
    info$len[drop$parent] <- info$len[drop$parent] + drop$len
    info <- info[seq_len(nrow(info)-2),]
  }

  attr(info, "t") <- t
  attr(info, "lineages_thru_time") <- n.taxa
  attr(info, "ages") <- ages
  info
}


# Function to evolve the tree under a given regime
# modified to accept diversity dependent lambda and mu
# from package diversitree.
run.until.change <- function(lineages, info, k, lambda, lambda_d, mu, mu_d, char,
                             max.t) {
  i <- 1
  time <- 0
  n.extant <- length(lineages)
  p.change <- 1/k
  niter <- 1
  repeat {
    state <- info$state[lineages]
    lx <- lambda(state) + lambda_d(length(lineages))
    mx <- mu(state) + mu_d(length(lineages))
    r <- sum(lx + mx)
    dt <- 1/(r*k)
    
    if ( runif(1) < p.change ) {
      if ( runif(1) < sum(lx)/r ) { # speciation
        i <- sample(n.extant, 1, prob=lx)
        info <- speciate(info, lineages[i])
        lineages <- c(lineages[-i], c(-1,0) + nrow(info))
      } else {
        i <- sample(n.extant, 1, prob=mx)
        info$extinct[lineages[i]] <- TRUE
        lineages <- lineages[-i]
      }
      info$len[lineages]   <- info$len[lineages] + dt
      info$state[lineages] <- char(info$state[lineages], dt)
      time <- time + dt
      break
    }

    info$len[lineages]   <- info$len[lineages] + dt
    info$state[lineages] <- char(info$state[lineages], dt)
    niter <- niter + 1
    time <- time + dt

    if ( time > max.t )
      break
  }

  list(lineages, info, length(lineages) - n.extant, time)
}


# This function rejigs the 'lineages' structure when speciation
# happens, creating new species.
# from package diversitree
speciate <- function(info, i) {
  j <- 1:2 + nrow(info)
  info[j,"idx"] <- j
  info[j,-1] <- list(len=0, parent=i,
                     state=info$state[i],
                     extinct=FALSE, split=FALSE)
  info$split[i] <- TRUE
  info
}


# Transform the lineages structure produced by sim.tree into an ape
# phylogeny.
# from package diversitree
me.to.ape.quasse <- function(info) {
  if ( nrow(info) == 0 )
    return(NULL)
  Nnode <- sum(!info$split) - 1
  n.tips <- sum(!info$split)

  info$idx2 <- NA
  info$idx2[!info$split] <- 1:n.tips
  info$idx2[ info$split] <- order(info$idx[info$split]) + n.tips + 1

  i <- match(info$parent, info$idx)
  info$parent2 <- info$idx2[i]
  info$parent2[is.na(info$parent2)] <- n.tips + 1

  tip.label <- ifelse(subset(info, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)

  info$name <- NA
  info$name[!info$split] <- tip.label

  tip.state <- info$state[match(1:n.tips, info$idx2)]
  names(tip.state) <- tip.label
  phy <- reorder(structure(list(edge=cbind(info$parent2, info$idx2),
                               Nnode=Nnode,
                               tip.label=tip.label,
                               tip.state=tip.state,
                               node.label=node.label,
                               edge.length=info$len,
                               orig=info),
                          class="phylo"))

  phy$edge.state <- info$state[match(phy$edge[,2], info$idx2)]
  phy
}


# Modified traitgram function to plot given simulated and estimated ancestral trait values.
# If method="sim" then plot black trait lines.
# If method="est" add to an existing plot the estimated trait values with red dashed lines.
# Modified from package picante.
traitgram_given <- function (x, internal_node_values, phy, finishing_time, color="red", xaxt = "s", underscore = FALSE, show.names = TRUE,
    show.xaxis.values = TRUE, method = c("sim", "est", "ML", "pic"), ...)
{
    method <- match.arg(method)
    Ntaxa = length(phy$tip.label)
    Ntot = Ntaxa + phy$Nnode
    phy = node.age(phy)
    ages = phy$ages[match(1:Ntot, phy$edge[, 2])]
    ages[Ntaxa + 1] = 0
    if (class(x) %in% c("matrix", "array")) {
        xx = as.numeric(x)
        names(xx) = row.names(x)
    }
    else xx = x
    if (!is.null(names(xx))) {
        umar = 0.1
        if (!all(names(xx) %in% phy$tip.label)) {
            print("trait and phy names do not match")
            return()
        }
        xx = xx[match(phy$tip.label, names(xx))]
    }
    else umar = 0.1
    lmar = 0.2
    if (xaxt == "s")
        if (show.xaxis.values)
            lmar = 1
        else lmar = 0.5
    if (method == "sim" || method == "est")
        xanc <- internal_node_values[2:length(internal_node_values)] # skip the root state
    else
        xanc <- ace(xx, phy, method = method)$ace
    xall = c(xx, xanc)
    a0 = ages[phy$edge[, 1]]
    a1 = ages[phy$edge[, 2]]
    x0 = xall[phy$edge[, 1]]
    x1 = xall[phy$edge[, 2]]
    if (method == "est" || method == "sim") {
        # add the root state and age
        offset <- finishing_time - max(a1)
        a1 <- c(offset + a0[1], offset + a1)
        a0 <- c(0, offset + a0)
        x1 <- c(x0[1], x1)
        x0 <- c(internal_node_values[1], x0)
    }
    tg = par(bty = "n", mai = c(lmar, 0.1, umar, 0.1))
    if (show.names) {
        #maxNameLength = max(nchar(names(xx)))
        ylim = c(0, finishing_time + 2) #* (1 + maxNameLength/50))
        if (!underscore)
            names(xx) = gsub("_", " ", names(xx))
    }
    else ylim = c(0, finishing_time)
    if (method != "est") {
        plot(range(c(x0, x1)), range(c(a0, a1)), type = "n", xaxt = "n",
            yaxt = "n", xlab = "", ylab = "", bty = "n", ylim = ylim, xlim=c(-10,10),
            cex.axis = 0.8)
        if (xaxt == "s")
            if (show.xaxis.values)
                axis(1, labels = TRUE)
            else axis(1, labels = FALSE)
        segments(x0, a0, x1, a1)
    } else {
        segments(x0, a0, x1, a1, col=color, lty=2)
    }
    if (show.names) {
        text(sort(xx), finishing_time, labels = names(xx)[order(xx)],
            adj = -0, srt = 90)
    }
    on.exit(par(tg))
}

