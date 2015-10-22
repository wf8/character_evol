
library(diversitree)
library(phytools)
library(picante)


# Modified tree.quasse function to accept regimes. 
# lambda, mu, and char are each input as a list of functions where
# each function corresponds to a regime. The regime times are
# input as the vector 'regimes'.
tree.quasse.regimes <- function(pars, regimes=NA, max.taxa=Inf, max.t=Inf,
                        include.extinct=FALSE, x0=NA,
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
  else
    prune(phy)
}


# helper function to build simulate quasse trees with regimes
# modified from package diversitree
make.tree.quasse.regimes <- function(pars, regimes, max.taxa=Inf, max.t=Inf, x0,
                             single.lineage=TRUE,
                             verbose=FALSE, k=500, ...) {
  lambda <- pars[[1]]
  mu     <- pars[[2]]
  char   <- pars[[3]]
  
  if ( single.lineage ) {
    info <- data.frame(idx=1, len=1e-8, parent=0, state=x0,
                       extinct=FALSE, split=FALSE)
  } else {
    info <- data.frame(idx=1:2, len=1e-8, parent=0, state=x0,
                       extinct=FALSE, split=FALSE)
  }

  lineages <- which(!info$extinct & !info$split)
  n.taxa <- length(lineages)
  t <- 0
  t.left_regime <- regimes[1]
  t.left_total <- max.t
  i <- 1

  while ( n.taxa <= max.taxa && n.taxa > 0 && t.left_total > 0 ) {
      verbose = TRUE
      while ( t.left_regime > 0 && t.left_total > 0 ) {
        x <- run.until.change(lineages, info, k, lambda[[i]], mu[[i]], char[[i]], t.left_regime)
        lineages <- x[[1]]
        info <- x[[2]]
        n.taxa <- length(lineages)
        t <- t + x[[4]]
        t.left_total <- t.left_total - x[[4]]
        t.left_regime <- t.left_regime - x[[4]]
        if ( verbose )
          cat(sprintf("%s: %d [%2.3f]\n",
                      c("-", " ", "+")[sign(x[[3]])+2], n.taxa, t))
      }
      i <- i + 1
      t.left_regime <- regimes[i]
  }

  if ( n.taxa > max.taxa ) {
    ## Drop final speciation event.
    drop <- info[nrow(info)-1,]
    info$split[drop$parent] <- FALSE
    info$state[drop$parent] <- drop$state
    info$len[drop$parent] <- info$len[drop$parent] + drop$len
    info <- info[seq_len(nrow(info)-2),]
  }

  attr(info, "t") <- t
  info
}


# function to evolve the tree under a give regime
# from package diversitree
run.until.change <- function(lineages, info, k, lambda, mu, char,
                             max.t) {
  i <- 1
  time <- 0
  n.extant <- length(lineages)
  p.change <- 1/k
  niter <- 1
  repeat {
    state <- info$state[lineages]
    lx <- lambda(state)
    mx <- mu(state)
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


# modified traitgram function to plot simulated ancestral values
traitgram_sim <- function (x, internal_node_values, phy, xaxt = "s", underscore = FALSE, show.names = TRUE,
    show.xaxis.values = TRUE, method = c("ML", "pic"), ...)
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
    #xanc <- ace(xx, phy, method = method)$ace
    xanc <- internal_node_values
    xall = c(xx, xanc)
    a0 = ages[phy$edge[, 1]]
    a1 = ages[phy$edge[, 2]]
    x0 = xall[phy$edge[, 1]]
    x1 = xall[phy$edge[, 2]]
    tg = par(bty = "n", mai = c(lmar, 0.1, umar, 0.1))
    if (show.names) {
        maxNameLength = max(nchar(names(xx)))
        ylim = c(min(ages), max(ages) * (1 + maxNameLength/50))
        if (!underscore)
            names(xx) = gsub("_", " ", names(xx))
    }
    else ylim = range(ages)
    plot(range(c(x0, x1)), range(c(a0, a1)), type = "n", xaxt = "n",
        yaxt = "n", xlab = "", ylab = "", bty = "n", ylim = ylim, xlim=c(-10,10),
        cex.axis = 0.8)
    if (xaxt == "s")
        if (show.xaxis.values)
            axis(1, labels = TRUE)
        else axis(1, labels = FALSE)
    segments(x0, a0, x1, a1)
    if (show.names) {
        text(sort(xx), max(ages), labels = names(xx)[order(xx)],
            adj = -0, srt = 90)
    }
    on.exit(par(tg))
}


# modified traitgram function to add estimated ancestral values to existing plot
traitgram_est <- function (x, phy, xaxt = "s", underscore = FALSE, show.names = TRUE,
    show.xaxis.values = TRUE, method = c("ML", "pic"), ...)
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
    xanc <- ace(xx, phy, method = method)$ace
    xall = c(xx, xanc)
    a0 = ages[phy$edge[, 1]]
    a1 = ages[phy$edge[, 2]]
    x0 = xall[phy$edge[, 1]]
    x1 = xall[phy$edge[, 2]]
    tg = par(bty = "n", mai = c(lmar, 0.1, umar, 0.1))
    if (show.names) {
        maxNameLength = max(nchar(names(xx)))
        ylim = c(min(ages), max(ages) * (1 + maxNameLength/50))
        if (!underscore)
            names(xx) = gsub("_", " ", names(xx))
    }
    else ylim = range(ages)
    #plot(range(c(x0, x1)), range(c(a0, a1)), type = "n", xaxt = "n",
    #    yaxt = "n", xlab = "", ylab = "", bty = "n", ylim = ylim,
    #    cex.axis = 0.8)
    #if (xaxt == "s")
    #    if (show.xaxis.values)
    #        axis(1, labels = TRUE)
    #    else axis(1, labels = FALSE)
    segments(x0, a0, x1, a1, col="red", lty=2)
    if (show.names) {
        text(sort(xx), max(ages), labels = names(xx)[order(xx)],
            adj = -0, srt = 90)
    }
    on.exit(par(tg))
}

# get start time
start_time <- Sys.time()

# a list of sigmoidal speciation functions for the regime shifts
#lambda <- list(function(x) sigmoid.x(x, y0=0, y1=8, xmid=2, r=2), function(x) sigmoid.x(x, y0=0, y1=8, xmid=1, r=2), function(x) sigmoid.x(x, y0=0, y1=8, xmid=0, r=2))

# a list of normal optimal speciation functions for the regime shifts
#lambda <- list(function(x) noroptimal.x(x, y0=0, y1=30, xmid=0, s2=10))
lambda <- list(function(x) constant.x(x, 0.3))

# a list of contant extinction functions for the regime shifts
mu <- list(function(x) constant.x(x, 0.0))

# character evolving with brownian motion
char <- list(make.brownian.with.drift(0, 1.0))

phy <- list()
tip_states <- list()
anc_states <- list()
branch_times <- list()
estimated_anc_states <- list()
simulated_anc_states <- list()
rescaled_branch_time <- list()
differences <- list()
root_differences <- vector()
tree_depths <- vector()

replicates <- 1
for (i in 1:replicates) {

    # simulate tree and character
    #phy[[i]] <- tree.quasse(c(lambda, mu, char), max.taxa=50, x0=0, single.lineage=FALSE)
    phy[[i]] <- tree.quasse.regimes(list(lambda, mu, char), regimes=c(10), max.t=10, x0=0, single.lineage=FALSE)

    # get tip states
    tip_states[[i]] <- phy[[i]]$tip.state

    # infer ancestral states
    anc_states[[i]] <- ace(tip_states[[i]], phy[[i]], type="continuous")

    #root_index <- 51

    branch_times[[i]] <- as.vector(branching.times(phy[[i]]))
    estimated_anc_states[[i]] <- as.vector(anc_states[[i]]$ace)

    # reorder and extract the simulated ancestral states 
    simulated_anc_states[[i]] <- 0.0
    for (j in length(phy[[i]]$tip.state)+1:max(phy[[i]]$orig$idx2)) {

       simulated_anc_states[[i]] <- c(simulated_anc_states[[i]], phy[[i]]$orig$state[ phy[[i]]$orig$idx2==j ] )

    }

    tree_depths[i] <- max(branch_times[[i]])
    rescaled_branch_time[[i]] <- branch_times[[i]] / max(branch_times[[i]])
    differences[[i]] <- simulated_anc_states[[i]] - estimated_anc_states[[i]]
    root_differences[i] <- simulated_anc_states[[i]][1] - estimated_anc_states[[i]][1]
}

processing_time <- Sys.time() - start_time

#hist(tree_depths, breaks=20)
#hist(root_differences)

# plot rescaled branch times against difference between simulated and estimated states
#plot(rescaled_branch_time, simulated_anc_states - estimated_anc_states)

# plot the two traitgrams on top of one another
traitgram_sim(tip_states[[1]], simulated_anc_states[[1]], phy[[1]])
traitgram_est(tip_states[[1]], phy[[1]])


# plot the ancestral states as circles on the nodes
#plot(phy, label.offset=.05)
#nodelabels()
#tiplabels()
#tiplabels(pch=21, cex=tip_states, bg="yellow")
#nodelabels(pch=21, cex= anc_states$ace, bg="yellow")
