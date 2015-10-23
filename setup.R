#
# this file loads the libraries and custom functions necessary to
# run the simulations
#

library(diversitree)
library(phytools)
library(picante)


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
  else
    prune(phy)
}


# helper function to build simulate quasse trees with regimes
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
  while ( n.taxa <= max.taxa && n.taxa > 0 && t.left_total > 0 ) {
      verbose = TRUE
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

  if ( n.taxa > max.taxa ) {
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


# function to evolve the tree under a give regime
# modified to accept diversity dependent lambda and mu
# from package diversitree
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


#prune <- function (phy, to.drop = NULL)
#{
#    if (is.null(to.drop))
#        to.drop <- subset(phy$orig, !split)$extinct
#    if (sum(!to.drop) < 2) {
#        NULL
#    }
#    else if (any(to.drop)) {
#        #orig_max <- max(branching.times(phy))
#        phy2 <- drop.tip.fixed(phy, phy$tip.label[to.drop])
#        phy2$orig <- phy2$orig[!phy2$orig$extinct, ]
#        phy2$tip.state <- phy2$tip.state[!to.drop]
#        phy2$node.state <- phy2$node.state[phy2$node.label]
#        phy2$hist <- prune.hist(phy, phy2)
#        #new_max <- max(branching.times(phy2))
#        #phy2$edge.length[1] <- phy2$edge.length[1] + (orig_max - new_max)
#    }
#    else {
#        phy
#    }
#}


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
traitgram_est <- function (x, internal_node_values, phy, xaxt = "s", underscore = FALSE, show.names = TRUE,
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


### This function aims to covert the "hist" object.  This is fairly
### complicated and possibly can be streamlined a bit.  The big issue
### here is that when extinct species are removed from the tree, it
### leaves unbranched nodes - the history along a branch with such a
### node needs to be joined.
#prune.hist <- function(phy, phy2) {
#  hist <- phy$hist
#  if ( is.null(hist) || nrow(hist) == 0 )
#    return(hist)
#
#  ## More interesting is to collect up all of the names and look at the
#  ## branches that terminate
#  phy.names <- c(phy$tip.label, phy$node.label)
#  phy2.names <- c(phy2$tip.label, phy2$node.label)
#
#  ## Next, check what the parent of the nodes is in the new tree, using
#  ## the standard names (parent-offspring)
#  ## First, for phy2
#  po.phy <- cbind(from=phy.names[phy$edge[,1]],
#                  to=phy.names[phy$edge[,2]])
#  po.phy2 <- cbind(from=phy2.names[phy2$edge[,1]],
#                   to=phy2.names[phy2$edge[,2]])
#
#  ## Then find out where the parent/offspring relationship changed:
#  ## i <- match(po.phy2[,2], po.phy[,2])
#  j <- which(po.phy[match(po.phy2[,2], po.phy[,2]),1] != po.phy2[,1])
#
#  for ( idx in j ) {
#    to <- po.phy2[idx,2]
#    from <- po.phy2[idx,1]
#    ans <- to
#    offset <- 0
#    repeat {
#      to <- po.phy[po.phy[,2] == to,1]
#      ans <- c(to, ans)
#      if ( is.na(to) )
#        stop("Horrible error")
#      if ( to == from )
#        break
#    }
#    
#    if ( any(ans[-1] %in% hist$name2) ) {
#      k <- hist$name2 %in% ans[-1]
#      offset <- cumsum(phy$edge.length[match(ans[-1], po.phy[,2])])
#      offset <- c(0, offset[-length(offset)])
#      hist$x0[k] <- hist$x0[k] - offset[match(hist$name2[k], ans[-1])]
#      hist$tc[k] <- hist$t[k] - hist$x0[k]
#      hist$name2[k] <- ans[length(ans)]
#    }
#  }
#
#  ## Prune out the extinct species and nodes that lead to them.  Note
#  ## that the root must be excluded as history objects that lead to
#  ## the new root (if it has changed) should not be allowed.
#  phy2.names.noroot <- phy2.names[phy2.names != phy2$node.label[1]]
#  hist <- hist[hist$name2 %in% phy2.names.noroot,]
#  #hist <- hist[hist$name2 %in% phy2.names,]
#
#  ## Remake idx2 to point at the new tree.
#  hist$idx2 <- match(hist$name2, phy2.names)
#
#  hist[order(hist$idx2, hist$t),]
#}
#
#
### This is a patched version of drop.tip that will keep the nodes in
### the correct order.  Otherwise it is exactly the same as the
### drop.tip in ape version 2.4-1 (right down to the use of dim(x)[1]
### instead of nrow(x)).  See REMOVE/REPLACE/DONE at the end for changes.
#drop.tip.fixed <- function(phy, tip, trim.internal = TRUE, subtree =
#                           FALSE, root.edge = 0, rooted = is.rooted(phy)) {
#  if (!inherits(phy, "phylo")) 
#    stop("object \"phy\" is not of class \"phylo\"")
#  Ntip <- length(phy$tip.label)
#  if (is.character(tip)) 
#    tip <- which(phy$tip.label %in% tip)
#  if (!rooted && subtree) {
#    phy <- root(phy, (1:Ntip)[-tip][1])
#    root.edge <- 0
#  }
#  phy <- reorder(phy)
#  NEWROOT <- ROOT <- Ntip + 1
#  Nnode <- phy$Nnode
#  Nedge <- dim(phy$edge)[1]
#  if (subtree) {
#    trim.internal <- TRUE
#    tr <- reorder(phy, "pruningwise")
#    N <- node.depth(phy)
#  }
#  wbl <- !is.null(phy$edge.length)
#  edge1 <- phy$edge[, 1]
#  edge2 <- phy$edge[, 2]
#  keep <- !logical(Nedge)
#  if (is.character(tip)) 
#    tip <- which(phy$tip.label %in% tip)
#  if (!rooted && subtree) {
#    phy <- root(phy, (1:Ntip)[-tip][1])
#    root.edge <- 0
#  }
#  keep[match(tip, edge2)] <- FALSE
#  if (trim.internal) {
#    ints <- edge2 > Ntip
#    repeat {
#      sel <- !(edge2 %in% edge1[keep]) & ints & keep
#      if (!sum(sel)) 
#        break
#      keep[sel] <- FALSE
#    }
#    if (subtree) {
#      subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
#      keep[subt] <- TRUE
#    }
#    if (root.edge && wbl) {
#      degree <- tabulate(edge1[keep])
#      if (degree[ROOT] == 1) {
#        j <- integer(0)
#        repeat {
#          i <- which(edge1 == NEWROOT & keep)
#          j <- c(i, j)
#          NEWROOT <- edge2[i]
#          degree <- tabulate(edge1[keep])
#          if (degree[NEWROOT] > 1) 
#            break
#        }
#        keep[j] <- FALSE
#        if (length(j) > root.edge) 
#          j <- 1:root.edge
#        NewRootEdge <- sum(phy$edge.length[j])
#        if (length(j) < root.edge && !is.null(phy$root.edge)) 
#          NewRootEdge <- NewRootEdge + phy$root.edge
#        phy$root.edge <- NewRootEdge
#      }
#    }
#  }
#  if (!root.edge) 
#    phy$root.edge <- NULL
#  phy$edge <- phy$edge[keep, ]
#  if (wbl) 
#    phy$edge.length <- phy$edge.length[keep]
#  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
#  oldNo.ofNewTips <- phy$edge[TERMS, 2]
#  n <- length(oldNo.ofNewTips)
#  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
#  if (subtree || !trim.internal) {
#    tips.kept <- oldNo.ofNewTips <= Ntip & !(oldNo.ofNewTips %in% 
#                   tip)
#    new.tip.label <- character(n)
#    new.tip.label[tips.kept] <- phy$tip.label[-tip]
#    node2tip <- oldNo.ofNewTips[!tips.kept]
#    new.tip.label[!tips.kept] <- if (subtree) {
#      paste("[", N[node2tip], "_tips]", sep = "")
#    }
#    else {
#      if (is.null(phy$node.label)) 
#        rep("NA", length(node2tip))
#      else phy$node.label[node2tip - Ntip]
#    }
#    if (!is.null(phy$node.label)) 
#      phy$node.label <- phy$node.label[-(node2tip - Ntip)]
#    phy$tip.label <- new.tip.label
#  }
#  else phy$tip.label <- phy$tip.label[-tip]
#  if (!is.null(phy$node.label)) 
#    phy$node.label <- phy$node.label[sort(unique(phy$edge[, 
#                                                          1])) - Ntip]
#  phy$Nnode <- dim(phy$edge)[1] - n + 1L
#
#  ## REMOVE:
#  ##     newNb <- integer(n + phy$Nnode)
#  ##     newNb[NEWROOT] <- n + 1L
#  ##     sndcol <- phy$edge[, 2] > n
#  ##     phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- (n + 
#  ##         2):(n + phy$Nnode)
#  ##     phy$edge[, 1] <- newNb[phy$edge[, 1]]
#  ## REPLACE:
#  i <- phy$edge > n
#  phy$edge[i] <- match(phy$edge[i], sort(unique(phy$edge[i]))) + n
#  ## DONE:
#  
#  storage.mode(phy$edge) <- "integer"
#  collapse.singles(phy)
#}

