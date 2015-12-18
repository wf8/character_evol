#
# Script to create a dataframe summarizing the simulation results
# Skewness and kurtosis computed with library e1071.
# Tree imbalance (colless' and sackin's indices) computed with library apTreeshape.
# Gamma statistic of Pybus and Harvey and Faith's PD computed using ape.
#

file_path = "output/"
summary_path = "results/"
sim_types = c("uniform",
              "static-peak",
              "jumping-peak",
              "moving-peak",
              "moving-peak-density05",
              "moving-peak-density1")

library(e1071)
library(apTreeshape)


# function to return a vector of the number of nodes
# from each tip to the root
get_node_dist = function(tree) {
    node_dist = vector()
    n_tips = length( tree$tip.label )
    t = reorder(tree, "postorder")
    for (i in 1:n_tips) {
        next_node = i
        distance = 0
        for (j in 1:length(t$edge[,1])) {
            begin = t$edge[j, 1]
            end = t$edge[j, 2]
            if (end == next_node) {
                distance = distance + 1
                next_node = begin
            }
        }
        node_dist = c(node_dist, distance)
    }
    return( node_dist )
}

# loop through results files and build dataframe 
for ( i in 1:length(sim_types) ) {

    load( paste(file_path, sim_types[i], ".RData", sep="") )

    # simulation parameters
    sim_type = list()
    birth = list()
    death = list()
    beta = list()
    ppeak_birth = list()
    sigma_birth = list()
    jump = list()
    slide = list()
    density_dep = list()
    sim_time = list()

    d = data.frame()

    # loop through each simulation scenario
    for ( j in 1:length(results) ) {

        d[j, "sim_type"] = sim_types[i]
        d[j, "sim_hours"] = as.numeric( results[j][[1]]$processing_time, units="hours" )
        d[j, "n_replicates"] = length( results[j][[1]]$data )

        # tree shape stats
        n_too_many = 0
        n_extinct = 0
        n_error = 0
        num_survivors = vector()
        node_ages = vector()
        colless_indices = vector()
        sackins_indices = vector()
        pd = vector()
        gamma = vector()

        # character simulation stats
        tip_ranges = vector()
        tip_skews = vector()
        tip_kurtosis = vector()
        corr_tip_values_to_nodal_distance = vector()

        # ancestral state estimates
        bm_root_diff = vector()
        bm_mean_node_diff = vector()
        bm_root_ci_min = vector()
        bm_root_ci_max = vector()
        bm_root_diff_div_tip_range = vector()
        lp_root_diff = vector()
        lp_mean_node_diff = vector()
        lp_root_diff_div_tip_range = vector()

        # loop through each replicate
        for ( k in 1:length( results[j][[1]]$data ) ) {
           
            rep = results[j][[1]]$data[k][[1]]
            
            if (rep == "Lineage threshold of 250 exceeded") {
                
                n_too_many = n_too_many + 1

            } else if (rep == "extinct") { 

                n_extinct = n_extinct + 1

            } else if (rep == "error") {
            
                n_error = n_error + 1

            } else {

                # tree shape stats
                num_survivors = c(num_survivors, length( rep$tip_states ) )
                node_ages = c(node_ages, rep$branch_times )        
                treeshape = as.treeshape( rep$sim_tree )
                if (length(rep$tip_states) > 2) {
                    colless_indices = c(colless_indices, colless( treeshape ) )
                    sackins_indices = c(sackins_indices, sackin( treeshape ) )
                }
                gamma = c(gamma, gammaStat( rep$sim_tree ) )
                pd = c(pd, sum(rep$sim_tree$edge.length) + rep$sim_tree$root.edge)

                # character simulation stats
                tip_ranges = c(tip_ranges, max( rep$tip_states ) - min( rep$tip_states ) )
                tip_skews = c(tip_skews, skewness( rep$tip_states ) )
                tip_kurtosis = c(tip_kurtosis, kurtosis( rep$tip_states ) )
                node_dist = get_node_dist( rep$sim_tree )
                corr_tip_values_to_nodal_distance = c(corr_tip_values_to_nodal_distance, cor(node_dist, rep$tip_states))

                # ancestral state estimates
                bm_root_diff = c(bm_root_diff, rep$root_difference_bm)
                bm_mean_node_diff = c(bm_mean_node_diff, mean(rep$node_differences_bm))
                bm_root_ci_min = c(bm_root_ci_min, min( rep$est_data_bm$CI95[ length(rep$est_data_bm$CI95[,2]),] ) )
                bm_root_ci_max = c(bm_root_ci_max, max( rep$est_data_bm$CI95[ length(rep$est_data_bm$CI95[,2]),] ) )
                bm_root_diff_div_tip_range = c(bm_root_diff_div_tip_range, rep$root_difference_bm / max( rep$tip_states ) - min( rep$tip_states ) )
                lp_root_diff = c(lp_root_diff, rep$root_difference_lp)
                lp_mean_node_diff = c(lp_mean_node_diff, mean(rep$node_differences_lp))
                lp_root_diff_div_tip_range = c(lp_root_diff_div_tip_range, rep$root_difference_lp / max( rep$tip_states ) - min( rep$tip_states ) )

           }
        } # end looping thru replicates

        # add data to dataframe

        # tree shape stats
        d[j, "num_rep_over_250_lineages"] = n_too_many
        d[j, "num_rep_extinct"] = n_extinct
        d[j, "num_rep_error"] = n_error
        d[j, "mean_num_survivor_lineages"] = mean(num_survivors, na.rm=TRUE)
        d[j, "median_node_age"] = median(node_ages)
        d[j, "mean_tree_imbalance_colless"] = mean( colless_indices, na.rm=TRUE )
        d[j, "mean_tree_imbalance_sackin"] = mean( sackins_indices, na.rm=TRUE )
        d[j, "mean_gamma_stat"] = mean( gamma, na.rm=TRUE )
        d[j, "mean_phylo_diversity"] = mean( pd, na.rm=TRUE )

        # character simulation stats
        d[j, "mean_tip_range"] = mean( tip_ranges, na.rm=TRUE )
        d[j, "mean_tip_skewness"] = mean( tip_skews, na.rm=TRUE )
        d[j, "mean_tip_kutosis"] = mean( tip_kurtosis, na.rm=TRUE )
        d[j, "mean_corr_tip_values_to_nodal_distance"] = mean( corr_tip_values_to_nodal_distance, na.rm=TRUE )

        # ancestral state estimates
        d[j, "mean_bm_root_diff"] = mean( bm_root_diff, na.rm=TRUE )
        d[j, "mean_bm_node_diff"] = mean( bm_mean_node_diff, na.rm=TRUE )
        d[j, "mean_bm_root_ci_min"] = mean( bm_root_ci_min, na.rm=TRUE )
        d[j, "mean_bm_root_ci_max"] = mean( bm_root_ci_max, na.rm=TRUE )
        d[j, "mean_bm_root_diff_div_tip_range"] = mean( bm_root_diff_div_tip_range, na.rm=TRUE )
        d[j, "mean_lp_root_diff"] = mean( lp_root_diff, na.rm=TRUE )
        d[j, "mean_lp_mean_node_diff"] = mean( lp_mean_node_diff, na.rm=TRUE )
        d[j, "mean_lp_root_diff_div_tip_range"] = mean( lp_root_diff_div_tip_range, na.rm=TRUE )

    } # end looping through all simulation scenarios
    print(paste("Writing file: ", summary_path, sim_types[i], ".csv", sep=""))
    write.table(d, sep=",", row.names=FALSE, file=paste(summary_path, sim_types[i], ".csv", sep=""))
}
print("Done.")
