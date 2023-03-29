# Sample from the uniform distribution over partitions of size n,
# using the method of Stam 1983
# First samples the number of clusters of the partition, u, and then
# allocates the n units into u bins. 
runiformpartition <- function(n, initial_factor = 10){
  initial <- initial_factor * n
  initial_seq <- 1:initial
  log_bell <- salso::lbell(n)
  log_u_probs <- n * log(initial_seq) - lfactorial(initial_seq) - 1 - log_bell
  u_probs <- exp(log_u_probs)
  sum_u_probs <- sum(u_probs)
  
  u_samp <- sample(x = c(initial_seq, initial + 1), size = 1,
                   prob = c(u_probs, 1 - sum_u_probs))
  while(u_samp == (initial + 1)){
    initial_factor <- initial_factor + 1
    initial_seq <- (initial + 1):(initial_factor * n)
    initial <- initial_factor * n
    
    log_u_probs <- n * log(initial_seq) - lfactorial(initial_seq) - 1 - log_bell
    u_probs <- exp(log_u_probs)
    sum_u_probs_temp <- sum(u_probs)
    sum_u_probs <- sum_u_probs +sum_u_probs_temp
    
    u_samp <- sample(x = c(initial_seq, initial + 1), size = 1,
                     prob = c(u_probs, 1 - sum_u_probs) / 
                       (sum_u_probs_temp + 1 - sum_u_probs))
  }
  
  
  
  label_samp <- sample(1:u_samp, size = n, replace = T,
                       prob = rep(1 / u_samp, u_samp))
  return(as.numeric(as.factor(label_samp)))
}

# Given an indexing scheme, find the connected components defined by the 
# indexing scheme, and sample a uniform partition over each connected
# component. 
uniform_index_sample <- function(comparison_list, 
                                 pairs_to_keep,
                                 num_samples = 4, 
                                 duplicates,
                                 seed = 43){
  
  set.seed(seed)
  
  turn_F <- function(x){return(rep(F, length(x)))}
  comparison_list_initial <- comparison_list
  comparison_list_initial$comparisons <- 
    comparison_list_initial$comparisons %>% as.data.frame() %>%
    mutate_all(turn_F) %>% as.matrix()
  
  reduced_comparison_list_initial <- 
    reduce_comparison_data(comparison_list_initial, pairs_to_keep, cc = 1)
  
  pairs_to_keep <- reduced_comparison_list_initial$pairs_to_keep
  
  r <- sum(comparison_list$file_sizes)
  blocks <- matrix(0, nrow = r, ncol = r)
  which_to_keep <- which(pairs_to_keep)
  for(i in 1:length(which_to_keep)){
    temp_pair <- comparison_list$record_pairs[which_to_keep[i], ]
    blocks[temp_pair$i, temp_pair$j] <- 1
    blocks[temp_pair$j, temp_pair$i] <- 1
  }
  
  blocks <- 
    blocks[reduced_comparison_list_initial$labels[order(reduced_comparison_list_initial$labels[, 2]), 1],
           reduced_comparison_list_initial$labels[order(reduced_comparison_list_initial$labels[, 2]), 1]]
  
  g <- igraph::graph_from_adjacency_matrix(blocks, mode = "undirected")
  ccs <- igraph::components(g)
  max_block_size <- max(ccs$csize)
  
  print("Sampling initializations")
  Z_init_temp <- matrix(0, 
                        nrow = sum(reduced_comparison_list_initial$file_sizes),
                        ncol = num_samples)
  partitions <- vector(mode = "list", length = max_block_size)
  partition_sizes <- rep(0, max_block_size)
  # Can only use enumerated partitions for connectected components of size < 13
  for(i in 1:min(max_block_size, 12)){
    partitions[[i]] <- partitions::setparts(i)
    partition_sizes[i] <- ncol(partitions[[i]])
  }
  
  for(i in 1:ccs$no){
    records <- which(ccs$membership == i)
    files <- reduced_comparison_list_initial$file_labels[records]
    block_size <- ccs$csize[i]
    part_size <- partition_sizes[block_size]
    for(j in 1:num_samples){
      samp <- NA
      while(is.na(sum(samp))){
        if(block_size < 13){
          temp <- which(rmultinom(1, size = 1, 
                                  prob = rep(1 / part_size, part_size)) == 1)
          part <- partitions[[block_size]][, temp]
        }
        else{
          part <- runiformpartition(block_size)
        }
        
        flag <- TRUE
        for(k in 1:length(unique(files))){
          temp_part <- part[which(files == k)]
          if((duplicates[k] == 0) & 
             (length(temp_part) != length(unique(temp_part)))){
            flag <- FALSE
          }
        }
        if(flag){
          samp <- part
        }
      }
      
      current <- max(Z_init_temp[, j])
      samp <- samp + current
      Z_init_temp[records, j] <- samp
    }
  }
  
  Z_init <- matrix(NA,
                   nrow = sum(comparison_list_initial$file_sizes), 
                   ncol = num_samples)
  for(i in 1:num_samples){
    Z_init[, i] <- relabel_bayes_estimate(reduced_comparison_list_initial,
                                          Z_init_temp[, i])$link_id %>%
      as.factor %>% as.numeric
  }
  
  
  return(Z_init)
  
  
}

# Function to get ESS at different time points
get_ess <- function(dat, num_times = 20){
  min_samp <- min(dat$`.iteration`)
  max_samp <- max(dat$`.iteration`)
  
  times <- floor(seq(min_samp, max_samp, length.out = num_times)[-1])
  
  output <- data.frame(variable = NA, ess = NA, time = NA)
  colnames_dat <- colnames(dat %>% select(-c(".chain", ".iteration", ".draw")))
  
  for(i in 1:length(times)){
    temp_dat <- dat %>% filter(.iteration <= times[i])
    for(j in 1:length(colnames_dat)){
      temp_ess <- ess_bulk(extract_variable_matrix(temp_dat, colnames_dat[j]))
      output <- rbind(output,
                      data.frame(variable = colnames_dat[j],
                                 ess = temp_ess, 
                                 time = times[i]))
    }
  }
  
  return(output)
  
}