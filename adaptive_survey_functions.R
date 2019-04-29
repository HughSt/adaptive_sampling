
## Adaptive surveys functions
exceedance_prob <- function(gam_model, prediction_data, n_sims, threshold){

  Cg <- predict(gam_model, prediction_data, type = "lpmatrix")
  sims <- rmvn(n_sims, mu = coef(gam_model), V = vcov(gam_model))
  fits <- Cg %*% t(sims)
  fits_prev <- exp(fits) / (1 + exp(fits))
  exceedance_prob <- apply(fits_prev, 1, function(x) {sum(x>threshold, na.rm=T)/n_sims})
  
  # If exceedance = 1, make it 0.999
  exceedance_prob[exceedance_prob==1] <- 0.999
  exceedance_prob[exceedance_prob==0] <- 0.001
  exceedance_prob
}

optimal_range <- function(min_dist, max_dist, length.out = 20, model_data){
  
  REML <- r <- seq(min_dist, max_dist, length.out = length.out)
  for (i in seq_along(r)) {
    m <- mgcv::gam(cbind(npos, nneg) ~ s(bioclim1, k=-1) +
                     s(bioclim4, k=-1) +
                     s(bioclim12, k=-1) +
                     s(bioclim15, k=-1) +
                     s(lng, lat, k=-1, bs="gp", m = c(3, r[i])),
                   data = init_sample,
                   family = "binomial",
                   method = "REML")
    REML[i] <- m$gcv.ubre
  }
  return(list(REML = REML,
              best_m = r[which.min(REML)]))
}

sample_entropy <- function(candidates, entropy, batch_size){
  
  ncol <- ncol(candidates)
  candidates$entropy <- entropy
  candidates$entropy_prob <- entropy / sum(entropy)
  
  
  in_sample <- sample(1:nrow(candidates), batch_size, prob = candidates$entropy_prob)
  samp_prob <- candidates$entropy_prob[in_sample]
  samp_entropy <- candidates$entropy[in_sample]
  
  return(list(sample = candidates[in_sample,1:ncol],
              probs = samp_prob,
              entropy = samp_entropy))
  
}
  

sample_penalized_entropy <- function(candidates, entropy, batch_size){

  ncol <- ncol(candidates)
  candidates$entropy <- entropy
  candidates$entropy_prob <- entropy / sum(entropy)
  
  # Optional - you can transform probability to make it even more likely to sample high entropy areas
  #candidates$entropy_prob <- candidates$entropy_prob^8 / sum(candidates$entropy_prob^8)
  
  candidates$id <- 1:nrow(candidates)
  
  in_sample <- sample(1:nrow(candidates), 1, prob = candidates$entropy_prob)
  samp_prob <- candidates$entropy_prob[in_sample]
  samp_entropy <- candidates$entropy[in_sample]
  
  if(batch_size > 1){
  for(i in 1:(batch_size-1)){
    
    candiates_in_sample <- candidates[in_sample,]
    candiates_not_in_sample <- candidates[-in_sample,]
    
    # First calc distance between the in_sample and the rest
    nn <- nn2(candiates_in_sample[,c("lng", "lat")], candiates_not_in_sample[,c("lng", "lat")])
    
    # convert distances to selection probability
    min_dist_to_other_points <- apply(nn$nn.dists, 1, min)
    min_dist_to_other_points <- min_dist_to_other_points / sum(min_dist_to_other_points)
    
    # Calc mean distance between each potential site and its 10 nearest neighbours
    # mean_dist_to_nearest_neighbours <- nn2(candiates_not_in_sample[,c("lng", "lat")], candiates_not_in_sample[,c("lng", "lat")])$nn.dist[,2:10]
    # mean_dist_to_nearest_neighbours_inv <- 1 / apply(mean_dist_to_nearest_neighbours, 1, mean)
    # mean_dist_to_nearest_neighbours_inv <- mean_dist_to_nearest_neighbours_inv / sum(mean_dist_to_nearest_neighbours_inv)
    
    # Multiply by entropy
    candiates_not_in_sample$pen_entropy <- candiates_not_in_sample$entropy_prob * min_dist_to_other_points #* mean_dist_to_nearest_neighbours_inv
    candiates_not_in_sample$pen_entropy_prob <- candiates_not_in_sample$pen_entropy / sum(candiates_not_in_sample$pen_entropy)
    
    # Optional weight to add to probability sample
    #candiates_not_in_sample$pen_entropy_prob <- candiates_not_in_sample$pen_entropy_prob^3 / 
    #  sum(candiates_not_in_sample$pen_entropy_prob^3)
    
    # Sample
    entropy_sample <- sample(1:nrow(candiates_not_in_sample), 1, prob = candiates_not_in_sample$pen_entropy_prob)
    next_site <- candiates_not_in_sample$id[entropy_sample]
    samp_prob <- c(samp_prob, candiates_not_in_sample$pen_entropy_prob[entropy_sample])
    samp_entropy <- c(samp_entropy, candiates_not_in_sample$entropy[entropy_sample])
    in_sample <- c(in_sample,next_site)
    }
  }
  return(list(sample = candidates[in_sample,1:ncol],
              probs = samp_prob,
              entropy = samp_entropy))
}

