library(mgcv)
library(leaflet)
library(oro.nifti)
library(RANN)

rm(list=ls())
source("https://raw.githubusercontent.com/HughSt/adaptive_sampling/master/adaptive_survey_functions.R")
# load up data
prev_data <- read.csv("/Users/sturrockh/Downloads/hti_data.csv")

# Plot
col_pal <- colorNumeric(tim.colors(), c(0, max(prev_data$theta)))
leaflet() %>% addProviderTiles("CartoDB.Positron") %>% 
  addCircleMarkers(prev_data$lng, prev_data$lat, radius = 3, col = col_pal(prev_data$theta)) %>%
  leaflet::addLegend(pal=col_pal, values=prev_data$theta, title="Prevalence")

# Loop the experiment n times

n_ind <- 100
initial_sample_size <- 100
batch_size <- 10
n_step <- 2
threshold <- 0.1
#exceedance_threshold <- 0.2
prop_corr_rand <- NULL
prop_corr_pe <- NULL
sens_rand <- NULL
sens_pe <- NULL
mean_entropy_pe <- NULL
mean_entropy_rand <- NULL

for(i in 1:10){
  
      # first take a binomial sample at each site using prevalence
      prev_data$npos <- rbinom(nrow(prev_data), n_ind, prev_data$theta)
      prev_data$nneg <- n_ind - prev_data$npos

      # Now override the original initial sample and take a new one
      prev_data$sample_0 <- 0
      prev_data$sample_0[sample(1:nrow(prev_data), initial_sample_size)] <- 1
      init_sample <- prev_data[prev_data$sample_0==1,]
      
      # Fit model
      # First calc optimal range param for gp
      REML_estimates <- optimal_range(min_dist = 0.01, 
                                      max_dist = max(dist(init_sample[,c("lng", "lat")])),
                                      model_data = model_data)
      
          # Then fit model
          gam_mod <- mgcv::gam(cbind(npos, nneg) ~ s(bioclim1, k=-1) +
                               s(bioclim4, k=-1) +
                               s(bioclim12, k=-1) +
                               s(bioclim15, k=-1) +
                                s(dist_to_water_m, k=-1) +
                                 s(elev_m, k=-1) +
                                 s(lng, lat, k=-1, bs="gp", m = c(3, REML_estimates$best_m)),
                               method = "REML",
                               data = init_sample,
                               family = "binomial")
      
      # Predict to all other locations and compare results
      #leaflet() %>% addTiles() %>% addCircleMarkers(prev_data$lng, prev_data$lat, radius = 1, col = col_pal(predict(gam_mod, prev_data, type="response")))
      #leaflet() %>% addTiles() %>% addCircleMarkers(prev_data$lng, prev_data$lat, radius = 1, col = col_pal(exceedance_prob(gam_mod, prev_data, 500, threshold)))    
      #plot(predict(gam_mod, prev_data, type="response"), prev_data$theta)
      
      # Generate exceedance probabilities
      candidate_locations <- prev_data[prev_data$sample_0==0,]
      exceeds <- exceedance_prob(gam_mod, candidate_locations, 500, threshold)
      entropy <- 0.5 - abs(0.5 - exceeds)
      #entropy = - exceeds * log(exceeds) - (1-exceeds) * log (1-exceeds)
      
      
      # Take a random sample of batch_size * n_step, add to sample and update model
      batch_sample_random <- candidate_locations[sample(1:nrow(candidate_locations), batch_size * n_step),]
      new_data_random <- rbind(init_sample, batch_sample_random)
      
      # Update model 
      # First calc optimal range param for gp
      gam_mod_rand <- mgcv::gam(cbind(npos, nneg) ~ s(bioclim1, k=-1) +
                             s(bioclim4, k=-1) +
                             s(bioclim12, k=-1) +
                             s(bioclim15, k=-1) +
                               s(dist_to_water_m, k=-1) +
                               s(elev_m, k=-1) +
                             s(lng, lat, k=-1, bs="gp", m = c(3, REML_estimates$best_m)),
                           data = new_data_random,
                           method = "REML",
                           family = "binomial")
      #plot(predict(gam_mod_rand, prev_data, type="response"), prev_data$theta)
      
      ## Adaptive sample ##
      new_data_pen_entropy <- init_sample
      weights <- 1/(rep(1/nrow(prev_data), initial_sample_size))
      
      for(step in 1:n_step){
        
              # Now take sample using penalized entropy
              batch_sample_penalized_entropy <- sample_penalized_entropy(candidate_locations, entropy, batch_size)
              new_data_pen_entropy <- rbind(new_data_pen_entropy, batch_sample_penalized_entropy$sample)
              
              # Calc weights
              weights = c(weights, 1/batch_sample_penalized_entropy$probs)
              
              # Scale to sum to the correct sample size
              corr_weights <- length(weights) * (weights / sum(weights))
              
              # Update model and 
              gam_mod_pen_entropy <- mgcv::gam(cbind(npos, nneg) ~ s(bioclim1, k=-1) +
                                          s(bioclim4, k=-1) +
                                          s(bioclim12, k=-1) +
                                          s(bioclim15, k=-1) +
                                            s(dist_to_water_m, k=-1) +
                                            s(elev_m, k=-1) +
                                          s(lng, lat, k=-1, bs="gp", m = c(3, REML_estimates$best_m)),
                                        data = new_data_pen_entropy,
                                        method = "REML",
                                        weights = corr_weights,
                                        family = "binomial")
              
              # PRedict and get entropy
              candidate_locations <- candidate_locations[-which(candidate_locations$index %in% batch_sample_penalized_entropy$sample$index),]
              exceeds <- exceedance_prob(gam_mod_pen_entropy, candidate_locations, 500, threshold)
              entropy <- 0.5 - abs(0.5 - exceeds)
              #entropy = - exceeds * log(exceeds) - (1-exceeds) * log (1-exceeds)
              }
      
      # Compare results. Classify based on exceedance probabilities
      exceeds_rand <- exceedance_prob(gam_mod_rand, prev_data, 500, threshold)
      # predictions_random <- as.numeric(exceeds_rand >= 0.5)
      # predictions_pen_entropy <- as.numeric(exceedance_prob(gam_mod_pen_entropy, prev_data, 500, threshold) >= 0.5)
      
      # classify based on point estimates
      predictions_random <- as.numeric(predict(gam_mod_rand, prev_data, type="response") > threshold)
      predictions_pen_entropy <- as.numeric(predict(gam_mod_pen_entropy, prev_data, type="response") > threshold)
      
      # Use observed prevalence values at sampled sites to classify site
      predictions_random[which(prev_data$index %in% new_data_random$index)] <- 
        as.numeric((new_data_random$npos[order(new_data_random$index)] / n_ind) > threshold)
      predictions_pen_entropy[which(prev_data$index %in% new_data_pen_entropy$index)] <- 
        as.numeric((new_data_pen_entropy$npos[order(new_data_pen_entropy$index)] / n_ind) > threshold)
      
      prop_corr_rand <- c(prop_corr_rand, mean(predictions_random == as.numeric(prev_data$theta > threshold)))
      sens_rand <- c(sens_rand, sum(predictions_random == 1 & prev_data$theta > threshold) / sum(prev_data$theta > threshold))
      prop_corr_pe <- c(prop_corr_pe, mean(predictions_pen_entropy == as.numeric(prev_data$theta > threshold)))
      sens_pe <- c(sens_pe, sum(predictions_pen_entropy == 1 & prev_data$theta > threshold) / sum(prev_data$theta > threshold))
      mean_entropy_pe <- c(mean_entropy_pe, mean(entropy))
      mean_entropy_rand <- c(mean_entropy_rand, mean(0.5 - abs(0.5 - exceeds_rand)))
      
}

# Compare 
mean(prop_corr_rand)
mean(prop_corr_pe)
mean(sens_rand)
mean(sens_pe)
mean(mean_entropy_rand)
mean(mean_entropy_pe)



