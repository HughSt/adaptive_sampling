
## Script to identify characteristics of optimal 
# survey sites 
library(mgcv)
library(leaflet)
library(oro.nifti)
library(FactoMineR)

# load up data
prev_data <- read.csv("/Users/sturrockh/Downloads/hti_data.csv")

# Plot
col_pal <- colorNumeric(tim.colors(), c(0, max(prev_data$theta)))
leaflet() %>% addTiles() %>% addCircleMarkers(prev_data$lng, prev_data$lat, radius = 1, col = col_pal(prev_data$theta))

# Take binomial sample at survey locations (should probably be done each time)
n_ind <- 200


# Set params
initial_sample_size <- 75
batch_size <- 50
n_step <- 1
threshold <- 0.02
prop_corr <- NULL
mean_dist_nearest <- NULL
mean_dist_nearest_n_neighbours <- NULL
mean_dist_nearest_obs <- NULL
mean_entropy <- NULL
mean_dist_nearest_bioclim1_obs <- NULL
mean_dist_nearest_covar_obs <- NULL
mean_dist_nearest_covar <- NULL


for(i in 1:50){
  
        # Take initial sample
        prev_data$npos <- rbinom(nrow(prev_data), n_ind, prev_data$theta)
        prev_data$nneg <- n_ind - prev_data$npos
        prev_data$sample_0 <- 0
        prev_data$sample_0[sample(1:nrow(prev_data), initial_sample_size)] <- 1
        init_sample <- prev_data[prev_data$sample_0==1,]
        
        # Fit model
        gam_mod <- mgcv::gam(cbind(npos, nneg) ~ s(bioclim1, k=-1) +
                               s(bioclim4, k=-1) +
                               s(bioclim12, k=-1) +
                               s(bioclim15, k=-1) +
                               s(lng, lat, k=-1, bs="gp"),
                             data = init_sample,
                             family = "binomial")
        
        # Id candidate locations
        candidate_locations <- prev_data[prev_data$sample_0==0,]
        exceeds <- exceedance_prob(gam_mod, candidate_locations, 500, threshold)
        entropy <- 0.5 - abs(0.5 - exceeds)
        
        # Take a random sample of batch_size * n_step, add to sample and update model
        sample_rows <- sample(1:nrow(candidate_locations), batch_size * n_step)
        batch_sample_random <- candidate_locations[sample_rows,]
        new_data_random <- rbind(init_sample, batch_sample_random)
        
        # Update model 
        gam_mod_rand <- mgcv::gam(cbind(npos, nneg) ~ s(bioclim1, k=-1) +
                                    s(bioclim4, k=-1) +
                                    s(bioclim12, k=-1) +
                                    s(bioclim15, k=-1) +
                                    s(lng, lat, k=-1, bs="gp"),
                                  #select=TRUE,
                                  data = new_data_random,
                                  family = "binomial")
        
        # Calculate proportion correct
        hotspot_prediction <- as.numeric(predict(gam_mod_rand, prev_data, type="response") > threshold)
        prop_corr <- c(prop_corr, mean(hotspot_prediction == as.numeric(prev_data$theta > threshold)))
        
        # get characteristics of the random sample
        mean_dist_nearest <- c(mean_dist_nearest, 
                               mean(nn2(batch_sample_random[,c("lng", "lat")],
                                 batch_sample_random[,c("lng", "lat")])$nn.dist[,2]))
        
        mean_dist_nearest_n_neighbours <- c(mean_dist_nearest_n_neighbours,
                                            mean(nn2(batch_sample_random[,c("lng", "lat")],
                                                     batch_sample_random[,c("lng", "lat")])$nn.dist[,2:10]))
        
        mean_dist_nearest_pred <- c(mean_dist_nearest, 
                                    mean(nn2(batch_sample_random[,c("lng", "lat")],
                                             batch_sample_random[,c("lng", "lat")])$nn.dist[,2]))
        
        mean_dist_nearest_obs <- c(mean_dist_nearest_obs,
                                   mean(nn2(init_sample[,c("lng", "lat")],
                                      batch_sample_random[,c("lng", "lat")])$nn.dist[,1]))
        
        mean_entropy <- c(mean_entropy,
                          mean(entropy[sample_rows]))
        mean_dist_nearest_bioclim1_obs <- c(mean_dist_nearest_bioclim1_obs,
                                            mean(nn2(init_sample$bioclim1, batch_sample_random$bioclim1)$nn.dist[,1]))
        
        mean_dist_nearest_covar <-c(mean_dist_nearest_covar, 
                                        mean(nn2(batch_sample_random[,c("bioclim1", "bioclim4", "bioclim12", "bioclim15")], 
                                                 batch_sample_random[,c("bioclim1", "bioclim4", "bioclim12", "bioclim15")])$nn.dist[,1]))
        
        
        mean_dist_nearest_covar_obs <-c(mean_dist_nearest_covar_obs, 
                                        mean(nn2(init_sample[,c("bioclim1", "bioclim4", "bioclim12", "bioclim15")], 
                                           batch_sample_random[,c("bioclim1", "bioclim4", "bioclim12", "bioclim15")])$nn.dist[,1]))

        print(i)
}

# Look for associations
par(mfrow = c(2,3))
scatter.smooth(mean_dist_nearest, prop_corr)
scatter.smooth(mean_dist_nearest_obs, prop_corr)
scatter.smooth(mean_entropy, prop_corr)
scatter.smooth(mean_dist_nearest_bioclim1_obs, prop_corr)
scatter.smooth(mean_dist_nearest_covar_obs, prop_corr)
mod<-gam(prop_corr ~ mean_dist_nearest + 
           mean_dist_nearest_obs + 
           mean_entropy + 
           mean_dist_nearest_bioclim1_obs + 
           mean_dist_nearest_covar_obs,
         select=TRUE)
summary(mod)
  