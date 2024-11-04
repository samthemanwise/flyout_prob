library(baseballr)
library(tidyverse)

# Likelihood model of catching fly ball

# out probability function
out_prob <- function(t, x_ball, y_ball) {
  
  # average xy positions of defense
  fielders <- data.frame(
    position = c("P", "C", "1B", "2B", "3B", "SS", "LF", "CF", "RF"),
    x = c(0, 0, 62.07, 30.77, -63.06, -33.07, -134.84, 0, 133.47),
    y = c(60.5, 0, 92.02, 144.77, 100.94, 143.23, 264.63, 322, 261.96)
  )
  
  # find closest fielder
  closest <- function(fielders, x_ball, y_ball) {
    distances <- sqrt((fielders$x - x_ball)^2 + (fielders$y - y_ball)^2)
    closest_index <- which.min(distances)
    closest <- fielders$position[closest_index]
    return(closest)
  }
  
  closest_fielder <- closest(fielders, x_ball, y_ball)
  
  # xy coords of closest fielder
  x_closest <- fielders$x[fielders$position == closest_fielder]
  y_closest <- fielders$y[fielders$position == closest_fielder]
  
  # distance from closest fielder to ball landing spot
  d = sqrt((x_ball-x_closest)^2+(y_ball-y_closest)^2)
  
  # depths of ball / fielder for behind variable
  ball_depth = sqrt(x_ball^2 + y_ball^2)
  fielder_depth = sqrt(x_closest^2 + y_closest^2)
  
  # accel_speed <- 10 # 7.8 OG Player's average accel in feet per second^2
  # accel_sd <- 1 # Player's standard deviation of speed in feet per second^2
  
  # Lets try jump metric (feet covered in 3 seconds) instead of acceleration
  jump_avg <- 33.3 # ft, statcast average
  jump_sd <- 3 # ft
  
  set.seed(123)
  num_sims <- 10000
  outs <- 0
  
  for (i in 1:num_sims) {
    
    # accel <- rnorm(1, mean = accel_speed, sd = accel_sd)
    
    jump <- rnorm(1, mean = jump_avg, sd = jump_sd)
    
    # slower speed when running backwards
    if (ball_depth > fielder_depth) {
      max_speed <- 27 # ft/sec
      behind = 1
    } else {
      max_speed <- 29 #ft/sec, league average
      behind = 0
    }

    if (t <= 3) {
      a <- (2 * jump) / 9  # calculated acceleration
      d_covered <- 0.5 * a * t^2
      if (d_covered >= d) {
        outs <- outs + 1
      }
    } else {
      if (jump >= d) {
        outs <- outs + 1
      } else {
        remaining_d <- d - jump
        time_at_max <- remaining_d / max_speed
        total_time <- 3 + time_at_max
        if (total_time <= t) {
          outs <- outs + 1
        }
      }
    }
    
    # if (t <= jump) {
    #   a <- max_speed/jump
    #   d_covered <- 0.5 * a * t^2
    #   if (d_covered >= d) {
    #     outs <- outs + 1
    #   }
    # } else {
    #   a <- max_speed/jump
    #   d_covered <- 0.5 * a * jump^2
    #   if (d_covered >= d) {
    #     outs <- outs + 1
    #   } else {
    #     remaining_d <- d - d_covered
    #     time_at_max <- remaining_d / max_speed
    #     total_time <- jump + time_at_max
    #     if (total_time <= t) {
    #       outs <- outs + 1
    #     }
    #   }
    # }
  
    # time_to_max <- max_speed / accel
    # 
    # accel_d <- 0.5 * accel * time_to_max^2
    # 
    # if (d <= accel_d) {
    #   total_time <- sqrt(2 * d/accel)
    # } else {
    #   remaining_d <- d - accel_d
    #   time_at_max <- remaining_d / max_speed
    #   total_time <- time_to_max + time_at_max
    # }
    
    # if (total_time <= t) {
    #   outs <- outs + 1
    # }
    
  }
  
  out_probability <- outs / num_sims
  
  out <- rbinom(1, size = 1, prob = out_probability)
  
  # time_diff <- (t - (d / 27))
  
  return_matrix <- cbind(d, behind, out)
  return(return_matrix)
}

# probability of out using Rays home season batted balls

rays_home_new <- rays_home_new %>% filter(t > 1.5)
rays_home_no_75 <- rays_home_no_75  %>% filter(t > 1.5)

# apply function to data set
new_rays_d_probs <- apply(rays_home_no_75, 1, function(row) {
  out_prob(row[1], row[2], row[3])
})

# put probabilities in df with original data

output_df <- data.frame(rays_d_probs)
output_df <- t(data.frame(output_df))
output_df <- cbind(rays_home_new, output_df)
output_df <- data.frame(output_df)
colnames(output_df) <- c("t", "x", "y", "d", "behind", "time_diff", "out")

new_output_df <- data.frame(new_rays_d_probs)
new_output_df <- t(data.frame(new_output_df))
new_output_df <- cbind(rays_home_no_75, new_output_df)
new_output_df <- data.frame(new_output_df)
colnames(new_output_df) <- c("t", "x", "y", "d", "behind", "out")

new_output_df %>%
  ggplot() + 
  geom_point(mapping = aes(x = x, y = y, color = out)) +
  scale_color_gradient(low = "blue", high = "red")

write.csv(new_output_df, "new_rays_fly_outs_jump.csv", row.names = FALSE)

count_result <- table(new_output_df$out)
count_value <- count_result[1]




