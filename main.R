Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8")
library(ggplot2)

dat_raw <- read.csv("sat_modis_proc_l4/slush-limit_output_table.csv")
dat_raw$date <- as.POSIXct(dat_raw$date, format = "%Y-%m-%d")

years   <- sort(unique(dat_raw$year))
stripes <- sort(unique(dat_raw$stripe), decreasing = TRUE)


for (year_id in 1:length(years)) {
  
  year_removed_X <- NULL # For the plot of all stripes in a year.
  
  for (stripe_id in 1:length(stripes)) {

    year_cur <- years[year_id]
    stripe_cur <- stripes[stripe_id]
    
    points_cur_ids <- which((dat_raw$year == year_cur) & (dat_raw$stripe == stripe_cur))
    
    if (length(points_cur_ids) > 0) {
    
      points_cur <- dat_raw[points_cur_ids,]
      
      
      # Find conflicts: a data point within a stripe conflicts
      # with all other points which are both higher and earlier,
      # as well as all which are both lower and later.
      # We make a (symmetric) square table where a cell is TRUE if the
      # two points at the corresponding row and column are in conflict.
      npoints_cur <- length(points_cur[,1])
      conflicts_cur <- matrix(data = NA, nrow = npoints_cur, ncol = npoints_cur)
      points_cur$n_conflicts <- NA
      for (point_id in 1:npoints_cur) {
        conflicts_cur[,point_id] <- ((points_cur$date < points_cur$date[point_id]) & (points_cur$SL > points_cur$SL[point_id])) | ((points_cur$date > points_cur$date[point_id]) & (points_cur$SL < points_cur$SL[point_id]))
        points_cur$n_conflicts[point_id] <- length(which(conflicts_cur[,point_id]))
      }
      
      # Create pruned data frame which we will strip of conflicting points.
      points_cur_pruned <- points_cur
      points_removed_X <- NULL
      
      # Iterative removal of the worst conflicting point.
      # The worst conflicting point is not simply the one
      # with the largest number of conflicts.
      # We also consider the rate of change of the slush limit
      # compared to the neighbors.
      while (any(points_cur_pruned$n_conflicts) > 0) {
        
        cat("Still", sum(points_cur_pruned$n_conflicts), "conflict(s) to solve...\n")
        
        # Analysis of the rate of rise of the slush limit.
        points_cur_pruned_diff_sl <- diff(points_cur_pruned$SL)
        points_cur_pruned_diff_date <- as.numeric(diff(points_cur_pruned$date))
        points_cur_pruned_rate_with_next <- points_cur_pruned_diff_sl / points_cur_pruned_diff_date
        points_cur_pruned_rate_with_next_abs <- abs(points_cur_pruned_rate_with_next)
        points_cur_pruned_rate_mean <- (c(points_cur_pruned_rate_with_next_abs[1],points_cur_pruned_rate_with_next_abs) + c(points_cur_pruned_rate_with_next_abs,points_cur_pruned_rate_with_next_abs[length(points_cur_pruned_rate_with_next_abs)]))/2
        
        points_cur_pruned$score <- points_cur_pruned$n_conflicts * (mean(points_cur_pruned_rate_mean) + points_cur_pruned_rate_mean) # We add the mean to decrease the effect of the rate of rise (else a 3-element cluster outlier with zero rate of rise could remain there forever).
        
        id_worst <- which.max(points_cur_pruned$score)
        points_removed_X <- append(points_removed_X, points_cur_pruned$X[id_worst]) # To keep track of the removed points.
        points_cur_pruned <- points_cur_pruned[-id_worst,]
        npoints_cur_pruned <- length(points_cur_pruned[,1])
        
        conflicts_cur <- matrix(data = NA, nrow = npoints_cur_pruned, ncol = npoints_cur_pruned)
        points_cur_pruned$n_conflicts <- NA
        for (point_id in 1:npoints_cur_pruned) {
          conflicts_cur[,point_id] <- ((points_cur_pruned$date < points_cur_pruned$date[point_id]) & (points_cur_pruned$SL > points_cur_pruned$SL[point_id])) | ((points_cur_pruned$date > points_cur_pruned$date[point_id]) & (points_cur_pruned$SL < points_cur_pruned$SL[point_id]))
          points_cur_pruned$n_conflicts[point_id] <- length(which(conflicts_cur[,point_id]))
        }
        
        # ggplot(points_cur_pruned) +
        #   geom_point(aes(x = date, y = SL, fill = n_conflicts), shape = 21, stroke = 0, size = 3) +
        #   geom_line(aes(x = date, y = SL), size = 0.1) +
        #   geom_text(aes(x = date, y = SL + 50, label = n_conflicts), size = 2) +
        #   scale_fill_gradientn(colors = c("#888888", "#DE26D7", "#2121CF", "#44AAF0", "#86F6C6", "#F9F961", "#D42E3B", "#000000")) +
        #   scale_x_datetime(date_labels = "%Y-%m-%d") +
        #   theme_bw()
      }
      
      npoints_removed <- length(points_cur[,1]) - length(points_cur_pruned[,1])
      cat("All conflicts solved! Removed", npoints_removed, "point(s).\n")
      
      point_removed_ids <- pmatch(points_removed_X, points_cur$X)
      points_removed <- points_cur[point_removed_ids,]
      
      year_removed_X <- append(year_removed_X, points_removed_X) # Keep annual track of the removed Xs.
      
      
      ggplot() +
        geom_point(data = points_cur, aes(x = date, y = SL, fill = n_conflicts), shape = 21, stroke = 0, size = 1.7) +
        geom_line(data = points_cur_pruned, aes(x = date, y = SL), size = 0.1) +
        geom_point(data = points_removed, aes(x = date, y = SL), shape = 4, size = 0.9, stroke = 0.3, color = "#FFFFFF") +
        geom_text(data = points_cur, aes(x = date, y = SL + 50, label = n_conflicts), size = 1.2) +
        xlab("Date") +
        ylab("Slush limit [m a.s.l.]") +
        ggtitle(paste("Slush limits year ", year_cur, " - centre latitude ", stripe_cur, sep = "")) +
        ylim(400,2200) +
        scale_fill_gradientn(name = "Initial #\nof conflicts",
                             colors = c("#888888", "#DE26D7", "#2121CF", "#44AAF0", "#86F6C6", "#D9D951", "#D42E3B", "#000000")) +
        scale_x_datetime(date_labels = "%Y-%m-%d", date_breaks = "1 month", limits = as.POSIXct(paste(year_cur, c(5,10), 1), format = "%Y %m %d")) +
        theme_bw(base_size = 8) +
        theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1),
              plot.title = element_text(hjust = 0.5))
      ggsave(filename = paste("modis_filter2/aSL_", stripe_cur, "_", year_cur, ".png", sep=""), width = 5, height = 3)
      
      
      

      # 
      # 
      # ggplot(df_rate) +
      #   geom_point(aes(x = date, y = rate)) +
      #   geom_line(aes(x = date, y = rate), size = 0.1) +
      #   scale_x_datetime() +
      #   theme_bw()
    
    }
  }
  
  dat_year <- dat_raw[which(dat_raw$year == year_cur),]
  year_removed_ids <- pmatch(year_removed_X, dat_year$X)
  dat_year_pruned <- dat_year[-year_removed_ids,]
  
  # ggplot() +
  #   geom_point(data = points_cur, aes(x = date, y = SL, fill = n_conflicts), shape = 21, stroke = 0, size = 2) +
  #   geom_text(data = points_cur, aes(x = date, y = SL + 50, label = n_conflicts), size = 2) +
  #   geom_line(data = points_cur_pruned, aes(x = date, y = SL), size = 0.1) +
  #   geom_point(data = points_removed, aes(x = date, y = SL), shape = 4, size = 1, stroke = 0.5, color = "#FFFFFF") +
  #   xlab("Date") +
  #   ylab("Slush limit [m]") +
  #   ylim(400,2200) +
  #   scale_fill_gradientn(name = "Initial #\nof conflicts",
  #                        colors = c("#888888", "#DE26D7", "#2121CF", "#44AAF0", "#86F6C6", "#F9F961", "#D42E3B", "#000000")) +
  #   scale_x_datetime(date_labels = "%Y-%m-%d", date_breaks = "1 month") +
  #   theme_bw()
  # ggsave(filename = paste("modis_filter1/aSL_", stripe_cur, "_", year_cur, ".png", sep=""), width = 5, height = 3)
  # 
  
}




