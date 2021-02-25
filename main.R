# NOTE: you need to install the used packages first.
# install.packages("ggplot2")
# install.packages("openxlsx")
# install.packages("scales")

#Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8") # To have month names in English.
library(ggplot2)
library(openxlsx)
library(scales)


# Read the raw data.
# dat_raw <- read.csv("sat_modis_proc_l4/slush-limit_output_table.csv")
dat_raw <- openxlsx::read.xlsx("D:/MODIS/sat_modis_proc_SW_l4/_slush-limit_output_table.xlsx", sheet = "Sheet1", colNames = TRUE, detectDates = TRUE)
dat_raw$date <- as.POSIXct(dat_raw$date, format = "%Y-%m-%d")

# Setup the loop: which years and stripes are available?
years   <- sort(unique(dat_raw$year))
stripes <- sort(unique(dat_raw$stripe), decreasing = TRUE)


# Iterate over the years and then (inner loop) over the stripes.
for (year_id in 1:length(years)) {
  
  year_removed_X <- NULL # For the plot of all stripes in a year.
  
  for (stripe_id in 1:length(stripes)) {

    year_cur <- years[year_id]
    stripe_cur <- stripes[stripe_id]
    
    points_cur_ids <- which((dat_raw$year == year_cur) & (dat_raw$stripe == stripe_cur))
    
    # There is one stripe with no retrievals!
    if (length(points_cur_ids) > 0) {
    
      points_cur <- dat_raw[points_cur_ids,]
    
      # Find conflicts: a data point within a stripe conflicts
      # with all other points which are both (significantly) higher and earlier,
      # as well as all which are both (significantly) lower and later.
      # (Significantly) means that we actually have a tolerance, since the
      # assignment of a point to an elevation band can be slightly inaccurate.
      # We make a square table with one column per retrieval, and also one column
      # per retrieval (compared with the one of the column).
      # The final number of conflicts for a point is the sum of its COLUMN.
      # The table is no longer necessarily symmetric (as in an early version):
      # a point with a good quality has some elevation tolerance before being assigned
      # a conflict (to its column), while a point with low quality has less to no tolerance,
      # so that in a conflicting pair the conflict could also be assigned to only one
      # of the two points, depending on their quality.
      # I.e., a good point can be a bit lower than an earlier one, but a bad point cannot.
      # For points with quality index > 11 (top 5 %), we give 85 m elevation tolerance (4 bands).
      # For points with quality index in [8,11] (75th-95th percentiles), we give 65 m elevation tolerance (3 bands).
      # For points with quality index in [5,8] (25th-75th percentiles), we give 45 m elevation tolerance (2 bands).
      # For points with quality index in [3,5] (5th-25th percentiles), we tive 25 m elevation tolerance (1 band).
      # For points with quality index < 3 (bottom 5 %), we give 5 m tolerance (i.e. no bands, but safety in case of small numerical inaccuracies).
      # NOTE: tolerance could also be computed using the quality index RELATIVE to the other points.
      points_cur_tolerances <- c(5,25,45,65,85)[as.integer(cut(points_cur$SL_qindex, c(0,3,5,8,11,Inf)))]
      
      # Now find conflicts.
      npoints_cur <- length(points_cur[,1])
      conflicts_cur <- matrix(data = NA, nrow = npoints_cur, ncol = npoints_cur)
      points_cur$n_conflicts <- NA
      for (point_id in 1:npoints_cur) {
        
        conflicts_cur[,point_id] <- ((points_cur$date < points_cur$date[point_id]) &
                                     (points_cur$SL - points_cur_tolerances[point_id] > points_cur$SL[point_id])) |
                                    ((points_cur$date > points_cur$date[point_id]) &
                                     (points_cur$SL + points_cur_tolerances[point_id] < points_cur$SL[point_id]))
        points_cur$n_conflicts[point_id] <- length(which(conflicts_cur[,point_id]))
      }
      
      # Create "pruned" data frame which we will strip of conflicting points.
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
        points_cur_pruned_rate_with_next <- points_cur_pruned_diff_sl / points_cur_pruned_diff_date # Rate = rise amount / duration.
        points_cur_pruned_rate_with_next_abs <- abs(points_cur_pruned_rate_with_next)
        points_cur_pruned_rate_mean <- (c(points_cur_pruned_rate_with_next_abs[1],points_cur_pruned_rate_with_next_abs) + c(points_cur_pruned_rate_with_next_abs,points_cur_pruned_rate_with_next_abs[length(points_cur_pruned_rate_with_next_abs)]))/2 # We compute the mean of the rates w.r.t. the points immediately before and after.
        
        # Compute the final point score. It corresponds to the number of conflicts (effect
        # moderated by addition of a constant, if there are conflicts),
        # corrected by the rate of rise (we favor slow rise, probably over-simplified).
        # We also add the *global mean* of the rate of rise to reduce the effect of each
        # individual rate of rise.
        # New correction: we also consider the point quality index.
        # This was already considered while establishing tolerances for conflict assessment,
        # but we re-introduce it here to help distinguish between points with e.g. a same
        # number of conflicts. We moderate its effect on the final score by summing a
        # constant to it.
        # points_cur_pruned$score <- points_cur_pruned$n_conflicts * (mean(points_cur_pruned_rate_mean) + points_cur_pruned_rate_mean) 
        # points_cur_pruned$score <- (points_cur_pruned$n_conflicts > 0) * (2 + points_cur_pruned$n_conflicts) * (mean(points_cur_pruned_rate_mean) + points_cur_pruned_rate_mean) * (20 - points_cur_pruned$SL_qindex)
        points_cur_pruned$score <- (points_cur_pruned$n_conflicts > 0) * (2 + points_cur_pruned$n_conflicts) * (20 - points_cur_pruned$SL_qindex)
        
        
        id_worst <- which.max(points_cur_pruned$score)
        points_removed_X <- append(points_removed_X, points_cur_pruned$X[id_worst]) # To keep track of the removed points.
        points_cur_pruned <- points_cur_pruned[-id_worst,]
        npoints_cur_pruned <- length(points_cur_pruned[,1])
        points_cur_pruned_tolerances <- c(5,25,45,65,85)[as.integer(cut(points_cur_pruned$SL_qindex, c(0,3,5,8,11,Inf)))]
        
        conflicts_cur <- matrix(data = NA, nrow = npoints_cur_pruned, ncol = npoints_cur_pruned)
        points_cur_pruned$n_conflicts <- NA
        for (point_id in 1:npoints_cur_pruned) {
          conflicts_cur[,point_id] <- ((points_cur_pruned$date < points_cur_pruned$date[point_id]) &
                                       (points_cur_pruned$SL - points_cur_pruned_tolerances[point_id] > points_cur_pruned$SL[point_id])) |
                                      ((points_cur_pruned$date > points_cur_pruned$date[point_id]) &
                                       (points_cur_pruned$SL + points_cur_pruned_tolerances[point_id] < points_cur_pruned$SL[point_id]))
          points_cur_pruned$n_conflicts[point_id] <- length(which(conflicts_cur[,point_id]))
        }
        
        # Intermediate (debug) plot after each point removal.
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
      #cat("IDs of points removed", point_removed_ids, ".\n")
      
      year_removed_X <- append(year_removed_X, points_removed_X) # Keep annual track of the removed Xs.
      
      # Plot the strip.
      mult <- 2
      ggplot() +
        geom_point(data = points_cur, aes(x = date, y = SL, fill = SL_qindex), shape = 21, stroke = 0.1 * mult, size = 1.7 * mult) +
        geom_line(data = points_cur_pruned, aes(x = date, y = SL), size = 0.1 * mult) +
        geom_point(data = points_removed, aes(x = date, y = SL), shape = 4, size = 0.8 * mult, stroke = 0.15 * mult, color = "#000000") +
        geom_text(data = points_cur, aes(x = date, y = SL + 50, label = n_conflicts), size = 1.2 * mult) +
        xlab("Date") +
        ylab("Slush limit [m a.s.l.]") +
        ggtitle(paste("Slush limits year ", year_cur, " - centre latitude ", stripe_cur, sep = "")) +
        scale_y_continuous(limits = c(400,2200), breaks = seq(400, 2200, 200),
                           minor_breaks = seq(400, 2200, 40)) +
        scale_fill_fermenter(name = "Quality\nindex", type = "div", palette = "RdYlGn",
                             direction = 1, limits = c(0,15), breaks = c(2,3,4,5,6,7,8,9,10,11)) +
                             # labels = c("1.0", "3.0","5", "6", "7", "8", "9", "11.0", "13.0")) +
        # scale_fill_gradientn(name = "Initial #\nof conflicts",
                             # colors = c("#888888", "#DE26D7", "#2121CF", "#44AAF0", "#86F6C6", "#D9D951", "#D42E3B", "#000000")) +
        scale_x_datetime(date_labels = "%Y-%m-%d", date_breaks = "1 month", limits = as.POSIXct(paste(year_cur, c(5,10), 1), format = "%Y %m %d")) +
        theme_bw(base_size = 8 * mult) +
        theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1),
              plot.title = element_text(hjust = 0.5),
              panel.grid.minor = element_line(color = "#000000", size = 0.05 * mult),
              panel.grid.major = element_line(color = "#000000", size = 0.1 * mult),
              legend.key.height = unit(0.3 * mult, "in"))
      ggsave(filename = paste("D:/MODIS/sat_modis_proc_SW_l4/_modis_filter3/aSL_", stripe_cur, "_", year_cur, ".png", sep=""), width = 5*mult, height = 3*mult)
    
    }
  }
  
  dat_year <- dat_raw[which(dat_raw$year == year_cur),]
  year_removed_ids <- pmatch(year_removed_X, dat_year$X)
  dat_year_pruned <- dat_year[-year_removed_ids,]
  
}
