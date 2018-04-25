# plot a single-subject design

plot_MBD <- function(yplot, names, dv, lims, filename) {
  require(ggplot2)
  library(gridExtra)
  n.cases <- nrow(yplot)/2 #each case has two rows, baseline and int phase
  counter <- 1
  plots <- list()
  for (l in 1:n.cases){
    y <- yplot[counter:(counter+1), ]
    P <- nrow(y)  
    T <- apply(y, 1, function(x) max(which(!is.na(x))))
    
    observations <- rep(NA, sum(T))
    k = 1
    for (i in 1:P) {
      for (j in 1:T[i]) {
        observations[k] = y[i, j]
        k = k + 1
      }
    }
    changepoints <- cumsum(T)
    
    DF <- data.frame(
      phase = rep(seq(P), times = T),
      time = seq(sum(T)),
      value = observations
    )
    
    plots[[l]] <- ggplot(DF, aes(x = time, y = value)) +
      geom_point() +
      geom_line(aes(group = phase)) +
      geom_vline(xintercept = head(changepoints, -1) + 0.5, linetype = "solid") +
      theme_bw() +
      labs(y = dv) +
      labs(title = names[l]) + xlim(1, lims) + ylim(min(yplot, na.rm=T),
                                                    max(yplot, na.rm=T))
    counter <- counter + 2
  }
  jpeg(paste0(filename,".jpg"))
  do.call(grid.arrange,  plots)
  dev.off()
}
