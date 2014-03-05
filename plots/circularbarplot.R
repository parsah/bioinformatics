library("plyr")
library("ecolitk")

prepCircularBarChart <- function(m, k=12) {
  # Prepares a user-provided matrix so that it can be mapped
  # and visualized as a circular-stacked plot. This matrix
  # is first clustered using k-means clustering, with intent
  # of each cluster being represented as its own unique color.
  #
  # Args:
  #   m: Matrix of dimensions i * j.
  #   k: Number of user-provided k-means clusters.
  #
  # Returns:
  #   data-frame of details used for generating circular plots.
  
  km <- kmeans(m, centers=k) # perform k-means clustering on matrix
  df <- data.frame('cluster' = km$cluster, m) # add a new clusters column
  
  # enumerate how many data-points are in each cluster.
  clus.sizes <- ddply(df, .(cluster), function(x) { length(x$cluster)}) 
  colnames(clus.sizes) <- c('cluster', 'n.items')
  position <- as.vector(unlist(apply(as.matrix(clus.sizes$n.items), 1, seq))) # set integer positions per cluster.
  df <- as.matrix(df[order(df$cluster, decreasing=F), ])
  df <- cbind(position, df)
  return(as.data.frame(df))
}

plotCircularBarChart = function(df, color=rainbow(ncol(df)-2), circ.height=0.8, spacing=50) {
  # Takes a data-frame and ultimately represents its columnar contents as a 
  # circular barchart. This data-frame must contain two mandatory columns: position
  # and cluster, respectively. The cluster column is indicative of which cluster each
  # data-point maps to. There are N data-points per cluster, with each integer in the
  # range of N indicative of its respective position. Upon completion, a visual
  # graph will be created made of rings as-per cell-line.
  #
  # Args:
  #   df: Data-frame with required 'position' and 'cluster' column, and numeric columns thereof.
  #   color: Color vector for visualizing individual clusters.
  #   circ.height: Degree of bar-chart height; reflective of each ring in the circle.
  #   spacing: Amount of space between each cluster.
  
  plot.new() # create a new plotting window and define sizes.
  plot.window(c(-10, 10), c(-10, 10))
  
  # break data-frame into the chr column and find maximum position value in each.
  lengths = ddply(df, .(cluster), function(x) { max(x$position) })  
  n.chroms = nrow(lengths)
  
  # off-set is with respect to the position of each data-point in the column.
  offsets = cumsum(c(0, lengths[, 2])) + cumsum(c(0, rep(spacing, n.chroms)))
  tot.length = offsets[length(offsets)] + spacing
  scales = circ.height / apply(abs(df[, 3:ncol(df)]), 2, max) # normalize values given largest point.
  
  # create a data-structure for holding colors per cluster
  col.matrix <- matrix(NA, nrow=length(unique(df$cluster)), ncol=2) # col 1 => cluster, col 2 => color
  colnames(col.matrix) <- c('Cluster', 'Color')
  
  for (i in 1:nrow(df)) {
    clus.num <- df$cluster[i]
    col.matrix[clus.num, 1] <- clus.num # set cluster number
    col.matrix[clus.num, 2] <- color[clus.num] # set cluster color
    cat('Row', i, '/', nrow(df), 
        '; cluster', clus.num, ';', color[clus.num], '\n')
    for (j in 3:ncol(df)) {
      # get beginning of each cluster and all positions referencing this region.
      start = offsets[df[i, 2]] + df[i, 1]
      polygonChrom(begin = start, end = start + 1, len.chrom = tot.length, 
                   radius.in = j - 2, 
                   radius.out = (df[i, j] * scales[j - 2]) + j - 2, 
                   col=color[clus.num])
    }
  }
  return(as.data.frame(col.matrix))
}
