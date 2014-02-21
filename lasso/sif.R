library("reshape")

parseCounts <- function(f) {
  # Parses a user-provided CSV file that is representative of a 
  # count matrix.
  #
  # Args:
  #   f: Input CSV file.
  #
  # Returns:
  #   x: Count matrix of dimensions i * j.
  
  df <- read.csv(file=f, header=T) # read matrix
  seqs <- df$Sequence
  df$Length <- NULL # remove length vector
  df$Target <- NULL # remove target vector when extracted
  df$X <- NULL # an integer column solely for indexing
  df$Sequence <- NULL # first column not needed
  x <- as.matrix(df) # everything else are PWM counts
  rownames(x) <- seqs
  cor.x <- cor(x) # derive the matrix correlation
  cor.x[lower.tri(cor.x)] <- 0 # remove the lower-triangular segment of the matrix
  return(list('x'=x, 'corr'=cor.x)) # return matrix and its correlation
}

filterMatrix <- function(x, threshold=0.01) {
  # Linearizes a matrix so that each specific individual element can be easily manipulated.
  #
  # Args:
  #   x: Matrix of dimensions i * i.
  #
  # Returns:
  #   v: Vector of length i * j.
  
  melt.m <- melt(x)
  # pull-out positively-correlated or negatively correlated entries.
  melt.m <- melt.m[melt.m$value >= threshold | melt.m$value <= -threshold, ]
  melt.m <- melt.m[melt.m$X1 != melt.m$X2, ] # remove paired correlations, eg. A and A.
  return(melt.m)
}

toSIF <- function(x, out='./network.sif') {
  # Writes contents of a melted matrix as an SIF file.
  #
  # Args:
  #   x: Matrix (melted) of dimensions i * 3.
  #   out: Output filename to save SIF as.
  
  stopifnot(dim(x)[2] == 3) # melted matrices must be 3x columns wide; 1, 2 => IDs, 3 => value
  sif.data <- paste(x[[1]], x[[2]], sep=' pp ')
  write.table(sif.data, out, quote=F, row.names=F, col.names=F)
}

toEdgeAttributes <- function(x, out='./edge-attribs.csv') {
  # Writes contents of a melted matrix such that matrix correlations
  # represent edge attributes.
  #
  # Args:
  #   x: Matrix (melted) of dimensions i * 3.
  #   out: Output filename to save attributes to.
  
  stopifnot(dim(x)[2] == 3) # melted matrices must be 3x columns wide; 1, 2 => IDs, 3 => value
  sif.data <- data.frame('Interaction'= paste(x[[1]], x[[2]], sep=' (pp) '))
  sif.data$Correlation <- x[[3]] # save the corresponding correlations.
  write.csv(sif.data, out, quote=F, row.names=F)
}

parseReports <- function(...) {
#  args <- list(...)
#  for (idx in 1:length(args) {
#    cat(args[idx],'\n')
#  }
}
