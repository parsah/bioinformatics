source("tableFunctions.R")

toCytoscapeNodeAttributes <- function(x, rank.cutoff=200, out='./node-attribs.csv') {
  # Writes contents of a matrix such that the number of times a PWMs
  # ordered rank is below the cutoff is enumerated.
  #
  # Args:
  #   x: Matrix of dimensions i * j.
  #   rank.cutoff: Weights equal or below this rank are enumated.
  #   out: Output filename to save attributes to.
  
  x.ordered <- orderMatrix(x) # order the matrix
  counts <- as.matrix(rowSums(x.ordered <= rank.cutoff))
  counts[counts == 0] <- 0
  df <- data.frame('Seq'=rownames(counts), 'Num.lines'=counts)
  write.csv(df, out, row.names=F)  
}

toGephiNodeAttributes <- function(x, rank.cutoff=200, out='./node-attribs.csv') {
  # Writes contents of a matrix such that the number of times a PWMs
  # ordered rank is below the cutoff is enumerated.
  #
  # Args:
  #   x: Matrix of dimensions i * j.
  #   rank.cutoff: Weights equal or below this rank are enumated.
  #   out: Output filename to save attributes to.
  
  x.ordered <- orderMatrix(x) # order the matrix
  counts <- as.matrix(rowSums(x.ordered <= rank.cutoff))
  counts[counts == 0] <- 0
  df <- data.frame('Seq'=rownames(counts), 'Num.lines'=counts)
  write.csv(df, out, row.names=F)  
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

toCytoscapeEdgeAttributes <- function(x, out='./edge-attribs.csv') {
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

toGephiEdgeAttributes <- function(x, out='./edge-attribs.csv') {
  # Writes contents of a melted matrix such that matrix correlations
  # represent edge attributes.
  #
  # Args:
  #   x: Matrix (melted) of dimensions i * 3.
  #   out: Output filename to save attributes to.
  
  stopifnot(dim(x)[2] == 3) # melted matrices must be 3x columns wide; 1, 2 => IDs, 3 => value
  sif.data <- data.frame('Source'=x[[1]], 'Target'=x[[2]], 'Correlation'=x[[3]])
  write.csv(sif.data, out, quote=F, row.names=F)
}
