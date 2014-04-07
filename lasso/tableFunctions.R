library("affy")
library("plyr")
library("reshape")
library("RColorBrewer")

orderMatrix <- function(x) {
  # Iteratively ranks a matrix given its respective columns.
  #
  # Args:
  #   x: Matrix of dimensions i * j.
  #   ordered.mat: Ordered matrix of dimensions i * j.
  # Returns:
  #   ordered.mat: Integer-ranked matrix of dimensions i * j
  
  n.row <- dim(x)[1]
  n.col <- dim(x)[2]
  ordered.mat <- matrix(nrow=n.row, ncol=n.col) # create matrix of NAs
  for (col.num in 1:n.col) { # iterate over each column.
    col <- x[, col.num]
    ordered.mat[, col.num] <- rank(-col) # rank column and set to matrix
  }
  colnames(ordered.mat) <- colnames(x)
  rownames(ordered.mat) <- rownames(x)
  return(ordered.mat) # return the ordered matrix
}

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

parseReports <- function(...) {
  # Iteratively ranks a matrix given its respective columns.
  #
  # Args:
  #   ... : Reports generated from the generateReport(...) LASSO function.
  # Returns:
  #   merged.m: Matrix of how many unique rows and columns are in all reports.
  
  files <- list(...)
  list.files <- list() # stored processed files.
  col.names <- c()
  for (idx in 1:length(files)) {
    col.names <- cbind(col.names, basename(files[[idx]])) # save column names
    df <- read.csv(files[[idx]]) # read-in the respective report file.
    df$Num.Query[df$Num.Query == 0] <- 1 # replace 0 with 1; avoid div-0 errors
    df$Num.Control[df$Num.Control == 0] <- 1
    df$Ratio <- df$Num.Query / df$Num.Control
    df$Importance <- df$Weight * df$Ratio
    df <- data.frame('Seq'=df$X, df$Weight) # pull-out important columns
    names(df)[2] <- paste(basename(files[[idx]]))
    list.files[[idx]] <- df # add data-frame to list for sorting    
  }
  merged.df <- join_all(dfs=list.files, by='Seq', match='first', type='full')
  rownames(merged.df) <- merged.df$Seq # set Sequence (Seq) i.e. PWM as row-names
  merged.df$Seq <- NULL
  merged.m <- na.omit(as.matrix(merged.df)) # return ordered matrix as well as merged weights
  return(merged.m)
}

pairElements <- function(x, threshold=0.01) {
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

getTopObservations <- function(x, top.n=100) {
  # Identifies the top N ranked values in a matrix. Values with
  # a rank closer to 1 in any matrix column are saved, however
  # all other values are set to 0. In doing so, any row with
  # a row-sum of 0 is removed since no values across all its
  # columns were ranked in the top N.
  #
  # Args:
  #   x: Matrix of dimensions i * j
  #   top.n: integer referencing the top N ranked values to keep.
  #
  # Returns:
  #   list of top N values and their respective values.
  
  ordered.mat <- orderMatrix(x) # order and filter ranks
  ordered.mat[ordered.mat > top.n] <- 0
  values <- x[rowSums(ordered.mat) != 0, ]
  colnames(values) <- seq(1, ncol(values))
  return(list('values'=values, 'ordered.values'=orderMatrix(values)))
}

normalizeMatrix <- function(x, constants) {
  # Normalizes an i * j matrix given a column-vector of j constants.
  #
  # Args:
  #   x: Matrix of dimensions i * j
  #   constants: Vector of constants of length j.
  #
  # Returns:
  #   norm.mat: Normalized matrix of dimensions i * j.
  
  n.cols <- ncol(x)
  norm.mat <- matrix(data=NA, ncol=ncol(x), nrow=nrow(x))
  rownames(norm.mat) <- row.names(x)
  colnames(x) <- seq(1, n.cols) # set column names
  for (i in 1: n.cols) { # normalize each column.
    norm.mat[, i] <- normalize.constant(x[, i], refconstant=constants[i])
  }
  return(norm.mat)
}

drawHeatmap <- function(x, k) {
  # Draws a traditional 2-dimensional heatmap and generates
  # a corresponding side-bar that helps cluster rows into
  # respective groups.
  #
  # Args:
  #   x: Matrix of dimensions i * j
  #   k: Integer referencing the number of k-means clusters.
  
  distances <- dist(x, method='euclidean')
  hclus <- hclust(distances) # cluster observations
  clusters <- cutree(hclus, k)

  # color rows into k clusters
  cols <- colorRampPalette(brewer.pal(k,"Set3"))
  clustcols = cols(k);
  
  # draw heatmap
  heatmap.2(x, lhei=c(1, 6),RowSideColors=clustcols[clusters], cexRow=0.5,
            col=colorRampPalette(c('orange', 'white', 'blue'))(100),
            trace=NULL, Colv='none', tracecol=NA, density.info='none')
}
