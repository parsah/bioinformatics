library("reshape")
library("plyr")

orderMatrix <- function(x) {
  # Iteratively ranks a matrix given its respective columns.
  #
  # Args:
  #   x: Matrix of dimensions i * j.
  #   ordered.mat: Ordered matrix of dimensions i * j.
  # Returns:
  #   ordered.mat: Integer-ranked matrix of dimensions i * j
  
  n.row <- dim(x)[1]; n.col <- dim(x)[2]
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
