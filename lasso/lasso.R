library("glmnet")

parseCSV <- function(f) {
  # Parses a user-provided CSV file that is representative of a 
  # count matrix and returns its target vector and counts.
  # Args:
  #   f: Input CSV file.
  #
  # Returns:
  #   list of x (matrix) and y (vector) given the input CSV file.
  
  df <- read.csv(file=f, header=T) # read matrix
  y <- as.matrix(df$Target) # extract target vector
  colnames(y) <- c('Target')
  seqs <- df$Sequence
  df$Length <- NULL # remove length vector
  df$Target <- NULL # remove target vector when extracted
  df$X <- NULL # an integer column solely for indexing
  df$Sequence <- NULL # first column not needed
  x <- as.matrix(df) # everything else are PWM counts
  rownames(x) <- seqs
  return(list('x'=x, 'y'=y)) # parse matrix and target vector
}

buildLASSOClassifier <- function(x, y, nfold=3) {
  # Builds a LASSO classifier given a data matrix and target vector.
  # Args:
  #   x: Count matrix of dimensions i * j
  #   y: Target vector of dimensions i * 1
  #   nfold: Number of cross-validation (CV) folds; default = 3.
  # Returns:
  #   fit.cv: LASSO classifier.
  
  stopifnot(length(table(y)) == 2) # y vector must have 0 or 1 values only 
  fit.cv <- cv.glmnet(x, y, nfolds=nfold, type.measure='auc', family='binomial')
  return(fit.cv)
}

getRatios <- function(x, y) {
  # Derives count-ratios per attribute within the count matrix. Resultant 
  # ratios capture attribute abundance within the target vector and serve
  # to shed light on attribute over-representation.
  # Args:
  #   x: Count matrix of dimensions i * j
  #   y: Target vector of dimensions i * 1
  #
  # Returns:
  #   m: Matrix containing attribute ratios and count enumerations.
  
  stopifnot(length(table(y)) == 2) # y vector must have 0 or 1 values only 
  cols <- colnames(x)
  m <- matrix(nrow=dim(x)[2], ncol=3) # columns: ratio, #/query, #/control
  for (i in 1: length(cols)) {
    col.data <- x[, i] # get matrix data for respective column
    num.query <- 0
    num.control <- 0
    for (j in 1: length(col.data)) {
      if (y[j] == 1) { # if data-point target is 1, increment num.query
        num.query <- num.query + col.data[j]
      }
      else {
        num.control <- num.control + col.data[j] # increment control
      }
    }
    ratio <- round(num.query / num.control, 4) # derive ratio
    m[i, ] <- c(ratio, num.query, num.control) # add results to matrix
  }
  colnames(m) <- c('Ratio', 'Num.Query', 'Num.Control')
  rownames(m) <- cols # set matrix columns
  m <- m[order(m[, 1]), ] # sort the columns by their ratio
  return(m)
}

homogenize <- function(x, y, iter=1, nfold=5, threshold=0.5) {
  # Homogenization of the count-matrix and corresponding target vector
  # (see Marlikar et. al., Genome Research, 2010) facilitates iterative 
  # n-fold cross-validation and derivation of an AUC per iteration.
  # Args:
  #   x: Count matrix of dimensions i * j
  #   y: Target vector of dimensions i * 1
  #   iter: Number of iterations to perform
  #   nfold: Number of cross-validations to perform
  
  stopifnot(length(table(y)) == 2) # y vector must have 0 or 1 values only 
  matrix.control <- x[which(y == 0), ]
  matrix.query <- x[which(y == 1), ]
  perc.train <- (nfold -1) / nfold # fraction as-to how much of matrix is for testing
  all.preds <- matrix(nrow=nrow(x), ncol=iter) # create matrix of all LASSO predictions
  rownames(all.preds) <- rownames(x) # save row names
  
  for (i in 1: iter) {
    # extract testing, training data from the shuffled control matrix
    #split.control <- splitMatrix(matrix.control, perc.train)
    #split.query <- splitMatrix(matrix.query, perc.train)
    
    query.train <-   matrix.query[sample.int(nrow(matrix.query), size=floor(nrow(matrix.query) * perc.train)), ]
    control.train <- matrix.control[sample.int(nrow(matrix.control), size=floor(nrow(matrix.control) * perc.train)), ]
    sample.x <- rbind(query.train, control.train) # join two matrices
    sample.y <- as.matrix(c(rep(1, nrow(query.train)),  
                            rep(0, nrow(control.train)))) # add target vector
    cat("[ Iteration ", i, '/', iter, ' ] ... building classifier.\n')
    sample.fit.cv <- buildLASSOClassifier(x=sample.x, y=sample.y, nfold=nfold)
    preds <- predict(sample.fit.cv, x, type='response', s=c('lambda.min')) # derive predictions
    all.preds[, i] <- preds # save predictions produced from iteration to column i
  }
  all.preds <- all.preds >= threshold # assign boolean for each iteration's prediction
  seq.results <- rowSums(all.preds) / iter # enumerate all predictions using probability
  only.true.preds <- seq.results >= threshold # determine if the predictions passed threshold
  indices.true.preds <- which(only.true.preds == TRUE) # get indices of valid predictions
  return(list('x'=x[indices.true.preds, ], 'y'=y[indices.true.preds])) # return valid data-points
}

splitMatrix <- function(x, perc) {
  # Splits a matrix into a training a test data-set given a training-set fraction.
  # For example: if the fraction is 0.8, then 80% of the matrix is designated as the
  # training set, while the remaining 20% is designated as the testing set.
  # Args:
  #   x: Matrix of dimensions i * j
  #   perc: Fraction representing percentage of matrix rows, i, being training set.
  # 
  # Returns:
  #   l: list representing training and testing data-sets.
  
  shuffled.matrix <- x[sample.int(nrow(x), replace=F), ] # shuffle provided matrix
  idx <- floor(perc * nrow(shuffled.matrix)) # find index representing percentile
  train.matrix <- x[1: idx, ] # derive training matrix
  test.matrix <- x[(idx+1): nrow(x), ] # derive smaller testing matrix
  l <- (list('test'=test.matrix, 'train'=train.matrix))
  return(l)
}

getWeights <- function(fit.cv) {
  # Yields attribute weights (coefficients) given LASSO classifier.
  # Args:
  #   fit.cv: LASSO cross-validated (CV) classifier.
  #
  # Returns:
  #   weights: Dataframe of attributes and their respective weights.

  m <- as.matrix(coef(fit.cv))
  m <- as.matrix(m[order(m), ]) # sort counts by their weight
  df <- data.frame('PWM' = rownames(m), 'Weight' = m[, 1])
  rownames(df) <- NULL
  return(df)
}

getAUC <- function(fit.cv) {
  # Identifies the Area Under Curve (AUC) of a LASSO classifier.
  # Args:
  #   fit.cv: LASSO cross-validated (CV) classifier.
  #
  # Returns:
  #   auc: Classifier AUC.
  
  min.lambda <- log(fit.cv$lambda.min) # lambda of maximal AUC; i.e. minimal CV error
  pos.auc <- which(log(fit.cv$lambda) == log(fit.cv$lambda.min), arr.ind=T) # get min lambda position
  auc <- fit.cv$cvm[pos.auc]
  return(auc)
}