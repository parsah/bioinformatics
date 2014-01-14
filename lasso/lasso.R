parseMatrix <- function(f) {
  # Parses a user-provided CSV file that is representative of a 
  # count matrix and returns its target vector and counts.
  # Args:
  #   f: Input CSV file.
  #
  # Returns:
  #   list of x (matrix) and y (vector) given the input CSV file.
  
  df <- read.csv(file=f, header=T) # read matrix
  y <- as.matrix(df$Target) # extract target vector
  df$Length <- NULL # remove length vector
  df$Sequence <- NULL # first column not needed
  df$Target <- NULL # remove target vector when extracted
  df$X <- NULL # an integer column solely for indexing
  x <- as.matrix(df) # everything else are PWM counts
  return(list('x'=x, 'y'=y)) # parse matrix and target vector
}

buildLASSOClassifier <- function(x, y) {
  # Builds a LASSO classifier given a data matrix and target vector.
  # Args:
  #   x: Count matrix of dimensions i * j
  #   y: Target vector of dimensions i * 1
  #
  # Returns:
  #   fit.cv: LASSO classifier.
  
  library("glmnet")
  fit.cv <- cv.glmnet(x, y, type.measure='auc', family='binomial')
  return(fit.cv)
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