library("glmnet")

parseClassifiable <- function(f) {
  # Parses a user-provided CSV file that is representative of a 
  # count matrix and returns its target vector and counts.
  #
  # Args:
  #   f: Input CSV file.
  #
  # Returns:
  #   x: Count matrix of dimensions i * j.
  #   y: Target vector of dimensions i * 1.
  
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
  return(list('x'=x, 'y'=y)) # return matrix and target vector
}

buildLASSOClassifier <- function(x, y, nfold=3) {
  # Builds a LASSO classifier given a data matrix and target vector.
  #
  # Args:
  #   x: Count matrix of dimensions i * j.
  #   y: Target vector of dimensions i * 1.
  #   nfold: Number of cross-validation (CV) folds; default = 3.
  #
  # Returns:
  #   LASSO classifier.
  
  fit.cv <- cv.glmnet(x, y, nfolds=nfold, type.measure='auc', family='binomial')
  return(fit.cv)
}

getRatios <- function(x, y) {
  # Derives count-ratios per attribute within the count matrix. Resultant 
  # ratios capture attribute abundance within the target vector and serve
  # to shed light on attribute over-representation.
  #
  # Args:
  #   x: Count matrix of dimensions i * j
  #   y: Target vector of dimensions i * 1
  #
  # Returns:
  #   Matrix containing attribute ratios and count enumerations.
  
  stopifnot(length(table(y)) == 2) # y vector must have 0 or 1 values only 
  sum.query <- colSums(x[which(y == 1), ]) # sum rows referencing query targets
  sum.control <- colSums(x[which(y == 0), ]) # repeat for control
  ratio <- round(sum.query / sum.control, 2) # derive ratio
  m <- as.matrix(cbind(ratio, sum.query, sum.control)) # gather results in matrix
  m <- as.matrix(m[order(row.names(m)), ]) # order by row-name
  colnames(m) <- c('Ratio', 'Num.Query', 'Num.Control')
  return(m)
}

getPValues <- function(x, y, adj.method="BH") {
  # Derives adjusted p-values for each count matrix attribute based on 
  # abundance of the observation.
  #
  # Args:
  #   x: Count matrix of dimensions i * j
  #   y: Target vector of dimensions i * 1
  #   method: P-value adjustment method. See p.adjust(method,...)
  #
  # Returns:
  #   P-value (adjusted) vector.
  
  ratios <- getRatios(x, y) # contains feature counts
  total.query <- sum(ratios[, 2]) # sum query column
  total.control <- sum(ratios[, 3]) # sum control column
  p.vals <- matrix(ncol=1, nrow=nrow(ratios)) # for storing p-values
  rownames(p.vals) <- rownames(ratios)

  for (i in 1: nrow(ratios)) { # iteratively compute matrix p-value
    val.query <- ratios[i, 2]
    val.control <- ratios[i, 3]
    m <- matrix(c(val.query, val.control, total.query-val.query, total.control-val.control),
                ncol=2, byrow=T)
    p.vals[i, 1] <- fisher.test(m)$p.value # save p-value
  }
  p.vals[, 1] <- p.adjust(p.vals, method=adj.method) # adjust p-values
  colnames(p.vals) <- c('P.value')
  return(p.vals)
}

toPredictionVector <- function(x, y, iter=1, nfold=5) {
  # A prediction matrix is the result of iteratively (and randomly) building 
  # classifiers and deriving predictions of each observation.
  # Upon completion, each observation references predictions derived per
  # iteration, useful in possible dimensionality reduction.
  #
  # Args:
  #   x: Count matrix of dimensions i * j.
  #   y: Target vector of dimensions i * 1.
  #   iter: Number of iterations to perform.
  #   nfold: Number of cross-validations to perform.
  #
  # Returns:
  #   Vector of predictions of dimensions i * 1. 
  
  stopifnot(length(table(y)) == 2) # y vector must have 0 or 1 values only 
  matrix.control <- x[which(y == 0), ]
  matrix.query <- x[which(y == 1), ]
  perc.train <- (nfold -1) / nfold # fraction as-to how much of matrix is for testing
  list.aucs <- rep(0, iter) # generate vector of AUCs representative of each iteration
  all.preds <- matrix(nrow=nrow(x), ncol=iter) # create matrix of all LASSO predictions
  rownames(all.preds) <- rownames(x) # save row names
  
  for (i in 1: iter) {
    time.start <- Sys.time() # get start time
    query.train <-   matrix.query[sample.int(nrow(matrix.query), size=floor(nrow(matrix.query) * perc.train)), ]
    control.train <- matrix.control[sample.int(nrow(matrix.control), size=floor(nrow(matrix.control) * perc.train)), ]
    sample.x <- rbind(query.train, control.train) # join two matrices
    sample.y <- as.matrix(c(rep(1, nrow(query.train)),  
                            rep(0, nrow(control.train)))) # add target vector
    cat("[ Iteration ", i, '/', iter, ' ] ... ')
    sample.fit.cv <- buildLASSOClassifier(x=sample.x, y=sample.y, nfold=nfold)
    iter.preds <- doPrediction(sample.fit.cv, x) # derive predictions
    all.preds[, i] <- iter.preds # save predictions produced from iteration to column i
    list.aucs[i] <- getAUC(sample.fit.cv) # set the desired AUC per iteration    
    time.iter <- format(.POSIXct(difftime(Sys.time(), time.start, units=c("secs")),tz="GMT"), "%H:%M:%S")
    cat(time.iter, 'h:m:s \n') # print on existing line.
  }
  all.preds <- rowSums(all.preds) / iter # average predictions for each observation
  return(list('preds'=all.preds, 'aucs'=list.aucs)) # return observation-specific predictions
}

generateReport <- function(weights, ratios, out='./report.csv') {
  # Generates a report given a LASSO fit, PWM ratios, and their p-values.
  #
  # Args:
  #   weights: LASSO weights derived using getWeights.
  #   ratios: Feature ratios derived using getRatios
  #   out: Output file to save results to.
  #
  # Returns:
  #   Output results file defined by the out argument.
  
  imp <- weights * ratios[, 1] # importance => weight * ratio
  colnames(imp) <- c('Importance')
  report.mat <- cbind(weights, ratios, imp)

  # test that row-names are shared between weights and ratios
  stopifnot(identical(rownames(weights), rownames(ratios)))
  write.csv(report.mat, out) # write merged martix to a file
}

homogenize <- function(x, y, preds, threshold = 0.5) {
  # Given a count matrix and a target vector, each row is enumerated against 
  # to determine if its predictions cumulatively exceed a user-provided threshold.
  # For query data-points to exceed, their cumulative prediction must exceed the
  # threshold. Query data-points not exceeding this threshold and removed along
  # with their corresponding control.
  # Implementation is based on Narlikar et. al., Genome Research, 2010.
  #
  # Args:
  #   x: Count matrix of dimensions i * j
  #   y: Target vector of dimensions i * 1
  #   preds:  Matrix of observation predictions of dimensions i * k
  #   threshold: prediction cutoff.
  #
  # Returns:
  #   List referencing homogenized count matrix and target vector.
  
  control.counts <- x[which(y == 0), ] # get matrix representing control and query data
  query.counts <- x[which(y == 1), ]
  query.preds <- preds[which(y == 1)] #>= threshold # uninterested in control predictions.
  
  # iterate over each query and check if it has a corresponding control; a
  # test purely determined if the query substring is within the control.
  new.query <- c() # create vector so rows can be added without matrix usage.
  new.control <- c() 
  names.query <-c() # references row-names as it is unknown which observations are/aren't valid.
  names.control <- c()
  rownames.control <- row.names(control.counts)
  
  for (i in 1: nrow(query.counts)) {
    query.name <- row.names(query.counts)[i]
    query.pred <- query.preds[i] # whether the query passes threshold
    match.pos <- grep(query.name, rownames.control, fixed=T) # vector of indices

    # if query prediction is true, add control and query data-points.
    if (query.pred >= threshold) {
      new.query <- rbind(new.query, query.counts[i, ]) # add passed queries.
      names.query <- rbind(names.query, query.name) # add query name
      
      # in the case of small datasets, a query may not have a corresponding
      # control, thus handle instances whereby no match is found and add if so.
      if (length(match.pos) > 0) { # if query maps to control, get first hit.
        control.name <- rownames.control[match.pos[1]] # fetch first hit only
        #cat( '   ',control.name, '\n')
        new.control <- rbind(new.control, control.counts[match.pos, ]) # add its respective control.
        names.control <- rbind(names.control, control.name) # add control name
      }
    }
  }
  
  rownames(new.query) <- names.query # add row-names to each data-frame.
  rownames(new.control) <- names.control
  new.x <- rbind(new.query, new.control) # merge counts into new matrix
  new.y <- as.matrix(c(rep(1, nrow(new.query)), # build target vector
                       rep(0, nrow(new.control))))
  return(list('x'=new.x, 'y'=new.y))
}

saveJob <- function(proj.name) {
  # Saves the current workspace and corresponding history.
  #
  # Args:
  #   proj.name: Project name to save workspace and history as.
  
  save.image(file=paste(proj.name, '.RData', sep=''))
  savehistory(file=paste(proj.name,'.RHistory', sep=''))
  cat(paste('Workspace, history saved to', getwd(), '\n' ,sep=' ')) 
}

doPrediction <- function(fit.cv, x) {
  # Performs predictions given a count-matrix and a corresponding classifier.
  # Args:
  #   fit.cv: LASSO cross-validated (CV) classifier.
  #   x: Matrix of dimensions i * j.
  #
  # Returns:
  #   preds: vector of size i * 1 representing predictions per matrix row.
  
  preds <- predict(fit.cv, x, type='response', s=c('lambda.min')) # derive predictions
  return(preds)
}

getWeights <- function(fit.cv, scale.weights=T) {
  # Yields attribute weights (coefficients) given LASSO classifier.
  # Weights can be scaled based on whether they are positive or
  # negative. Such scaling is based off Taher, et. al., 2013.
  #
  # Args:
  #   fit.cv: LASSO cross-validated (CV) classifier.
  #   scale.weights: Z-normalize LASSO weights.
  #
  # Returns:
  #   weights: Dataframe of attributes and their respective weights.

  m <- as.matrix(coef(fit.cv, s=min(fit.cv$lambda)))
  if (scale.weights) { # standardize LASSO weights, if sought.
    weights.min <- min(m)
    weights.max <- max(m)
    for (i in 1: length(m)) {
      val <- m[i, 1] # retrieve the current value
      if (val < 0) { # normalize weights less than zero
        m[i, 1] <- -1 * ( 1 - ( ( val - weights.min ) / -weights.min ) )
      }
      else { # normalize weights greater-than or equal to zero
        m[i, 1] <- val / weights.max # maximum weight is set to 1.0
      }
    }
  }
  
  m <- as.matrix(m[order(row.names(m)), ]) # order weights by row-name
  m <- as.matrix(m[row.names(m) != '(Intercept)', ]) # intercept weight not needed
  colnames(m) <- c('Weight')
  return(m)
}

getAUC <- function(fit.cv) {
  # Identifies the Area Under Curve (AUC) of a LASSO classifier.
  #
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
