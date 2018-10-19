library(data.table)
library(caret) # classification
library(pROC) #ROC curve
library(AppliedPredictiveModeling)

target.prioritization.pipeline <- function() {
    ############################################################
    ### Prepare data
    ############################################################
    targets <- as.character(unlist(read.table(targets.path)))

    # Import node scores
    analysis <- read.csv(features.path, fileEncoding = "UTF-8", sep = "\t")
    # add labels (target or not)
    analysis$IsTarget <- analysis$GeneID %in% targets
    analysis[which(analysis$IsTarget == T), length(analysis)] <- "yes"
    analysis[which(analysis$IsTarget == F), length(analysis)] <- "no"
    analysis$IsTarget <- as.factor(analysis$IsTarget)
    # pre-format
    analysis <- na.omit(as.data.frame(analysis))
    rownames(analysis) <- analysis$GeneID
    return(analysis[ , 2:6])


    ############################################################
    ### Classify targets
    ############################################################
    # assess the performance of classification method
    auc.mat <- matrix(0, nrow = 5, ncol = 10)
    for (i in 1 : 10) {
        aucs <- cross.validate(analysis)
        auc.mat[, i] <- aucs
    }

    write.table(colMeans(auc.mat),
                file = auc.output.path,
                fileEncoding = "UTF-8",
                sep = "\t",
                row.names = FALSE
    )
}

cross.validate <- function(analysis, k.fold=5) {
  # splitting
  folds <- createFolds(analysis$IsTarget, k = k.fold)
  aucs = vector(mode="numeric", length = k.fold)
  for (i in 1:k.fold) {
    pre.train <- analysis[-folds[[i]], ]
    pre.test <- analysis[folds[[i]], ]
    result = classify.fold(pre.train, pre.test)
    aucs[i] <- as.numeric(result$auc)
  }
  return(aucs)
}

classify.fold <- function(training, test.set){
  training.feat <- training[, 1:4]
  test.set.feat <- test.set[, 1:4]

  # training
  train_control <- trainControl(method="boot", number=100)
  model <- train(x          = as.matrix(training.feat),
                 y          = training$IsTarget,
                 trControl  = train_control,
                 method     = "glm",
                 family     = "binomial")
  # calculate the likelihoods of the test set and evaluate the model
  test.set.likelihood <- predict(model, newdata=test.set.feat, type = "prob")
  roc.curve   <- roc(response = test.set$IsTarget,
                     predictor = test.set.likelihood$yes,
                     levels = rev(levels(test.set$IsTarget)))

  return(list(likelihoods = test.set.likelihood, auc = roc.curve$auc))
}

target.prioritization.pipeline()

