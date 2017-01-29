#!/usr/bin/Rscript

library(caret)
library(ROCR)

library(xgboost)

library(doMC)
registerDoMC(48)

args <- commandArgs(trailingOnly=TRUE)


# datasets
positive <- rownames(read.csv("./data/Dang.tsv", sep="\t", head=TRUE, row.names=1))
negative <- as.character(read.csv("./data/1kg_LoFT.tsv", head=FALSE)$V1)

print(paste("positive (raw): ", length(positive), sep=""))
print(paste("negative (raw): ", length(negative), sep=""))
print(" ")

# remove duplicate(s)
dup <- intersect(positive, negative)
positive <- setdiff(positive, dup)
negative <- setdiff(negative, dup)

print(paste("dup: ", length(dup), sep=""))
print(paste("positive (no dup): ", length(positive), sep=""))
print(paste("negative (no dup): ", length(negative), sep=""))
print(" ")

# fetch annotation data
df <- read.csv(gzfile(args[1]), sep=",", head=TRUE, row.names=1)

X <- df[intersect(rownames(df), unique(c(positive, negative))), ]
Y <- ifelse(rownames(X) %in% positive, 1, 0)

print(table(Y))

# remove zero-/near zero-Variance predictors
nzv <- nearZeroVar(X)
X <- X[, -nzv]

# store column names for later
fname <- colnames(X)

# normalize data
normalizer <- preProcess(X, method=c("center", "scale"))
X <- predict(normalizer, X)

# by default, caret does not report all the performance statistics 
# we need, so we need to pass a custom function to report these metrics ...
metrics <- function (data, lev=NULL, model=NULL) 
{
  # adapted from caret::twoClassSummary
  xtab <- as.vector(table(data[, 1], data[, 2])) # ... predicted, observed
  names(xtab) <- c("TN", "FP", "FN", "TP")
  
  Accuracy <- (xtab["TP"] + xtab["TN"]) / sum(xtab)
  Precision <- xtab["TP"] / (xtab['TP'] + xtab['FP'])   # aka. PPV
  Sensitivity <- xtab["TP"] / (xtab['TP'] + xtab['FN']) # aka. Recall
  Specificity <- xtab["TN"] / (xtab['TN'] + xtab['FP'])
  NPV <- xtab["TN"] / (xtab['TN'] + xtab['FN'])
  MCC <- ((xtab["TP"] * xtab["TN"]) - (xtab["FP"] + xtab["FN"])) / sqrt((xtab["TP"] + xtab["FP"]) * (xtab["TP"] + xtab["FN"]) * (xtab["TN"] + xtab["FP"]) * (xtab["TN"] + xtab["FN"]))
  
  # Note: use ROCR package for AUC
  pred <- prediction(data[, 4], data[, 2])       # ... predicted (prob), observed
  AUC  <- performance(pred, "auc")@y.values[[1]]
  
  perf <- c(Accuracy, Precision, Sensitivity, Specificity, NPV, MCC, AUC)
  names(perf) <- c("Accuracy", "Precision", "Sensitivity", "Specificity", "NPV", "MCC", "AUC")
  perf
}

# set up the grid search ...
grid <- expand.grid(
  nrounds=10,           # maximum number of iterations
  max_depth=c(1:10),    # maximum depth of the tree
  eta=c(0.1, 0.2, 0.3), # learning rate
  gamma=0,
  colsample_bytree=1,
  min_child_weight=1
)

# ... and pack the cross-validation parameters
ctrl <- trainControl(
  method="repeatedcv",
  number=10,               # k-fold
  repeats=30,              # repeated CV
  summaryFunction=metrics, # use our own defined metric function ...
  classProbs=TRUE,         # return probabilities so that we can derive AUCs
  allowParallel=TRUE,
  verboseIter=TRUE
)

# now, perform the grid search/cross validation
set.seed(42)

xgbcv <- train(
  x=X,
  y=make.names(Y), # need to use make.names ... why?!
  trControl=ctrl,
  tuneGrid=grid,
  method="xgbTree"
)

# dump incremental performance
write.csv(xgbcv$results, file=args[2], sep=",", quote=FALSE, row.names=FALSE)

# pick the best performing model
pick <- best(xgbcv$results, metric="AUC", maximize=TRUE)

mod <- xgboost(
  data=as.matrix(X), label=as.matrix(Y),
  objective="binary:logistic",
  nrounds=xgbcv$results[pick,]$nrounds,
  max_depth=xgbcv$results[pick,]$max_depth,
  eta=xgbcv$results[pick,]$eta,
  gamma=xgbcv$results[pick,]$gamma,          
  colsample_bytree=xgbcv$results[pick,]$colsample_bytree,
  min_child_weight=xgbcv$results[pick,]$min_child_weight,
  nthread=detectCores()
)

# perform feature selection
imp <- xgb.importance(model=mod)

sink(args[3])
for (idx in rev(c(1: nrow(imp)))) {
  p <- as.integer(imp[idx]$Feature)
  q <- as.numeric(imp[idx]$Gain)
  
  print(paste("Feature:", p, fname[p], "Gain:", q, sep=" "))
}
sink()

imp <- imp[1:10,]
imp$Feature <- fname[as.numeric(imp$Feature)]

