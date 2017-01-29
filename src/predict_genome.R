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
df <- read.csv("./data/Annotations/median.csv.gz", sep=",", head=TRUE, row.names=1)

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

mod <- xgboost(
  data=as.matrix(X),
  label=as.matrix(Y),
  objective="binary:logistic",
  nrounds=10,
  max_depth=2,
  eta=0.2,
  gamma=0,          
  colsample_bytree=1,
  min_child_weight=1,
  nthread=detectCores()
)


# now assign predictions to the genome ...

data <- as.matrix(predict(normalizer, df[fname]))
prob <- predict(mod, data)

pred <- cbind(rownames(df), prob)
colnames(pred) <- c("Gene", "Prob")

write.table(pred, file="../HIPred.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

