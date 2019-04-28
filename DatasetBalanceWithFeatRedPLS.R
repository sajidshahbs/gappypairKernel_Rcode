library(kebabs)
library(caret)
library(caretEnsemble)
library(plsgenomics)

dtype<-"b" # dataset type; b for balanced dataset imb for imblanced.

if (dtype == "b")
{
  dataset<-readDNAStringSet("D:/Teaching/MS/Altaf/CodeResult&Graphs/CodeResult&Graphs/dnadatasetbalance.fa")
  p<-787 # positive samples
  n<-787 # Negative samples
} else {
  dataset<-readDNAStringSet("D:/Teaching/MS/Altaf/CodeResult&Graphs/CodeResult&Graphs/dnadataset.fa")
  p<-787 # positive samples
  n<-1639 # negative samples 
  
  }



gappy<- gappyPairKernel(k=7, m=23)

dnamethkernal <- getKernelMatrix(gappy, dataset)
#dim(dnamethkernal)
dnamethkernal<-as.matrix(dnamethkernal)
#dim(dnamethkernal)

#exrep <- getExRep(dataset, gappy, sparse=FALSE)
#km<-as.matrix(exrep)

label<-c(rep(1,p),rep(0,n))
#print(label)

#print(transformy(label))

output.pls <- pls.regression(cbind(dnamethkernal),transformy(label), ncomp=110)                             
#dnamethkernal<-as.matrix(dnamethkernal)
data.learn <- scale(cbind(dnamethkernal), scale=TRUE, center=output.pls$meanX)%*%output.pls$R
dim(data.learn)

#label<-as.factor(label)
#dnaFrame<-data.frame(label,data.learn)

#algorithmList <- c('lda', 'rf', 'nnet', 'knn', 'svmRadial')
#algorithmList <- c('lda', 'rf', 'nnet')

algorithmList <- c('lda','nnet')

#trainControl = trainControl(method = "cv",number=5)
trainControl = trainControl(method = "LOOCV", classProbs=TRUE)
label<-c(rep('p',p),rep('n',n))
dnaFrame<-data.frame(label,data.learn)

#models <- caretList(label~., data=dnaFrame,trControl=trainControl(method="cv",number=5,classProbs=TRUE), methodList=algorithmList)
models <- caretList(label~., data=dnaFrame,trControl=trainControl(method = "cv",number=5, classProbs=TRUE), methodList=algorithmList)
models

results <- resamples(models)
summary(results)
#dotplot(results)


# correlation between results
#modelCor(results)
#splom(results)

seed<-999

# stack using glm
#stackControl <- trainControl(method="repeatedcv", number=5, repeats=3, savePredictions=TRUE, classProbs=TRUE)
stackControl <- trainControl(method="LOOCV", savePredictions=TRUE, classProbs=TRUE)
#set.seed(seed)
#stack.glm <- caretStack(models, method="glm", metric="Accuracy", trControl=stackControl)
#print(stack.glm)

# stack using svmRadial
#set.seed(seed)
#stack.svmRadial <- caretStack(models, method="svmRadial", metric="Accuracy", trControl=stackControl)
#print(stack.svmRadial)

# stack using random forest
#set.seed(seed)
#stack.rf <- caretStack(models, method="rf", metric="Accuracy", trControl=stackControl)
#print(stack.rf)

# stack using nnet
#set.seed(seed)
#stack.nnet <- caretStack(models, method="nnet", metric="Accuracy", trControl=stackControl)
#print(stack.nnet)

# stack using rpart
#set.seed(seed)
#stack.rpart <- caretStack(models, method="rpart", metric="Accuracy", trControl=stackControl)
#print(stack.rpart)

# stack using LDA
set.seed(seed)
stack.lda <- caretStack(models, method="lda", metric="Accuracy", trControl=stackControl)
print(stack.lda)

#Stack using knn 
#set.seed(seed)
#stack.knn <- caretStack(models, method="knn", metric="Accuracy", trControl=stackControl)
#print(stack.knn)


pred <- predict(models$lda, newdata=dnaFrame)
confusionMatrix(data = pred, reference = dnaFrame$label)

#evaluatePrediction(pred,dnaFrame$label, allLabels=unique(dnaFrame$label))