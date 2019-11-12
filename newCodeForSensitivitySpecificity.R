library(kebabs)
library(caret)
library(caretEnsemble)
library(plsgenomics)

dtype<-"imb" # dataset type; b for balanced dataset imb for imblanced.

if (dtype == "mouse")
{
  dataset<-readDNAStringSet("D:/Teaching/MS/Altaf/CodeResult&Graphs/CodeResult&Graphs/mouse_m6A.txt")
  p<-1934 # positive samples
  n<-1934 # Negative samples
}

if (dtype == "b")
{
  dataset<-readDNAStringSet("D:/Teaching/MS/Altaf/CodeResult&Graphs/CodeResult&Graphs/dnadatasetbalance.fa")
  p<-787 # positive samples
  n<-787 # Negative samples
}
if (dtype == "imb")
{
  dataset<-readDNAStringSet("D:/Teaching/MS/Altaf/CodeResult&Graphs/CodeResult&Graphs/dnadataset.fa")
  p<-787 # positive samples
  n<-1639 # negative samples 
  
}

#sink(paste("D:/Teaching/MS/Altaf/CodeResult&Graphs/CodeResult&Graphs/",'newResults_MouseNComponentsFifty2_99_','_.txt',sep=''))
for (km in 8:8)
{
  for (ncom in 72:72)
  {
    mm=23
    
    cat('k = ', km, ',', 'm = ',mm, 'ncomp = ',ncom, 'dType: ',dtype, '\n')
    gappy<- gappyPairKernel(k=km, m=mm)
    
    dnamethkernal <- getKernelMatrix(gappy, dataset)
    #dim(dnamethkernal)
    dnamethkernal<-as.matrix(dnamethkernal)
    #dim(dnamethkernal)
    
    #exrep <- getExRep(dataset, gappy, sparse=FALSE)
    #km<-as.matrix(exrep)
    
    label<-c(rep(1,p),rep(0,n))
    #print(label)
    
    #print(transformy(label))
    
    output.pls <- pls.regression(cbind(dnamethkernal),transformy(label), ncomp=ncom)                             
    #dnamethkernal<-as.matrix(dnamethkernal)
    data.learn <- scale(cbind(dnamethkernal), scale=TRUE, center=output.pls$meanX)%*%output.pls$R
    dim(data.learn)
    
    #label<-as.factor(label)
    #dnaFrame<-data.frame(label,data.learn)
    
    #algorithmList <- c('lda', 'rf', 'nnet', 'knn', 'svmRadial')
    #algorithmList <- c('lda', 'rf', 'nnet')
    
    algorithmList <- c('knn','rf')
    
    fitControl <- trainControl(method = "LOOCV",number=5,classProbs=TRUE, summaryFunction = twoClassSummary )
    #trainControl = trainControl(method = "LOOCV", classProbs=TRUE)
    label<-c(rep('p',p),rep('n',n))
    dnaFrame<-data.frame(label,data.learn)
    
    #models <- caretList(label~., data=dnaFrame,trControl=trainControl(method="cv",number=5,classProbs=TRUE), methodList=algorithmList)
    #models <- caretList(label~., data=dnaFrame,trControl=trainControl(method = "cv",number=5, classProbs=TRUE), methodList=algorithmList)
    
    rfFit3 <- train(label~., data = dnaFrame, 
                     method = "svmRadial", 
                     trControl = fitControl, 
                     verbose = FALSE, 
                     ## Specify which metric to optimize
                     metric = "ROC")
    print(rfFit3)
    
    #results <- resamples(models)
    #print(results)
    
    #print(summary(results))
    #dotplot(results)
    
    
    # correlation between results
    #modelCor(results)
    #splom(results)
    
    #pred <- predict(models$rf, newdata=dnaFrame)
    #print(confusionMatrix(data = pred, reference = dnaFrame$label))
    #evaluatePrediction(pred,dnaFrame$label, allLabels=unique(dnaFrame$label))
    print('\n##################################################################################\n')
    
    #remove(models)
    #remove(results)
    
  }
  
  
}

#sink()
