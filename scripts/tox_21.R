
# Jun 27 2020
# This script was written to run classification models on tox data...
#...for CAMDA 2020

# Temi

library(rstudioapi)
currPath <- getActiveDocumentContext()$path
setwd(dirname(currPath))

# load libraries
source('../scripts/packages.R')

# read data
tox.data <- read.csv('../data/p9-tox21-camda2020.csv',
                     header=T, stringsAsFactors=F, row.names = 'CAM_ID')

#dim(tox.data)
# target/label/training
target.data <- read.csv('../data/targets-camda2020.csv',
                        header=T, stringsAsFactors = F, check.names = F)

# tox.data has all the data
# select the training set based on what is available in target.data
training.set <- subset(target.data, Training_Validation=='Training Set')
tox.train <- tox.data[rownames(tox.data) %in% training.set$CAM_ID, ] 

complete.columns <- as.data.frame(apply(tox.train, 1, function(x){
    sum(is.na(x))
}))
names(complete.columns)[1] <- 'No.of.NAs'
complete.col.names <- rownames(subset(complete.columns, No.of.NAs==0))


#table(complete.columns$No.of.NAs)

# remove NAs
comp.tox.train <- tox.train[complete.cases(tox.train), ]
highly.corr <- findCorrelation(cor(comp.tox.train), cutoff=0.82)
comp.tox.train <- comp.tox.train[, -highly.corr]


#corrplot::corrplot(cor(comp.tox.train[, -ncol(comp.tox.train)]))

#View(cor(comp.tox.train))

# target data for dili 1 to 6
# target.labels.1 <- as.factor(comp.tox.1$Class)
target.labels.1 <- training.set[training.set$CAM_ID %in% row.names(comp.tox.train), ]$DILI1
target.labels.1 <- as.factor(mapvalues(target.labels.1, 
                                       from=c(0, 1), to=c('Positive', 'Negative')))
target.labels.3 <- training.set[training.set$CAM_ID %in% row.names(comp.tox.train), ]$DILI3
target.labels.3 <- as.factor(mapvalues(target.labels.3, 
                                       from=c(0, 1), to=c('Positive', 'Negative')))
target.labels.5 <- training.set[training.set$CAM_ID %in% row.names(comp.tox.train), ]$DILI5
target.labels.5 <- as.factor(mapvalues(target.labels.5, 
                                       from=c(0, 1), to=c('Positive', 'Negative')))
target.labels.6 <- training.set[training.set$CAM_ID %in% row.names(comp.tox.train), ]$DILI6
target.labels.6 <- as.factor(mapvalues(target.labels.6, 
                                       from=c(0, 1), to=c('Positive', 'Negative')))

# normalize the dataset to z-scores
#scaled.tox <- data.frame(scale(comp.tox.train, center = T, scale=T))

# bind with a target label
bound.tox.train.1 <- data.frame(cbind(comp.tox.train, Class=as.factor(target.labels.1)))
bound.tox.train.3 <- data.frame(cbind(comp.tox.train, Class=as.factor(target.labels.3)))
bound.tox.train.5 <- data.frame(cbind(comp.tox.train, Class=as.factor(target.labels.5)))
bound.tox.train.6 <- data.frame(cbind(comp.tox.train, Class=as.factor(target.labels.6)))

# smote dataset
set.seed(1)
smote1.tox <- bimba::SMOTE(data=bound.tox.train.1, perc_min = 50, k = 5)
set.seed(1)
smote3.tox <- bimba::SMOTE(data=bound.tox.train.3, perc_min = 50, k = 5)
set.seed(1)
smote5.tox <- bimba::SMOTE(data=bound.tox.train.5, perc_min = 50, k = 5)
set.seed(1)
smote6.tox <- caret::downSample(x=bound.tox.train.6[, -ncol(bound.tox.train.6)],
                                y=as.factor(bound.tox.train.6$Class))

# rose dataset
rose1.tox <- ROSE(Class ~., data=bound.tox.train.1, seed=1, N=max(table(bound.tox.train.1$Class))*2)$data
rose3.tox <- ROSE(Class ~., data=bound.tox.train.3, seed=1, N=max(table(bound.tox.train.3$Class))*2)$data
rose5.tox <- ROSE(Class ~., data=bound.tox.train.5, seed=1, N=max(table(bound.tox.train.5$Class))*2)$data
rose6.tox <- ROSE(Class ~., data=bound.tox.train.6, seed=1, N=max(table(bound.tox.train.6$Class))*2)$data


# remove these objects
rm('bound.tox.train.1','bound.tox.train.3','bound.tox.train.5','bound.tox.train.6')

# data mining ======
# define a general control
general.ctrl <- trainControl(method='repeatedcv',
                             repeats=100, number=5,
                             savePredictions=T,
                             classProbs = T,
                             summaryFunction = twoClassSummary,
                             allowParallel = T)

# cluster
#cl <- makeCluster(detectCores() - 56)
registerDoParallel(clusters)

print('Started training tox models...')

# ============= DILI1 ==================================
# glm ======
set.seed(1)
glm1.scaled.model <- train(x=comp.tox.train, y=target.labels.1,
                           maxit=200, method='glm',
                           family=binomial(link='logit'),
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
glm1.rose.model <- train(Class ~., data=rose1.tox,
                         maxit=200, method='glm',
                         family=binomial(link='logit'),
                         trControl = general.ctrl, metric='ROC')

set.seed(1)
glm1.smote.model <- train(Class ~., data=smote1.tox,
                          maxit=200, method='glm',
                          family=binomial(link='logit'),
                          trControl = general.ctrl, metric='ROC')

glm1.res <- resamples(list(glm1=glm1.scaled.model,
                           glm1.rose=glm1.rose.model,
                           glm1.smote=glm1.smote.model))

summary(glm1.res)

# qda ==========
set.seed(1)
qda1.scaled.model <- train(x=comp.tox.train, y=target.labels.1, 
                           method='qda', trControl = general.ctrl, metric='ROC')

set.seed(1)
qda1.rose.model <- train(Class ~., data=rose1.tox, method='qda',
                         trControl = general.ctrl, metric='ROC')

set.seed(1)
qda1.smote.model <- train(Class ~., data=smote1.tox, method='qda',
                          trControl = general.ctrl, metric='ROC')

qda1.res <- resamples(list(qda1=qda1.scaled.model,
                           qda1.rose=qda1.rose.model,
                           qda1.smote=qda1.smote.model))

summary(qda1.res)

# random forest ======
set.seed(1)
rf1.scaled.model <- train(x=comp.tox.train, y=target.labels.1, 
                          method='rf', trControl = general.ctrl, metric='ROC')

set.seed(1)
rf1.rose.model <- train(Class ~., data=rose1.tox, method='rf',
                        trControl = general.ctrl, metric='ROC')

set.seed(1)
rf1.smote.model <- train(Class ~., data=smote1.tox, method='rf',
                         trControl = general.ctrl, metric='ROC')

rf1.res <- resamples(list(rf1=rf1.scaled.model,
                          rf1.rose=rf1.rose.model,
                          rf1.smote=rf1.smote.model))
summary(rf1.res)


# rpart=========
set.seed(1)
rpart1.scaled.model <- train(x=comp.tox.train, y=target.labels.1, 
                             method='rpart', trControl = general.ctrl, metric='ROC')

set.seed(1)
rpart1.rose.model <- train(Class ~., data=rose1.tox, method='rpart',
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
rpart1.smote.model <- train(Class ~., data=smote1.tox, method='rpart',
                            trControl = general.ctrl, metric='ROC')

rpart1.res <- resamples(list(rpart1=rpart1.scaled.model,
                             rpart1.rose=rpart1.rose.model,
                             rpart1.smote=rpart1.smote.model))
summary(rpart1.res)

# svmRadial==========
set.seed(1)
tr_data <- as.data.frame(cbind(Class=factor(target.labels.1, levels = c('Positive', 'Negative')), comp.tox.train))
svmRadial1.scaled.model <- train(Class~., data=tr_data,
                                 method='svmRadial', trControl = general.ctrl, metric='ROC',
                                 tuneLength=50)

# svmRadial1.scaled.model <- train(x=comp.tox.train, y=target.labels.1,
#                                  method='svmRadial', trControl = general.ctrl, metric='ROC',
#                                  tuneLength=50)

set.seed(1)
svmRadial1.rose.model <- train(Class ~., data=rose1.tox,
                               method='svmRadial', trControl = general.ctrl, metric='ROC',
                               tuneLength=50)

set.seed(1)
svmRadial1.smote.model <- train(Class ~., data=smote1.tox,
                                method='svmRadial', trControl = general.ctrl, metric='ROC',
                                tuneLength=50)

svmRadial1.res <- resamples(list(svmRadial1=svmRadial1.scaled.model,
                                 svmRadial1.rose=svmRadial1.rose.model,
                                 svmRadial1.smote=svmRadial1.smote.model))
summary(svmRadial1.res)


# svmPoly==========
set.seed(1)
tr_data <- as.data.frame(cbind(Class=factor(target.labels.1, levels = c('Positive', 'Negative')), comp.tox.train))
svmPoly1.scaled.model <- train(Class~., data=tr_data,
                               method='svmPoly', trControl = general.ctrl, metric='ROC')

set.seed(1)
svmPoly1.rose.model <- train(Class ~., data=rose1.tox,
                             method='svmPoly', trControl = general.ctrl, metric='ROC')

set.seed(1)
svmPoly1.smote.model <- train(Class ~., data=smote1.tox,
                              method='svmPoly', trControl = general.ctrl, metric='ROC')

svmPoly1.res <- resamples(list(svmPoly1=svmPoly1.scaled.model,
                               svmPoly1.rose=svmPoly1.rose.model,
                               svmPoly1.smote=svmPoly1.smote.model))
summary(svmPoly1.res)

# nnet==========
set.seed(1)
nnet1.scaled.model <- train(x=comp.tox.train, y=target.labels.1,
                            method='nnet', trControl = general.ctrl, metric='ROC',
                            maxit=200)

set.seed(1)
nnet1.rose.model <- train(Class ~., data=rose1.tox,
                          method='nnet', trControl = general.ctrl, metric='ROC',
                          maxit=200)

set.seed(1)
nnet1.smote.model <- train(Class ~., data=smote1.tox,
                           method='nnet', trControl = general.ctrl, metric='ROC',
                           maxit=200)

nnet1.res <- resamples(list(nnet1=nnet1.scaled.model,
                            nnet1.rose=nnet1.rose.model,
                            nnet1.smote=nnet1.smote.model))
summary(nnet1.res)

# naive bayes==========
set.seed(1)
nb1.grid <- expand.grid(laplace=c(0, 1),
                        usekernel=c('TRUE', 'FALSE'),
                        adjust=c(1:5))

nb1.scaled.model <- train(x=comp.tox.train, y=target.labels.1,
                          method='naive_bayes', trControl = general.ctrl,
                          metric='ROC', tuneGrid = nb1.grid)

set.seed(1)
nb1.rose.model <- train(Class ~., data=rose1.tox,
                        method='naive_bayes', trControl = general.ctrl,
                        metric='ROC', tuneGrid = nb1.grid)

set.seed(1)
nb1.smote.model <- train(Class ~., data=smote1.tox,
                         method='naive_bayes', trControl = general.ctrl, 
                         tuneGrid = nb1.grid, metric='ROC')

nb1.res <- resamples(list(nb1=nb1.scaled.model,
                          nb1.rose=nb1.rose.model,
                          nb1.smote=nb1.smote.model))
summary(nb1.res)


# =============== DILI3 ========================
# glm ======
set.seed(1)
glm3.scaled.model <- train(x=comp.tox.train, y=target.labels.3,
                           maxit=200, method='glm',
                           family=binomial(link='logit'),
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
glm3.rose.model <- train(Class ~., data=rose3.tox,
                         maxit=200, method='glm',
                         family=binomial(link='logit'),
                         trControl = general.ctrl, metric='ROC')

set.seed(1)
glm3.smote.model <- train(Class ~., data=smote3.tox,
                          maxit=200, method='glm',
                          family=binomial(link='logit'),
                          trControl = general.ctrl, metric='ROC')

glm3.res <- resamples(list(glm3=glm3.scaled.model,
                           glm3.rose=glm3.rose.model,
                           glm3.smote=glm3.smote.model))

summary(glm3.res)

# qda===========
set.seed(1)
qda3.scaled.model <- train(x=comp.tox.train, y=target.labels.3, method='qda',
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
qda3.rose.model <- train(Class ~., data=rose3.tox, method='qda',
                         trControl = general.ctrl, metric='ROC')

set.seed(1)
qda3.smote.model <- train(Class ~., data=smote3.tox, method='qda',
                          trControl = general.ctrl, metric='ROC')

qda3.res <- resamples(list(qda3=qda3.scaled.model,
                           qda3.rose=qda3.rose.model,
                           qda3.smote=qda3.smote.model))

summary(qda3.res)

# random forest ======
set.seed(1)
rf3.scaled.model <- train(x=comp.tox.train, y=target.labels.3, method='rf',
                          trControl = general.ctrl, metric='ROC')

set.seed(1)
rf3.rose.model <- train(Class ~., data=rose3.tox, method='rf',
                        trControl = general.ctrl, metric='ROC')

set.seed(1)
rf3.smote.model <- train(Class ~., data=smote3.tox, method='rf',
                         trControl = general.ctrl, metric='ROC')

rf3.res <- resamples(list(rf3=rf3.scaled.model,
                          rf3.rose=rf3.rose.model,
                          rf3.smote=rf3.smote.model))
summary(rf3.res)

# rpart=========
set.seed(1)
rpart3.scaled.model <- train(x=comp.tox.train, y=target.labels.3, method='rpart',
                             trControl = general.ctrl, metric='ROC')

set.seed(1)
rpart3.rose.model <- train(Class ~., data=rose3.tox, method='rpart',
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
rpart3.smote.model <- train(Class ~., data=smote3.tox, method='rpart',
                            trControl = general.ctrl, metric='ROC')

rpart3.res <- resamples(list(rpart3=rpart3.scaled.model,
                             rpart3.rose=rpart3.rose.model,
                             rpart3.smote=rpart3.smote.model))
summary(rpart3.res)

# svmRadial==========
set.seed(1)
tr_data <- as.data.frame(cbind(Class=factor(target.labels.3, levels = c('Positive', 'Negative')), comp.tox.train))
svmRadial3.scaled.model <- train(Class~., data=tr_data,
                                 method='svmRadial', trControl = general.ctrl, metric='ROC',
                                 tuneLength=50)

set.seed(1)
svmRadial3.rose.model <- train(Class ~., data=rose3.tox,
                               method='svmRadial', trControl = general.ctrl, metric='ROC',
                               tuneLength=50)

set.seed(1)
svmRadial3.smote.model <- train(Class ~., data=smote3.tox,
                                method='svmRadial', trControl = general.ctrl, metric='ROC',
                                tuneLength=50)

svmRadial3.res <- resamples(list(svmRadial3=svmRadial3.scaled.model,
                                 svmRadial3.rose=svmRadial3.rose.model,
                                 svmRadial3.smote=svmRadial3.smote.model))
summary(svmRadial3.res)

# svmPoly==========
set.seed(1)
tr_data <- as.data.frame(cbind(Class=factor(target.labels.3, levels = c('Positive', 'Negative')), comp.tox.train))
svmPoly3.scaled.model <- train(Class~., data = tr_data,
                               method='svmPoly', trControl = general.ctrl, metric='ROC')

set.seed(1)
svmPoly3.rose.model <- train(Class ~., data=rose3.tox,
                             method='svmPoly', trControl = general.ctrl, metric='ROC')

set.seed(1)
svmPoly3.smote.model <- train(Class ~., data=smote3.tox,
                              method='svmPoly', trControl = general.ctrl, metric='ROC')

svmPoly3.res <- resamples(list(svmPoly3=svmPoly3.scaled.model,
                               svmPoly3.rose=svmPoly3.rose.model,
                               svmPoly3.smote=svmPoly3.smote.model))
summary(svmPoly3.res)

# nnet==========
set.seed(1)
nnet3.scaled.model <- train(x=comp.tox.train, y=target.labels.3,
                            method='nnet', trControl = general.ctrl, metric='ROC',
                            maxit=200)

set.seed(1)
nnet3.rose.model <- train(Class ~., data=rose3.tox,
                          method='nnet', trControl = general.ctrl, metric='ROC',
                          maxit=200)

set.seed(1)
nnet3.smote.model <- train(Class ~., data=smote3.tox,
                           method='nnet', trControl = general.ctrl, metric='ROC',
                           maxit=200)

nnet3.res <- resamples(list(nnet3=nnet3.scaled.model,
                            nnet3.rose=nnet3.rose.model,
                            nnet3.smote=nnet3.smote.model))
summary(nnet3.res)

# naive bayes==========
set.seed(1)
nb3.grid <- expand.grid(laplace=c(0, 1),
                        usekernel=c('TRUE', 'FALSE'),
                        adjust=c(1:5))

nb3.scaled.model <- train(x=comp.tox.train, y=target.labels.3,
                          method='naive_bayes', trControl = general.ctrl,
                          metric='ROC', tuneGrid = nb3.grid)

set.seed(1)
nb3.rose.model <- train(Class ~., data=rose3.tox,
                        method='naive_bayes', trControl = general.ctrl, 
                        tuneGrid = nb3.grid, metric='ROC')

set.seed(1)
nb3.smote.model <- train(Class ~., data=smote3.tox,
                         method='naive_bayes', trControl = general.ctrl, 
                         tuneGrid = nb3.grid, metric='ROC')

nb3.res <- resamples(list(nb3=nb3.scaled.model,
                          nb3.rose=nb3.rose.model,
                          nb3.smote=nb3.smote.model))
summary(nb3.res)


# ========== DILI5 =======================
# glm ======
set.seed(1)
glm5.scaled.model <- train(x=comp.tox.train, y=target.labels.5,
                           maxit=200, method='glm',
                           family=binomial(link='logit'),
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
glm5.rose.model <- train(Class ~., data=rose5.tox,
                         maxit=200, method='glm',
                         family=binomial(link='logit'),
                         trControl = general.ctrl, metric='ROC')

set.seed(1)
glm5.smote.model <- train(Class ~., data=smote5.tox,
                          maxit=200, method='glm',
                          family=binomial(link='logit'),
                          trControl = general.ctrl, metric='ROC')

glm5.res <- resamples(list(glm5=glm5.scaled.model,
                           glm5.rose=glm5.rose.model,
                           glm5.smote=glm5.smote.model))

summary(glm5.res)

# qda===========
set.seed(1)
qda5.scaled.model <- train(x=comp.tox.train, y=target.labels.5, method='qda',
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
qda5.rose.model <- train(Class ~., data=rose5.tox, method='qda',
                         trControl = general.ctrl, metric='ROC')

set.seed(1)
qda5.smote.model <- train(Class ~., data=smote5.tox, method='qda',
                          trControl = general.ctrl, metric='ROC')

qda5.res <- resamples(list(qda5=qda5.scaled.model,
                           qda5.rose=qda5.rose.model,
                           qda5.smote=qda5.smote.model))

summary(qda5.res)

# random forest ======
set.seed(1)
rf5.scaled.model <- train(x=comp.tox.train, y=target.labels.5, method='rf',
                          trControl = general.ctrl, metric='ROC')

set.seed(1)
rf5.rose.model <- train(Class ~., data=rose5.tox, method='rf',
                        trControl = general.ctrl, metric='ROC')

set.seed(1)
rf5.smote.model <- train(Class ~., data=smote5.tox, method='rf',
                         trControl = general.ctrl, metric='ROC')

rf5.res <- resamples(list(rf5=rf5.scaled.model,
                          rf5.rose=rf5.rose.model,
                          rf5.smote=rf5.smote.model))
summary(rf5.res)

# rpart=========
set.seed(1)
rpart5.scaled.model <- train(x=comp.tox.train, y=target.labels.5, method='rpart',
                             trControl = general.ctrl, metric='ROC')

set.seed(1)
rpart5.rose.model <- train(Class ~., data=rose5.tox, method='rpart',
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
rpart5.smote.model <- train(Class ~., data=smote5.tox, method='rpart',
                            trControl = general.ctrl, metric='ROC')

rpart5.res <- resamples(list(rpart5=rpart5.scaled.model,
                             rpart5.rose=rpart5.rose.model,
                             rpart5.smote=rpart5.smote.model))
summary(rpart5.res)

# svmRadial==========
set.seed(1)
tr_data <- as.data.frame(cbind(Class=factor(target.labels.5, levels = c('Positive', 'Negative')), comp.tox.train))
svmRadial5.scaled.model <- train(Class ~ ., data=tr_data,
                                 method='svmRadial', trControl = general.ctrl, metric='ROC',
                                 tuneLength=50)
# 
# svmRadial5.scaled.model <- train(x=comp.tox.train, y=target.labels.5,
#                                  method='svmRadial', trControl = general.ctrl, metric='ROC',
#                                  tuneLength=50)

set.seed(1)
svmRadial5.rose.model <- train(Class ~., data=rose5.tox,
                               method='svmRadial', trControl = general.ctrl, metric='ROC',
                               tuneLength=50)

set.seed(1)
svmRadial5.smote.model <- train(Class ~., data=smote5.tox,
                                method='svmRadial', trControl = general.ctrl, metric='ROC',
                                tuneLength=50)

svmRadial5.res <- resamples(list(svmRadial5=svmRadial5.scaled.model,
                                 svmRadial5.rose=svmRadial5.rose.model,
                                 svmRadial5.smote=svmRadial5.smote.model))
summary(svmRadial5.res)

# svmPoly==========
set.seed(1)
tr_data <- as.data.frame(cbind(Class=factor(target.labels.5, levels = c('Positive', 'Negative')), comp.tox.train))
svmPoly5.scaled.model <- train(Class~., data=tr_data,
                               method='svmPoly', trControl = general.ctrl, metric='ROC')

set.seed(1)
svmPoly5.rose.model <- train(Class ~., data=rose5.tox,
                             method='svmPoly', trControl = general.ctrl, metric='ROC')

set.seed(1)
svmPoly5.smote.model <- train(Class ~., data=smote5.tox,
                              method='svmPoly', trControl = general.ctrl, metric='ROC')

svmPoly5.res <- resamples(list(svmPoly5=svmPoly5.scaled.model,
                               svmPoly5.rose=svmPoly5.rose.model,
                               svmPoly5.smote=svmPoly5.smote.model))
summary(svmPoly5.res)

# nnet==========
set.seed(1)
nnet5.scaled.model <- train(x=comp.tox.train, y=target.labels.5,
                            method='nnet', trControl = general.ctrl, metric='ROC',
                            maxit=200)

set.seed(1)
nnet5.rose.model <- train(Class ~., data=rose5.tox,
                          method='nnet', trControl = general.ctrl, metric='ROC',
                          maxit=200)

set.seed(1)
nnet5.smote.model <- train(Class ~., data=smote5.tox,
                           method='nnet', trControl = general.ctrl, metric='ROC',
                           maxit=200)

nnet5.res <- resamples(list(nnet5=nnet5.scaled.model,
                            nnet5.rose=nnet5.rose.model,
                            nnet5.smote=nnet5.smote.model))
summary(nnet5.res)

# naive bayes==========
set.seed(1)
nb5.grid <- expand.grid(laplace=c(0, 1),
                        usekernel=c('TRUE', 'FALSE'),
                        adjust=c(1:5))

nb5.scaled.model <- train(x=comp.tox.train, y=target.labels.5,
                          method='naive_bayes', trControl = general.ctrl,
                          metric='ROC', tuneGrid = nb5.grid)

set.seed(1)
nb5.rose.model <- train(Class ~., data=rose5.tox,
                        method='naive_bayes', trControl = general.ctrl, 
                        tuneGrid = nb5.grid, metric='ROC')

set.seed(1)
nb5.smote.model <- train(Class ~., data=smote5.tox,
                         method='naive_bayes', trControl = general.ctrl, 
                         tuneGrid = nb5.grid, metric='ROC')

nb5.res <- resamples(list(nb5=nb5.scaled.model,
                          nb5.rose=nb5.rose.model,
                          nb5.smote=nb5.smote.model))
summary(nb5.res)


# ================ DILI6 ====================
# glm ======
set.seed(1)
glm6.scaled.model <- train(x=comp.tox.train, y=target.labels.6,
                           maxit=200, method='glm',
                           family=binomial(link='logit'),
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
glm6.rose.model <- train(Class ~., data=rose6.tox,
                         maxit=200, method='glm',
                         family=binomial(link='logit'),
                         trControl = general.ctrl, metric='ROC')

set.seed(1)
glm6.smote.model <- train(Class ~., data=smote6.tox,
                          maxit=200, method='glm',
                          family=binomial(link='logit'),
                          trControl = general.ctrl, metric='ROC')

glm6.res <- resamples(list(glm6=glm6.scaled.model,
                           glm6.rose=glm6.rose.model,
                           glm6.smote=glm6.smote.model))

summary(glm6.res)

# qda===========
set.seed(1)
qda6.scaled.model <- train(x=comp.tox.train, y=target.labels.6, method='qda',
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
qda6.rose.model <- train(Class ~., data=rose6.tox, method='qda',
                         trControl = general.ctrl, metric='ROC')

set.seed(1)
qda6.smote.model <- train(Class ~., data=smote6.tox, method='qda',
                          trControl = general.ctrl, metric='ROC')

qda6.res <- resamples(list(qda6=qda6.scaled.model,
                           qda6.rose=qda6.rose.model,
                           qda6.smote=qda6.smote.model))

summary(qda6.res)

# random forest ======
set.seed(1)
rf6.scaled.model <- train(x=comp.tox.train, y=target.labels.6, method='rf',
                          trControl = general.ctrl, metric='ROC')

set.seed(1)
rf6.rose.model <- train(Class ~., data=rose6.tox, method='rf',
                        trControl = general.ctrl, metric='ROC')

set.seed(1)
rf6.smote.model <- train(Class ~., data=smote6.tox, method='rf',
                         trControl = general.ctrl, metric='ROC')

rf6.res <- resamples(list(rf6=rf6.scaled.model,
                          rf6.rose=rf6.rose.model,
                          rf6.smote=rf6.smote.model))
summary(rf6.res)

# rpart=========
set.seed(1)
rpart6.scaled.model <- train(x=comp.tox.train, y=target.labels.6, method='rpart',
                             trControl = general.ctrl, metric='ROC')

set.seed(1)
rpart6.rose.model <- train(Class ~., data=rose6.tox, method='rpart',
                           trControl = general.ctrl, metric='ROC')

set.seed(1)
rpart6.smote.model <- train(Class ~., data=smote6.tox, method='rpart',
                            trControl = general.ctrl, metric='ROC')

rpart6.res <- resamples(list(rpart6=rpart6.scaled.model,
                             rpart6.rose=rpart6.rose.model,
                             rpart6.smote=rpart6.smote.model))
summary(rpart6.res)

# svmRadial==========
set.seed(1)
tr_data <- as.data.frame(cbind(Class=factor(target.labels.6, levels = c('Positive', 'Negative')), comp.tox.train))
svmRadial6.scaled.model <- train(Class~., data=tr_data,
                                 method='svmRadial', trControl = general.ctrl, metric='ROC',
                                 tuneLength=50)

set.seed(1)
svmRadial6.rose.model <- train(Class ~., data=rose6.tox,
                               method='svmRadial', trControl = general.ctrl, metric='ROC',
                               tuneLength=50)

set.seed(1)
svmRadial6.smote.model <- train(Class ~., data=smote6.tox,
                                method='svmRadial', trControl = general.ctrl, metric='ROC',
                                tuneLength=50)

svmRadial6.res <- resamples(list(svmRadial6=svmRadial6.scaled.model,
                                 svmRadial6.rose=svmRadial6.rose.model,
                                 svmRadial6.smote=svmRadial6.smote.model))
summary(svmRadial6.res)

# svmPoly==========
set.seed(1)
tr_data <- as.data.frame(cbind(Class=factor(target.labels.6, levels = c('Positive', 'Negative')), comp.tox.train))
svmPoly6.scaled.model <- train(Class~., data=tr_data,
                               method='svmPoly', trControl = general.ctrl, metric='ROC')

set.seed(1)
svmPoly6.rose.model <- train(Class ~., data=rose6.tox,
                             method='svmPoly', trControl = general.ctrl, metric='ROC')

set.seed(1)
svmPoly6.smote.model <- train(Class ~., data=smote6.tox,
                              method='svmPoly', trControl = general.ctrl, metric='ROC')

svmPoly6.res <- resamples(list(svmPoly6=svmPoly6.scaled.model,
                               svmPoly6.rose=svmPoly6.rose.model,
                               svmPoly6.smote=svmPoly6.smote.model))
summary(svmPoly6.res)

# nnet==========
set.seed(1)
nnet6.scaled.model <- train(x=comp.tox.train, y=target.labels.6,
                            method='nnet', trControl = general.ctrl, metric='ROC',
                            maxit=200)

set.seed(1)
nnet6.rose.model <- train(Class ~., data=rose6.tox,
                          method='nnet', trControl = general.ctrl, metric='ROC',
                          maxit=200)

set.seed(1)
nnet6.smote.model <- train(Class ~., data=smote6.tox,
                           method='nnet', trControl = general.ctrl, metric='ROC',
                           maxit=200)

nnet6.res <- resamples(list(nnet6=nnet6.scaled.model,
                            nnet6.rose=nnet6.rose.model,
                            nnet6.smote=nnet6.smote.model))
summary(nnet6.res)

# naive bayes==========
set.seed(1)
nb6.grid <- expand.grid(laplace=c(0, 1),
                        usekernel=c('TRUE', 'FALSE'),
                        adjust=c(1:5))

nb6.scaled.model <- train(x=comp.tox.train, y=target.labels.6,
                          method='naive_bayes', trControl = general.ctrl,
                          metric='ROC', tuneGrid = nb6.grid)

set.seed(1)
nb6.rose.model <- train(Class ~., data=rose6.tox,
                        method='naive_bayes', trControl = general.ctrl, 
                        tuneGrid = nb6.grid, metric='ROC')

set.seed(1)
nb6.smote.model <- train(Class ~., data=smote6.tox,
                         method='naive_bayes', trControl = general.ctrl, 
                         tuneGrid = nb6.grid, metric='ROC')

nb6.res <- resamples(list(nb6=nb6.scaled.model,
                          nb6.rose=nb6.rose.model,
                          nb6.smote=nb6.smote.model))
summary(nb6.res)

print('Finished training tox models...')

stopCluster(clusters)
# 



# =========== Evaluate ===================

summary.list.1 <- list(summary(glm1.res),
                       summary(qda1.res),
                       summary(nb1.res),
                       summary(nnet1.res),
                       summary(rpart1.res),
                       summary(rf1.res),
                       summary(svmPoly1.res),
                       summary(svmRadial1.res))

names(summary.list.1) <- paste0(c('glm', 'qda', 'nb', 'nnet', 'rpart', 'rf', 'svmPoly', 'svmRadial'), 1)

model.list.1 <- list(glm1=glm1.scaled.model, glm1.rose=glm1.rose.model, glm1.smote=glm1.smote.model, 
                     rf1=rf1.scaled.model, rf1.rose=rf1.rose.model, rf1.smote=rf1.smote.model, 
                     svmRadial1=svmRadial1.scaled.model, svmRadial1.rose=svmRadial1.rose.model, svmRadial1.smote=glm1.smote.model, 
                     svmPoly1=svmPoly1.scaled.model, svmPoly1.rose=svmPoly1.rose.model, svmPoly1.smote=svmPoly1.smote.model, 
                     qda1=qda1.scaled.model, qda1.rose=qda1.rose.model, qda1.smote=qda1.smote.model, 
                     rpart1=rpart1.scaled.model, rpart1.rose=rpart1.rose.model, rpart1.smote=glm1.smote.model, 
                     nnet1=nnet1.scaled.model, nnet1.rose=nnet1.rose.model, nnet1.smote=nnet1.smote.model, 
                     nb1=nb1.scaled.model, nb1.rose=nb1.rose.model, nb1.smote=nb1.smote.model)

summary.list.3 <- list(summary(glm3.res),
                       summary(qda3.res),
                       summary(nb3.res),
                       summary(nnet3.res),
                       summary(rpart3.res),
                       summary(rf3.res),
                       summary(svmPoly3.res),
                       summary(svmRadial3.res))

names(summary.list.3) <- paste0(c('glm', 'qda', 'nb', 'nnet', 'rpart', 'rf', 'svmPoly', 'svmRadial'), 3)

model.list.3 <- list(glm3=glm3.scaled.model, glm3.rose=glm3.rose.model, glm3.smote=glm3.smote.model, 
                     rf3=rf3.scaled.model, rf3.rose=rf3.rose.model, rf3.smote=rf3.smote.model, 
                     svmRadial3=svmRadial3.scaled.model, svmRadial3.rose=svmRadial3.rose.model, svmRadial3.smote=glm3.smote.model, 
                     svmPoly3=svmPoly3.scaled.model, svmPoly3.rose=svmPoly3.rose.model, svmPoly3.smote=svmPoly3.smote.model, 
                     qda3=qda3.scaled.model, qda3.rose=qda3.rose.model, qda3.smote=qda3.smote.model, 
                     rpart3=rpart3.scaled.model, rpart3.rose=rpart3.rose.model, rpart3.smote=glm3.smote.model, 
                     nnet3=nnet3.scaled.model, nnet3.rose=nnet3.rose.model, nnet3.smote=nnet3.smote.model, 
                     nb3=nb3.scaled.model, nb3.rose=nb3.rose.model, nb3.smote=nb3.smote.model)

summary.list.5 <- list(summary(glm5.res),
                       summary(qda5.res),
                       summary(nb5.res),
                       summary(nnet5.res),
                       summary(rpart5.res),
                       summary(rf5.res),
                       summary(svmPoly5.res),
                       summary(svmRadial5.res))

names(summary.list.5) <- paste0(c('glm', 'qda', 'nb', 'nnet', 'rpart', 'rf', 'svmPoly', 'svmRadial'), 5)

model.list.5 <- list(glm5=glm5.scaled.model, glm5.rose=glm5.rose.model, glm5.smote=glm5.smote.model, 
                     rf5=rf5.scaled.model, rf5.rose=rf5.rose.model, rf5.smote=rf5.smote.model, 
                     svmRadial5=svmRadial5.scaled.model, svmRadial5.rose=svmRadial5.rose.model, svmRadial5.smote=glm5.smote.model, 
                     svmPoly5=svmPoly5.scaled.model, svmPoly5.rose=svmPoly5.rose.model, svmPoly5.smote=svmPoly5.smote.model, 
                     qda5=qda5.scaled.model, qda5.rose=qda5.rose.model, qda5.smote=qda5.smote.model, 
                     rpart5=rpart5.scaled.model, rpart5.rose=rpart5.rose.model, rpart5.smote=glm5.smote.model, 
                     nnet5=nnet5.scaled.model, nnet5.rose=nnet5.rose.model, nnet5.smote=nnet5.smote.model, 
                     nb5=nb5.scaled.model, nb5.rose=nb5.rose.model, nb5.smote=nb5.smote.model)

summary.list.6 <- list(summary(glm6.res),
                       summary(qda6.res),
                       summary(nb6.res),
                       summary(nnet6.res),
                       summary(rpart6.res),
                       summary(rf6.res),
                       summary(svmPoly6.res),
                       summary(svmRadial6.res))

names(summary.list.6) <- paste0(c('glm', 'qda', 'nb', 'nnet', 'rpart', 'rf', 'svmPoly', 'svmRadial'), 6)

model.list.6 <- list(glm6=glm6.scaled.model, glm6.rose=glm6.rose.model, glm6.smote=glm6.smote.model, 
                     rf6=rf6.scaled.model, rf6.rose=rf6.rose.model, rf6.smote=rf6.smote.model, 
                     svmRadial6=svmRadial6.scaled.model, svmRadial6.rose=svmRadial6.rose.model, svmRadial6.smote=glm6.smote.model, 
                     svmPoly6=svmPoly6.scaled.model, svmPoly6.rose=svmPoly6.rose.model, svmPoly6.smote=svmPoly6.smote.model, 
                     qda6=qda6.scaled.model, qda6.rose=qda6.rose.model, qda6.smote=qda6.smote.model, 
                     rpart6=rpart6.scaled.model, rpart6.rose=rpart6.rose.model, rpart6.smote=glm6.smote.model, 
                     nnet6=nnet6.scaled.model, nnet6.rose=nnet6.rose.model, nnet6.smote=nnet6.smote.model, 
                     nb6=nb6.scaled.model, nb6.rose=nb6.rose.model, nb6.smote=nb6.smote.model)

# save the models into a list
tox_models <- list(tox_dili1=model.list.1, tox_dili3=model.list.3, tox_dili5=model.list.5, tox_dili6=model.list.6)
save(tox_models, file = '../models/tox_models.RData')

rm(list=ls())

# tox_models$tox_dili1$svmRadial1 <- svmRadial1.scaled.model
# tox_models$tox_dili1$svmPoly1 <- svmPoly1.scaled.model
# 
# tox_models$tox_dili3$svmRadial3 <- svmRadial3.scaled.model
# tox_models$tox_dili3$svmPoly3 <- svmPoly3.scaled.model
# 
# tox_models$tox_dili5$svmRadial5 <- svmRadial5.scaled.model
# tox_models$tox_dili5$svmPoly5 <- svmPoly5.scaled.model
# 
# tox_models$tox_dili6$svmRadial6 <- svmRadial6.scaled.model
# tox_models$tox_dili6$svmPoly6 <- svmPoly6.scaled.model

#save.image(file='/extData/NGS/hurlab/temi/projects/camda/tox21/latest_tox21_data.RData')

# using the top 3 models in each category ==================
# b.model.list.1 <- list(rf1.smote=rf1.smote.model, 
#                        svmRadial1.smote=svmRadial1.smote.model,
#                        nnet1.smote=nnet1.smote.model,
#                        svmRadial1=svmRadial1.scaled.model,
#                        nb1=nb1.scaled.model,
#                        svmPoly1=svmPoly1.scaled.model)
# b.model.list.3 <- list(rf3.smote=rf3.smote.model, 
#                        svmRadial3.rose=svmRadial3.rose.model,
#                        rf3.rose=rf3.rose.model,
#                        svmPoly3=svmPoly3.scaled.model,
#                        nnet3=nnet3.scaled.model,
#                        rf3=rf3.scaled.model)
# b.model.list.5 <- list(rf5.rose=rf5.rose.model, 
#                        nb5.rose=nb5.rose.model,
#                        qda5.rose=qda5.rose.model,
#                        glm5=glm5.scaled.model,
#                        qda5=qda5.scaled.model,
#                        nb5=nb5.scaled.model)
# b.model.list.6 <- list(rf6.rose=rf6.rose.model, 
#                        svmRadial6.rose=svmRadial6.rose.model,
#                        qda6.rose=qda6.rose.model,
#                        rf6=rf6.scaled.model,
#                        svmPoly6=svmPoly6.scaled.model,
#                        qda6=qda6.scaled.model)
# 
# # ROC plots
# dili1.roc <- lapply(b.model.list.1, function(x){
#     roc(predictor=x$pred$Positive, 
#         response=x$pred$obs)
# })
# dili3.roc <- lapply(b.model.list.3, function(x){
#     roc(predictor=x$pred$Positive, 
#         response=x$pred$obs)
# })
# dili5.roc <- lapply(b.model.list.5, function(x){
#     roc(predictor=x$pred$Positive, 
#         response=x$pred$obs)
# })
# dili6.roc <- lapply(b.model.list.6, function(x){
#     roc(predictor=x$pred$Positive, 
#         response=x$pred$obs)
# })
# 
# collect <- function(roc.object){
#     
#     FPR <- rev(1-roc.object[['specificities']])
#     TPR <- rev(roc.object[['sensitivities']])
#     bound <- cbind(TPR, FPR)
#     return(bound)
# }
# 
# dili1.out <- list()
# for (i in 1:length(dili1.roc)){
#     
#     values <- data.frame(collect(dili1.roc[[i]]))
#     model <- rep(paste(names(dili1.roc[i]), round(dili1.roc[[i]]$auc, 3), sep=': '),
#                  nrow(values))
#     dili1.out[[i]] <- data.frame(cbind(values, model))
# }
# comb1 <- do.call(rbind, dili1.out)
# 
# dili3.out <- list()
# for (i in 1:length(dili3.roc)){
#     
#     values <- data.frame(collect(dili3.roc[[i]]))
#     model <- rep(paste(names(dili3.roc[i]), round(dili3.roc[[i]]$auc, 3), sep=': '),
#                  nrow(values))
#     dili3.out[[i]] <- data.frame(cbind(values, model))
# }
# comb3 <- do.call(rbind, dili3.out)
# 
# dili5.out <- list()
# for (i in 1:length(dili5.roc)){
#     
#     values <- data.frame(collect(dili5.roc[[i]]))
#     model <- rep(paste(names(dili5.roc[i]), round(dili5.roc[[i]]$auc, 3), sep=': '),
#                  nrow(values))
#     dili5.out[[i]] <- data.frame(cbind(values, model))
# }
# comb5 <- do.call(rbind, dili5.out)
# 
# dili6.out <- list()
# for (i in 1:length(dili6.roc)){
#     
#     values <- data.frame(collect(dili6.roc[[i]]))
#     model <- rep(paste(names(dili6.roc[i]), round(dili6.roc[[i]]$auc, 3), sep=': '),
#                  nrow(values))
#     dili6.out[[i]] <- data.frame(cbind(values, model))
# }
# comb6 <- do.call(rbind, dili6.out)
# 
# roc.plot <- ggplot(comb1, aes(x=as.numeric(FPR), y=as.numeric(TPR), col=model)) + geom_line() + theme_bw() + 
#     geom_abline(intercept = 0, slope = 1, alpha=0.8, col='grey') + ylim(0, 1) + xlim(0, 1) +
#     theme(legend.text=element_text(size=14), 
#           legend.title = element_text(color='black', size=15),
#           axis.title.x=element_text(size=15, face='bold'), 
#           axis.title.y=element_text(size=15, face='bold'),
#           plot.title = element_text(face='bold', size=20),
#           axis.line = element_line(colour ='black', size=1)) + 
#     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#     labs(col='Model and AUC value', x='False Positive Rate', y='True Positive Rate', title='Tox21 Dili1 Model Performance: AUC-ROC Curve')
# 
# png(file="../tox_results/dili1_roc_plot.png", width = 800, height=600)
# roc.plot
# dev.off()
# 
# roc.plot <- ggplot(comb3, aes(x=as.numeric(FPR), y=as.numeric(TPR), col=model)) + geom_line() + theme_bw() + 
#     geom_abline(intercept = 0, slope = 1, alpha=0.8, col='grey') + ylim(0, 1) + xlim(0, 1) +
#     theme(legend.text=element_text(size=14), 
#           legend.title = element_text(color='black', size=15),
#           axis.title.x=element_text(size=15, face='bold'), 
#           axis.title.y=element_text(size=15, face='bold'),
#           plot.title = element_text(face='bold', size=20),
#           axis.line = element_line(colour ='black', size=1)) + 
#     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#     labs(col='Model and AUC value', x='False Positive Rate', y='True Positive Rate', title='Tox21 Dili3 Model Performance: AUC-ROC Curve')
# 
# png(file="../tox_results/dili3_roc_plot.png", width = 800, height=600)
# roc.plot
# dev.off()
# 
# roc.plot <- ggplot(comb5, aes(x=as.numeric(FPR), y=as.numeric(TPR), col=model)) + geom_line() + theme_bw() + 
#     geom_abline(intercept = 0, slope = 1, alpha=0.8, col='grey') + ylim(0, 1) + xlim(0, 1) +
#     theme(legend.text=element_text(size=14), 
#           legend.title = element_text(color='black', size=15),
#           axis.title.x=element_text(size=15, face='bold'), 
#           axis.title.y=element_text(size=15, face='bold'),
#           plot.title = element_text(face='bold', size=20),
#           axis.line = element_line(colour ='black', size=1)) + 
#     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#     labs(col='Model and AUC value', x='False Positive Rate', y='True Positive Rate', title='Tox21 Dili5 Model Performance: AUC-ROC Curve')
# 
# png(file="../tox_results/dili5_roc_plot.png", width = 800, height=600)
# roc.plot
# dev.off()
# 
# roc.plot <- ggplot(comb6, aes(x=as.numeric(FPR), y=as.numeric(TPR), col=model)) + geom_line() + theme_bw() + 
#     geom_abline(intercept = 0, slope = 1, alpha=0.8, col='grey') + ylim(0, 1) + xlim(0, 1) +
#     theme(legend.text=element_text(size=14), 
#           legend.title = element_text(color='black', size=15),
#           axis.title.x=element_text(size=15, face='bold'), 
#           axis.title.y=element_text(size=15, face='bold'),
#           plot.title = element_text(face='bold', size=20),
#           axis.line = element_line(colour ='black', size=1)) + 
#     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#     labs(col='Model and AUC value', x='False Positive Rate', y='True Positive Rate', title='Tox21 Dili6 Model Performance: AUC-ROC Curve')
# 
# png(file="../tox_results/dili6_roc_plot.png", width = 800, height=600)
# roc.plot
# dev.off()
# 
# 
# # mcc ==============
# library(mccr)
# library(tidyr)
# library(tidyverse)
# 
# rpart.model.list <- list(rpart1.scaled.model$pred,
#                          rpart3.scaled.model$pred,
#                          rpart5.scaled.model$pred,
#                          rpart6.scaled.model$pred)
# glm.model.list <- list(glm1.scaled.model$pred,
#                        glm3.scaled.model$pred,
#                        glm5.scaled.model$pred,
#                        glm6.scaled.model$pred)
# nb.model.list <- list(nb1.scaled.model$pred,
#                       nb3.scaled.model$pred,
#                       nb5.scaled.model$pred,
#                       nb6.scaled.model$pred)
# rf.model.list <- list(rf1.scaled.model$pred,
#                       rf3.scaled.model$pred,
#                       rf5.scaled.model$pred,
#                       rf6.scaled.model$pred)
# svmPoly.model.list <- list(svmPoly1.scaled.model$pred,
#                            svmPoly3.scaled.model$pred,
#                            svmPoly5.scaled.model$pred,
#                            svmPoly6.scaled.model$pred)
# svmRadial.model.list <- list(svmRadial1.scaled.model$pred,
#                              svmRadial3.scaled.model$pred,
#                              svmRadial5.scaled.model$pred,
#                              svmRadial6.scaled.model$pred)
# qda.model.list <- list(qda1.scaled.model$pred,
#                        qda3.scaled.model$pred,
#                        qda5.scaled.model$pred,
#                        qda6.scaled.model$pred)
# nnet.model.list <- list(nnet1.scaled.model$pred,
#                         nnet3.scaled.model$pred,
#                         nnet5.scaled.model$pred,
#                         nnet6.scaled.model$pred)
# 
# collect.model.mcc <- function(x.model.list, model.name=''){
#     
#     out <- list()
#     
#     for(i in 1:length(x.model.list)){
#         x_mcc <- as.data.frame(x.model.list[[i]]$pred)
#         x_mcc <- x_mcc[complete.cases(x_mcc), ]
#         x_mcc[,1] <- as.character(x_mcc[,1])
#         x_mcc[,2] <- as.character(x_mcc[,2])
#         x_mcc[x_mcc=="Positive"] <- 0
#         x_mcc[x_mcc=="Negative"] <- 1
#         
#         x.details <- x_mcc %>%
#             separate(Resample,c("cv","rep"),sep="\\.") %>%
#             group_by(rep) %>%
#             dplyr::summarise(mcc=mccr::mccr(obs, pred))
#         
#         names(x.details) <- c("5_fold_CV","DILI1")
#         out[[i]] <- x.details
#     }
#     output <- as.data.frame(do.call(cbind, out))
#     `5_fold_CV` <- paste0('Run', seq(1:100))
#     output[names(output) %in% c('5_fold_CV')] <- NULL
#     output <- as.data.frame(cbind(`5_fold_CV`, output))
#     names(output)[2:ncol(output)] <- c('DILI1', 'DILI3', 'DILI5', 'DILI6')
#     write.csv(output, file = paste('../output_files/p9.', model.name, '-crossvalidation-camda2020-UND', '.csv', sep=''), row.names=F)
# }
# 
# top1 <- list(DILI1=nb1.scaled.model,
#              DILI3=nnet3.scaled.model,
#              DILI5=glm5.scaled.model,
#              DILI6=rf6.scaled.model)
# top2 <- list(DILI1=nnet1.scaled.model,
#              DILI3=rf3.scaled.model,
#              DILI5=qda5.scaled.model,
#              DILI6=qda6.scaled.model)
# top3 <- list(DILI1=rpart1.scaled.model,
#              DILI3=glm3.scaled.model,
#              DILI5=nb5.scaled.model,
#              DILI6=glm6.scaled.model)
# 
# top1.res <- list(DILI1=rf1.smote.model,
#                  DILI3=rf3.smote.model,
#                  DILI5=rf5.rose.model,
#                  DILI6=rf6.rose.model)
# top2.res <- list(DILI1=svmRadial1.smote.model,
#                  DILI3=nb3.rose.model,
#                  DILI5=nb5.rose.model,
#                  DILI6=svmRadial6.rose.model)
# top3.res <- list(DILI1=nnet1.smote.model,
#                  DILI3=svmRadial3.rose.model,
#                  DILI5=qda5.rose.model,
#                  DILI6=qda6.rose.model)
# 
# collect.model.mcc(top1, model.name = 'top1')
# collect.model.mcc(top2, model.name = 'top2')
# collect.model.mcc(top3, model.name = 'top3')
# 
# collect.model.mcc(top1.res, model.name = 'resampled-top1')
# collect.model.mcc(top2.res, model.name = 'resampled-top2')
# collect.model.mcc(top3.res, model.name = 'resampled-top3')
# 
# collect.model.mcc(rpart.model.list, model.name = 'rpart')
# collect.model.mcc(nnet.model.list, model.name = 'nnet')
# collect.model.mcc(rf.model.list, model.name ='rf')
# collect.model.mcc(svmRadial.model.list, model.name ='svmRadial')
# collect.model.mcc(nb.model.list, model.name ='naive.bayes')
# collect.model.mcc(glm.model.list, model.name ='glm')
# collect.model.mcc(svmPoly.model.list, model.name ='svmPoly')
# collect.model.mcc(qda.model.list, model.name ='qda')
# 
# 
# stopCluster(cl)
# 
# 
# # # a function to calculate the matthew's correlation co-efficient ====
# # mcc <- function(x.model){
# #     cm <- confusionMatrix(x.model)[['table']]
# #     true.pos <- cm[1,1]
# #     false.pos <- cm[1,2]
# #     false.neg <- cm[2,1]
# #     true.neg <- cm[2,2]
# #
# #     above <- true.pos*true.neg - false.pos*false.neg
# #     below <- as.double(true.pos+false.pos)*as.double(true.pos+false.neg)*as.double(true.neg+false.pos)*as.double(true.neg+false.neg)
# #
# #     mcc <- above/sqrt(below)
# #     mcc
# # }
# #
# # collate.metrics <- function(sum.list, x.list, dili.name=''){
# #
# #     out.roc <- list()
# #     out.sens <- list()
# #     out.spec <- list()
# #     out.mcc <- list()
# #
# #     for(i in 1:length(sum.list)){
# #         roc <- as.data.frame(sum.list[[i]]$statistics$ROC[, 'Mean'])
# #         roc$dataset <- row.names(roc)
# #         row.names(roc) <- NULL
# #         out.roc[[i]] <- roc
# #
# #         sens <- as.data.frame(sum.list[[i]]$statistics$Sens[, 'Mean'])
# #         sens$dataset <- row.names(sens)
# #         row.names(sens) <- NULL
# #         out.sens[[i]] <- sens
# #
# #         spec <- as.data.frame(sum.list[[i]]$statistics$Spec[, 'Mean'])
# #         spec$dataset <- row.names(spec)
# #         row.names(spec) <- NULL
# #         out.spec[[i]] <- spec
# #     }
# #
# #     for(i in 1:length(x.list)){
# #         mcc.i <- data.frame(dataset=names(x.list)[i], mcc=mcc(x.list[[i]]))
# #         out.mcc[[i]] <- mcc.i
# #     }
# #
# #     roc <- do.call(rbind, out.roc)
# #     names(roc)[1] <- 'ROC'
# #     sens <- do.call(rbind, out.sens)
# #     names(sens)[1] <- 'Sens'
# #     spec <- do.call(rbind, out.spec)
# #     names(spec)[1] <- 'Spec'
# #     mcc <- do.call(rbind, out.mcc)
# #
# #     first <- merge(roc, sens, by='dataset')
# #     mid <- merge(first, spec, by='dataset')
# #     last <- merge(mid, mcc, by='dataset')
# #     last <- last %>% arrange(ROC)
# #
# #     write.csv(last, file = paste('../tox21/tox_results/tox_', dili.name, '.csv', sep=''), row.names=F)
# #     #print(last)
# # }
# #
# # collate.metrics(sum.list=summary.list.1, x.list=model.list.1, dili.name = 'dili1')
# # collate.metrics(sum.list=summary.list.3, x.list=model.list.3, dili.name = 'dili3')
# # collate.metrics(sum.list=summary.list.5, x.list=model.list.5, dili.name = 'dili5')
# # collate.metrics(sum.list=summary.list.6, x.list=model.list.6, dili.name = 'dili6')
# #
# #
# # # ============== Prediction and Validation ==================
# # # use the best models to predict
# # # none-resampled datasets
# #
# # non.res.best <- resamples(list(nb1=nb1.scaled.model,
# #                            nnet3=nnet3.scaled.model,
# #                            glm5=glm5.scaled.model,
# #                            rf6=rf6.scaled.model))
# # res.best <- resamples(list(rf1.smote=rf1.smote.model,
# #                            rf3.smote=rf3.smote.model,
# #                            rf5.rose=rf5.rose.model,
# #                            rf6.rose=rf6.rose.model))
# # best.non.res.model.list <- list(nb1=nb1.scaled.model,
# #                                 nnet3=nnet3.scaled.model,
# #                                 glm5=glm5.scaled.model,
# #                                 rf6=rf6.scaled.model)
# # best.res.model.list <- list(rf1.smote=rf1.smote.model,
# #                              rf3.smote=rf3.smote.model,
# #                              rf5.rose=rf5.rose.model,
# #                              rf6.rose=rf6.rose.model)
# #
# #
# # There are 190 validation samples
# # 18 of them have 24 NA columns >>> taking them out
# 
# 
# 
# # validation set
# validation.set <- subset(target.data, Training_Validation=='') # 195 obs
# tox.validate <- tox.data[rownames(tox.data) %in% validation.set$CAM_ID, ] # 190 obs
# 
# df.validation.dataset <- data.frame(sort(apply(tox.validate, 1, function(x){sum(is.na(x))}), decreasing = T))
# names(df.validation.dataset) <- 'No.of.NAs'
# tox.validate.names <- row.names(subset(df.validation.dataset, No.of.NAs==0))
# tox.validate <- tox.validate[row.names(tox.validate) %in% tox.validate.names, ]
# 
# # run predictions
# # run predictions
# run.predictions <- function(model.list, new.data, model.name='', resampled=F, test.data=T){
#     
#     output <- list()
#     for (i in 1:length(model.list)){
#         predicted.values <- as.vector(predict(model.list[[i]], newdata=new.data))
#         predicted.values[predicted.values == 'Positive'] <- 0
#         predicted.values[predicted.values == 'Negative'] <- 1
#         output[[paste(names(model.list)[i], sep='')]] <- as.numeric(predicted.values)
#     }
#     
#     # print(output)
#     # str(do.call(cbind, output))
#     result <- data.frame(cbind(CAM_ID=row.names(new.data), as.data.frame(do.call(cbind, output))))
#     #row.names(result) <- row.names(new.data)
#     #names(result)[1] <- 'CAM_ID'
#     
#     if (resampled==F) {
#         if (test.data==T) {
#             write.csv(result, file=paste('../output_files/p9.', model.name, '-predictions-camda2020-UND.csv', sep=''), row.names=F)
#         } else {
#             write.csv(result, file='../output_files/train_tox2_predictions.csv')
#         }
#     } else if (resampled==T) {
#         if (test.data==T) {
#             write.csv(result, file=paste('../output_files/p9.resampled-', model.name, '-predictions-camda2020-UND.csv', sep=''), row.names=F)
#         } else {
#             write.csv(result, file='../output_files/train_tox2_resampled_predictions.csv')
#         }
#     }
#     result
# }
# 
# run.predictions(top1, tox.validate, model.name = 'top1', resampled=F, test.data = T)
# run.predictions(top2, tox.validate, model.name = 'top2', resampled=F, test.data = T)
# run.predictions(top3, tox.validate, model.name = 'top3', resampled=F, test.data = T)
# 
# run.predictions(top1.res, tox.validate, model.name = 'top1', resampled=T, test.data = T)
# run.predictions(top2.res, tox.validate, model.name = 'top2', resampled=T, test.data = T)
# run.predictions(top3.res, tox.validate, model.name = 'top3', resampled=T, test.data = T)
# 
# run.predictions(best.non.res.model.list, tox.validate, resampled=F, test.data = T)
# run.predictions(best.res.model.list, tox.validate, resampled=T, test.data = T)
# 
# # run predictions on the training set ===============
# run.predictions(best.non.res.model.list, comp.tox.train, resampled=F, test.data = F)
# run.predictions(best.res.model.list, comp.tox.train, resampled=T, test.data = F)
# #
# #
# #
# #
# #
# #
# # # collect mcc ==============
# # library(mccr)
# # library(tidyr)
# # library(tidyverse)
# #
# # #rf is my random forest model for dili1
# #
# # collect.model.mcc <- function(x.model.list, model.name=''){
# #
# #     out <- list()
# #
# #     for(i in 1:length(x.model.list)){
# #         x_mcc <- as.data.frame(x.model.list[[i]]$pred)
# #         x_mcc[,1] <- as.character(x_mcc[,1])
# #         x_mcc[,2] <- as.character(x_mcc[,2])
# #         x_mcc[x_mcc=="Positive"] <- 0
# #         x_mcc[x_mcc=="Negative"] <- 1
# #
# #         x.details <- x_mcc %>%
# #             separate(Resample,c("cv","rep"),sep="\\.") %>%
# #             group_by(rep) %>%
# #             summarise(mcc=mccr(obs, pred))
# #
# #         names(x.details) <- c("5_fold_CV",paste("DILI", i, sep=''))
# #         out[[i]] <- x.details
# #     }
# #     output <- as.data.frame(do.call(cbind, out))
# #     `5_fold_CV` <- paste0('Run', seq(1:100))
# #     output[names(output) %in% c('5_fold_CV')] <- NULL
# #     output <- as.data.frame(cbind(`5_fold_CV`, output))
# #     names(output)[2:ncol(output)] <- c('DILI1', 'DILI3', 'DILI5', 'DILI6')
# #     print(output)
# #     #write.csv(output, file = paste('../a_results/to_dropbox/p1.', model.name, '-crossvalidation-camda2020-UND', '.csv', sep=''), row.names=F)
# # }
# #
# # collect.model.mcc(best.non.res.model.list, model.name = 'rf')
# #
# #
# #
# # x_mcc <- as.data.frame(nb1.scaled.model$pred)
# # x_mcc[,1] <- as.character(x_mcc[,1])
# # x_mcc[,2] <- as.character(x_mcc[,2])
# # x_mcc[x_mcc=="Positive"] <- 0
# # x_mcc[x_mcc=="Negative"] <- 1
# # x.details <- x_mcc %>%
# #     separate(Resample,c("cv","rep"),sep="\\.") %>%
# #     group_by(rep) %>%
# #     summarise(mcc=mccr(obs, pred))
# 
# 
# 
# 
# 
# 
# 
# 
# save()
# 
# 
# 
# 
# 
# 
# 
# 
# # mcc
# #
# # collect.model.mcc <- function(x.model.list, model.name=''){
# #
# #     out <- list()
# #
# #     for(i in 1:length(x.model.list)){
# #         x_mcc <- as.data.frame(x.model.list[[i]])
# #         x_mcc[,1] <- as.character(x_mcc[,1])
# #         x_mcc[,2] <- as.character(x_mcc[,2])
# #         x_mcc[x_mcc=="Positive"] <- 0
# #         x_mcc[x_mcc=="Negative"] <- 1
# #
# #         x.details <- x_mcc %>%
# #             separate(Resample,c("cv","rep"),sep="\\.") %>%
# #             group_by(rep) %>%
# #             summarise(mcc=mccr(obs, pred))
# #
# #         names(x.details) <- c("5_fold_CV","DILI1")
# #         out[[i]] <- x.details
# #     }
# #     output <- as.data.frame(do.call(cbind, out))
# #     `5_fold_CV` <- paste0('Run', seq(1:100))
# #     output[names(output) %in% c('5_fold_CV')] <- NULL
# #     output <- as.data.frame(cbind(`5_fold_CV`, output))
# #     names(output)[2:ncol(output)] <- c('DILI1', 'DILI3', 'DILI5', 'DILI6')
# #     write.csv(output, file = paste('../a_results/to_dropbox/p1.', model.name, '-crossvalidation-camda2020-UND', '.csv', sep=''), row.names=F)
# # }
