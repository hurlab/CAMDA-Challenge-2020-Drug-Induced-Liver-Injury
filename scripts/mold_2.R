
# Jun 27 2020
# This script was written to run classification models on mold2 data...
#...for CAMDA 2020

# Temi

library(rstudioapi)
currPath <- getActiveDocumentContext()$path
setwd(dirname(currPath))

# load libraries
source('../scripts/packages.R')

#setwd('/extData/NGS/hurlab/temi/projects/camda/mold2/')

# Read data
mold2.data <- read.csv('../data/p1-mold2-camda2020.csv',
                       header=T, stringsAsFactors = F, check.names = F, row.names='CAM_ID')
# target/label/training
target.data <- read.csv('../data/targets-camda2020.csv',
                        header=T, stringsAsFactors = F, check.names = F)

# mold2.data has all the data
# select the training set based on what is available in target.data
training.set <- subset(target.data, Training_Validation=='Training Set') # 422 obs
mold2.train.set <- mold2.data[rownames(mold2.data) %in% training.set$CAM_ID, ] # 422 obs
# >>> balanced

as.data.frame(apply(mold2.train.set, 1, function(x){
    sum(is.na(x))
}))

# create dili1 to dili6 datasets
dili1.train <- cbind(mold2.train.set, Class=as.factor(training.set$DILI1))
dili1.train$Class <- as.factor(mapvalues(dili1.train$Class, 
                                         from=c(0, 1), to=c('Positive', 'Negative')))

dili3.train <- cbind(mold2.train.set, Class=as.factor(training.set$DILI3))
dili3.train$Class <- as.factor(mapvalues(dili3.train$Class, 
                                         from=c(0, 1), to=c('Positive', 'Negative')))

dili5.train <- cbind(mold2.train.set, Class=as.factor(training.set$DILI5))
dili5.train$Class <- as.factor(mapvalues(dili5.train$Class, 
                                         from=c(0, 1), to=c('Positive', 'Negative')))

dili6.train <- cbind(mold2.train.set, Class=as.factor(training.set$DILI6))
dili6.train$Class <- as.factor(mapvalues(dili6.train$Class, 
                                         from=c(0, 1), to=c('Positive', 'Negative')))

# upsampled dataset
set.seed(1)
up.dili1.train <- upSample(x=dili1.train[, -ncol(dili1.train)],
                           y=as.factor(dili1.train$Class),
                           yname='Class')
set.seed(1)
up.dili3.train <- upSample(x=dili3.train[, -ncol(dili3.train)],
                           y=as.factor(dili3.train$Class),
                           yname='Class')
set.seed(1)
up.dili5.train <- upSample(x=dili5.train[, -ncol(dili5.train)],
                           y=as.factor(dili5.train$Class),
                           yname='Class')
set.seed(1)
up.dili6.train <- upSample(x=dili6.train[, -ncol(dili6.train)],
                           y=as.factor(dili6.train$Class),
                           yname='Class')
# smote dataset
set.seed(1)
smote.dili1.train <- bimba::SMOTE(data=dili1.train, perc_min = 50, k = 5)
set.seed(1)
smote.dili3.train <- bimba::SMOTE(data=dili3.train, perc_min = 50, k = 5)
set.seed(1)
smote.dili5.train <- bimba::SMOTE(data=dili5.train, perc_min = 50, k = 5)
set.seed(1)
smote.dili6.train <- bimba::SMOTE(data=dili6.train, perc_min = 50, k = 5)

# rose dataset
set.seed(1)
rose.dili1.train <- ROSE::ovun.sample(Class ~., data=dili1.train, 
                                      N=max(table(dili1.train$Class))*2, method='over', seed=1)$data
set.seed(1)
rose.dili3.train <- ROSE::ovun.sample(Class ~., data=dili3.train, 
                                      N=max(table(dili3.train$Class))*2, method='over', seed=1)$data
set.seed(1)
rose.dili5.train <- ROSE::ovun.sample(Class ~., data=dili5.train, 
                                      N=max(table(dili5.train$Class))*2, method='over', seed=1)$data
set.seed(1)
rose.dili6.train <- ROSE::ovun.sample(Class ~., data=dili6.train, 
                                      N=max(table(dili6.train$Class))*2, method='over', seed=1)$data

cl <- makeCluster(detectCores() - 56)
registerDoParallel(cl)

# use glm to remove zero variance variables
reduce.using.glm <- function(df){
    set.seed(1)
    fe.glm.fit <- glm(Class~., data=df, maxit=200, family=binomial)
    
    # eliminate variables with NA coefficients
    na.df <- data.frame(fe.glm.fit$coefficients)
    not.na <- row.names(na.df)[!is.na(na.df)][-1]
    
    df[, names(df) %in% not.na]
}


# dili1 dataset
dili1.train <- cbind(reduce.using.glm(dili1.train), as.factor(dili1.train$Class))
names(dili1.train)[ncol(dili1.train)] <- 'Class'

up.dili1.train <- cbind(reduce.using.glm(up.dili1.train), as.factor(up.dili1.train$Class))
names(up.dili1.train)[ncol(up.dili1.train)] <- 'Class'

smote.dili1.train <- cbind(reduce.using.glm(smote.dili1.train), as.factor(smote.dili1.train$Class))
names(smote.dili1.train)[ncol(smote.dili1.train)] <- 'Class'

rose.dili1.train <- cbind(reduce.using.glm(rose.dili1.train), as.factor(rose.dili1.train$Class))
names(rose.dili1.train)[ncol(rose.dili1.train)] <- 'Class'

# dili3 dataset
dili3.train <- cbind(reduce.using.glm(dili3.train), as.factor(dili3.train$Class))
names(dili3.train)[ncol(dili3.train)] <- 'Class'

up.dili3.train <- cbind(reduce.using.glm(up.dili3.train), as.factor(up.dili3.train$Class))
names(up.dili3.train)[ncol(up.dili3.train)] <- 'Class'

smote.dili3.train <- cbind(reduce.using.glm(smote.dili3.train), as.factor(smote.dili3.train$Class))
names(smote.dili3.train)[ncol(smote.dili3.train)] <- 'Class'

rose.dili3.train <- cbind(reduce.using.glm(rose.dili3.train), as.factor(rose.dili3.train$Class))
names(rose.dili3.train)[ncol(rose.dili3.train)] <- 'Class'

# dili5 dataset
dili5.train <- cbind(reduce.using.glm(dili5.train), as.factor(dili5.train$Class))
names(dili5.train)[ncol(dili5.train)] <- 'Class'

up.dili5.train <- cbind(reduce.using.glm(up.dili5.train), as.factor(up.dili5.train$Class))
names(up.dili5.train)[ncol(up.dili5.train)] <- 'Class'

smote.dili5.train <- cbind(reduce.using.glm(smote.dili5.train), as.factor(smote.dili5.train$Class))
names(smote.dili5.train)[ncol(smote.dili5.train)] <- 'Class'

rose.dili5.train <- cbind(reduce.using.glm(rose.dili5.train), as.factor(rose.dili5.train$Class))
names(rose.dili5.train)[ncol(rose.dili5.train)] <- 'Class'

# dili6 dataset
dili6.train <- cbind(reduce.using.glm(dili6.train), as.factor(dili6.train$Class))
names(dili6.train)[ncol(dili6.train)] <- 'Class'

up.dili6.train <- cbind(reduce.using.glm(up.dili6.train), as.factor(up.dili6.train$Class))
names(up.dili6.train)[ncol(up.dili6.train)] <- 'Class'

smote.dili6.train <- cbind(reduce.using.glm(smote.dili6.train), as.factor(smote.dili6.train$Class))
names(smote.dili6.train)[ncol(smote.dili6.train)] <- 'Class'

rose.dili6.train <- cbind(reduce.using.glm(rose.dili6.train), as.factor(rose.dili6.train$Class))
names(rose.dili6.train)[ncol(rose.dili6.train)] <- 'Class'


# data mining ======
# define a general control
general.ctrl <- trainControl(method='repeatedcv',
                             repeats=100, number=5,
                             savePredictions=T,
                             classProbs = T,
                             summaryFunction = twoClassSummary,
                             allowParallel = T)

print('Started training mold models...')

# ============ DILI1 ========================
# glm ============
set.seed(1)
glm1.model <- train(x=dili1.train[, -ncol(dili1.train)],
                    y=as.factor(dili1.train$Class),
                    family=binomial(link = 'logit'),
                    maxit=200, method='glm',
                    metric='ROC', trControl=general.ctrl,
                    preProcess=c('scale', 'center'))

set.seed(1)
up.glm1.model <- train(x=up.dili1.train[, -ncol(up.dili1.train)],
                       y=as.factor(up.dili1.train$Class),
                       family=binomial(link = 'logit'),
                       maxit=200, method='glm',
                       metric='ROC',trControl = general.ctrl,
                       preProcess=c('scale', 'center'))

set.seed(1)
smote.glm1.model <- train(x=smote.dili1.train[, -ncol(smote.dili1.train)],
                          y=as.factor(smote.dili1.train$Class),
                          family=binomial(link = 'logit'),
                          maxit=200, method='glm',
                          metric='ROC',trControl = general.ctrl,
                          preProcess=c('scale', 'center'))

set.seed(1)
rose.glm1.model <- train(x=rose.dili1.train[, -ncol(rose.dili1.train)],
                         y=as.factor(rose.dili1.train$Class),
                         family=binomial(link = 'logit'),
                         maxit=200, method='glm',
                         metric='ROC',trControl = general.ctrl,
                         preProcess=c('scale', 'center'))

# rpart =============
set.seed(1)
rpart1.model <- train(Class ~ ., data=dili1.train, method='rpart',
                      metric='ROC', trControl=general.ctrl,
                      preProcess=c('scale', 'center'))
set.seed(1)
up.rpart1.model <- train(Class ~ ., data=up.dili1.train, method='rpart',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))
set.seed(1)
smote.rpart1.model <- train(Class ~ ., data=smote.dili1.train, method='rpart',
                            metric='ROC', trControl=general.ctrl,
                            preProcess=c('scale', 'center'))
set.seed(1)
rose.rpart1.model <- train(Class ~ ., data=rose.dili1.train, method='rpart',
                           metric='ROC', trControl=general.ctrl,
                           preProcess=c('scale', 'center'))

# lda =======
set.seed(1)
lda1.model <- train(Class ~ ., data=dili1.train, method='lda',
                    metric='ROC', trControl=general.ctrl,
                    preProcess=c('scale', 'center'))
set.seed(1)
up.lda1.model <- train(Class ~ ., data=up.dili1.train, method='lda',
                       metric='ROC', trControl=general.ctrl,
                       preProcess=c('scale', 'center'))
set.seed(1)
smote.lda1.model <- train(Class ~ ., data=smote.dili1.train, method='lda',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
rose.lda1.model <- train(Class ~ ., data=rose.dili1.train, method='lda',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))

# svmLinear =====
set.seed(1)
svmLinear1.model <- train(Class ~ ., data=dili1.train, method='svmLinear',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
up.svmLinear1.model <- train(Class ~ ., data=up.dili1.train, method='svmLinear',
                             metric='ROC', trControl=general.ctrl,
                             preProcess=c('scale', 'center'))
set.seed(1)
smote.svmLinear1.model <- train(Class ~ ., data=smote.dili1.train, method='svmLinear',
                                metric='ROC', trControl=general.ctrl,
                                preProcess=c('scale', 'center'))
set.seed(1)
rose.svmLinear1.model <- train(Class ~ ., data=rose.dili1.train, method='svmLinear',
                               metric='ROC', trControl=general.ctrl,
                               preProcess=c('scale', 'center'))

# svmRadial ============
set.seed(1)
svmRadial1.model <- train(Class ~ ., data=dili1.train, method='svmRadial',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
up.svmRadial1.model <- train(Class ~ ., data=up.dili1.train, method='svmRadial',
                             metric='ROC', trControl=general.ctrl,
                             preProcess=c('scale', 'center'))
set.seed(1)
smote.svmRadial1.model <- train(Class ~ ., data=smote.dili1.train, method='svmRadial',
                                metric='ROC', trControl=general.ctrl,
                                preProcess=c('scale', 'center'))
set.seed(1)
rose.svmRadial1.model <- train(Class ~ ., data=rose.dili1.train, method='svmRadial',
                               metric='ROC', trControl=general.ctrl,
                               preProcess=c('scale', 'center'))

# naive bayes ============
set.seed(1)
nb1.model <- train(Class ~ ., data=dili1.train, method='naive_bayes',
                   metric='ROC', trControl=general.ctrl,
                   preProcess=c('scale', 'center'))
set.seed(1)
up.nb1.model <- train(Class ~ ., data=up.dili1.train, method='naive_bayes',
                      metric='ROC', trControl=general.ctrl,
                      preProcess=c('scale', 'center'))
set.seed(1)
smote.nb1.model <- train(Class ~ ., data=smote.dili1.train, method='naive_bayes',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))
set.seed(1)
rose.nb1.model <- train(Class ~ ., data=rose.dili1.train, method='naive_bayes',
                        metric='ROC', trControl=general.ctrl,
                        preProcess=c('scale', 'center'))

# nnet ============
set.seed(1)
nnet1.model <- train(Class ~ ., data=dili1.train, method='nnet',
                     metric='ROC', trControl=general.ctrl,
                     preProcess=c('scale', 'center'))
set.seed(1)
up.nnet1.model <- train(Class ~ ., data=up.dili1.train, method='nnet',
                        metric='ROC', trControl=general.ctrl,
                        preProcess=c('scale', 'center'))
set.seed(1)
smote.nnet1.model <- train(Class ~ ., data=smote.dili1.train, method='nnet',
                           metric='ROC', trControl=general.ctrl,
                           preProcess=c('scale', 'center'))
set.seed(1)
rose.nnet1.model <- train(Class ~ ., data=rose.dili1.train, method='nnet',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))

# svmPoly ============
set.seed(1)
svmPoly1.grid <- expand.grid(degree=c(1:3),
                             scale=seq(0.01, 0.02, length=5),
                             C=seq(0.1, 0.3, length=5))

set.seed(1)
svmPoly1.model <- train(Class ~ ., data=dili1.train, method='svmPoly',
                        metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly1.grid,
                        preProcess=c('scale', 'center'))
set.seed(1)
up.svmPoly1.model <- train(Class ~ ., data=up.dili1.train, method='svmPoly',
                           metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly1.grid,
                           preProcess=c('scale', 'center'))
set.seed(1)
smote.svmPoly1.model <- train(Class ~ ., data=smote.dili1.train, method='svmPoly',
                              metric='ROC', trControl=general.ctrl,tuneGrid=svmPoly1.grid,
                              preProcess=c('scale', 'center'))
set.seed(1)
rose.svmPoly1.model <- train(Class ~ ., data=rose.dili1.train, method='svmPoly',
                             metric='ROC', trControl=general.ctrl,tuneGrid=svmPoly1.grid,
                             preProcess=c('scale', 'center'))

# =============== DILI3 ================
# glm ============
set.seed(1)
glm3.model <- train(x=dili3.train[, -ncol(dili3.train)],
                    y=as.factor(dili3.train$Class),
                    family=binomial(link = 'logit'),
                    maxit=200, method='glm',
                    metric='ROC', trControl=general.ctrl,
                    preProcess=c('scale', 'center'))
set.seed(1)
up.glm3.model <- train(x=up.dili3.train[, -ncol(up.dili3.train)],
                       y=as.factor(up.dili3.train$Class),
                       family=binomial(link = 'logit'),
                       maxit=200, method='glm',
                       metric='ROC',trControl = general.ctrl,
                       preProcess=c('scale', 'center'))
set.seed(1)
smote.glm3.model <- train(x=smote.dili3.train[, -ncol(smote.dili3.train)],
                          y=as.factor(smote.dili3.train$Class),
                          family=binomial(link = 'logit'),
                          maxit=200, method='glm',
                          metric='ROC',trControl = general.ctrl,
                          preProcess=c('scale', 'center'))

set.seed(1)
rose.glm3.model <- train(x=rose.dili3.train[, -ncol(rose.dili3.train)],
                         y=as.factor(rose.dili3.train$Class),
                         family=binomial(link = 'logit'),
                         maxit=200, method='glm',
                         metric='ROC',trControl = general.ctrl,
                         preProcess=c('scale', 'center'))

# rpart =============
set.seed(1)
rpart3.model <- train(Class ~ ., data=dili3.train, method='rpart',
                      metric='ROC', trControl=general.ctrl,
                      preProcess=c('scale', 'center'))
set.seed(1)
up.rpart3.model <- train(Class ~ ., data=up.dili3.train, method='rpart',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))
set.seed(1)
smote.rpart3.model <- train(Class ~ ., data=smote.dili3.train, method='rpart',
                            metric='ROC', trControl=general.ctrl,
                            preProcess=c('scale', 'center'))
set.seed(1)
rose.rpart3.model <- train(Class ~ ., data=rose.dili3.train, method='rpart',
                           metric='ROC', trControl=general.ctrl,
                           preProcess=c('scale', 'center'))

# lda =======
set.seed(1)
lda3.model <- train(Class ~ ., data=dili3.train, method='lda',
                    metric='ROC', trControl=general.ctrl,
                    preProcess=c('scale', 'center'))
set.seed(1)
up.lda3.model <- train(Class ~ ., data=up.dili3.train, method='lda',
                       metric='ROC', trControl=general.ctrl,
                       preProcess=c('scale', 'center'))
set.seed(1)
smote.lda3.model <- train(Class ~ ., data=smote.dili3.train, method='lda',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
rose.lda3.model <- train(Class ~ ., data=rose.dili3.train, method='lda',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))

# svmLinear =====
set.seed(1)
svmLinear3.model <- train(Class ~ ., data=dili3.train, method='svmLinear',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
up.svmLinear3.model <- train(Class ~ ., data=up.dili3.train, method='svmLinear',
                             metric='ROC', trControl=general.ctrl,
                             preProcess=c('scale', 'center'))
set.seed(1)
smote.svmLinear3.model <- train(Class ~ ., data=smote.dili3.train, method='svmLinear',
                                metric='ROC', trControl=general.ctrl,
                                preProcess=c('scale', 'center'))
set.seed(1)
rose.svmLinear3.model <- train(Class ~ ., data=rose.dili3.train, method='svmLinear',
                               metric='ROC', trControl=general.ctrl,
                               preProcess=c('scale', 'center'))

# svmRadial ============
set.seed(1)
svmRadial3.model <- train(Class ~ ., data=dili3.train, method='svmRadial',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
up.svmRadial3.model <- train(Class ~ ., data=up.dili3.train, method='svmRadial',
                             metric='ROC', trControl=general.ctrl,
                             preProcess=c('scale', 'center'))
set.seed(1)
smote.svmRadial3.model <- train(Class ~ ., data=smote.dili3.train, method='svmRadial',
                                metric='ROC', trControl=general.ctrl,
                                preProcess=c('scale', 'center'))
set.seed(1)
rose.svmRadial3.model <- train(Class ~ ., data=rose.dili3.train, method='svmRadial',
                               metric='ROC', trControl=general.ctrl,
                               preProcess=c('scale', 'center'))

# naive bayes ============
set.seed(1)
nb3.model <- train(Class ~ ., data=dili3.train, method='naive_bayes',
                   metric='ROC', trControl=general.ctrl,
                   preProcess=c('scale', 'center'))
set.seed(1)
up.nb3.model <- train(Class ~ ., data=up.dili3.train, method='naive_bayes',
                      metric='ROC', trControl=general.ctrl,
                      preProcess=c('scale', 'center'))
set.seed(1)
smote.nb3.model <- train(Class ~ ., data=smote.dili3.train, method='naive_bayes',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))
set.seed(1)
rose.nb3.model <- train(Class ~ ., data=rose.dili3.train, method='naive_bayes',
                        metric='ROC', trControl=general.ctrl,
                        preProcess=c('scale', 'center'))

# nnet ============
set.seed(1)
nnet3.model <- train(Class ~ ., data=dili3.train, method='nnet',
                     metric='ROC', trControl=general.ctrl,
                     preProcess=c('scale', 'center'))
set.seed(1)
up.nnet3.model <- train(Class ~ ., data=up.dili3.train, method='nnet',
                        metric='ROC', trControl=general.ctrl,
                        preProcess=c('scale', 'center'))
set.seed(1)
smote.nnet3.model <- train(Class ~ ., data=smote.dili3.train, method='nnet',
                           metric='ROC', trControl=general.ctrl,
                           preProcess=c('scale', 'center'))
set.seed(1)
rose.nnet3.model <- train(Class ~ ., data=rose.dili3.train, method='nnet',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))

# svmPoly ============
set.seed(1)
svmPoly3.grid <- expand.grid(degree=c(1:3),
                             scale=seq(0.01, 0.02, length=5),
                             C=seq(0.1, 0.3, length=5))

set.seed(1)
svmPoly3.model <- train(Class ~ ., data=dili3.train, method='svmPoly',
                        metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly3.grid,
                        preProcess=c('scale', 'center'))
set.seed(1)
up.svmPoly3.model <- train(Class ~ ., data=up.dili3.train, method='svmPoly',
                           metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly3.grid,
                           preProcess=c('scale', 'center'))
set.seed(1)
smote.svmPoly3.model <- train(Class ~ ., data=smote.dili3.train, method='svmPoly',
                              metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly3.grid,
                              preProcess=c('scale', 'center'))
set.seed(1)
rose.svmPoly3.model <- train(Class ~ ., data=rose.dili3.train, method='svmPoly',
                             metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly3.grid,
                             preProcess=c('scale', 'center'))

# ============ DILI5 =================
# glm ============
set.seed(1)
glm5.model <- train(x=dili5.train[, -ncol(dili5.train)],
                    y=as.factor(dili5.train$Class),
                    family=binomial(link = 'logit'),
                    maxit=200, method='glm',
                    metric='ROC', trControl=general.ctrl,
                    preProcess=c('scale', 'center'))

set.seed(1)
up.glm5.model <- train(x=up.dili5.train[, -ncol(up.dili5.train)],
                       y=as.factor(up.dili5.train$Class),
                       family=binomial(link = 'logit'),
                       maxit=200, method='glm',
                       metric='ROC',trControl = general.ctrl,
                       preProcess=c('scale', 'center'))

set.seed(1)
smote.glm5.model <- train(x=smote.dili5.train[, -ncol(smote.dili5.train)],
                          y=as.factor(smote.dili5.train$Class),
                          family=binomial(link = 'logit'),
                          maxit=200, method='glm',
                          metric='ROC',trControl = general.ctrl,
                          preProcess=c('scale', 'center'))

set.seed(1)
rose.glm5.model <- train(x=rose.dili5.train[, -ncol(rose.dili5.train)],
                         y=as.factor(rose.dili5.train$Class),
                         family=binomial(link = 'logit'),
                         maxit=200, method='glm',
                         metric='ROC',trControl = general.ctrl,
                         preProcess=c('scale', 'center'))

# rpart =============
set.seed(1)
rpart5.model <- train(Class ~ ., data=dili5.train, method='rpart',
                      metric='ROC', trControl=general.ctrl,
                      preProcess=c('scale', 'center'))
set.seed(1)
up.rpart5.model <- train(Class ~ ., data=up.dili5.train, method='rpart',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))
set.seed(1)
smote.rpart5.model <- train(Class ~ ., data=smote.dili5.train, method='rpart',
                            metric='ROC', trControl=general.ctrl,
                            preProcess=c('scale', 'center'))
set.seed(1)
rose.rpart5.model <- train(Class ~ ., data=rose.dili5.train, method='rpart',
                           metric='ROC', trControl=general.ctrl,
                           preProcess=c('scale', 'center'))

# lda =======
set.seed(1)
lda5.model <- train(Class ~ ., data=dili5.train, method='lda',
                    metric='ROC', trControl=general.ctrl,
                    preProcess=c('scale', 'center'))
set.seed(1)
up.lda5.model <- train(Class ~ ., data=up.dili5.train, method='lda',
                       metric='ROC', trControl=general.ctrl,
                       preProcess=c('scale', 'center'))
set.seed(1)
smote.lda5.model <- train(Class ~ ., data=smote.dili5.train, method='lda',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
rose.lda5.model <- train(Class ~ ., data=rose.dili5.train, method='lda',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))

# svmLinear =====
set.seed(1)
svmLinear5.model <- train(Class ~ ., data=dili5.train, method='svmLinear',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
up.svmLinear5.model <- train(Class ~ ., data=up.dili5.train, method='svmLinear',
                             metric='ROC', trControl=general.ctrl,
                             preProcess=c('scale', 'center'))
set.seed(1)
smote.svmLinear5.model <- train(Class ~ ., data=smote.dili5.train, method='svmLinear',
                                metric='ROC', trControl=general.ctrl,
                                preProcess=c('scale', 'center'))
set.seed(1)
rose.svmLinear5.model <- train(Class ~ ., data=rose.dili5.train, method='svmLinear',
                               metric='ROC', trControl=general.ctrl,
                               preProcess=c('scale', 'center'))

# svmRadial ============
set.seed(1)
svmRadial5.model <- train(Class ~ ., data=dili5.train, method='svmRadial',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
up.svmRadial5.model <- train(Class ~ ., data=up.dili5.train, method='svmRadial',
                             metric='ROC', trControl=general.ctrl,
                             preProcess=c('scale', 'center'))
set.seed(1)
smote.svmRadial5.model <- train(Class ~ ., data=smote.dili5.train, method='svmRadial',
                                metric='ROC', trControl=general.ctrl,
                                preProcess=c('scale', 'center'))
set.seed(1)
rose.svmRadial5.model <- train(Class ~ ., data=rose.dili5.train, method='svmRadial',
                               metric='ROC', trControl=general.ctrl,
                               preProcess=c('scale', 'center'))

# naive bayes ============
set.seed(1)
nb5.model <- train(Class ~ ., data=dili5.train, method='naive_bayes',
                   metric='ROC', trControl=general.ctrl,
                   preProcess=c('scale', 'center'))
set.seed(1)
up.nb5.model <- train(Class ~ ., data=up.dili5.train, method='naive_bayes',
                      metric='ROC', trControl=general.ctrl,
                      preProcess=c('scale', 'center'))
set.seed(1)
smote.nb5.model <- train(Class ~ ., data=smote.dili5.train, method='naive_bayes',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))
set.seed(1)
rose.nb5.model <- train(Class ~ ., data=rose.dili5.train, method='naive_bayes',
                        metric='ROC', trControl=general.ctrl,
                        preProcess=c('scale', 'center'))

# nnet ============
set.seed(1)
nnet5.model <- train(Class ~ ., data=dili5.train, method='nnet',
                     metric='ROC', trControl=general.ctrl,
                     preProcess=c('scale', 'center'))
set.seed(1)
up.nnet5.model <- train(Class ~ ., data=up.dili5.train, method='nnet',
                        metric='ROC', trControl=general.ctrl,
                        preProcess=c('scale', 'center'))
set.seed(1)
smote.nnet5.model <- train(Class ~ ., data=smote.dili5.train, method='nnet',
                           metric='ROC', trControl=general.ctrl,
                           preProcess=c('scale', 'center'))
set.seed(1)
rose.nnet5.model <- train(Class ~ ., data=rose.dili5.train, method='nnet',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))

# svmPoly ============
set.seed(1)
svmPoly5.grid <- expand.grid(degree=c(1:3),
                             scale=seq(0.01, 0.02, length=5),
                             C=seq(0.1, 0.3, length=5))
set.seed(1)
svmPoly5.model <- train(Class ~ ., data=dili5.train, method='svmPoly',
                        metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly5.grid,
                        preProcess=c('scale', 'center'))
set.seed(1)
up.svmPoly5.model <- train(Class ~ ., data=up.dili5.train, method='svmPoly',
                           metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly5.grid,
                           preProcess=c('scale', 'center'))
set.seed(1)
smote.svmPoly5.model <- train(Class ~ ., data=smote.dili5.train, method='svmPoly',
                              metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly5.grid,
                              preProcess=c('scale', 'center'))
set.seed(1)
rose.svmPoly5.model <- train(Class ~ ., data=rose.dili5.train, method='svmPoly',
                             metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly5.grid,
                             preProcess=c('scale', 'center'))

# ========== DILI6 =========================
# glm ============
set.seed(1)
glm6.model <- train(x=dili6.train[, -ncol(dili6.train)],
                    y=as.factor(dili6.train$Class),
                    family=binomial(link = 'logit'),
                    maxit=200, method='glm',
                    metric='ROC', trControl=general.ctrl,
                    preProcess=c('scale', 'center'))

set.seed(1)
up.glm6.model <- train(x=up.dili6.train[, -ncol(up.dili6.train)],
                       y=as.factor(up.dili6.train$Class),
                       family=binomial(link = 'logit'),
                       maxit=200, method='glm',
                       metric='ROC',trControl = general.ctrl,
                       preProcess=c('scale', 'center'))

set.seed(1)
smote.glm6.model <- train(x=smote.dili6.train[, -ncol(smote.dili6.train)],
                          y=as.factor(smote.dili6.train$Class),
                          family=binomial(link = 'logit'),
                          maxit=200, method='glm',
                          metric='ROC',trControl = general.ctrl,
                          preProcess=c('scale', 'center'))

set.seed(1)
rose.glm6.model <- train(x=rose.dili6.train[, -ncol(rose.dili6.train)],
                         y=as.factor(rose.dili6.train$Class),
                         family=binomial(link = 'logit'),
                         maxit=200, method='glm',
                         metric='ROC',trControl = general.ctrl,
                         preProcess=c('scale', 'center'))

# rpart =============
set.seed(1)
rpart6.model <- train(Class ~ ., data=dili6.train, method='rpart',
                      metric='ROC', trControl=general.ctrl,
                      preProcess=c('scale', 'center'))
set.seed(1)
up.rpart6.model <- train(Class ~ ., data=up.dili6.train, method='rpart',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))
set.seed(1)
smote.rpart6.model <- train(Class ~ ., data=smote.dili6.train, method='rpart',
                            metric='ROC', trControl=general.ctrl,
                            preProcess=c('scale', 'center'))
set.seed(1)
rose.rpart6.model <- train(Class ~ ., data=rose.dili6.train, method='rpart',
                           metric='ROC', trControl=general.ctrl,
                           preProcess=c('scale', 'center'))

# lda =======
set.seed(1)
lda6.model <- train(Class ~ ., data=dili6.train, method='lda',
                    metric='ROC', trControl=general.ctrl,
                    preProcess=c('scale', 'center'))
set.seed(1)
up.lda6.model <- train(Class ~ ., data=up.dili6.train, method='lda',
                       metric='ROC', trControl=general.ctrl,
                       preProcess=c('scale', 'center'))
set.seed(1)
smote.lda6.model <- train(Class ~ ., data=smote.dili6.train, method='lda',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
rose.lda6.model <- train(Class ~ ., data=rose.dili6.train, method='lda',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))

# svmLinear =====
set.seed(1)
svmLinear6.model <- train(Class ~ ., data=dili6.train, method='svmLinear',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
up.svmLinear6.model <- train(Class ~ ., data=up.dili6.train, method='svmLinear',
                             metric='ROC', trControl=general.ctrl,
                             preProcess=c('scale', 'center'))
set.seed(1)
smote.svmLinear6.model <- train(Class ~ ., data=smote.dili6.train, method='svmLinear',
                                metric='ROC', trControl=general.ctrl,
                                preProcess=c('scale', 'center'))
set.seed(1)
rose.svmLinear6.model <- train(Class ~ ., data=rose.dili6.train, method='svmLinear',
                               metric='ROC', trControl=general.ctrl,
                               preProcess=c('scale', 'center'))

# svmRadial ============
set.seed(1)
svmRadial6.model <- train(Class ~ ., data=dili6.train, method='svmRadial',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))
set.seed(1)
up.svmRadial6.model <- train(Class ~ ., data=up.dili6.train, method='svmRadial',
                             metric='ROC', trControl=general.ctrl,
                             preProcess=c('scale', 'center'))
set.seed(1)
smote.svmRadial6.model <- train(Class ~ ., data=smote.dili6.train, method='svmRadial',
                                metric='ROC', trControl=general.ctrl,
                                preProcess=c('scale', 'center'))
set.seed(1)
rose.svmRadial6.model <- train(Class ~ ., data=rose.dili6.train, method='svmRadial',
                               metric='ROC', trControl=general.ctrl,
                               preProcess=c('scale', 'center'))

# naive bayes ============
set.seed(1)
nb6.model <- train(Class ~ ., data=dili6.train, method='naive_bayes',
                   metric='ROC', trControl=general.ctrl,
                   preProcess=c('scale', 'center'))
set.seed(1)
up.nb6.model <- train(Class ~ ., data=up.dili6.train, method='naive_bayes',
                      metric='ROC', trControl=general.ctrl,
                      preProcess=c('scale', 'center'))
set.seed(1)
smote.nb6.model <- train(Class ~ ., data=smote.dili6.train, method='naive_bayes',
                         metric='ROC', trControl=general.ctrl,
                         preProcess=c('scale', 'center'))
set.seed(1)
rose.nb6.model <- train(Class ~ ., data=rose.dili6.train, method='naive_bayes',
                        metric='ROC', trControl=general.ctrl,
                        preProcess=c('scale', 'center'))

# nnet ============
set.seed(1)
nnet6.model <- train(Class ~ ., data=dili6.train, method='nnet',
                     metric='ROC', trControl=general.ctrl,
                     preProcess=c('scale', 'center'))
set.seed(1)
up.nnet6.model <- train(Class ~ ., data=up.dili6.train, method='nnet',
                        metric='ROC', trControl=general.ctrl,
                        preProcess=c('scale', 'center'))
set.seed(1)
smote.nnet6.model <- train(Class ~ ., data=smote.dili6.train, method='nnet',
                           metric='ROC', trControl=general.ctrl,
                           preProcess=c('scale', 'center'))
set.seed(1)
rose.nnet6.model <- train(Class ~ ., data=rose.dili6.train, method='nnet',
                          metric='ROC', trControl=general.ctrl,
                          preProcess=c('scale', 'center'))

# svmPoly ============
set.seed(1)
svmPoly6.grid <- expand.grid(degree=c(1:3),
                             scale=seq(0.01, 0.02, length=5),
                             C=seq(0.1, 0.3, length=5))

set.seed(1)
svmPoly6.model <- train(Class ~ ., data=dili6.train, method='svmPoly',
                        metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly6.grid,
                        preProcess=c('scale', 'center'))
set.seed(1)
up.svmPoly6.model <- train(Class ~ ., data=up.dili6.train, method='svmPoly',
                           metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly6.grid,
                           preProcess=c('scale', 'center'))
set.seed(1)
smote.svmPoly6.model <- train(Class ~ ., data=smote.dili6.train, method='svmPoly',
                              metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly6.grid,
                              preProcess=c('scale', 'center'))
set.seed(1)
rose.svmPoly6.model <- train(Class ~ ., data=rose.dili6.train, method='svmPoly',
                             metric='ROC', trControl=general.ctrl, tuneGrid=svmPoly6.grid,
                             preProcess=c('scale', 'center'))

stopCluster(cl)

print('Finished training mold models...')


# ============= Evaluation =====================
# dili1 ===================================
res.glm1 <- resamples(list(glm1=glm1.model,
                           up.glm1=up.glm1.model,
                           smote.glm1=smote.glm1.model,
                           rose.glm1=rose.glm1.model))
sum.glm1 <- summary(res.glm1)

res.rpart1 <- resamples(list(rpart1=rpart1.model,
                             up.rpart1=up.rpart1.model,
                             smote.rpart1=smote.rpart1.model,
                             rose.rpart1=rose.rpart1.model))
sum.rpart1 <- summary(res.rpart1)

res.lda1 <- resamples(list(lda1=lda1.model,
                           up.lda1=up.lda1.model,
                           smote.lda1=smote.lda1.model,
                           rose.lda1=rose.lda1.model))
sum.lda1 <- summary(res.lda1)

res.svmLinear1 <- resamples(list(svmLinear1=svmLinear1.model,
                                 up.svmLinear1=up.svmLinear1.model,
                                 smote.svmLinear1=smote.svmLinear1.model,
                                 rose.svmLinear1=rose.svmLinear1.model))
sum.svmLinear1 <- summary(res.svmLinear1)

res.svmRadial1 <- resamples(list(svmRadial1=svmRadial1.model,
                                 up.svmRadial1=up.svmRadial1.model,
                                 smote.svmRadial1=smote.svmRadial1.model,
                                 rose.svmRadial1=rose.svmRadial1.model))
sum.svmRadial1 <- summary(res.svmRadial1)

res.nb1 <- resamples(list(nb1=nb1.model,
                          up.nb1=up.nb1.model,
                          smote.nb1=smote.nb1.model,
                          rose.nb1=rose.nb1.model))
sum.nb1 <- summary(res.nb1)

res.nnet1 <- resamples(list(nnet1=nnet1.model,
                            up.nnet1=up.nnet1.model,
                            smote.nnet1=smote.nnet1.model,
                            rose.nnet1=rose.nnet1.model))
sum.nnet1 <- summary(res.nnet1)

res.svmPoly1 <- resamples(list(svmPoly1=svmPoly1.model,
                               up.svmPoly1=up.svmPoly1.model,
                               smote.svmPoly1=smote.svmPoly1.model,
                               rose.svmPoly1=rose.svmPoly1.model))
sum.svmPoly1 <- summary(res.svmPoly1)

summary.list.1 <- list(sum.glm1,
                       sum.lda1,
                       sum.nb1,
                       sum.nnet1,
                       sum.rpart1,
                       sum.svmLinear1,
                       sum.svmPoly1,
                       sum.svmRadial1)

names(summary.list.1) <- paste0(c('glm', 'lda', 'nb', 'nnet', 'rpart', 'svmLinear', 'svmPoly', 'svmRadial'), 1)

model.list.1 <- list(glm1=glm1.model, up.glm1=up.glm1.model, rose.glm1=rose.glm1.model, smote.glm1=smote.glm1.model, 
                     nb1=nb1.model, rose.nb1=rose.nb1.model, smote.nb1=smote.nb1.model, up.nb1=up.nb1.model,
                     lda1=lda1.model, up.lda1=up.lda1.model, smote.lda1=smote.lda1.model, rose.lda1=rose.lda1.model,
                     rpart1=rpart1.model, smote.rpart1=smote.rpart1.model, up.rpart1=up.rpart1.model, rose.rpart1=rose.rpart1.model,
                     nnet1=nnet1.model, up.nnet1=up.nnet1.model, rose.nnet1=rose.nnet1.model, smote.nnet1=smote.nnet1.model, 
                     svmPoly1=svmPoly1.model, up.svmPoly1=up.svmPoly1.model, rose.svmPoly1=rose.svmPoly1.model, smote.svmPoly1=smote.svmPoly1.model,
                     up.svmLinear1=up.svmLinear1.model, rose.svmLinear1=rose.svmLinear1.model, smote.svmLinear1=smote.svmLinear1.model, svmLinear1=svmLinear1.model,
                     up.svmRadial1=up.svmRadial1.model, smote.svmRadial1=smote.svmRadial1.model, rose.svmRadial1=rose.svmRadial1.model, svmRadial1=svmRadial1.model)


# dili3 =========================================
res.glm3 <- resamples(list(glm3=glm3.model,
                           up.glm3=up.glm3.model,
                           smote.glm3=smote.glm3.model,
                           rose.glm3=rose.glm3.model))
sum.glm3 <- summary(res.glm3)

res.rpart3 <- resamples(list(rpart3=rpart3.model,
                             up.rpart3=up.rpart3.model,
                             smote.rpart3=smote.rpart3.model,
                             rose.rpart3=rose.rpart3.model))
sum.rpart3 <- summary(res.rpart3)

res.lda3 <- resamples(list(lda3=lda3.model,
                           up.lda3=up.lda3.model,
                           smote.lda3=smote.lda3.model,
                           rose.lda3=rose.lda3.model))
sum.lda3 <- summary(res.lda3)

res.svmLinear3 <- resamples(list(svmLinear3=svmLinear3.model,
                                 up.svmLinear3=up.svmLinear3.model,
                                 smote.svmLinear3=smote.svmLinear3.model,
                                 rose.svmLinear3=rose.svmLinear3.model))
sum.svmLinear3 <- summary(res.svmLinear3)

res.svmRadial3 <- resamples(list(svmRadial3=svmRadial3.model,
                                 up.svmRadial3=up.svmRadial3.model,
                                 smote.svmRadial3=smote.svmRadial3.model,
                                 rose.svmRadial3=rose.svmRadial3.model))
sum.svmRadial3 <- summary(res.svmRadial3)

res.nb3 <- resamples(list(nb3=nb3.model,
                          up.nb3=up.nb3.model,
                          smote.nb3=smote.nb3.model,
                          rose.nb3=rose.nb3.model))
sum.nb3 <- summary(res.nb3)

res.nnet3 <- resamples(list(nnet3=nnet3.model,
                            up.nnet3=up.nnet3.model,
                            smote.nnet3=smote.nnet3.model,
                            rose.nnet3=rose.nnet3.model))
sum.nnet3 <- summary(res.nnet3)

res.svmPoly3 <- resamples(list(svmPoly3=svmPoly3.model,
                               up.svmPoly3=up.svmPoly3.model,
                               smote.svmPoly3=smote.svmPoly3.model,
                               rose.svmPoly3=rose.svmPoly3.model))
sum.svmPoly3 <- summary(res.svmPoly3)

summary.list.3 <- list(sum.glm3,
                       sum.lda3,
                       sum.nb3,
                       sum.nnet3,
                       sum.rpart3,
                       sum.svmLinear3,
                       sum.svmPoly3,
                       sum.svmRadial3)

names(summary.list.3) <- paste0(c('glm', 'lda', 'nb', 'nnet', 'rpart', 'svmLinear', 'svmPoly', 'svmRadial'), 3)

model.list.3 <- list(glm3=glm3.model, up.glm3=up.glm3.model, rose.glm3=rose.glm3.model, smote.glm3=smote.glm3.model, 
                     nb3=nb3.model, rose.nb3=rose.nb3.model, smote.nb3=smote.nb3.model, up.nb3=up.nb3.model,
                     lda3=lda3.model, up.lda3=up.lda3.model, smote.lda3=smote.lda3.model, rose.lda3=rose.lda3.model,
                     rpart3=rpart3.model, smote.rpart3=smote.rpart3.model, up.rpart3=up.rpart3.model, rose.rpart3=rose.rpart3.model,
                     nnet3=nnet3.model, up.nnet3=up.nnet3.model, rose.nnet3=rose.nnet3.model, smote.nnet3=smote.nnet3.model, 
                     svmPoly3=svmPoly3.model, up.svmPoly3=up.svmPoly3.model, rose.svmPoly3=rose.svmPoly3.model, smote.svmPoly3=smote.svmPoly3.model,
                     up.svmLinear3=up.svmLinear3.model, rose.svmLinear3=rose.svmLinear3.model, smote.svmLinear3=smote.svmLinear3.model, svmLinear3=svmLinear3.model,
                     up.svmRadial3=up.svmRadial3.model, smote.svmRadial3=smote.svmRadial3.model, rose.svmRadial3=rose.svmRadial3.model, svmRadial3=svmRadial3.model)

# dili5 ======================================================
res.glm5 <- resamples(list(glm5=glm5.model,
                           up.glm5=up.glm5.model,
                           smote.glm5=smote.glm5.model,
                           rose.glm5=rose.glm5.model))
sum.glm5 <- summary(res.glm5)

res.rpart5 <- resamples(list(rpart5=rpart5.model,
                             up.rpart5=up.rpart5.model,
                             smote.rpart5=smote.rpart5.model,
                             rose.rpart5=rose.rpart5.model))
sum.rpart5 <- summary(res.rpart5)

res.lda5 <- resamples(list(lda5=lda5.model,
                           up.lda5=up.lda5.model,
                           smote.lda5=smote.lda5.model,
                           rose.lda5=rose.lda5.model))
sum.lda5 <- summary(res.lda5)

res.svmLinear5 <- resamples(list(svmLinear5=svmLinear5.model,
                                 up.svmLinear5=up.svmLinear5.model,
                                 smote.svmLinear5=smote.svmLinear5.model,
                                 rose.svmLinear5=rose.svmLinear5.model))
sum.svmLinear5 <- summary(res.svmLinear5)

res.svmRadial5 <- resamples(list(svmRadial5=svmRadial5.model,
                                 up.svmRadial5=up.svmRadial5.model,
                                 smote.svmRadial5=smote.svmRadial5.model,
                                 rose.svmRadial5=rose.svmRadial5.model))
sum.svmRadial5 <- summary(res.svmRadial5)

res.nb5 <- resamples(list(nb5=nb5.model,
                          up.nb5=up.nb5.model,
                          smote.nb5=smote.nb5.model,
                          rose.nb5=rose.nb5.model))
sum.nb5 <- summary(res.nb5)

res.nnet5 <- resamples(list(nnet5=nnet5.model,
                            up.nnet5=up.nnet5.model,
                            smote.nnet5=smote.nnet5.model,
                            rose.nnet5=rose.nnet5.model))
sum.nnet5 <- summary(res.nnet5)

res.svmPoly5 <- resamples(list(svmPoly5=svmPoly5.model,
                               up.svmPoly5=up.svmPoly5.model,
                               smote.svmPoly5=smote.svmPoly5.model,
                               rose.svmPoly5=rose.svmPoly5.model))
sum.svmPoly5 <- summary(res.svmPoly5)

summary.list.5 <- list(sum.glm5,
                       sum.lda5,
                       sum.nb5,
                       sum.nnet5,
                       sum.rpart5,
                       sum.svmLinear5,
                       sum.svmPoly5,
                       sum.svmRadial5)

names(summary.list.5) <- paste0(c('glm', 'lda', 'nb', 'nnet', 'rpart', 'svmLinear', 'svmPoly', 'svmRadial'), 5)

model.list.5 <- list(glm5=glm5.model, up.glm5=up.glm5.model, rose.glm5=rose.glm5.model, smote.glm5=smote.glm5.model, 
                     nb5=nb5.model, rose.nb5=rose.nb5.model, smote.nb5=smote.nb5.model, up.nb5=up.nb5.model,
                     lda5=lda5.model, up.lda5=up.lda5.model, smote.lda5=smote.lda5.model, rose.lda5=rose.lda5.model,
                     rpart5=rpart5.model, smote.rpart5=smote.rpart5.model, up.rpart5=up.rpart5.model, rose.rpart5=rose.rpart5.model,
                     nnet5=nnet5.model, up.nnet5=up.nnet5.model, rose.nnet5=rose.nnet5.model, smote.nnet5=smote.nnet5.model, 
                     svmPoly5=svmPoly5.model, up.svmPoly5=up.svmPoly5.model, rose.svmPoly5=rose.svmPoly5.model, smote.svmPoly5=smote.svmPoly5.model,
                     up.svmLinear5=up.svmLinear5.model, rose.svmLinear5=rose.svmLinear5.model, smote.svmLinear5=smote.svmLinear5.model, svmLinear5=svmLinear5.model,
                     up.svmRadial5=up.svmRadial5.model, smote.svmRadial5=smote.svmRadial5.model, rose.svmRadial5=rose.svmRadial5.model, svmRadial5=svmRadial5.model)


# dili6 =========================================
res.glm6 <- resamples(list(glm6=glm6.model,
                           up.glm6=up.glm6.model,
                           smote.glm6=smote.glm6.model,
                           rose.glm6=rose.glm6.model))
sum.glm6 <- summary(res.glm6)

res.rpart6 <- resamples(list(rpart6=rpart6.model,
                             up.rpart6=up.rpart6.model,
                             smote.rpart6=smote.rpart6.model,
                             rose.rpart6=rose.rpart6.model))
sum.rpart6 <- summary(res.rpart6)

res.lda6 <- resamples(list(lda6=lda6.model,
                           up.lda6=up.lda6.model,
                           smote.lda6=smote.lda6.model,
                           rose.lda6=rose.lda6.model))
sum.lda6 <- summary(res.lda6)

res.svmLinear6 <- resamples(list(svmLinear6=svmLinear6.model,
                                 up.svmLinear6=up.svmLinear6.model,
                                 smote.svmLinear6=smote.svmLinear6.model,
                                 rose.svmLinear6=rose.svmLinear6.model))
sum.svmLinear6 <- summary(res.svmLinear6)

res.svmRadial6 <- resamples(list(svmRadial6=svmRadial6.model,
                                 up.svmRadial6=up.svmRadial6.model,
                                 smote.svmRadial6=smote.svmRadial6.model,
                                 rose.svmRadial6=rose.svmRadial6.model))
sum.svmRadial6 <- summary(res.svmRadial6)

res.nb6 <- resamples(list(nb6=nb6.model,
                          up.nb6=up.nb6.model,
                          smote.nb6=smote.nb6.model,
                          rose.nb6=rose.nb6.model))
sum.nb6 <- summary(res.nb6)

res.nnet6 <- resamples(list(nnet6=nnet6.model,
                            up.nnet6=up.nnet6.model,
                            smote.nnet6=smote.nnet6.model,
                            rose.nnet6=rose.nnet6.model))
sum.nnet6 <- summary(res.nnet6)

res.svmPoly6 <- resamples(list(svmPoly6=svmPoly6.model,
                               up.svmPoly6=up.svmPoly6.model,
                               smote.svmPoly6=smote.svmPoly6.model,
                               rose.svmPoly6=rose.svmPoly6.model))
sum.svmPoly6 <- summary(res.svmPoly6)

summary.list.6 <- list(sum.glm6,
                       sum.lda6,
                       sum.nb6,
                       sum.nnet6,
                       sum.rpart6,
                       sum.svmLinear6,
                       sum.svmPoly6,
                       sum.svmRadial6)

names(summary.list.6) <- paste0(c('glm', 'lda', 'nb', 'nnet', 'rpart', 'svmLinear', 'svmPoly', 'svmRadial'), 6)

model.list.6 <- list(glm6=glm6.model, up.glm6=up.glm6.model, rose.glm6=rose.glm6.model, smote.glm6=smote.glm6.model, 
                     nb6=nb6.model, rose.nb6=rose.nb6.model, smote.nb6=smote.nb6.model, up.nb6=up.nb6.model,
                     lda6=lda6.model, up.lda6=up.lda6.model, smote.lda6=smote.lda6.model, rose.lda6=rose.lda6.model,
                     rpart6=rpart6.model, smote.rpart6=smote.rpart6.model, up.rpart6=up.rpart6.model, rose.rpart6=rose.rpart6.model,
                     nnet6=nnet6.model, up.nnet6=up.nnet6.model, rose.nnet6=rose.nnet6.model, smote.nnet6=smote.nnet6.model, 
                     svmPoly6=svmPoly6.model, up.svmPoly6=up.svmPoly6.model, rose.svmPoly6=rose.svmPoly6.model, smote.svmPoly6=smote.svmPoly6.model,
                     up.svmLinear6=up.svmLinear6.model, rose.svmLinear6=rose.svmLinear6.model, smote.svmLinear6=smote.svmLinear6.model, svmLinear6=svmLinear6.model,
                     up.svmRadial6=up.svmRadial6.model, smote.svmRadial6=smote.svmRadial6.model, rose.svmRadial6=rose.svmRadial6.model, svmRadial6=svmRadial6.model)

# save the models into a list

mold_models <- list(mold_dili1=model.list.1, mold_dili3=model.list.3, mold_dili5=model.list.5, mold_dili6=model.list.6)

save(mold_models, file = '../models/mold_models.RData')

rm(list=ls())


# Automatically select the top 3 models









# 
# 
# 
# # using the top 3 models in each category ==================
# b.model.list.1 <- list(smote.svmPoly1=smote.svmPoly1.model, 
#                        rose.svmPoly1=rose.svmPoly1.model,
#                        up.svmPoly1=up.svmPoly1.model,
#                        svmPoly1=svmPoly1.model,
#                        lda1=lda1.model,
#                        rpart1=rpart1.model)
# b.model.list.3 <- list(smote.svmPoly3=smote.svmPoly3.model, 
#                        rose.svmPoly3=rose.svmPoly3.model,
#                        up.svmPoly3=up.svmPoly3.model,
#                        svmPoly3=svmPoly3.model,
#                        svmRadial3=svmRadial3.model,
#                        svmLinear3=svmLinear3.model)
# b.model.list.5 <- list(rose.svmLinear5=rose.svmLinear5.model, 
#                        rose.svmPoly5=rose.svmPoly5.model,
#                        rose.svmRadial5=rose.svmRadial5.model,
#                        nb5=nb5.model,
#                        nnet5=nnet5.model,
#                        svmRadial5=svmRadial5.model)
# b.model.list.6 <- list(rose.rpart6=rose.rpart6.model, 
#                        up.rpart6=up.rpart6.model,
#                        smote.rpart6=smote.rpart6.model,
#                        rpart6=rpart6.model,
#                        nnet6=nnet6.model,
#                        svmRadial6=svmRadial6.model)
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
# roc.plot <- ggplot(comb3, aes(x=as.numeric(FPR), y=as.numeric(TPR), col=model)) + geom_line() + theme_bw() + 
#     geom_abline(intercept = 0, slope = 1, alpha=0.8, col='grey') + ylim(0, 1) + xlim(0, 1) +
#     theme(legend.text=element_text(size=14), 
#           legend.title = element_text(color='black', size=15),
#           axis.title.x=element_text(size=15, face='bold'), 
#           axis.title.y=element_text(size=15, face='bold'),
#           plot.title = element_text(face='bold', size=20),
#           axis.line = element_line(colour ='black', size=1)) + 
#     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#     labs(col='Model and AUC value', x='False Positive Rate', y='True Positive Rate', title='Mold2 Dili3 Model Performance: AUC-ROC Curve')
# 
# png(file="../output_plots/mold_results_dili3_roc_plot.png", width = 800, height=600)
# roc.plot
# dev.off()
# 
# 
# # mcc ==============
# 
# #rf is my random forest model for dili1 
# rpart.model.list <- list(rpart1.model$pred,
#                          rpart3.model$pred,
#                          rpart5.model$pred,
#                          rpart6.model$pred)
# glm.model.list <- list(glm1.model$pred,
#                        glm3.model$pred,
#                        glm5.model$pred,
#                        glm6.model$pred)
# nb.model.list <- list(nb1.model$pred,
#                       nb3.model$pred,
#                       nb5.model$pred,
#                       nb6.model$pred)
# svmLinear.model.list <- list(svmLinear1.model$pred,
#                              svmLinear3.model$pred,
#                              svmLinear5.model$pred,
#                              svmLinear6.model$pred)
# svmPoly.model.list <- list(svmPoly1.model$pred,
#                            svmPoly3.model$pred,
#                            svmPoly5.model$pred,
#                            svmPoly6.model$pred)
# svmRadial.model.list <- list(svmRadial1.model$pred,
#                              svmRadial3.model$pred,
#                              svmRadial5.model$pred,
#                              svmRadial6.model$pred)
# lda.model.list <- list(lda1.model$pred,
#                        lda3.model$pred,
#                        lda5.model$pred,
#                        lda6.model$pred)
# nnet.model.list <- list(nnet1.model$pred,
#                         nnet3.model$pred,
#                         nnet5.model$pred,
#                         nnet6.model$pred)
# 
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
#     write.csv(output, file = paste('../output_files/p1.', model.name, '-crossvalidation-camda2020-UND', '.csv', sep=''), row.names=F)
#     #output
# }
# 
# top1 <- list(DILI1=svmPoly1.model,
#              DILI3=svmPoly3.model,
#              DILI5=nb5.model,
#              DILI6=rpart6.model)
# top2 <- list(DILI1=lda1.model,
#              DILI3=nnet3.model,
#              DILI5=nnet5.model,
#              DILI6=nnet6.model)
# top3 <- list(DILI1=rpart1.model,
#              DILI3=lda3.model,
#              DILI5=rpart5.model,
#              DILI6=svmRadial6.model)
# 
# top1.res <- list(DILI1=smote.svmPoly1.model,
#                  DILI3=rose.svmPoly3.model,
#                  DILI5=rose.svmPoly5.model,
#                  DILI6=rose.rpart6.model)
# top2.res <- list(DILI1=rose.svmPoly1.model,
#                  DILI3=smote.svmPoly3.model,
#                  DILI5=rose.svmRadial5.model,
#                  DILI6=up.rpart6.model)
# top3.res <- list(DILI1=up.svmPoly1.model,
#                  DILI3=up.svmPoly3.model,
#                  DILI5=rose.svmLinear5.model,
#                  DILI6=smote.rpart6.model)
# 
# collect.model.mcc(top1, model.name = 'top1')
# collect.model.mcc(top2, model.name = 'top2')
# collect.model.mcc(top3, model.name = 'top3')
# 
# collect.model.mcc(top1.res, model.name = 'resampled-top1')
# collect.model.mcc(top2.res, model.name = 'resampled-top2')
# collect.model.mcc(top3.res, model.name = 'resampled-top3')
# 
# 
# 
# # collect.model.mcc(rpart.model.list, model.name = 'rpart')
# # collect.model.mcc(nnet.model.list, model.name = 'nnet')
# # collect.model.mcc(svmLinear.model.list, model.name ='svmLinear')
# # collect.model.mcc(svmRadial.model.list, model.name ='svmRadial')
# # collect.model.mcc(nb.model.list, model.name ='naive.bayes')
# # collect.model.mcc(glm.model.list, model.name ='glm')
# # collect.model.mcc(svmPoly.model.list, model.name ='svmPoly')
# # collect.model.mcc(lda.model.list, model.name ='lda')
# 
# 
# 
# # 
# # # a function to calculate the matthew's correlation co-efficient ====
# # mcc <- function(x.model){
# #     cm <- confusionMatrix(x.model)[['table']]
# #     
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
# #     write.csv(last, file = paste('../mold2/mold_results/', dili.name, 'metrics.csv', sep=''), row.names=T)
# # }
# # 
# # collate.metrics(sum.list=summary.list.1, x.list=model.list.1, dili.name = 'dili1')
# # collate.metrics(sum.list=summary.list.3, x.list=model.list.3, dili.name = 'dili3')
# # collate.metrics(sum.list=summary.list.5, x.list=model.list.5, dili.name = 'dili5')
# # collate.metrics(sum.list=summary.list.6, x.list=model.list.6, dili.name = 'dili6')
# # 
# # 
# ============== Prediction and Validation ==================
# use the best models to predict
# none-resampled datasets

non.res.best <- resamples(list(rpart6=rpart6.model,
                               nb5=nb5.model,
                               svmPoly3=svmPoly3.model,
                               svmPoly1=svmPoly1.model))
res.best <- resamples(list(rose.svmPoly5=rose.svmPoly5.model,
                           rose.rpart6=rose.rpart6.model,
                           smote.svmPoly1=smote.svmPoly1.model,
                           rose.svmPoly3=rose.svmPoly3.model))

best.non.res.model.list <- list(rpart6=rpart6.model,
                                nb5=nb5.model,
                                svmPoly3=svmPoly3.model,
                                svmPoly1=svmPoly1.model)
best.res.model.list <- list(rose.svmPoly5=rose.svmPoly5.model,
                            rose.rpart6=rose.rpart6.model,
                            smote.svmPoly1=smote.svmPoly1.model,
                            rose.svmPoly3=rose.svmPoly3.model)


# There are 190 validation samples
# 18 of them have 24 NA columns >>> taking them out


#mold2.train <- subset(mold2.train, all_gender_all != 0)

# validation set
validation.set <- subset(target.data, Training_Validation=='') # 195 obs
mold2.validate <- mold2.data[rownames(mold2.data) %in% validation.set$CAM_ID, ] # 190 obs

df.validation.dataset <- data.frame(sort(apply(mold2.validate, 1, function(x){sum(is.na(x))}), decreasing = T))
names(df.validation.dataset) <- 'No.of.NAs'
mold2.validate.names <- row.names(subset(df.validation.dataset, No.of.NAs==0))
mold2.validate <- mold2.validate[row.names(mold2.validate) %in% mold2.validate.names, ]

# v.names <- row.names(mold2.validate)
# mold2.validate <- mold2.validate %>%
#     mutate(ratio_dili_all=dili_gender_all/all_gender_all,
#            male_rate=(dili_gender_all*dili_gender_male_percentage)/(all_gender_male_percentage*all_gender_all),
#            female_rate=(dili_gender_all*dili_gender_female_percentage)/(all_gender_female_percentage*all_gender_all))
# row.names(mold2.validate) <- v.names

# model list
# rpart.model.list <- list(DILI1=rpart1.model,
#                          DILI3=rpart3.model,
#                          DILI5=rpart5.model,
#                          DILI6=rpart6.model)




# glm.model.list <- list(DILI1=glm1.model,
#                        DILI3=glm3.model,
#                        DILI5=glm5.model,
#                        DILI6=glm6.model)
# nb.model.list <- list(DILI1=nb1.model,
#                       DILI3=nb3.model,
#                       DILI5=nb5.model,
#                       DILI6=nb6.model)
# svmLinear.model.list <- list(DILI1=svmLinear1.model,
#                              DILI3=svmLinear3.model,
#                              DILI5=svmLinear5.model,
#                              DILI6=svmLinear6.model)
# svmPoly.model.list <- list(DILI1=svmPoly1.model,
#                            DILI3=svmPoly3.model,
#                            DILI5=svmPoly5.model,
#                            DILI6=svmPoly6.model)
# svmRadial.model.list <- list(DILI1=svmRadial1.model,
#                              DILI3=svmRadial3.model,
#                              DILI5=svmRadial5.model,
#                              DILI6=svmRadial6.model)
# lda.model.list <- list(DILI1=lda1.model,
#                        DILI3=lda3.model,
#                        DILI5=lda5.model,
#                        DILI6=lda6.model)
# nnet.model.list <- list(DILI1=nnet1.model,
#                         DILI3=nnet3.model,
#                         DILI5=nnet5.model,
#                         DILI6=nnet6.model)

# run predictions
run.predictions <- function(model.list, new.data, model.name='', resampled=F, test.data=T){

    output <- list()
    for (i in 1:length(model.list)){
        predicted.values <- as.vector(predict(model.list[[i]], newdata=new.data))
        predicted.values[predicted.values == 'Positive'] <- 0
        predicted.values[predicted.values == 'Negative'] <- 1
        output[[paste(names(model.list)[i], sep='')]] <- as.numeric(predicted.values)
    }

    # print(output)
    # str(do.call(cbind, output))
    result <- data.frame(cbind(CAM_ID=row.names(new.data), as.data.frame(do.call(cbind, output))))
    #row.names(result) <- row.names(new.data)
    #names(result)[1] <- 'CAM_ID'

    if (resampled==F) {
        if (test.data==T) {
            write.csv(result, file=paste('../output_files/p1.', model.name, '-predictions-camda2020-UND.csv', sep=''), row.names = F)
        } else {
            write.csv(result, file='../output_files/mold_results_train_mold2_predictions.csv')
        }
    } else if (resampled==T) {
        if (test.data==T) {
            write.csv(result, file=paste('../output_files/p1.resampled-', model.name, '-predictions-camda2020-UND.csv', sep=''), row.names = F)
        } else {
            write.csv(result, file='../output_files/mold_results_train_mold2_resampled_predictions.csv')
        }
    }
    result
}

# run.predictions(rpart.model.list, mold2.validate, model.name = 'rpart', resampled=F, test.data = T)
# run.predictions(glm.model.list, mold2.validate, model.name = 'glm', resampled=F, test.data = T)
# run.predictions(svmLinear.model.list, mold2.validate, model.name = 'svmLinear', resampled=F, test.data = T)
# run.predictions(nb.model.list, mold2.validate, model.name = 'nb', resampled=F, test.data = T)
# run.predictions(svmPoly.model.list, mold2.validate, model.name = 'svmPoly', resampled=F, test.data = T)
# run.predictions(svmRadial.model.list, mold2.validate, model.name = 'svmRadial', resampled=F, test.data = T)
# run.predictions(lda.model.list, mold2.validate, model.name = 'lda', resampled=F, test.data = T)
# run.predictions(nnet.model.list, mold2.validate, model.name = 'nnet', resampled=F, test.data = T)

run.predictions(top1, mold2.validate, model.name = 'top1', resampled=F, test.data = T)
run.predictions(top2, mold2.validate, model.name = 'top2', resampled=F, test.data = T)
run.predictions(top3, mold2.validate, model.name = 'top3', resampled=F, test.data = T)

run.predictions(top1.res, mold2.validate, model.name = 'top1', resampled=T, test.data = T)
run.predictions(top2.res, mold2.validate, model.name = 'top2', resampled=T, test.data = T)
run.predictions(top3.res, mold2.validate, model.name = 'top3', resampled=T, test.data = T)



#run.predictions(best.res.model.list, mold2.validate, resampled=T, test.data = T)

# run predictions on the training set ===============

run.predictions(best.non.res.model.list, mold2.train.set, resampled=F, test.data = F)
run.predictions(best.res.model.list, mold2.train.set, resampled=T, test.data = F)



