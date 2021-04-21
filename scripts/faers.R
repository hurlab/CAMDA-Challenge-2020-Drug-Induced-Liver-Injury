
# Jun 27 2020
# This script was written to run classification models on faers data...
#...for CAMDA 2020

# Temi


library(rstudioapi)
currPath <- getActiveDocumentContext()$path
setwd(dirname(currPath))

# load libraries
source('../scripts/packages.R')

# read data
faers.data <- read.csv('../data/p2-faers-camda2020.csv',
                       header=T, stringsAsFactors=F, row.names = 'CAM_ID')
# target/label/training
target.data <- read.csv('../data/targets-camda2020.csv',
                        header=T, stringsAsFactors = F, check.names = F)


# faers.data has all the data
# select the training set based on what is available in target.data
training.set <- subset(target.data, Training_Validation=='Training Set')
faers.train <- faers.data[rownames(faers.data) %in% training.set$CAM_ID, ] 


complete.columns <- as.data.frame(apply(faers.train, 1, function(x){
    sum(is.na(x))
}))
names(complete.columns)[1] <- 'No.of.NAs'
complete.col.names <- rownames(subset(complete.columns, No.of.NAs==0))


table(complete.columns$No.of.NAs)

# processing this dataset
faers.train <- subset(faers.train, all_gender_all != 0)

faers.train <- faers.train %>% 
    mutate(ratio_dili_all=dili_gender_all/all_gender_all, 
           male_rate=(dili_gender_all*dili_gender_male_percentage)/(all_gender_male_percentage*all_gender_all),
           female_rate=(dili_gender_all*dili_gender_female_percentage)/(all_gender_female_percentage*all_gender_all))  

# %>%
#     dplyr::select(ratio_dili_all, male_rate, female_rate, dili_age_neonate_percentage,
#                   dili_age_infant_percentage, dili_age_child_percentage,
#                   dili_age_adolescent_percentage, dili_age_adult_percentage,
#                   dili_age_elderly_percentage, all_age_neonate_percentage )

faers.train <- faers.train[-which(is.nan(faers.train$male_rate)), ]
cor.faers <- cor(faers.train)
highly.cor <- findCorrelation(cor.faers, cutoff=0.80)
faers.train <- faers.train[, -highly.cor]

# dili1 to dili6 targets dataset
target.labels.1 <- training.set[training.set$CAM_ID %in% row.names(faers.train), ]$DILI1
target.labels.1 <- as.factor(mapvalues(target.labels.1, 
                                       from=c(0, 1), to=c('Positive', 'Negative')))
target.labels.3 <- training.set[training.set$CAM_ID %in% row.names(faers.train), ]$DILI3
target.labels.3 <- as.factor(mapvalues(target.labels.3, 
                                       from=c(0, 1), to=c('Positive', 'Negative')))
target.labels.5 <- training.set[training.set$CAM_ID %in% row.names(faers.train), ]$DILI5
target.labels.5 <- as.factor(mapvalues(target.labels.5, 
                                       from=c(0, 1), to=c('Positive', 'Negative')))
target.labels.6 <- training.set[training.set$CAM_ID %in% row.names(faers.train), ]$DILI6
target.labels.6 <- as.factor(mapvalues(target.labels.6, 
                                       from=c(0, 1), to=c('Positive', 'Negative')))

# upsampled dataset
up.faers1.train <- ovun.sample(Class~., data=data.frame(cbind(faers.train, Class=target.labels.1)),
                               seed=1, N=max(table(target.labels.1))*2)$data
up.faers3.train <- ovun.sample(Class~., data=data.frame(cbind(faers.train, Class=target.labels.3)),
                               seed=1, N=max(table(target.labels.3))*2)$data
up.faers5.train <- ovun.sample(Class~., data=data.frame(cbind(faers.train, Class=target.labels.5)),
                               seed=1, N=max(table(target.labels.5))*2)$data
up.faers6.train <- ovun.sample(Class~., data=data.frame(cbind(faers.train, Class=target.labels.6)),
                               seed=1, N=max(table(target.labels.6))*2)$data


# some columns have a lot of zeros. Take them out, please.
# apply(up.faers1.train, 2, function(x){
#     sum(x==0)
# })

# cluster
cl <- makeCluster(100)
registerDoParallel(cl)

# feature selection ============
rfe.ctrl <- rfeControl(method = "repeatedcv", repeats = 5, number=100, 
                       functions = rfFuncs)
set.seed(1)
rfe.faers1.model <- rfe(x=faers.train, y=target.labels.1, 
                        sizes = c(2:20), rfeControl = rfe.ctrl, preProcess=c('center', 'scale'))
#rfe.faers1.model
#ggplot(rfe.faers1.model) + theme_bw()
rfe1 <- rfe.faers1.model$optVariables

set.seed(1)
rfe.faers3.model <- rfe(x=faers.train, y=target.labels.3, 
                        sizes = c(2:20), rfeControl = rfe.ctrl, preProcess=c('center', 'scale'))
#rfe.faers3.model
#ggplot(rfe.faers3.model) + theme_bw()
rfe3 <- rfe.faers3.model$optVariables

set.seed(1)
rfe.faers5.model <- rfe(x=faers.train, y=target.labels.5, 
                        sizes = c(2:20), rfeControl = rfe.ctrl, preProcess=c('center', 'scale'))
# rfe.faers5.model
# ggplot(rfe.faers5.model) + theme_bw()
rfe5 <- rfe.faers5.model$optVariables

set.seed(1)
rfe.faers6.model <- rfe(x=faers.train, y=target.labels.6, 
                        sizes = c(2:20), rfeControl = rfe.ctrl, preProcess=c('center', 'scale'))
# rfe.faers6.model
# ggplot(rfe.faers6.model) + theme_bw()
rfe6 <- rfe.faers6.model$optVariables


sbf.ctrl <- sbfControl(functions = rfSBF, method = "repeatedcv", 
                       repeats = 5, number=100)
set.seed(1)
sbf.faers1.model <- sbf(x=faers.train, y=target.labels.1, 
                        sbfControl = sbf.ctrl, preProcess=c('center', 'scale'))
#sbf.faers1.model
sbf1 <- sbf.faers1.model$optVariables

set.seed(1)
sbf.faers3.model <- sbf(x=faers.train, y=target.labels.3, 
                        sbfControl = sbf.ctrl, preProcess=c('center', 'scale'))
#sbf.faers3.model
sbf3 <- sbf.faers3.model$optVariables

set.seed(1)
sbf.faers5.model <- sbf(x=faers.train, y=target.labels.5, 
                        sbfControl = sbf.ctrl, preProcess=c('center', 'scale'))
#sbf.faers5.model
sbf5 <- sbf.faers5.model$optVariables

set.seed(1)
sbf.faers6.model <- sbf(x=faers.train, y=target.labels.6, 
                        sbfControl = sbf.ctrl, preProcess=c('center', 'scale'))
#sbf.faers6.model
sbf6 <- sbf.faers6.model$optVariables

## data mining ======
# define a general control
general.ctrl <- trainControl(method='repeatedcv',
                             repeats=100, number=5,
                             savePredictions=T,
                             classProbs = T,
                             summaryFunction = twoClassSummary,
                             allowParallel = T)

# ============= DILI1 ==================================
# glm ======
set.seed(1)
glm1.model <- train(x=faers.train, y=target.labels.1,
                    maxit=200, method='glm',
                    family=binomial(link='logit'),
                    trControl = general.ctrl, metric='ROC')


set.seed(1)
glm1.up.model <- train(x=up.faers1.train[, -ncol(up.faers1.train)][, rfe1], 
                       y=as.factor(up.faers1.train$Class),
                       maxit=200, method='glm',
                       family=binomial(link='logit'),
                       trControl = general.ctrl, metric='ROC')


# rpart ================
set.seed(1)
rpart1.model <- train(x=faers.train, y=target.labels.1, 
                      method='rpart', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
rpart1.up.model <- train(x=up.faers1.train[, -ncol(up.faers1.train)][, rfe1], 
                         y=as.factor(up.faers1.train$Class), 
                         method='rpart', trControl = general.ctrl, metric='ROC',
                         preProcess = c('scale', 'center'), tuneLength = 20)


# svmPoly ================
set.seed(1)
svmPoly1.model <- train(x=faers.train, y=target.labels.1, 
                        method='svmPoly', trControl = general.ctrl, metric='ROC',
                        preProcess = c('scale', 'center'))


set.seed(1)
svmPoly1.up.model <- train(x=up.faers1.train[, -ncol(up.faers1.train)][, rfe1], 
                           y=as.factor(up.faers1.train$Class), 
                           method='svmPoly', trControl = general.ctrl, metric='ROC',
                           preProcess = c('scale', 'center'))


# nnet ================
set.seed(1)
nnet1.model <- train(x=faers.train, y=target.labels.1, 
                     method='nnet', trControl = general.ctrl, metric='ROC',
                     preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
nnet1.up.model <- train(x=up.faers1.train[, -ncol(up.faers1.train)][, rfe1], 
                        y=as.factor(up.faers1.train$Class), 
                        method='nnet', trControl = general.ctrl, metric='ROC',
                        preProcess = c('scale', 'center'), tuneLength = 20)


# lda ================
set.seed(1)
lda1.model <- train(x=faers.train, y=target.labels.1, 
                    method='lda', trControl = general.ctrl, metric='ROC',
                    preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
lda1.up.model <- train(x=up.faers1.train[, -ncol(up.faers1.train)][, rfe1], 
                       y=as.factor(up.faers1.train$Class), 
                       method='lda', trControl = general.ctrl, metric='ROC',
                       preProcess = c('scale', 'center'), tuneLength = 20)


# rf ================
set.seed(1)
rf1.model <- train(x=faers.train, y=target.labels.1, 
                   method='rf', trControl = general.ctrl, metric='ROC',
                   preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
rf1.up.model <- train(x=up.faers1.train[, -ncol(up.faers1.train)][, rfe1], 
                      y=as.factor(up.faers1.train$Class), 
                      method='rf', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


# nb ================
set.seed(1)
nb1.model <- train(x=faers.train, y=target.labels.1, 
                   method='naive_bayes', trControl = general.ctrl, metric='ROC',
                   preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
nb1.up.model <- train(x=up.faers1.train[, -ncol(up.faers1.train)][, rfe1], 
                      y=as.factor(up.faers1.train$Class), 
                      method='naive_bayes', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)



# ============= DILI3 ==================================
# glm ======
set.seed(1)
glm3.model <- train(x=faers.train, y=target.labels.3,
                    maxit=200, method='glm',
                    family=binomial(link='logit'),
                    trControl = general.ctrl, metric='ROC')


set.seed(1)
glm3.up.model <- train(x=up.faers3.train[, -ncol(up.faers3.train)][, rfe3], 
                       y=as.factor(up.faers3.train$Class),
                       maxit=200, method='glm',
                       family=binomial(link='logit'),
                       trControl = general.ctrl, metric='ROC')


# rpart ================
set.seed(1)
rpart3.model <- train(x=faers.train, y=target.labels.3, 
                      method='rpart', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
rpart3.up.model <- train(x=up.faers3.train[, -ncol(up.faers3.train)][, rfe3], 
                         y=as.factor(up.faers3.train$Class), 
                         method='rpart', trControl = general.ctrl, metric='ROC',
                         preProcess = c('scale', 'center'), tuneLength = 20)


# svmPoly ================
set.seed(1)
svmPoly3.model <- train(x=faers.train, y=target.labels.3, 
                        method='svmPoly', trControl = general.ctrl, metric='ROC',
                        preProcess = c('scale', 'center'))


set.seed(1)
svmPoly3.up.model <- train(x=up.faers3.train[, -ncol(up.faers3.train)][, rfe3], 
                           y=as.factor(up.faers3.train$Class), 
                           method='svmPoly', trControl = general.ctrl, metric='ROC',
                           preProcess = c('scale', 'center'))


# nnet ================
set.seed(1)
nnet3.model <- train(x=faers.train, y=target.labels.3, 
                     method='nnet', trControl = general.ctrl, metric='ROC',
                     preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
nnet3.up.model <- train(x=up.faers3.train[, -ncol(up.faers3.train)][, rfe3], 
                        y=as.factor(up.faers3.train$Class), 
                        method='nnet', trControl = general.ctrl, metric='ROC',
                        preProcess = c('scale', 'center'), tuneLength = 20)


# lda ================
set.seed(1)
lda3.model <- train(x=faers.train, y=target.labels.3, 
                    method='lda', trControl = general.ctrl, metric='ROC',
                    preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
lda3.up.model <- train(x=up.faers3.train[, -ncol(up.faers3.train)][, rfe3], 
                       y=as.factor(up.faers3.train$Class), 
                       method='lda', trControl = general.ctrl, metric='ROC',
                       preProcess = c('scale', 'center'), tuneLength = 20)


# rf ================
set.seed(1)
rf3.model <- train(x=faers.train, y=target.labels.3, 
                   method='rf', trControl = general.ctrl, metric='ROC',
                   preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
rf3.up.model <- train(x=up.faers3.train[, -ncol(up.faers3.train)][, rfe3], 
                      y=as.factor(up.faers3.train$Class), 
                      method='rf', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


# nb ================
set.seed(1)
nb3.model <- train(x=faers.train, y=target.labels.3, 
                   method='naive_bayes', trControl = general.ctrl, metric='ROC',
                   preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
nb3.up.model <- train(x=up.faers3.train[, -ncol(up.faers3.train)][, rfe3], 
                      y=as.factor(up.faers3.train$Class), 
                      method='naive_bayes', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)



# ============= DILI5 ==================================
# glm ======
set.seed(1)
glm5.model <- train(x=faers.train, y=target.labels.5,
                    maxit=200, method='glm',
                    family=binomial(link='logit'),
                    trControl = general.ctrl, metric='ROC')


set.seed(1)
glm5.up.model <- train(x=up.faers5.train[, -ncol(up.faers5.train)][, rfe5], 
                       y=as.factor(up.faers5.train$Class),
                       maxit=200, method='glm',
                       family=binomial(link='logit'),
                       trControl = general.ctrl, metric='ROC')


# rpart ================
set.seed(1)
rpart5.model <- train(x=faers.train, y=target.labels.5, 
                      method='rpart', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
rpart5.up.model <- train(x=up.faers5.train[, -ncol(up.faers5.train)][, rfe5], 
                         y=as.factor(up.faers5.train$Class), 
                         method='rpart', trControl = general.ctrl, metric='ROC',
                         preProcess = c('scale', 'center'), tuneLength = 20)


# svmPoly ================
set.seed(1)
svmPoly5.model <- train(x=faers.train, y=target.labels.5, 
                        method='svmPoly', trControl = general.ctrl, metric='ROC',
                        preProcess = c('scale', 'center'))


set.seed(1)
svmPoly5.up.model <- train(x=up.faers5.train[, -ncol(up.faers5.train)][, rfe5], 
                           y=as.factor(up.faers5.train$Class), 
                           method='svmPoly', trControl = general.ctrl, metric='ROC',
                           preProcess = c('scale', 'center'))


# nnet ================
set.seed(1)
nnet5.model <- train(x=faers.train, y=target.labels.5, 
                     method='nnet', trControl = general.ctrl, metric='ROC',
                     preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
nnet5.up.model <- train(x=up.faers5.train[, -ncol(up.faers5.train)][, rfe5], 
                        y=as.factor(up.faers5.train$Class), 
                        method='nnet', trControl = general.ctrl, metric='ROC',
                        preProcess = c('scale', 'center'), tuneLength = 20)


# lda ================
set.seed(1)
lda5.model <- train(x=faers.train, y=target.labels.5, 
                    method='lda', trControl = general.ctrl, metric='ROC',
                    preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
lda5.up.model <- train(x=up.faers5.train[, -ncol(up.faers5.train)][, rfe5], 
                       y=as.factor(up.faers5.train$Class), 
                       method='lda', trControl = general.ctrl, metric='ROC',
                       preProcess = c('scale', 'center'), tuneLength = 20)


# rf ================
set.seed(1)
rf5.model <- train(x=faers.train, y=target.labels.5, 
                   method='rf', trControl = general.ctrl, metric='ROC',
                   preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
rf5.up.model <- train(x=up.faers5.train[, -ncol(up.faers5.train)][, rfe5], 
                      y=as.factor(up.faers5.train$Class), 
                      method='rf', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


# nb ================
set.seed(1)
nb5.model <- train(x=faers.train, y=target.labels.5, 
                   method='naive_bayes', trControl = general.ctrl, metric='ROC',
                   preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
nb5.up.model <- train(x=up.faers5.train[, -ncol(up.faers5.train)][, rfe5], 
                      y=as.factor(up.faers5.train$Class), 
                      method='naive_bayes', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)



# ============= DILI6 ==================================
# glm ======
set.seed(1)
glm6.model <- train(x=faers.train, y=target.labels.6,
                    maxit=200, method='glm',
                    family=binomial(link='logit'),
                    trControl = general.ctrl, metric='ROC')


set.seed(1)
glm6.up.model <- train(x=up.faers6.train[, -ncol(up.faers6.train)][, rfe6], 
                       y=as.factor(up.faers6.train$Class),
                       maxit=200, method='glm',
                       family=binomial(link='logit'),
                       trControl = general.ctrl, metric='ROC')


# rpart ================
set.seed(1)
rpart6.model <- train(x=faers.train, y=target.labels.6, 
                      method='rpart', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


set.seed(1)
rpart6.up.model <- train(x=up.faers6.train[, -ncol(up.faers6.train)][, rfe6], 
                         y=as.factor(up.faers6.train$Class), 
                         method='rpart', trControl = general.ctrl, metric='ROC',
                         preProcess = c('scale', 'center'), tuneLength = 20)


# svmPoly ================
set.seed(1)
svmPoly6.model <- train(x=faers.train, y=target.labels.6, 
                        method='svmPoly', trControl = general.ctrl, metric='ROC',
                        preProcess = c('scale', 'center'))


set.seed(1)
svmPoly6.up.model <- train(x=up.faers6.train[, -ncol(up.faers6.train)][, rfe6], 
                           y=as.factor(up.faers6.train$Class), 
                           method='svmPoly', trControl = general.ctrl, metric='ROC',
                           preProcess = c('scale', 'center'))


# nnet ================
set.seed(1)
nnet6.model <- train(x=faers.train, y=target.labels.6, 
                     method='nnet', trControl = general.ctrl, metric='ROC',
                     preProcess = c('scale', 'center'), tuneLength = 20)

set.seed(1)
nnet6.up.model <- train(x=up.faers6.train[, -ncol(up.faers6.train)][, rfe6], 
                        y=as.factor(up.faers6.train$Class), 
                        method='nnet', trControl = general.ctrl, metric='ROC',
                        preProcess = c('scale', 'center'), tuneLength = 20)

# lda ================
set.seed(1)
lda6.model <- train(x=faers.train, y=target.labels.6, 
                    method='lda', trControl = general.ctrl, metric='ROC',
                    preProcess = c('scale', 'center'), tuneLength = 20)

set.seed(1)
lda6.up.model <- train(x=up.faers6.train[, -ncol(up.faers6.train)][, rfe6], 
                       y=as.factor(up.faers6.train$Class), 
                       method='lda', trControl = general.ctrl, metric='ROC',
                       preProcess = c('scale', 'center'), tuneLength = 20)

# rf ================
set.seed(1)
rf6.model <- train(x=faers.train, y=target.labels.6, 
                   method='rf', trControl = general.ctrl, metric='ROC',
                   preProcess = c('scale', 'center'), tuneLength = 20)

set.seed(1)
rf6.up.model <- train(x=up.faers6.train[, -ncol(up.faers6.train)][, rfe6], 
                      y=as.factor(up.faers6.train$Class), 
                      method='rf', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


# nb ================
set.seed(1)
nb6.model <- train(x=faers.train, y=target.labels.6, 
                   method='naive_bayes', trControl = general.ctrl, metric='ROC',
                   preProcess = c('scale', 'center'), tuneLength = 20)
set.seed(1)
nb6.up.model <- train(x=up.faers6.train[, -ncol(up.faers6.train)][, rfe6], 
                      y=as.factor(up.faers6.train$Class), 
                      method='naive_bayes', trControl = general.ctrl, metric='ROC',
                      preProcess = c('scale', 'center'), tuneLength = 20)


stopCluster(cl)

# Evaluation==================
# dili1==============
glm1.res <- resamples(list(glm1=glm1.model,
                           glm1.up=glm1.up.model))
summary(glm1.res)

rpart1.res <- resamples(list(rpart1=rpart1.model,
                             rpart1.up=rpart1.up.model))
summary(rpart1.res)

svmPoly1.res <- resamples(list(svmPoly1=svmPoly1.model,
                               svmPoly1.up=svmPoly1.up.model))
summary(svmPoly1.res)

nnet1.res <- resamples(list(nnet1=nnet1.model,
                            nnet1.up=nnet1.up.model))
summary(nnet1.res)

lda1.res <- resamples(list(lda1=lda1.model,
                           lda1.up=lda1.up.model))
summary(lda1.res)

rf1.res <- resamples(list(rf1=rf1.model,
                          rf1.up=rf1.up.model))
summary(rf1.res)

nb1.res <- resamples(list(nb1=nb1.model,
                          nb1.up=nb1.up.model))
summary(nb1.res)

# dili3==========
glm3.res <- resamples(list(glm3=glm3.model,
                           glm3.up=glm3.up.model))
summary(glm3.res)

rpart3.res <- resamples(list(rpart3=rpart3.model,
                             rpart3.up=rpart3.up.model))
summary(rpart3.res)

svmPoly3.res <- resamples(list(svmPoly3=svmPoly3.model,
                               svmPoly3.up=svmPoly3.up.model))
summary(svmPoly3.res)

nnet3.res <- resamples(list(nnet3=nnet3.model,
                            nnet3.up=nnet3.up.model))
summary(nnet3.res)

lda3.res <- resamples(list(lda3=lda3.model,
                           lda3.up=lda3.up.model))
summary(lda3.res)

rf3.res <- resamples(list(rf3=rf3.model,
                          rf3.up=rf3.up.model))
summary(rf3.res)

nb3.res <- resamples(list(nb3=nb3.model,
                          nb3.up=nb3.up.model))
summary(nb3.res)

# dili5==============
glm5.res <- resamples(list(glm5=glm5.model,
                           glm5.up=glm5.up.model))
summary(glm5.res)

rpart5.res <- resamples(list(rpart5=rpart5.model,
                             rpart5.up=rpart5.up.model))
summary(rpart5.res)

svmPoly5.res <- resamples(list(svmPoly5=svmPoly5.model,
                               svmPoly5.up=svmPoly5.up.model))
summary(svmPoly5.res)

nnet5.res <- resamples(list(nnet5=nnet5.model,
                            nnet5.up=nnet5.up.model))
summary(nnet5.res)

lda5.res <- resamples(list(lda5=lda5.model,
                           lda5.up=lda5.up.model))
summary(lda5.res)

rf5.res <- resamples(list(rf5=rf5.model,
                          rf5.up=rf5.up.model))
summary(rf5.res)

nb5.res <- resamples(list(nb5=nb5.model,
                          nb5.up=nb5.up.model))
summary(nb5.res)

# dili6==================
glm6.res <- resamples(list(glm6=glm6.model,
                           glm6.up=glm6.up.model))
summary(glm6.res)

rpart6.res <- resamples(list(rpart6=rpart6.model,
                             rpart6.up=rpart6.up.model))
summary(rpart6.res)

svmPoly6.res <- resamples(list(svmPoly6=svmPoly6.model,
                               svmPoly6.up=svmPoly6.up.model))
summary(svmPoly6.res)

nnet6.res <- resamples(list(nnet6=nnet6.model,
                            nnet6.up=nnet6.up.model))
summary(nnet6.res)

lda6.res <- resamples(list(lda6=lda6.model,
                           lda6.up=lda6.up.model))
summary(lda6.res)

rf6.res <- resamples(list(rf6=rf6.model,
                          rf6.up=rf6.up.model))
summary(rf6.res)

nb6.res <- resamples(list(nb6=nb6.model,
                          nb6.up=nb6.up.model))
summary(nb6.res)


summary.list.1 <- list(summary(glm1.res),
                       summary(rpart1.res),
                       summary(nnet1.res),
                       summary(svmPoly1.res),
                       summary(lda1.res),
                       summary(nb1.res),
                       summary(rf1.res))

names(summary.list.1) <- paste0(c('glm', 'rpart', 'nnet', 'svm', 'lda', 'nb', 'rf'), 1)

filter(as.data.frame(summary.list.1$glm1$values), `glm1.up~ROC` == max(`glm1.up~ROC`))

model.list.1 <- list(glm1=glm1.model, glm1.up=glm1.up.model, 
                     rf1=rf1.model, rf1.up=rf1.up.model, 
                     svmPoly1=svmPoly1.model, svmPoly1.up=svmPoly1.up.model, 
                     lda1=lda1.model, lda1.up=lda1.up.model, 
                     rpart1=rpart1.model, rpart1.up=rpart1.up.model, 
                     nnet1=nnet1.model, nnet1.up=nnet1.up.model,
                     nb1=nb1.model, nb1.up=nb1.up.model)

summary.list.3 <- list(summary(glm3.res),
                       summary(rpart3.res),
                       summary(nnet3.res),
                       summary(svmPoly3.res),
                       summary(lda3.res),
                       summary(nb3.res),
                       summary(rf3.res))

names(summary.list.3) <- paste0(c('glm', 'rpart', 'nnet', 'svm', 'lda', 'nb', 'rf'), 3)

model.list.3 <- list(glm3=glm3.model, glm3.up=glm3.up.model, 
                     rf3=rf3.model, rf3.up=rf3.up.model, 
                     svmPoly3=svmPoly3.model, svmPoly3.up=svmPoly3.up.model, 
                     lda3=lda3.model, lda3.up=lda3.up.model, 
                     rpart3=rpart3.model, rpart3.up=rpart3.up.model, 
                     nnet3=nnet3.model, nnet3.up=nnet3.up.model,
                     nb3=nb3.model, nb3.up=nb3.up.model)

summary.list.5 <- list(summary(glm5.res),
                       summary(rpart5.res),
                       summary(nnet5.res),
                       summary(svmPoly5.res),
                       summary(lda5.res),
                       summary(nb5.res),
                       summary(rf5.res))

names(summary.list.5) <- paste0(c('glm', 'rpart', 'nnet', 'svm', 'lda', 'nb', 'rf'), 5)

model.list.5 <- list(glm5=glm5.model, glm5.up=glm5.up.model, 
                     rf5=rf5.model, rf5.up=rf5.up.model, 
                     svmPoly5=svmPoly5.model, svmPoly5.up=svmPoly5.up.model, 
                     lda5=lda5.model, lda5.up=lda5.up.model, 
                     rpart5=rpart5.model, rpart5.up=rpart5.up.model, 
                     nnet5=nnet5.model, nnet5.up=nnet5.up.model,
                     nb5=nb5.model, nb5.up=nb5.up.model)

summary.list.6 <- list(summary(glm6.res),
                       summary(rpart6.res),
                       summary(nnet6.res),
                       summary(svmPoly6.res),
                       summary(lda6.res),
                       summary(nb6.res),
                       summary(rf6.res))

names(summary.list.6) <- paste0(c('glm', 'rpart', 'nnet', 'svm', 'lda', 'nb', 'rf'), 6)

model.list.6 <- list(glm6=glm6.model, glm6.up=glm6.up.model, 
                     rf6=rf6.model, rf6.up=rf6.up.model, 
                     svmPoly6=svmPoly6.model, svmPoly6.up=svmPoly6.up.model, 
                     lda6=lda6.model, lda6.up=lda6.up.model, 
                     rpart6=rpart6.model, rpart6.up=rpart6.up.model, 
                     nnet6=nnet6.model, nnet6.up=nnet6.up.model,
                     nb6=nb6.model, nb6.up=nb6.up.model)

# a function to calculate the matthew's correlation co-efficient ====
mcc <- function(x.model){
    cm <- confusionMatrix(x.model)[['table']]
    true.pos <- cm[1,1]
    false.pos <- cm[1,2]
    false.neg <- cm[2,1]
    true.neg <- cm[2,2]
    
    above <- true.pos*true.neg - false.pos*false.neg
    below <- as.double(true.pos+false.pos)*as.double(true.pos+false.neg)*as.double(true.neg+false.pos)*as.double(true.neg+false.neg)
    
    mcc <- above/sqrt(below)
    mcc
}

collate.metrics <- function(sum.list, x.list, dili.name='', write.file=F){
    
    # """ This function collates metrics for a list of models """
    
    out.roc <- list()
    out.sens <- list()
    out.spec <- list()
    out.mcc <- list()
    
    for(i in 1:length(sum.list)){
        roc <- as.data.frame(sum.list[[i]]$statistics$ROC[, 'Mean'])
        roc$dataset <- row.names(roc)
        row.names(roc) <- NULL
        out.roc[[i]] <- roc
        
        sens <- as.data.frame(sum.list[[i]]$statistics$Sens[, 'Mean'])
        sens$dataset <- row.names(sens)
        row.names(sens) <- NULL
        out.sens[[i]] <- sens
        
        spec <- as.data.frame(sum.list[[i]]$statistics$Spec[, 'Mean'])
        spec$dataset <- row.names(spec)
        row.names(spec) <- NULL
        out.spec[[i]] <- spec
    }
    
    for(i in 1:length(x.list)){
        mcc.i <- data.frame(dataset=names(x.list)[i], mcc=mcc(x.list[[i]]))
        out.mcc[[i]] <- mcc.i
    }
    
    roc <- do.call(rbind, out.roc)
    names(roc)[1] <- 'ROC'
    sens <- do.call(rbind, out.sens)
    names(sens)[1] <- 'Sens'
    spec <- do.call(rbind, out.spec)
    names(spec)[1] <- 'Spec'
    mcc <- do.call(rbind, out.mcc)
    
    first <- merge(roc, sens, by='dataset')
    mid <- merge(first, spec, by='dataset')
    last <- merge(mid, mcc, by='dataset')
    last <- last %>% arrange(ROC)
    
    if(write.file==T){
        write.csv(last, file = paste('../output_files/p2.', dili.name, '.csv', sep=''), row.names=F)
    } else {
        return(last)
    }
    
    write.csv(last, file = paste('../output_files/p2.', dili.name, '.csv', sep=''), row.names=F)
}

collate.metrics(sum.list=summary.list.1, x.list=model.list.1, dili.name = 'dili1', dt='p2', write.file = F)
collate.metrics(sum.list=summary.list.3, x.list=model.list.3, dili.name = 'dili3')
collate.metrics(sum.list=summary.list.5, x.list=model.list.5, dili.name = 'dili5')
collate.metrics(sum.list=summary.list.6, x.list=model.list.6, dili.name = 'dili6')


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
            write.csv(result, file=paste('../output_files/p2.', model.name, '-predictions-camda2020-UND.csv', sep=''), row.names = F)
        } else {
            write.csv(result, file='../output_files/faers_results_train_faers_predictions.csv')
        }
    } else if (resampled==T) {
        if (test.data==T) {
            write.csv(result, file=paste('../output_files/p1.resampled-', model.name, '-predictions-camda2020-UND.csv', sep=''), row.names = F)
        } else {
            write.csv(result, file='../output_files/faers_results_train_faers_resampled_predictions.csv')
        }
    }
    result
}

