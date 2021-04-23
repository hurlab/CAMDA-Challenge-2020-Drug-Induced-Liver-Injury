





curr <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(curr))

source('../scripts/packages.R')

load('../models/faers_models.RData')
load('../models/tox_models.RData')
load('../models/mold_models.RData')

#unlist(faers_models, recursive = F)

faers_summary <- lapply(faers_models, function(x){
    summary(caret::resamples(x))
})

tox_summary <- lapply(tox_models, function(x){
    summary(caret::resamples(x))
})

mold_summary <- lapply(mold_models, function(x){
    summary(caret::resamples(x))
})

# collate metrics 
dili_names <- paste0('dili', c(1,3,5,6))

faers_metrics <- lapply(seq_along(faers_summary), function(i, ss, xx, nn){
    
    metrics <- collate.metrics(sum.list = ss, x.list = xx[[i]], dili.name = nn[[i]], write.file = F)
    
    return(metrics)
    
}, ss=faers_summary, xx=faers_models, nn=dili_names)

tox_metrics <- lapply(seq_along(tox_summary), function(i, ss, xx, nn){
    
    metrics <- collate.metrics(sum.list = ss, x.list = xx[[i]], dili.name = nn[[i]], write.file = F)
    
    return(metrics)
    
}, ss=tox_summary, xx=tox_models, nn=dili_names)

mold_metrics <- lapply(seq_along(mold_summary), function(i, ss, xx, nn){
    
    metrics <- collate.metrics(sum.list = ss, x.list = xx[[i]], dili.name = nn[[i]], write.file = F)
    
    return(metrics)
    
}, ss=mold_summary, xx=mold_models, nn=dili_names)


# select top three models
faers_top <- lapply(seq_along(faers_metrics), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=faers_metrics, mods=faers_models)

tox_top <- lapply(seq_along(tox_metrics), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=tox_metrics, mods=tox_models)

mold_top <- lapply(seq_along(mold_metrics), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=mold_metrics, mods=mold_models)

names(faers_top) <- paste0('dili', c(1,3,5,6))
names(tox_top) <- paste0('dili', c(1,3,5,6))
names(mold_top) <- paste0('dili', c(1,3,5,6))

# set up predictions and predict on the test set

# tox, validation set

# target/label/training
target.data <- read.csv('../data/targets-camda2020.csv',
                        header=T, stringsAsFactors = F, check.names = F)

tox_data <- read.csv('../data/p9-tox21-camda2020.csv', row.names = 1)
faers.data <- read.csv('../data/p2-faers-camda2020.csv',
                       header=T, stringsAsFactors=F, row.names = 'CAM_ID')
mold2.data <- read.csv('../data/p1-mold2-camda2020.csv',
                       header=T, stringsAsFactors = F, check.names = F, row.names='CAM_ID')

# tox
validation.set <- subset(target.data, Training_Validation=='') # 195 obs
tox.validate <- tox_data[rownames(tox_data) %in% validation.set$CAM_ID, ] # 190 obs
df.validation.dataset <- data.frame(sort(apply(tox.validate, 1, function(x){sum(is.na(x))}), decreasing = T))
names(df.validation.dataset) <- 'No.of.NAs'
tox.validate.names <- row.names(subset(df.validation.dataset, No.of.NAs==0))
tox.validate <- tox.validate[row.names(tox.validate) %in% tox.validate.names, ]

# faers
faers.validate <- faers.data[rownames(faers.data) %in% validation.set$CAM_ID, ]
faers.validate <- faers.validate %>% 
    mutate(ratio_dili_all=dili_gender_all/all_gender_all, 
           male_rate=(dili_gender_all*dili_gender_male_percentage)/(all_gender_male_percentage*all_gender_all),
           female_rate=(dili_gender_all*dili_gender_female_percentage)/(all_gender_female_percentage*all_gender_all)) 

# mold
# validation set
#validation.set <- subset(target.data, Training_Validation=='') # 195 obs
mold2.validate <- mold2.data[rownames(mold2.data) %in% validation.set$CAM_ID, ] # 190 obs

df.validation.dataset <- data.frame(sort(apply(mold2.validate, 1, function(x){sum(is.na(x))}), decreasing = T))
names(df.validation.dataset) <- 'No.of.NAs'
mold2.validate.names <- row.names(subset(df.validation.dataset, No.of.NAs==0))
mold2.validate <- mold2.validate[row.names(mold2.validate) %in% mold2.validate.names, ]


# predict

faers_top_models <- lapply(faers_top, function(x){
    x[['top_models']]
})

tox_top_models <- lapply(tox_top, function(x){
    x[['top_models']]
})

mold_top_models <- lapply(mold_top, function(x){
    x[['top_models']]
})

run.predictions <- function(model.list, new.data){
    
    output <- list()
    
    out <- lapply(seq_along(model.list), function(i, mm, nd){
        
        d_type <- mm[[i]]
        
        ret <- lapply(d_type, function(mod){
            predicted.values <- as.vector(predict(mod, newdata=nd))
            predicted.values[predicted.values == 'Positive'] <- 0
            predicted.values[predicted.values == 'Negative'] <- 1
            
            return(as.numeric(predicted.values))
            
            #data.frame(cbind(CAM_ID=row.names(new.data), as.data.frame(do.call(cbind, output))))
        })
        
    }, mm=model.list, nd=new.data)
    
    
    names(out) <- names(model.list)
    
    result <- data.frame(cbind(CAM_ID=row.names(new.data), as.data.frame(do.call(data.frame, out))))
    
    result
    
    
}

tox_predictions <- run.predictions(tox_top_models, tox.validate)
faers_predictions <- run.predictions(faers_top_models, faers.validate)
mold_predictions <- run.predictions(mold_top_models, mold2.validate)


# =============
# run.predictions <- function(model.list, new.data){
#     
#     
#     output <- list()
#     out <- lapply(seq_along(model.list), function(i, mm, nd){
#         
#         d_type <- mm[[i]]
#         
#         for(j in 1:length(d_type)){
#             
#             d_type_name <- names(d_type)[j]
#             predicted.values <- as.vector(predict(d_type[[j]], newdata=nd))
#             predicted.values[predicted.values == 'Positive'] <- 0
#             predicted.values[predicted.values == 'Negative'] <- 1
#             #output[[paste(d_type_name, sep='')]] <- as.numeric(predicted.values)
#             
#             as.numeric(predicted.values)
#         }
#         
#     }, mm=model.list, nd=new.data)
#     
#     
#     #result <- data.frame(cbind(CAM_ID=row.names(new.data), as.data.frame(do.call(cbind, output))))
#     
#     print(out)
#     
#     # for (i in 1:length(model.list)){
#     #     predicted.values <- as.vector(predict(model.list[[i]], newdata=new.data))
#     #     predicted.values[predicted.values == 'Positive'] <- 0
#     #     predicted.values[predicted.values == 'Negative'] <- 1
#     #     output[[paste(names(model.list)[i], sep='')]] <- as.numeric(predicted.values)
#     # }
# 
#     # print(output)
#     # str(do.call(cbind, output))
#     #result <- data.frame(cbind(CAM_ID=row.names(new.data), as.data.frame(do.call(cbind, output))))
#     #row.names(result) <- row.names(new.data)
#     #names(result)[1] <- 'CAM_ID'
# 
#     # if (resampled==F) {
#     #     if (test.data==T) {
#     #         write.csv(result, file=paste('../output_files/p9.', model.name, '-predictions-camda2020-UND.csv', sep=''), row.names=F)
#     #     } else {
#     #         write.csv(result, file='../output_files/train_tox2_predictions.csv')
#     #     }
#     # } else if (resampled==T) {
#     #     if (test.data==T) {
#     #         write.csv(result, file=paste('../output_files/p9.resampled-', model.name, '-predictions-camda2020-UND.csv', sep=''), row.names=F)
#     #     } else {
#     #         write.csv(result, file='../output_files/train_tox2_resampled_predictions.csv')
#     #     }
#     # }
#     #result
# }















select.top.three(faers_metrics[[1]], faers_models$faers_dili1)

faers_metrics[[1]]


names(faers_models$faers_dili1)


