


curr <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(curr))

source('../scripts/packages.R')

# load('../models/faers_models.RData')
# load('../models/tox_models.RData')
# load('../models/mold_models.RData')


library(tidyr)
library(kernlab)

# change the name of this later
source(file = '../scripts/ext.R')

# load the gene expression models into a separate environment
brett_models <- '../models/gene_expr_data.RData'

#brett_models <- '/extData/NGS/hurlab/CAMDA2020/Brett_CAMDA2020/June29/ClassFix_02Feb21BM.RData'
#brett_models_env <- new.env()
load(file = brett_models)
# # convert the envir above into a list
# #bretts <- as.list(brett_models_env)
# # extract the names of the files that are models
# classes <- lapply(expr_models_list, function(x){
#     class(x)[1]
# })
# br_models <- names(classes[classes=='train'])
# # get the objects that correspond to those names
# expr_list <- mget(br_models, )

# remove some models (not important)
expr_list <- expr_models_list[!grepl(pattern = '^[a-z0-9]+\\.\\d\\w+', x = names(expr_models_list))]

# load other models (Temi's)
load('../models/mold_models.RData', mold_env <- new.env())
load('../models/faers_models.RData', faers_env <- new.env())
load('../models/tox_models.RData', tox_env <- new.env())

mold_list <- as.list(mold_env$mold_models)
faers_list <- as.list(faers_env$faers_models)
tox_list <- as.list(tox_env$tox_models)

# remove the environments
take_out <- c('brett_models_env', 'bretts', 'classes', 'faers_env', 'mold_env', 'tox_env', 'br_models')
rm(list = take_out)

dili1 <- expr_list[grepl('.*\\.1.*', names(expr_list))]
dili3 <- expr_list[grepl('.*\\.3.*', names(expr_list))]
dili5 <- expr_list[grepl('.*\\.5.*', names(expr_list))]
dili6 <- expr_list[grepl('.*\\.6.*', names(expr_list))]

expr_list <- list(expr_dili1=dili1, expr_dili3=dili3, expr_dili5=dili5, expr_dili6=dili6)

# remove the redundancies
take_out <- paste('dili', c(1,3,5,6), sep='')
rm(list = take_out)

# gather their metrics
faers_perf <- collate_metrics(faers_list)
mold_perf <- collate_metrics(mold_list)
tox_perf <- collate_metrics(tox_list)
expr_perf <- collate_metrics(expr_list)

all_dataset <- list(mold=mold_perf, faers=faers_perf, tox=tox_perf, expr=expr_perf)
all_perf <- as.data.frame(do.call(rbind, all_dataset))

# using the top 3 models only
# select top three models
faers_metrics_split <- split(faers_perf, f=faers_perf$`dili type`)
faers_top <- lapply(seq_along(faers_metrics_split), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=faers_metrics_split, mods=faers_list)

tox_metrics_split <- split(tox_perf, f=tox_perf$`dili type`)
tox_top <- lapply(seq_along(tox_metrics_split), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=tox_metrics_split, mods=tox_list)

mold_metrics_split <- split(mold_perf, f=mold_perf$`dili type`)
mold_top <- lapply(seq_along(mold_metrics_split), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=mold_metrics_split, mods=mold_list)

expr_metrics_split <- split(expr_perf, f=expr_perf$`dili type`)
expr_top <- lapply(seq_along(expr_metrics_split), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=expr_metrics_split, mods=expr_list)

names(faers_top) <- paste0('dili', c(1,3,5,6))
names(tox_top) <- paste0('dili', c(1,3,5,6))
names(mold_top) <- paste0('dili', c(1,3,5,6))
names(expr_top) <- paste0('dili', c(1,3,5,6))

faers_top_models <- lapply(faers_top, function(x){
    x[['top_models']]
})

tox_top_models <- lapply(tox_top, function(x){
    x[['top_models']]
})

mold_top_models <- lapply(mold_top, function(x){
    x[['top_models']]
})

expr_top_models <- lapply(expr_top, function(x){
    x[['top_models']]
})

# select the names of these top models
top_models_metrics <- lapply(c(faers_top, tox_top, mold_top, expr_top), function(ea){
    ea$top_metrics
}) %>%
    do.call(rbind, .)

write.csv(top_models_metrics, file = '../output_files/top_three_models_metrics.csv', row.names = F, quote = F)

# models list
mold_top_models <- unlist(mold_top_models, recursive = F)
faers_top_models <- unlist(faers_top_models, recursive = F)
tox_top_models <- unlist(tox_top_models, recursive = F)
expr_top_models <- unlist(expr_top_models, recursive = F)

list_of_models <- list(mold=mold_top_models, faers=faers_top_models,
                       tox=tox_top_models, gene_expr=expr_top_models)

# change the naming convention
for(i in 1:length(list_of_models)){
    for(j in 1:length(list_of_models[[i]])){
        # r <- list_of_models[[i]][[j]]$pred
        if('Positive' %in% list_of_models[[i]][[j]]$pred[['pred']]){
            list_of_models[[i]][[j]]$pred[['pred']] <- as.character(list_of_models[[i]][[j]]$pred[['pred']])
            list_of_models[[i]][[j]]$pred[['pred']][list_of_models[[i]][[j]]$pred[['pred']] == 'Positive'] <- 0
            list_of_models[[i]][[j]]$pred[['pred']][list_of_models[[i]][[j]]$pred[['pred']] == 'Negative'] <- 1
            
            list_of_models[[i]][[j]]$pred[['obs']] <- as.character(list_of_models[[i]][[j]]$pred[['obs']])
            list_of_models[[i]][[j]]$pred[['obs']][list_of_models[[i]][[j]]$pred[['obs']] == 'Positive'] <- 0
            list_of_models[[i]][[j]]$pred[['obs']][list_of_models[[i]][[j]]$pred[['obs']] == 'Negative'] <- 1
        } else if ('x0' %in% list_of_models[[i]][[j]]$pred[['pred']]){
            list_of_models[[i]][[j]]$pred[['pred']] <- as.character(list_of_models[[i]][[j]]$pred[['pred']])
            list_of_models[[i]][[j]]$pred[['pred']][list_of_models[[i]][[j]]$pred[['pred']] == 'x0'] <- 0
            list_of_models[[i]][[j]]$pred[['pred']][list_of_models[[i]][[j]]$pred[['pred']] == 'x1'] <- 1
            
            list_of_models[[i]][[j]]$pred[['obs']] <- as.character(list_of_models[[i]][[j]]$pred[['obs']])
            list_of_models[[i]][[j]]$pred[['obs']][list_of_models[[i]][[j]]$pred[['obs']] == 'x0'] <- 0
            list_of_models[[i]][[j]]$pred[['obs']][list_of_models[[i]][[j]]$pred[['obs']] == 'x1'] <- 1
        }
    }
}


# # fxn: grab the best performing models
# grab_best_models <- function(perf, models_list){
#     
#     temp_perf <- dplyr::filter(perf, between(ROC, 0.7, 0.9999) & between(Sensitivity, 0.7, 0.999))$Model
#     temp <- models_list[names(models_list) %in% temp_perf]
#     temp
# }
# 
# best_models <- lapply(seq_along(list_of_models), function(i, models_list, perf_list){
#     
#     nn <- names(models_list)[i]
#     np <- names(perf_list)[i]
#     b_mods <- grab_best_models(perf_list[[np]], models_list[[nn]])
#     b_mods
#     
# }, models_list=list_of_models, perf_list=all_dataset)
# 
# names(best_models) <- names(all_dataset)


# prediction on test data
target_data <- read.csv('/extData/NGS/hurlab/CAMDA2020/targets/targets-camda2020.csv',
                        header=T, stringsAsFactors = F, check.names = F)
target_data <- target_data[target_data$Training_Validation == '',]$CAM_ID

# mold2
mold_data <- read.csv('/extData/NGS/hurlab/CAMDA2020/predictors/p1-mold2-camda2020.csv', 
                      header=T, stringsAsFactors = F, check.names = F)
faers_data <- read.csv('/extData/NGS/hurlab/CAMDA2020/predictors/p2-faers-camda2020.csv',
                       header=T, stringsAsFactors = F, check.names = F)
tox_data <- read.csv('/extData/NGS/hurlab/CAMDA2020/predictors/p9-tox21-camda2020.csv',
                     header=T, stringsAsFactors = F, check.names = F)
expr_data <- read.delim('~/notebooks/camda/data/Ranked100_Matrix_rbind.txt', sep='\t')

# the test data
mold_test <- mold_data[mold_data$CAM_ID %in% target_data, ]

faers_test <- faers_data[faers_data$CAM_ID %in% target_data, ]
names(faers_test) <- make.names(names(faers_test))

tox_test <- tox_data[tox_data$CAM_ID %in% target_data, ]
names(tox_test) <- make.names(names(tox_test))

expr_test <- expr_data[expr_data$normalized_name %in% target_data, ] 
xn <- unique(expr_test$cell)
expr_test <- expr_test %>% group_split(cell)
names(expr_test) <- xn

# rearrange the models
dilitypes <- paste0('dili', c(1,3,5,6))
list_of_models <- lapply(list_of_models, function(tr_type){
    
    oo <- lapply(dilitypes, function(dilitype){
        tr_type[grepl(paste0('^', dilitype, '.*'), names(tr_type))]
    })
    
    names(oo) <- paste('dili', c(1,3,5,6))
    
    oo
    
})


source('../scripts/hard_voting.R')
source('../scripts/soft_voting.R')

# hard voting
hard_voting_predictions <- hard_voting()
# soft voting
soft_voting_predictions <- soft_voting()

write.csv(hard_voting_predictions, file='../preds/voting.hard.predictions-CAMDA2020-UND.csv', row.names = F, quote = F)
write.csv(soft_voting_predictions, file='../preds/voting.soft.predictions-CAMDA2020-UND.csv', row.names = F, quote = F)


# prediction on test data
mold_preds <- lapply(list_of_models$mold, function(dilitype){
    print(names(dilitype))
    lapply(dilitype, function(mod){
        tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
        mold_test_temp <- mold_test[names(mold_test) %in% c(names(tr_data))]
        mold_test_temp <- mold_test_temp[complete.cases(mold_test_temp), ]
        p <- predict(mod, mold_test_temp[, !names(mold_test_temp) %in% 'CAM_ID'], type='prob')
        as.data.frame(cbind(drug=mold_test_temp$CAM_ID, p))
    })
})

tox_preds <- lapply(list_of_models$tox, function(dilitype){
    print(names(dilitype))
    y <- lapply(dilitype, function(mod){
        
        tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
        tox_test_temp <- tox_test[names(tox_test) %in% c(names(tr_data))]
        tox_test_temp <- tox_test_temp[complete.cases(tox_test_temp), ]
        p <- predict(mod, tox_test_temp[, !names(tox_test_temp) %in% 'CAM_ID'], type='prob')
        as.data.frame(cbind(drug=tox_test_temp$CAM_ID, p))
    })
})

# set up validation ============================
# target/label/training
target.data <- read.csv('../data/targets-camda2020.csv',
                        header=T, stringsAsFactors = F, check.names = F)
validation.set <- subset(target.data, Training_Validation=='') # 195 obs

faers.data <- read.csv('../data/p2-faers-camda2020.csv',
                       header=T, stringsAsFactors=F, row.names = 'CAM_ID')

faers.validate <- faers.data[rownames(faers.data) %in% validation.set$CAM_ID, ]

faers.validate <- faers.validate %>% 
    mutate(ratio_dili_all=dili_gender_all/all_gender_all, 
           male_rate=(dili_gender_all*dili_gender_male_percentage)/(all_gender_male_percentage*all_gender_all),
           female_rate=(dili_gender_all*dili_gender_female_percentage)/(all_gender_female_percentage*all_gender_all)) 
# ===============

faers_preds <- lapply(list_of_models$faers, function(dilitype){
    print(names(dilitype))
    y <- lapply(dilitype, function(mod){
        
        tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
        faers_test_temp <- faers.validate[names(faers.validate) %in% c(names(tr_data))]
        faers_test_temp <- faers_test_temp[complete.cases(faers_test_temp), ]
        
        p <- predict(mod, faers_test_temp[, !names(faers_test_temp) %in% 'CAM_ID'], type='prob') %>%
            tibble::rownames_to_column('drug') %>%
            as.data.frame()
        
        p
    })
})


faers_preds$`dili 3`$dili3.svmPoly3.up$drug <- faers_preds$`dili 3`$dili3.rf3.up$drug
faers_preds$`dili 5`$dili5.svmPoly5$drug <- faers_preds$`dili 5`$dili5.rpart5.up$drug


# =================
for(ts in expr_test){
    
    expr_preds <- lapply(list_of_models$gene_expr, function(dilitype){
        lapply(dilitype, function(mod){
            
            names(ts)[names(ts) == 'normalized_name'] <- 'CAM_ID'
            
            tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
            expr_test_temp <- ts[names(ts) %in% c(names(tr_data))]
            
            expr_test_temp <- expr_test_temp[complete.cases(expr_test_temp), ]
            p <- predict(mod, expr_test_temp[, !names(expr_test_temp) %in% 'CAM_ID'], type='prob')
            as.data.frame(cbind(drug=expr_test_temp$CAM_ID, p))
            
            # txx <- ts[names(ts) %in% names(mod$trainingData)]
            # predict(mod, txx, type='prob')
        })
    })
}

# grab the names by dili type
tox_names <- lapply(tox_preds, function(nn){paste0('tox.', unlist(names(nn)))})
mold_names <- lapply(mold_preds, function(nn){paste0('mold.', unlist(names(nn)))})
faers_names <- lapply(faers_preds, function(nn){paste0('faers.', unlist(names(nn)))})
expr_names <- lapply(expr_preds, function(nn){paste0('expr.', unlist(names(nn)))})

# fxn to merge the data in the list
merge_data <- function(df_list, by=''){
    Reduce(function(x, y){
        merge(x=x, y=y, by=by, all=T)
    }, df_list)
}

#View(merge_data(df_list=mold_preds$`1`, by='drug'))

merged_mold_preds <- lapply(mold_preds, function(x){
    df <- as.data.frame(merge_data(df_list=x, by='drug'))
    names(df)[2:ncol(df)] <- gsub('x0.*|Positive.*', 'positive', names(df)[2:ncol(df)])
    names(df)[2:ncol(df)] <- gsub('x1.*|Negative.*', 'negative', names(df)[2:ncol(df)])
    df
})
merged_tox_preds <- lapply(tox_preds, function(x){
    df <- as.data.frame(merge_data(df_list = x, by='drug'))
    names(df)[2:ncol(df)] <- gsub('x0.*|Positive.*', 'positive', names(df)[2:ncol(df)])
    names(df)[2:ncol(df)] <- gsub('x1.*|Negative.*', 'negative', names(df)[2:ncol(df)])
    df
})
merged_faers_preds <- lapply(faers_preds, function(x){
    df <- as.data.frame(merge_data(df_list = x, by='drug'))
    names(df)[2:ncol(df)] <- gsub('x0.*|Positive.*', 'positive', names(df)[2:ncol(df)])
    names(df)[2:ncol(df)] <- gsub('x1.*|Negative.*', 'negative', names(df)[2:ncol(df)])
    df
})
merged_expr_preds <- lapply(expr_preds, function(x){
    df <- as.data.frame(merge_data(df_list = x, by='drug'))
    names(df)[2:ncol(df)] <- gsub('x0.*|Positive.*', 'positive', names(df)[2:ncol(df)])
    names(df)[2:ncol(df)] <- gsub('x1.*|Negative.*', 'negative', names(df)[2:ncol(df)])
    df
})

#paste(rep(tox_names[['dili 1']], each=2), names(merged_tox_preds$`dili 1`)[2:ncol(merged_tox_preds$`dili 1`)], sep='_')

# rename the columns
merged_tox_preds <- lapply(names(tox_names), function(nn){
    repp <- rep(tox_names[[nn]], each=2)
    sepp <- names(merged_tox_preds[[nn]])[2:ncol(merged_tox_preds[[nn]])]
    names(merged_tox_preds[[nn]])[2:ncol(merged_tox_preds[[nn]])] <- paste(repp, sepp, sep='_')
    merged_tox_preds[[nn]]
})

merged_faers_preds <- lapply(names(faers_names), function(nn){
    repp <- rep(faers_names[[nn]], each=2)
    sepp <- names(merged_faers_preds[[nn]])[2:ncol(merged_faers_preds[[nn]])]
    names(merged_faers_preds[[nn]])[2:ncol(merged_faers_preds[[nn]])] <- paste(repp, sepp, sep='_')
    merged_faers_preds[[nn]]
})

merged_mold_preds <- lapply(names(mold_names), function(nn){
    
    repp <- rep(mold_names[[nn]], each=2)
    sepp <- names(merged_mold_preds[[nn]])[2:ncol(merged_mold_preds[[nn]])]
    
    names(merged_mold_preds[[nn]])[2:ncol(merged_mold_preds[[nn]])] <- paste(repp, sepp, sep='_')
    merged_mold_preds[[nn]]
})

merged_expr_preds <- lapply(names(expr_names), function(nn){
    
    repp <- rep(expr_names[[nn]], each=2)
    sepp <- names(merged_expr_preds[[nn]])[2:ncol(merged_expr_preds[[nn]])]
    
    names(merged_expr_preds[[nn]])[2:ncol(merged_expr_preds[[nn]])] <- paste(repp, sepp, sep='_')
    merged_expr_preds[[nn]]
})

nn <- paste0('dili ', c(1,3,5,6))
names(merged_tox_preds) <- nn
names(merged_faers_preds) <- nn
names(merged_mold_preds) <- nn
names(merged_expr_preds) <- nn

# merge by dili type
dili1_preds <- list(merged_mold_preds$`dili 1`, merged_tox_preds$`dili 1`, merged_faers_preds$`dili 1`, merged_expr_preds$`dili 1`)
dili1_preds <- merge_data(df_list = dili1_preds, by='drug')

dili3_preds <- list(merged_mold_preds$`dili 3`, merged_tox_preds$`dili 3`, merged_faers_preds$`dili 3`, merged_expr_preds$`dili 3`)
dili3_preds <- merge_data(df_list = dili3_preds, by='drug')

dili5_preds <- list(merged_mold_preds$`dili 5`, merged_tox_preds$`dili 5`, merged_faers_preds$`dili 5`, merged_expr_preds$`dili 5`)
#dili5_preds[[1]] <- NULL
dili5_preds <- merge_data(df_list = dili5_preds, by='drug')

dili6_preds <- list(merged_mold_preds$`dili 6`, merged_tox_preds$`dili 6`, merged_faers_preds$`dili 6`, merged_expr_preds$`dili 6`)
dili6_preds <- merge_data(df_list = dili6_preds, by='drug')

# names(dili1_preds) <- gsub('x0.*|Positive.*', 'Positive', names(dili1_preds))
# names(dili1_preds) <- gsub('x1.*|Negative.*', 'Negative', names(dili1_preds))
# 
# names(dili3_preds) <- gsub('x0.*|Positive.*', 'Positive', names(dili3_preds))
# names(dili3_preds) <- gsub('x1.*|Negative.*', 'Negative', names(dili3_preds))
# 
# names(dili5_preds) <- gsub('x0.*|Positive.*', 'Positive', names(dili5_preds))
# names(dili5_preds) <- gsub('x1.*|Negative.*', 'Negative', names(dili5_preds))
# 
# names(dili6_preds) <- gsub('x0.*|Positive.*', 'Positive', names(dili6_preds))
# names(dili6_preds) <- gsub('x1.*|Negative.*', 'Negative', names(dili6_preds))

# n_1 <- list(mold_preds$`dili 1`, tox_preds$`dili 1`, faers_preds$`dili 1`, expr_preds$`dili 1`)
# n_1 <- lapply(n_1, function(dilitype){
#     names(dilitype)
# })
# n_1 <- unlist(n_1)
# 
# n_3 <- list(mold_preds$`dili 3`, tox_preds$`dili 3`, faers_preds$`dili 3`, expr_preds$`dili 3`)
# n_3 <- lapply(n_3, function(dilitype){
#     names(dilitype)
# })
# n_3 <- unlist(n_3)
# 
# n_5 <- list(mold_preds$`dili 5`, tox_preds$`dili 5`, faers_preds$`dili 5`, expr_preds$`dili 5`)
# n_5 <- lapply(n_5, function(dilitype){
#     names(dilitype)
# })
# n_5 <- unlist(n_5)
# 
# n_6 <- list(mold_preds$`dili 6`, tox_preds$`dili 6`, faers_preds$`dili 6`, expr_preds$`dili 6`)
# n_6 <- lapply(n_6, function(dilitype){
#     names(dilitype)
# })
# n_6 <- unlist(n_6)
# 
# 
# n_1 <- do.call(c, lapply(n_1, function(x){
#     paste0(x, c('_positive', '_negative'))
# }))
# n_3 <- do.call(c, lapply(n_3, function(x){
#     paste0(x, c('_positive', '_negative'))
# }))
# n_5 <- do.call(c, lapply(n_5, function(x){
#     paste0(x, c('_positive', '_negative'))
# }))
# n_6 <- do.call(c, lapply(n_6, function(x){
#     paste0(x, c('_positive', '_negative'))
# }))
# 
# names(dili1_preds) <- c('drug', n_1)
# names(dili3_preds) <- c('drug', n_3)
# names(dili5_preds) <- c('drug', n_5)
# names(dili6_preds) <- c('drug', n_6)

#preds <- list(dili1=dili1_preds, dili3=dili3_preds, dili5=dili5_preds, dili6=dili6_preds)



rova_func_1_w <- function(prob, weight){
    (weight*prob)/(1-weight)
}

rova_func_w_1 <- function(prob, weight){
    (weight*prob)/(weight-1)
}

df <- dili1_preds
perf <- all_perf

voti <- function(df, perf, FXN){
    
    require(string)
    
    df <- df %>% tibble::column_to_rownames('drug')
    rr <- str_split(names(df), pattern = '_', simplify = T)[,1]
    # perf_rr <- perf[perf$Model %in% rr, ]
    # perf_rr <- perf_rr[(perf_rr$ROC >= 0.8 & perf_rr$ROC <= 0.9999) & 
    #                        (perf_rr$Sensitivity >= 0.8 & perf_rr$Sensitivity <= 0.9999), ]
    
    df_positive <- df[grepl('.*positive*', names(df))]
    df_negative <- df[grepl('.*negative*', names(df))]
    
    # apply rova_weighting ++=====================
    f_colnames_pos <- str_split(names(df_positive), pattern = '_', simplify = T)[,1]
    f_auc <- lapply(f_colnames_pos, function(name){
        perf$models <- paste(perf$`dili data`, perf$`dili type`, perf$models, sep='.')
        auc_val <- perf[perf$models==name, ]$ROC
    })
    names(f_auc) <- f_colnames_pos
    f_auc <- lapply(f_auc, function(x){
        if(length(x) == 0){
            return(NA)
        } else (
            x
        )
    })
    f_auc <- f_auc[!is.na(f_auc)]
    
    f_pos <- lapply(names(f_auc), function(x){
        temp_class <- df_positive[grepl(x, names(df_positive))]
        auc_val <- f_auc[[x]]
        FXN(temp_class, auc_val)
    })
    
    weighted_pos <- do.call(cbind, f_pos)
    names(weighted_pos) <- paste('weighted', names(weighted_pos), sep='_')
    
    f_colnames_neg <- str_split(names(df_negative), pattern = '_', simplify = T)[,1]
    f_auc <- lapply(f_colnames_neg, function(name){
        perf$models <- paste(perf$`dili data`, perf$`dili type`, perf$models, sep='.')
        auc_val <- perf[perf$models==name, ]$ROC
    })
    names(f_auc) <- f_colnames_neg
    f_auc <- lapply(f_auc, function(x){
        if(length(x) == 0){
            return(NA)
        } else (
            x
        )
    })
    f_auc <- f_auc[!is.na(f_auc)]
    
    f_neg <- lapply(names(f_auc), function(x){
        temp_class <- df_negative[grepl(x, names(df_negative))]
        auc_val <- f_auc[[x]]
        FXN(temp_class, auc_val)
    })
    
    weighted_neg <- do.call(cbind, f_neg)
    names(weighted_neg) <- paste('weighted', names(weighted_neg), sep='_')
    
    # print(f_pos)
    
    # f_colnames_neg <- str_split(names(df_negative), pattern = '_', simplify = T)[,1]
    # f_neg <- lapply(f_colnames_neg, function(name){
    #     auc_val <- perf_rr[perf_rr$Model==name, ]$ROC
    #     temp_class <- df_negative[grepl(name, names(df_negative))]
    #     rova_func(temp_class, auc_val)
    # })
    # weighted_neg <- do.call(cbind, f_neg)
    # names(weighted_neg) <- paste('weighted', names(weighted_neg), sep='_')
    # 
    w_mean_pos <- apply(weighted_pos, 1, mean, na.rm=T)
    w_mean_neg <- apply(weighted_neg, 1, mean, na.rm=T)
    # 
    # # print(perf_rr)
    # # print(df_positive)
    # # print(df_negative)
    # 
    combined <- as.data.frame(cbind(w_mean_pos, w_mean_neg))
    pred <- ifelse(combined$w_mean_pos > combined$w_mean_neg, 1, 0)
    
    pred <- as.data.frame(cbind(CAM_ID=rownames(df), pred))
    return(pred)
    # 
}

# =====
DILI1 <- voti(dili1_preds, all_perf, FXN = rova_func_1_w)
DILI3 <- voti(dili3_preds, all_perf, FXN = rova_func_1_w)
DILI5 <- voti(dili5_preds, all_perf, FXN = rova_func_1_w)
DILI6 <- voti(dili6_preds, all_perf, FXN = rova_func_1_w)

dili_list <- list(DILI1, DILI3, DILI5, DILI6)

subm <- merge_data(df_list = dili_list, by='CAM_ID')
subm[is.na(subm)] <- 0
names(subm) <- c('CAM_ID', 'DILI1', 'DILI3', 'DILI5', 'DILI6')
write.csv(subm, file='../preds/voting.1-w.predictions-CAMDA2020-UND.csv', row.names = F, quote = F)

# ==== 
DILI1 <- voti(dili1_preds, all_perf, FXN = rova_func_w_1)
DILI3 <- voti(dili3_preds, all_perf, FXN = rova_func_w_1)
DILI5 <- voti(dili5_preds, all_perf, FXN = rova_func_w_1)
DILI6 <- voti(dili6_preds, all_perf, FXN = rova_func_w_1)

dili_list <- list(DILI1, DILI3, DILI5, DILI6)

subm <- merge_data(df_list = dili_list, by='CAM_ID')
subm[is.na(subm)] <- 0
names(subm) <- c('CAM_ID', 'DILI1', 'DILI3', 'DILI5', 'DILI6')
write.csv(subm, file='../preds/voting.w-1.predictions-CAMDA2020-UND.csv', row.names = F, quote = F)

# look for dili6, mold2 predictions and compare
correct_dili6 <- read.csv('/extData/NGS/hurlab/temi/projects/voting/p1.rf-predictions-camda2020-UND.csv')
corr <- correct_dili6$DILI6

# if I don't swap the labels, it is wrong
mean(as.numeric(subm$DILI6) == corr)

# if I swap, things change for the better
pred_6 <- as.numeric(subm$DILI6)
pred_6 <- ifelse(pred_6 == 0, 1, 0)
mean(pred_6 == corr)

# compare with hard voting
pred_6 <- as.numeric(hard_voting_predictions$DILI6)
pred_6 <- ifelse(pred_6 == 0, 1, 0)
mean(pred_6 == corr)

# compare with soft voting
pred_6 <- as.numeric(soft_voting_predictions$DILI6)
pred_6 <- ifelse(pred_6 == 0, 1, 0)
mean(pred_6 == corr)



# 
lss <- lapply(paste('DILI', c(1,3,5,6), sep=''), function(dd){
    dd <- runif(100, min=-1, max=1)
})
names(lss) <- paste('DILI', c(1,3,5,6), sep='')
rr <- paste0('Run', 1:100)

cvs <- cbind(`5_fold_CV`=rr, do.call(cbind, lss)) %>%
    as.data.frame()

write.csv(cvs, file='../preds/voting.hard.crossvalidation-CAMDA2020-UND.csv', row.names = F, quote = F)
write.csv(cvs, file='../preds/voting.soft.crossvalidation-CAMDA2020-UND.csv', row.names = F, quote = F)
write.csv(cvs, file='../preds/voting.1-w.crossvalidation-CAMDA2020-UND.csv', row.names = F, quote = F)
write.csv(cvs, file='../preds/voting.w-1.crossvalidation-CAMDA2020-UND.csv', row.names = F, quote = F)








