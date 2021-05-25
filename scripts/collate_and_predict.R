

# this script collects all the models, selects the top 3 models, and makes predictions on the test data

curr <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(curr))

source('../scripts/packages.R')

load('../models/faers_models.RData')
load('../models/tox_models.RData')
load('../models/mold_models.RData')

#unlist(faers_models, recursive = F)

# collate metrics 
faers_metrics <- collate_metrics(faers_models)
tox_metrics <- collate_metrics(tox_models)
mold_metrics <- collate_metrics(mold_models)

# select top five models
faers_metrics_split <- split(faers_metrics, f=faers_metrics$`dili type`)
faers_top <- lapply(seq_along(faers_metrics_split), function(i, mets, mods){

    select.top.three(mets[[i]], mods[[i]])
    
}, mets=faers_metrics_split, mods=faers_models)

tox_metrics_split <- split(tox_metrics, f=tox_metrics$`dili type`)
tox_top <- lapply(seq_along(tox_metrics_split), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=tox_metrics_split, mods=tox_models)

mold_metrics_split <- split(mold_metrics, f=mold_metrics$`dili type`)
mold_top <- lapply(seq_along(mold_metrics_split), function(i, mets, mods){
    
    select.top.three(mets[[i]], mods[[i]])
    
}, mets=mold_metrics_split, mods=mold_models)

names(faers_top) <- paste0('dili', c(1,3,5,6))
names(tox_top) <- paste0('dili', c(1,3,5,6))
names(mold_top) <- paste0('dili', c(1,3,5,6))

# save the metrics into one
faers_top_metrics <- lapply(faers_top, function(x){
    x[['top_metrics']]
})

tox_top_metrics <- lapply(tox_top, function(x){
    x[['top_metrics']]
})

mold_top_metrics <- lapply(mold_top, function(x){
    x[['top_metrics']]
})

all_dili_metrics <- as.data.frame(
    rbind(do.call(rbind, faers_top_metrics),
      do.call(rbind, tox_top_metrics),
      do.call(rbind, mold_top_metrics))
)

write.csv(x=all_dili_metrics, file = '../output_files/all_non_gene_expression_metrics.csv', row.names = F)

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

#tox_predictions[,grepl('.*dili3.*', names(tox_predictions))]

organize_predictions <- function(preds){
  
  temp_res <- preds[grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', colnames(preds))] 
  temp_non <- preds[!grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', colnames(preds))] %>%
    dplyr::select(-CAM_ID)
  
  temp_1 <- list(resampled=temp_res, non_resampled=temp_non)
  
  dilis <- c('dili1', 'dili3', 'dili5', 'dili6')
  
  temp_2 <- lapply(temp_1, function(df){
    tt <- lapply(dilis, function(x){
      df[, grepl(paste0('.*', x, '.*'), names(df))]
    })
    names(tt) <- dilis %>%
      toupper()
    
    tt
  })
  
  tops_names <- paste('top_', 1:3, sep='')
  
  every_col <- c(1,2,3)
  
  temp_3 <- lapply(temp_2, function(df){
    
    tt <- lapply(every_col, function(x){
      lapply(df, function(d){
        d[, x]
      })
    })
    
    names(tt) <- tops_names 
    
    tt
    
  })
  
  out <- lapply(temp_3, function(l1){
    lapply(l1, function(l2){
      tt <- do.call(data.frame, l2)
      
      as.data.frame(cbind(CAM_ID=preds$CAM_ID, tt))
    })
  })
  
  out
}

faers_preds <- organize_predictions(faers_predictions) %>%
  unlist(., recursive = F)
tox_preds <- organize_predictions(tox_predictions) %>%
  unlist(., recursive = F)
mold_preds <- organize_predictions(mold_predictions) %>%
  unlist(., recursive = F)


sapply(names(faers_preds), function(x){
  write.csv(faers_preds[[x]], file = paste0('../preds/p2.', x, '.predictions-CAMDA2020-UND.csv'), 
            quote = F, row.names = F)
})

sapply(names(mold_preds), function(x){
  write.csv(mold_preds[[x]], file = paste0('../preds/p1.', x, '.predictions-CAMDA2020-UND.csv'), 
            quote = F, row.names = F)
})

sapply(names(tox_preds), function(x){
  write.csv(tox_preds[[x]], file = paste0('../preds/p9.', x, '.predictions-CAMDA2020-UND.csv'), 
            quote = F, row.names = F)
})

# crossvalidations

mcc_cv <- function(xx){
  
  cm <- table(xx$obs, xx$pred)
  true.pos <- as.numeric(cm[1,1])
  false.pos <- as.numeric(cm[1,2])
  false.neg <- as.numeric(cm[2,1])
  true.neg <- as.numeric(cm[2,2])
  
  above <- true.pos*true.neg - false.pos*false.neg
  below <- as.double(true.pos+false.pos)*as.double(true.pos+false.neg)*as.double(true.neg+false.pos)*as.double(true.neg+false.neg)
  
  mcc <- above/sqrt(below)
  mcc
}

# temp.table <- table(x.model.list$nnet1.up$pred$pred, x.model.list$nnet1.up$pred$obs)
# 
# mcc_cv(temp.table)


collect_model_mcc <- function(x.model.list){
  
  vv <- lapply(x.model.list, function(xd){
    
    uu <- lapply(xd, function(xm){
      
      x_mcc <- as.data.frame(xm$pred)
      x_mcc <- x_mcc[complete.cases(x_mcc), ]
      x_mcc[,1] <- as.character(x_mcc[,1])
      x_mcc[,2] <- as.character(x_mcc[,2])
      # x_mcc[x_mcc=="Positive"] <- 0
      # x_mcc[x_mcc=="Negative"] <- 1
      
      x.details <- x_mcc %>%
        separate(Resample,c("cv","rep"),sep="\\.") %>%
        split(., .$rep) %>%
        lapply(., function(df){
          #tt <- table(df$obs, df$pred)
          mcc_cv(df)
        })
      
      x.details <- do.call(rbind, x.details) %>%
        as.data.frame() 
      
      names(x.details) <- c("mcc")
      
      x.details
      
    })
    
    # uu <- lapply(uu, function(u){
    #   u %>%
    #     tibble::column_to_rownames('5_fold_CV')
    # })
    
    yy <- as.data.frame(do.call(cbind, uu))
    colnames(yy) <- names(uu)
    rownames(yy) <- paste0('Run', seq(1:100))
    yy <- yy %>%
      tibble::rownames_to_column('5_fold_CV')
    
    yy
    
  })
  
  vv
  
}

faers_mcc <- collect_model_mcc(faers_top_models)
tox_mcc <- collect_model_mcc(tox_top_models)
mold_mcc <- collect_model_mcc(mold_top_models)


organize_crossvalidations <- function(cv_list){
  
  oo <- lapply(cv_list, function(cv){
    
    temp_res <- cv[grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', colnames(cv))] 
    temp_non <- cv[!grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', colnames(cv))] %>%
      dplyr::select(-`5_fold_CV`)
    
    temp_1 <- list(resampled=temp_res, non_resampled=temp_non)
    
    temp_1
    
  })
  
  tt <- unlist(oo, recursive = F)
  
  tt_res <- tt[!grepl('.*non_resampled.*', names(tt))]
  tt_non_res <- tt[grepl('.*non_resampled.*', names(tt))]
  
  tt_list <- list(resampled=tt_res, non_resampled=tt_non_res)
  
  tt_tops <- lapply(tt_list, function(ty){
    
    top_1 <- lapply(ty, function(t){
      t[, 1]
    })
    
    top_2 <- lapply(ty, function(t){
      t[, 2]
    })
    
    top_3 <- lapply(ty, function(t){
      t[, 3]
    })
    
    top_list <- list(top_1, top_2, top_3)
    top_list <- lapply(top_list, function(top){
      temp <- do.call(cbind, top) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('5_fold_CV') %>%
        mutate(`5_fold_CV`=paste0('Run', seq(1:100)))
      
      names(temp)[2:5] <- c('DILI1', 'DILI3', 'DILI5', 'DILI6')
      
      temp
    })
    
    names(top_list) <- c('top_1', 'top_2', 'top_3')
    
    top_list
    
    
  })
  
  tt_tops
  
  
  
  # top_list <- list(top_1, top_2, top_3)
  # top_list <- lapply(top_list, function(top){
  #   temp <- do.call(cbind, top) %>%
  #     as.data.frame() %>%
  #     tibble::rownames_to_column('5_fold_CV') %>%
  #     mutate(`5_fold_CV`=paste0('Run', seq(1:100)))
  #   
  #   names(temp)[2:5] <- c('DILI1', 'DILI3', 'DILI5', 'DILI6')
  #   
  #   temp
  # })
  # 
  # names(top_list) <- c('top_1', 'top_2', 'top_3')
  
  # dilis <- c('dili1', 'dili3', 'dili5', 'dili6')
  # 
  # temp_2 <- lapply(temp_1, function(df){
  #   tt <- lapply(dilis, function(x){
  #     df[, grepl(paste0('.*', x, '.*'), names(df))]
  #   })
  #   names(tt) <- dilis %>%
  #     toupper()
  #   
  #   tt
  # })
  # 
  # tops_names <- paste('top_', 1:3, sep='')
  # 
  # every_col <- c(1,2,3)
  # 
  # temp_3 <- lapply(temp_2, function(df){
  #   
  #   tt <- lapply(every_col, function(x){
  #     lapply(df, function(d){
  #       d[, x]
  #     })
  #   })
  #   
  #   names(tt) <- tops_names 
  #   
  #   tt
  #   
  # })
  # 
  # out <- lapply(temp_3, function(l1){
  #   lapply(l1, function(l2){
  #     tt <- do.call(data.frame, l2)
  #     
  #     as.data.frame(cbind(CAM_ID=preds$CAM_ID, tt))
  #   })
  # })
  # 
  # out
}

faers_cvs <- organize_crossvalidations(faers_mcc) %>%
  unlist(., recursive = F)
tox_cvs <- organize_crossvalidations(tox_mcc) %>%
  unlist(., recursive = F)
mold_cvs <- organize_crossvalidations(mold_mcc) %>%
  unlist(., recursive = F)

sapply(names(faers_cvs), function(x){
  write.csv(faers_cvs[[x]], file = paste0('../preds/p2.', x, '.crossvalidation-CAMDA2020-UND.csv'), 
            quote = F, row.names = F)
})

sapply(names(mold_cvs), function(x){
  write.csv(mold_cvs[[x]], file = paste0('../preds/p1.', x, '.crossvalidation-CAMDA2020-UND.csv'), 
            quote = F, row.names = F)
})

sapply(names(tox_cvs), function(x){
  write.csv(tox_cvs[[x]], file = paste0('../preds/p9.', x, '.crossvalidation-CAMDA2020-UND.csv'), 
            quote = F, row.names = F)
})




# temp_res <- mold_predictions[grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', colnames(mold_predictions))] 
# temp_non <- mold_predictions[!grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', colnames(mold_predictions))] %>%
#   dplyr::select(-CAM_ID)
# 
# temp_1 <- list(resampled=temp_res, non_resampled=temp_non)
# 
# dilis <- c('dili1', 'dili3', 'dili5', 'dili6')
# 
# temp_2 <- lapply(temp_1, function(df){
#   tt <- lapply(dilis, function(x){
#     df[, grepl(paste0('.*', x, '.*'), names(df))]
#   })
#   names(tt) <- dilis %>%
#     toupper()
#   
#   tt
# })
# 
# tops_names <- paste('top_', 1:3, sep='')
# 
# every_col <- c(1,2,3)
# 
# temp_3 <- lapply(temp_2, function(df){
#   
#   tt <- lapply(every_col, function(x){
#     lapply(df, function(d){
#       d[, x]
#     })
#   })
#   
#   names(tt) <- tops_names
#   
#   tt
#   
# })
# 
# 
# out <- lapply(temp_3, function(l1){
#   lapply(l1, function(l2){
#     tt <- do.call(data.frame, l2)
#     
#     as.data.frame(cbind(CAM_ID=tox_predictions$CAM_ID, tt))
#   })
# })
# 
# 
# 
# 
# 
# 
# 
# tox_predictions %>%
#   dplyr::select(grepl('.*dili3.*', names(tox_predictions)))
# 
# apply(tox_predictions, 2, function(x){
# })
# 
# 
# 
# # some predictions contained missing values
# tox_miss <- tox_predictions[, colSums(is.na(tox_predictions)) == nrow(tox_predictions)] %>%
#   colnames()
# mold_miss <- mold_predictions[, colSums(is.na(mold_predictions)) == nrow(mold_predictions)] %>%
#   colnames()
# faers_miss <- faers_predictions[, colSums(is.na(faers_predictions)) == nrow(faers_predictions)] %>%
#   colnames()
# 
# 
# 
# 
# 
# 
# 
# tox_predictions %>%
#   select_if(is.na)
# 
# apply(tox_predictions, 2, function(ytr){
#   is.na(ytr)
# }) %>%
#   tox_predictions[colSums(.) == 172]




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


















