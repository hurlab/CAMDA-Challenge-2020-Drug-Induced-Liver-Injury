

# predictions and evaluation 

curr <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(curr))

source('../scripts/packages.R')

create_validations <- function(d_type='', name=''){
    
    if(d_type == 'voting'){
        patt <- paste0('.*\\.', d_type, '\\..*\\.csv')
        # use the patterns to get the files paths
        files <- list.files('../evaluations/May2021', full.names = T, 
                            pattern = patt)
        # extract the names of the files
        nn <- lapply(files, function(f){
            stringr::str_match(f, pattern = paste('\\.*', d_type, '\\.(\\w+.\\w+)\\.predictions\\.*', sep = ''))[,2]
        })
        # read the files
        files_read <- lapply(files, function(file){
            df <- read.csv(file, header=T)
            df
        })
        # rename the files
        names(files_read) <- nn
    } else {
        patt <- paste0('.*\\.', d_type, '\\..*\\.csv')
        # use the patterns to get the files paths
        files <- list.files('../evaluations/May2021', full.names = T, 
                            pattern = patt)
        # extract the names of the files
        nn <- lapply(files, function(f){
            stringr::str_match(f, pattern = paste('\\.*', d_type, '\\.(\\w+.\\w+)\\.predictions\\.*', sep = ''))[,2]
        })
        # read the files
        files_read <- lapply(files, function(file){
            df <- read.csv(file, header=T)
            df
        })
        # rename the files
        names(files_read) <- nn
    }
    
    j <- lapply(seq_along(files_read), function(i){
        accuracy <- round(files_read[[i]]$acc, 3)
        dilitype <- files_read[[i]]$target
        recall <- round(files_read[[i]]$recall, 3)
        spec <- round(files_read[[i]]$spec, 3)
        mcc <- round(files_read[[i]]$mcc, 3)
        f1 <- round(files_read[[i]]$f1, 3)
        otu <- as.data.frame(cbind(accuracy, dilitype, sens=recall, spec, mcc, f1))
        otu$dataset <- name
        otu$model <- names(files_read)[i]
        otu
    })
    
    as.data.frame(do.call(rbind, j)) %>%
        dplyr::arrange(desc(accuracy))
}

mold2_validations <- create_validations('p1', 'mold2')
faers_validations <- create_validations('p2', 'faers')
tox21_validations <- create_validations('p9', 'tox21')
voting_validations <- create_validations('voting', 'voting') %>%
    dplyr::filter(model!='1-w')

all_validations <- rbind(mold2_validations, faers_validations, tox21_validations, voting_validations)
all_validations$set <- 'test set'

# top_3_models <- read.csv('../../validations/renewed_top_3_models.csv', header=T)
# 
# all_validations <- merge(all_validations, top_3_models, by=c('model', 'dataset', 'dilitype'))

top_models_metrics <- read.csv('../output_files/top_three_models_metrics.csv')
top_models_metrics$set <- 'training set'

top_models_metrics <- top_models_metrics %>%
    dplyr::group_by(dili.data, dili.type, training.type, desc(ROC)) %>%
    ungroup() %>%
    dplyr::select(-`desc(ROC)`)

top_models_metrics$ranking <- rep(paste('top_', c(1:3), sep = ''), nrow(top_models_metrics)/3)

all_validations <- all_validations %>%
    tidyr::separate(col=model, into=c('training.type', 'ranking'), sep='\\.', convert=TRUE)

all_validations$dilitype <- tolower(all_validations$dilitype)
all_validations$dataset <- dplyr::case_when(all_validations$dataset == 'mold2' ~ 'mold',
                                              all_validations$dataset == 'tox21' ~ 'tox',
                                              T ~ all_validations$dataset)
top_models_metrics$training.type <- gsub(pattern = '-', replacement = '_', x = top_models_metrics$training.type)


trupe <- top_models_metrics %>%
    dplyr::select(dili.data, dili.type, models, training.type, ranking)

test_set_metrics <- merge(all_validations, trupe, by.x=c('dataset', 'dilitype', 'training.type', 'ranking'), 
      by.y=c('dili.data', 'dili.type', 'training.type', 'ranking'), all.x=T)

names(test_set_metrics) <- c('dataset', 'dilitype', 'training_type', 'ranking', 'accuracy',
                             'sensitivity', 'specificity', 'mcc', 'f1', 'set', 'models')

test_set_metrics <- test_set_metrics %>%
  mutate(across(.cols = c(accuracy, sensitivity, specificity, mcc, f1), .fns = as.numeric)) %>% 
  mutate(balanced_accuracy = (sensitivity + specificity) / 2)

#test_set_metrics$binary_accuracy <- (as.numeric(test_set_metrics$sensitivity) + as.numeric(test_set_metrics$specificity))/2

training_set_metrics <- as.data.frame(top_models_metrics)
names(training_set_metrics) <- c('dataset', 'dilitype', 'models', 'training_type', 'AUC', 'sensitivity', 'specificity', 'mcc', 
                                 'balanced_accuracy', 'accuracy', 'precision', 'f1', 'set', 'ranking')



expr_metrics <- readxl::read_xlsx('../ModelSummary_v5.xlsx', sheet = 'FigTable', skip = 1)


expr_metrics <- expr_metrics %>%
  dplyr::select(-c(1))

names(expr_metrics) <- c('models', 'dilitype', 'algorithm', 'gene_size', 'AUC', 'sensitivity', 'specificity', 'predictors', 'MCC', 
                         'test_accuracy', 'test_precision', 'test_recall', 'test_specificity', 'test_F1', 'test_MCC')
expr_metrics$dilitype <- paste('DILI', expr_metrics$dilitype, sep='')
expr_metrics <- expr_metrics[!expr_metrics$dilitype=='DILINA', ]
expr_metrics$test_precision <- as.numeric(expr_metrics$test_precision)

# all_set_metrics
# 
# expr_metrics %>%
#   mutate(across(.cols = c(AUC, accuracy, sensitivity, specificity, MCC, f1, recall, test_accuracy, 
#                           test_precision, test_recall, test_specificity, test_F1, test_MCC), .fns = ~ round(.x, 2))) 

expr_metrics <- expr_metrics %>%
  dplyr::mutate(across(where(is.numeric), ~round(.x, 2)))

all_set_metrics <- dplyr::bind_rows(training_set_metrics, test_set_metrics)

test_expr_metrics <- expr_metrics %>%
  dplyr::filter(
    models %in% (all_set_metrics %>%
  dplyr::filter(set=='training set' & dataset=='expr') %>% .$models)
  ) 

test_expr_metrics <- test_expr_metrics %>%
  dplyr::select(-c(gene_size, AUC, sensitivity, specificity, predictors, MCC))

names(test_expr_metrics) <- c('models', 'dilitype', 'model_exact', 'accuracy', 'sensitivity', 'recall', 'specificity', 'f1', 'mcc')
test_expr_metrics$set <- 'test set'
#test_expr_metrics$dilitype <- paste('dili', test_expr_metrics$dilitype, sep = '')
test_expr_metrics$dataset <- 'expr'
test_expr_metrics$training_type <- 'non_resampled'
test_expr_metrics <- test_expr_metrics %>%
  mutate(across(.cols = c(accuracy, sensitivity, specificity, mcc, f1, recall), .fns = as.numeric)) %>% 
  mutate(balanced_accuracy = (sensitivity + specificity) / 2)

test_expr_metrics <- all_set_metrics %>%
  dplyr::filter(set=='training set' & dataset=='expr') %>%
  dplyr::select(models, ranking) %>%
  merge(., test_expr_metrics, by='models', all.y=T) %>%
  dplyr::arrange(dilitype, ranking)

# rbind with all_set_metrics

all_set_metrics <- dplyr::bind_rows(all_set_metrics, test_expr_metrics)

# plots ===========

# choose a color palette
rcartPalette <- colorRampPalette(rcartocolor::carto_pal(12, 'Safe'))

# define a function to extract legends
g_legend <- function(a.gplot){
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  
  legend <- tmp$grobs[[leg]]
  
  return(legend)
  
  
}

all_set_metrics <- all_set_metrics %>% mutate(model_exact=case_when(grepl(pattern = '.*rf.*', .$models) ~ 'rf',
                                                            grepl(pattern = '.*nb.*', .$models) ~ 'nb',
                                                            grepl(pattern = '.*svmRadial.*', .$models) ~ 'svmRadial',
                                                            grepl(pattern = '.*svmPoly.*', .$models) ~ 'svmPoly',
                                                            grepl(pattern = '.*nnet.*', .$models) ~ 'nnet',
                                                            grepl(pattern = '.*qda.*', .$models) ~ 'qda',
                                                            grepl(pattern = '.*glm.*', .$models) ~ 'glm',
                                                            grepl(pattern = '.*rpart.*', .$models) ~ 'rpart',
                                                            grepl(pattern = '.*lda.*', .$models) ~ 'lda',
                                                            grepl(pattern = '.*svmLinear.*', .$models) ~ 'svmLinear',
                                                            grepl(pattern = '.*svm.*', .$models) ~ 'svm',
                                                            TRUE ~ 'dunno'))

# === tox figures ===================
tox_top <- all_set_metrics[all_set_metrics$dataset=='tox', ]

# =====
tox_top <- tox_top %>%
  mutate(fill_val = case_when(
    (grepl('.*smote.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(SMOTE)',
    (grepl('.*rose.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(ROSE)',
    (grepl('.*up.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(upsampled)',
    (grepl('.*smote.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(SMOTE)',
    (grepl('.*rose.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(ROSE)',
    (grepl('.*up.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(upsampled)',
    (grepl('.*smote.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(SMOTE)',
    (grepl('.*rose.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(ROSE)',
    (grepl('.*up.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(upsampled)',
    (grepl('.*up.*', models) & grepl('.*rf.*', models)) ~ 'rf(upsampled)',
    (grepl('.*up.*', models) & grepl('.*nnet.*', models)) ~ 'nnet(upsampled)',
    (grepl('.*up.*', models) & grepl('.*rpart.*', models)) ~ 'rpart(upsampled)',
    (grepl('.*up.*', models) & grepl('.*qda.*', models)) ~ 'qda(upsampled)',
    (grepl('.*rose.*', models) & grepl('.*rf.*', models)) ~ 'rf(ROSE)',
    (grepl('.*rose.*', models) & grepl('.*nnet.*', models)) ~ 'nnet(ROSE)',
    (grepl('.*rose.*', models) & grepl('.*rpart.*', models)) ~ 'rpart(ROSE)',
    (grepl('.*rose.*', models) & grepl('.*qda.*', models)) ~ 'qda(ROSE)',
    (grepl('.*smote.*', models) & grepl('.*rf.*', models)) ~ 'rf(SMOTE)',
    (grepl('.*smote.*', models) & grepl('.*nnet.*', models)) ~ 'nnet(SMOTE)',
    (grepl('.*smote.*', models) & grepl('.*rpart.*', models)) ~ 'rpart(SMOTE)',
    (grepl('.*smote.*', models) & grepl('.*qda.*', models)) ~ 'qda(SMOTE)',
    (grepl('.*smote.*', models) & grepl('.*nb.*', models)) ~ 'nb(SMOTE)',
    (grepl('.*rose.*', models) & grepl('.*nb.*', models)) ~ 'nb(ROSE)',
    (grepl('.*up.*', models) & grepl('.*nb.*', models)) ~ 'nb(upsampled)',
    T ~ models)
  ) %>%
  mutate(fill_val = gsub(pattern = '\\d+', x = .$fill_val, replacement = ''))

#model_exact == 'svmPoly', tox_top$models, tox_top$model_exact)

# tox_top <- tox_top %>%
#   dplyr::filter(!is.na(AUC))

grader_melt <- melt(tox_top, id.vars = c('models', 'dataset', 'dilitype', 'model_exact', 'set', 'training_type', 'fill_val'), 
                    measure.vars = c('AUC', 'sensitivity', 'specificity', 'mcc', 'balanced_accuracy'))
grader_melt$set <- factor(grader_melt$set, levels=c('training set', 'test set'))

#grader_melt$value[is.na(grader_melt$value)] <- 0

grader_melt <- grader_melt %>%
  dplyr::filter(!is.na(value))
grader_melt$dilitype <- toupper(grader_melt$dilitype)
grader_melt$variable <- as.character(grader_melt$variable)
grader_melt$variable[grader_melt$variable == 'balanced_accuracy'] <- 'balanced accuracy'
grader_melt$variable[grader_melt$variable == 'mcc'] <- 'MCC'
grader_melt$variable <- as.factor(grader_melt$variable)

res_supply_labels <- c('Non-resampled', 'Resampled')
names(res_supply_labels) <- c('non_resampled', 'resampled')

set_supply_labels <- c('Training set', 'Test set')
names(set_supply_labels) <- c('training set', 'test set')

lb <- min(grader_melt$value, na.rm = T)

lbls <- grader_melt$fill_val %>% unique()

#grepl('.*smote.*', unique(tox_top$models))

p_tox <- ggplot(grader_melt, aes(variable, value, fill=fill_val)) +
  geom_bar(stat='identity', position = position_dodge(), color='black') +
  facet_grid(dilitype ~ set + training_type, scales = 'free_x', 
             labeller = labeller(training_type = res_supply_labels, set=set_supply_labels)) +
  scale_fill_manual('Model', values = rcartPalette(length(unique(grader_melt$fill_val))), labels=lbls) +
  theme_classic() +
  labs(title='Performance of algorithms using TOX21 data', x=NULL, y='Performance') +
  scale_y_continuous(limits=c(lb, 1), 
                     breaks=c(lb, seq(lb+(0-lb), 1, by = 0.2))) +
  theme(plot.title = element_text(face='bold'),
        axis.text.x = element_text(face='bold', size=10, angle = 90, hjust=0.95, vjust=0.2),
        axis.text.y = element_text(face ='bold', size=10),
        axis.title.x=element_text(face='bold', size=15),
        axis.title.y=element_text(face='bold', size=15),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.position="bottom")

p_tox

# === mold figures ===================
mold_top <- all_set_metrics[all_set_metrics$dataset=='mold', ]

mold_top <- mold_top %>%
  mutate(fill_val = case_when(
    (grepl('.*smote.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(SMOTE)',
                              (grepl('.*rose.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(ROSE)',
                              (grepl('.*up.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(upsampled)',
                              (grepl('.*smote.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(SMOTE)',
                              (grepl('.*rose.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(ROSE)',
                              (grepl('.*up.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(upsampled)',
                              (grepl('.*smote.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(SMOTE)',
                              (grepl('.*rose.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(ROSE)',
                              (grepl('.*up.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(upsampled)',
                              T ~ models)
    ) %>%
  mutate(fill_val = gsub(pattern = '\\d+', x = .$fill_val, replacement = ''))
                              
    #model_exact == 'svmPoly', mold_top$models, mold_top$model_exact)

grader_melt <- melt(mold_top, id.vars = c('models', 'dataset', 'dilitype', 'model_exact', 'set', 'training_type', 'fill_val'), 
                    measure.vars = c('AUC', 'sensitivity', 'specificity', 'mcc', 'balanced_accuracy'))
grader_melt$set <- factor(grader_melt$set, levels=c('training set', 'test set'))

grader_melt <- grader_melt %>%
  dplyr::filter(!is.na(value))
grader_melt$dilitype <- toupper(grader_melt$dilitype)
grader_melt$variable <- as.character(grader_melt$variable)
grader_melt$variable[grader_melt$variable == 'balanced_accuracy'] <- 'balanced accuracy'
grader_melt$variable[grader_melt$variable == 'mcc'] <- 'MCC'
grader_melt$variable <- as.factor(grader_melt$variable)

res_supply_labels <- c('Non-resampled', 'Resampled')
names(res_supply_labels) <- c('non_resampled', 'resampled')

set_supply_labels <- c('Training set', 'Test set')
names(set_supply_labels) <- c('training set', 'test set')

lb <- min(grader_melt$value, na.rm = T)

lbls <- grader_melt$fill_val %>% unique()

#grepl('.*smote.*', unique(mold_top$models))

p_mold <- ggplot(grader_melt, aes(variable, value, fill=fill_val)) +
  geom_bar(stat='identity', position = position_dodge(), color='black') +
  facet_grid(dilitype ~ set + training_type, scales = 'free_x', 
             labeller = labeller(training_type = res_supply_labels, set=set_supply_labels)) +
  scale_fill_manual('Model', values = rcartPalette(length(unique(grader_melt$fill_val))), labels=lbls) +
  theme_classic() +
  labs(title='Performance of algorithms using MOLD2 data', x=NULL, y='Performance') +
  scale_y_continuous(limits=c(lb, 1), 
                     breaks=c(lb, seq(lb+(0-lb), 1, by = 0.2))) +
  theme(plot.title = element_text(face='bold'),
        axis.text.x = element_text(face='bold', size=10, angle = 90, hjust=0.95, vjust=0.2),
        axis.text.y = element_text(face ='bold', size=10),
        axis.title.x=element_text(face='bold', size=15),
        axis.title.y=element_text(face='bold', size=15),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.position="bottom")

p_mold

# === faers figures ===================
faers_top <- all_set_metrics[all_set_metrics$dataset=='faers', ]

faers_top <- faers_top %>%
  mutate(fill_val = case_when(
    (grepl('.*smote.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(SMOTE)',
    (grepl('.*rose.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(ROSE)',
    (grepl('.*up.*', models) & grepl('.*svmPoly.*', models)) ~ 'svmPoly(upsampled)',
    (grepl('.*smote.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(SMOTE)',
    (grepl('.*rose.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(ROSE)',
    (grepl('.*up.*', models) & grepl('.*svmLinear.*', models)) ~ 'svmLinear(upsampled)',
    (grepl('.*smote.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(SMOTE)',
    (grepl('.*rose.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(ROSE)',
    (grepl('.*up.*', models) & grepl('.*svmRadial.*', models)) ~ 'svmRadial(upsampled)',
    (grepl('.*up.*', models) & grepl('.*rf.*', models)) ~ 'rf(upsampled)',
    (grepl('.*up.*', models) & grepl('.*nnet.*', models)) ~ 'nnet(upsampled)',
    (grepl('.*up.*', models) & grepl('.*rpart.*', models)) ~ 'rpart(upsampled)',
    T ~ models)
  ) %>%
  mutate(fill_val = gsub(pattern = '\\d+', x = .$fill_val, replacement = ''))

#model_exact == 'svmPoly', faers_top$models, faers_top$model_exact)

grader_melt <- melt(faers_top, id.vars = c('models', 'dataset', 'dilitype', 'model_exact', 'set', 'training_type', 'fill_val'), 
                    measure.vars = c('AUC', 'sensitivity', 'specificity', 'mcc', 'balanced_accuracy'))
grader_melt$set <- factor(grader_melt$set, levels=c('training set', 'test set'))

grader_melt <- grader_melt %>%
  dplyr::filter(!is.na(value))
grader_melt$dilitype <- toupper(grader_melt$dilitype)
grader_melt$variable <- as.character(grader_melt$variable)
grader_melt$variable[grader_melt$variable == 'balanced_accuracy'] <- 'balanced accuracy'
grader_melt$variable[grader_melt$variable == 'mcc'] <- 'MCC'
grader_melt$variable <- as.factor(grader_melt$variable)

res_supply_labels <- c('Non-resampled', 'Resampled')
names(res_supply_labels) <- c('non_resampled', 'resampled')

set_supply_labels <- c('Training set', 'Test set')
names(set_supply_labels) <- c('training set', 'test set')

lb <- min(grader_melt$value, na.rm = T)

lbls <- grader_melt$fill_val %>% unique()

#grepl('.*smote.*', unique(faers_top$models))

p_faers <- ggplot(grader_melt, aes(variable, value, fill=fill_val)) +
  geom_bar(stat='identity', position = position_dodge(), color='black') +
  facet_grid(dilitype ~ set + training_type, scales = 'free_x', 
             labeller = labeller(training_type = res_supply_labels, set=set_supply_labels)) +
  scale_fill_manual('Model', values = rcartPalette(length(unique(grader_melt$fill_val))), labels=lbls) +
  theme_classic() +
  labs(title='Performance of algorithms using FAERS data', x=NULL, y='Performance') +
  scale_y_continuous(limits=c(lb, 1), 
                     breaks=c(lb, seq(lb+(0-lb), 1, by = 0.2))) +
  theme(plot.title = element_text(face='bold'),
        axis.text.x = element_text(face='bold', size=10, angle=90, hjust=0.95, vjust=0.2),
        axis.text.y = element_text(face ='bold', size=10),
        axis.title.x=element_text(face='bold', size=15),
        axis.title.y=element_text(face='bold', size=15),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.position="bottom")

p_faers

# === expr figures ===================
expr_top <- all_set_metrics[all_set_metrics$dataset=='expr', ]

grader_melt <- melt(expr_top, id.vars = c('models', 'dataset', 'dilitype', 'model_exact', 'set', 'training_type'), 
                    measure.vars = c('AUC', 'sensitivity', 'specificity', 'mcc', 'balanced_accuracy'))
grader_melt$set <- factor(grader_melt$set, levels=c('training set', 'test set'))

grader_melt$model_exact <- gsub(pattern = '\\.\\d+', replacement = '', x = grader_melt$models)
grader_melt$value <- round(grader_melt$value, 3)

grader_melt <- grader_melt %>%
  dplyr::filter(!is.na(value))
grader_melt$dilitype <- toupper(grader_melt$dilitype)
grader_melt$variable <- as.character(grader_melt$variable)
grader_melt$variable[grader_melt$variable == 'balanced_accuracy'] <- 'balanced accuracy'
grader_melt$variable[grader_melt$variable == 'mcc'] <- 'MCC'
grader_melt$variable <- as.factor(grader_melt$variable)

res_supply_labels <- c('Non-resampled', 'Resampled')
names(res_supply_labels) <- c('non_resampled', 'resampled')

set_supply_labels <- c('Training set', 'Test set')
names(set_supply_labels) <- c('training set', 'test set')

lb <- min(grader_melt$value, na.rm = T)

p_expr <- ggplot(grader_melt, aes(variable, value, fill=model_exact)) +
  geom_bar(stat='identity', position = position_dodge(), color='black') +
  facet_grid(dilitype ~ set + training_type, scales = 'free_x', 
             labeller = labeller(training_type = res_supply_labels, set=set_supply_labels)) +
  scale_fill_manual('Model', values = rcartPalette(length(unique(grader_melt$model_exact)))) +
  theme_classic() +
  labs(title='Performance of algorithms using gene expression data', x=NULL, y='Performance') +
  scale_y_continuous(limits=c(lb, 1), 
                     breaks=c(lb, seq(lb+(0-lb), 1, by = 0.2))) +
  theme(plot.title = element_text(face='bold'),
        axis.text.x = element_text(face='bold', size=10, angle =90, hjust=0.95, vjust=0.2),
        axis.text.y = element_text(face ='bold', size=10),
        axis.title.x=element_text(face='bold', size=15),
        axis.title.y=element_text(face='bold', size=15),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.position="bottom")

p_expr

# === mold figures ===================
voting_top <- all_set_metrics[all_set_metrics$dataset=='voting', ]

grader_melt <- melt(voting_top, id.vars = c('dataset', 'dilitype', 'set', 'training_type'), 
                    measure.vars = c('sensitivity', 'specificity', 'balanced_accuracy'))
#grader_melt$set <- factor(grader_melt$set, levels=c('training set', 'test set'))

grader_melt$dilitype <- toupper(grader_melt$dilitype)
grader_melt$variable <- as.character(grader_melt$variable)
grader_melt$variable[grader_melt$variable == 'balanced_accuracy'] <- 'balanced accuracy'
grader_melt$variable <- as.factor(grader_melt$variable)

grader_melt$training_type <- as.character(grader_melt$training_type)
grader_melt$training_type[grader_melt$training_type == 'w-1'] <- 'weighted'
grader_melt$training_type <- as.factor(grader_melt$training_type)

supp_lab <- c('Test set')
names(supp_lab) <- c('test set')

p_voting <- ggplot(grader_melt, aes(variable, value, fill=training_type)) +
  geom_bar(stat='identity', position = position_dodge(), color='black') +
  facet_grid(dilitype ~ set, labeller = labeller(set=supp_lab)) +
  scale_fill_manual('Model', values = rcartPalette(length(unique(grader_melt$training_type)))) +
  theme_classic() +
  labs(title='Ensemble (Voting) performance on predicting DILI types', x=NULL, y='Performance') +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by = 0.2)) +
  theme(plot.title = element_text(face='bold'),
        axis.text.x = element_text(face='bold', size=10, angle = 90, hjust = 0.95, vjust = 0.2),
        axis.text.y = element_text(face ='bold', size=10),
        axis.title.x=element_text(face='bold', size=15),
        axis.title.y=element_text(face='bold', size=15),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.position="bottom")

p_voting


g_voting <- ggpubr::ggdotchart(grader_melt, x='variable', y='value', color = 'training_type',
                palette = c("#00AFBB", "#E7B800", "#FC4E07"), sorting = "descending",                      
                add = "segments", rotate = TRUE, dot.size = 3, 
                label = round(grader_melt$value, 2), repel = T,
                title = 'Ensemble (voting) performance on predicting DILI types', 
                xlab=F, facet.by = 'dilitype', legend.title='Ensemble (voting) type',
                ggtheme = theme_classic())


g_voting <- g_voting + theme(plot.title = element_text(face='bold'),
                 axis.text.x = element_text(face='bold', size=10, angle = 90, hjust = 0.95, vjust = 0.2),
                 axis.text.y = element_text(face ='bold', size=10),
                 axis.title.y=element_text(face='bold', size=15),
                 strip.text.x = element_text(size = 12, face = "bold"),
                 strip.text.y = element_text(size = 12, face = "bold")) +
  labs(x='Metric')

g_voting

# save the figures

figs <- list(faers=p_faers, tox=p_tox, mold=p_mold, voting=p_voting, voting_2=g_voting, expression=p_expr)

lapply(names(figs), function(x){
  file_name <- paste('../output_plots/', x, '_performance.jpeg', sep='')
  print(file_name)
  ggsave(filename = file_name, plot=figs[[x]], width = 8, height = 8)
})

lapply(names(figs), function(x){
  file_name <- paste('../output_plots/', x, '_performance.pdf', sep='')
  print(file_name)
  ggsave(filename = file_name, plot=figs[[x]], width = 8, height = 8, device = 'pdf')
})

dev.off()


# leg <- g_legend(p_voting)
# 
# pdf('/extData/NGS/hurlab/temi/projects/camda/journal_plots/tox_img.pdf', width = 7, height = 8)
# tiff("/extData/NGS/hurlab/temi/projects/camda/journal_plots/tox_img.tiff", units="in", width=7, height=8, res=300)
# gridExtra::grid.arrange(arrangeGrob(tox_training + theme(legend.position="none"),
#                                     tox_testing + theme(legend.position="none"),
#                                     nrow=1), leg, heights=c(10, 1))
# dev.off()

# properly formatting the summary tables

# expr_metrics <- readxl::read_xlsx('../ModelSummary_v5.xlsx', sheet = 'FigTable', skip = 1)
# 
# expr_metrics <- expr_metrics %>%
#   dplyr::select(-c(1))
# 
# names(expr_metrics) <- c('models', 'dilitype', 'algorithm', 'gene_size', 'AUC', 'sensitivity', 'specificity', 'predictors', 'MCC', 
#                          'test_accuracy', 'test_precision', 'test_recall', 'test_specificity', 'test_F1', 'test_MCC')
# expr_metrics$dilitype <- paste('DILI', expr_metrics$dilitype, sep='')
# expr_metrics <- expr_metrics[!expr_metrics$dilitype=='DILINA', ]
# expr_metrics$test_precision <- as.numeric(expr_metrics$test_precision)
# 
# # all_set_metrics
# # 
# # expr_metrics %>%
# #   mutate(across(.cols = c(AUC, accuracy, sensitivity, specificity, MCC, f1, recall, test_accuracy, 
# #                           test_precision, test_recall, test_specificity, test_F1, test_MCC), .fns = ~ round(.x, 2))) 
# 
# expr_metrics <- expr_metrics %>%
#   dplyr::mutate(across(where(is.numeric), ~round(.x, 2)))
# 
# 
# %>% 
#   mutate(test_balanced_accuracy = (test_recall + test_specificity) / 2)

