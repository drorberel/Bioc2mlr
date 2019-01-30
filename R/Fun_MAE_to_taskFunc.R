#' Convert MAE to mlr task
#'
#' @param MAE_obj MAE class
#' @param param.Y.name Vector of dependent variable name
#' @param param.covariates Vector of coaraiate variable(s) name
#' @param param_positive_y_level if ClassifTask, value (character or numeric) to be considered as the positive factor outcome
#' @return mlr's \code{ClassifTask} or \code{RegrTask}
#' @details In case of individual MAE assay (omic) with multiple sub-assays, only first sub-assay will be used.
#' MAE's helpers functions longFormat and wideFormat may not be best candidates, since mlr's functional data require
#' complete subject structure for all 'assays' ('functionals'). this can be achieved by either removing non-complete subjects,
#' or by creating dummy entities with NA. (which can later be imputed if required).
#'
#'Either ClassifTask or RegrTask will be returned, based on the type of the param.Y.name variable

#' @examples
#' data(miniACC, package = 'MultiAssayExperiment') # ExpressionSet
#' miniACC
#' Fun_MAE_to_taskFunc(miniACC, param.Y.name = 'vital_status', param.covariates = c('gender','days_to_death'), param_positive_y_level = '1')
#'
#' @export



Fun_MAE_to_taskFunc<-function(MAE_obj, param.Y.name, param.covariates, param_positive_y_level, task_type){
  # i=4; MAE = MAE_DF$MAE_selected[[i]]; param.Y.name = MAE_DF$param.Y.name[[i]]; param.covariates = MAE_DF$param.covariates[[i]]
  # MAE_obj<-miniACC; param.Y.name = 'vital_status'; param.covariates = c('gender','days_to_death'); param_positive_y_level = '1'
  # sampleMap(MAE) %>% data.frame %>% pull(primary) %>% table
  # task_type = 'classif'

  # MAE_obj = myMultiAssay; param.Y.name = 'sex'; param.covariates = NULL; param_positive_y_level = 'M'
  # i=4; MAE_obj = Meta_tibble$MAE[[i]]; param.Y.name = 'study_product_class_combination';  param.covariates = NULL; param_positive_y_level = Meta_tibble$param_positive_y_level_V[[i]]

  ## A. Extract data from SE: colData and assay(s)

  DF_ColData<-MAE_obj %>% colData %>% data.frame %>% rownames_to_column('primary')# %>% dplyr::slice(1:10)# %>% select(-!!param.Y.name)

  sample_map<-MAE_obj %>% sampleMap %>% data.frame %>%
    unite('Unique_sample_id', c('assay','primary','colname'), sep='_', remove = FALSE) %>%
    left_join(DF_ColData, by = 'primary')


  DF_functionals<-sample_map[,c('Unique_sample_id', 'assay', 'primary', 'colname', param.covariates, param.Y.name), drop = FALSE]
  # DF_functionals$target<-DF_ColData$cpep_model_decayrate

  ## all non-numeric coariates must be factors!
  DF_functionals %<>% mutate_if(is.character, as.factor)
  # DF_functionals %>% str

  DF_exprsS<-assays(MAE_obj) %>% as.list %>% map(t) ## assume single assay within each SE
  # str(DF_exprsS)
  # DF_exprsS %>% map(dim)
  # DF_exprsS %>% map(rownames) %>% gplots::venn()


  ## makeFunctionalData() require only complete cases.
  ## instead of removing non-overlap subjects, fill in dummy NAs: (rowbind (fill) dropouts with NAs)
  DF_exprsS_completeNA<-DF_exprsS %>% map(function(x){
    # x = DF_exprsS[[3]][,1:3]

    ## replace (assay) colname with 'Unique_sample_id'
    DF_exprsS_i<-x %>% data.frame %>% rownames_to_column('colname') %>%
      left_join(sample_map %>% select(colname, Unique_sample_id), by = 'colname') %>%
      select(-colname)

    ## add dummy NA columns for each assay (by add_row())
    if( (setdiff(DF_functionals$Unique_sample_id %>% as.character,
                 DF_exprsS_i$Unique_sample_id) %>% length) > 0){
      DF_exprsS_i<-DF_exprsS_i %>%
        add_row(Unique_sample_id = setdiff(DF_functionals$Unique_sample_id %>% as.character,
                                           .$Unique_sample_id))
    } else {DF_exprsS_i<-DF_exprsS_i}# end if
  })
  # DF_exprsS_completeNA %>% str
  # DF_exprsS_completeNA %>% map(dim)




  ## Add assays to DF_functionals (currently it has only the covariates)
  ## not sure how to do it alternatively to 'for'
  for (i in 1: length(MAE_obj %>% experiments) ){ # for each assay
    # i=1
    DF_functionals[,names(DF_exprsS_completeNA)[[i]]]<-DF_exprsS_completeNA[[i]] %>% column_to_rownames('Unique_sample_id') %>% as.matrix
  }
  # str(DF_functionals)


  ## B. Returned task. Functional

  ## Classif:
  if( task_type == 'classif'){
    y_levels<-DF_functionals[,param.Y.name] %>% factor %>% levels
    if( ( y_levels  %>% length )==1 ) cat('Error: target has a single level')

    DF_functionals[,param.Y.name]<-DF_functionals[,param.Y.name] %>% as.factor
    task_MAE_functionals<-makeClassifTask(data = DF_functionals, target = param.Y.name, positive = param_positive_y_level) # with covariates!
  }

  ## Regr
  if( task_type == 'regr'){
    DF_functionals[,param.Y.name] %<>% as.numeric
    task_MAE_functionals<-makeRegrTask(data = DF_functionals, target = param.Y.name) # with covariates!
  }

  return(task_MAE_functionals)
}
