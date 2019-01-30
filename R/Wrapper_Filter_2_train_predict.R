#' Customized train and predict functions for a customized preprocessing wrapper
#'
#' @param data data
#' @param target target (character)
#' @param args args (list)
#' @param control control (predict only)
#' @name Wrapper_Filter_2_train_predict
NULL



#' @return The sum of \code{x} and \code{y}
#' @details In
#' @examples
#' F_PreProc_1_Univ_Only_Train_MaG(data, target, args)
#' F_PreProc_3_UnivClust_Train_MaG(data, target, args)
#' F_PreProc13_BOTH_Predict_MaG(data, target, args, control)



#' @export
#' @rdname Wrapper_Filter_2_train_predict

# PrePrcocess 1: Univ only (Module/gene)
F_PreProc_1_Univ_Only_Train_MaG<-function(data, target,
                                          args = list(param.Univ.filt.top.n.features   = param.Univ.filt.top.n.features,
                                                      parame.gene.or.module            = parame.gene.or.module,
                                                      param.LASSO.n.features.arbitrary = param.LASSO.n.features.arbitrary)) {

  # data = data_to_bake; target = target_to_bake; args = args_UnivOnly
  ## task = task_i; data = task %>% getTaskData(); target = task %>% getTaskTargetNames(); args = args_UnivOnly
  ## args = list(param.Univ.filt.top.n.features = param.Univ.filt.top.n.features, parame.gene.or.module = parame.gene.or.module )
  ## args = lrn_1$par.vals


  ## Step_1: Decide which task type to create:
  # If Classif :
  if( data %>% pull(!!target) %>% levels %>% length ==2 ) task.int<-makeClassifTask(data = data, target = target, positive='1')
  # If Regr:
  if( data %>% pull(!!target) %>% levels %>% is.null)     task.int<-   makeRegrTask(data = data, target = target)


  ##############################################################################################################
  ## Step_2: Run Univaraite (covariates will be adjusted for, but will NOT be returned)
  ## module/gene level is handels internally, params is passed via args
  Features_and_DF_L<-Func_task_to_Passed_Features_DF_MaG(task.int, args) # Covariates are used to adjust univariate analysis, but not returned!
  ## Features_and_DF_L has two elements:
  ##  $Features_Master_all is used for feature ranking within clustering, as reference list
  ##  $DF_passed is used for the actual filtering
  ##############################################################################################################

  ## Step_3: features master (for cluster purposes):
  # Features_and_DF_L$Features_Master_all


  ## Step_4: dataOUT: (non-imputed!)
  dataOUT<-Features_and_DF_L$DF_passed
  # str(dataOUT)


  ### Step_5: PRE-Impute/PRE-clustering: (crucial for clustering)

  ## 5_A. PRE-imput: remove SUBJECTS with 100% NAs
  crit_1_subjects<-dataOUT %>% t %>% data.frame %>% map(~!all(is.na(.x))) %>% unlist
  table(crit_1_subjects)
  dataOUT<-dataOUT[crit_1_subjects,]

  ## 5_B. PRE-clustering: remove FEATURES with sd=0 (can be done as early as here. should not affect knn.impute)
  crit_2_features<-dataOUT %>% map(sd, na.rm=TRUE) %>% unlist!=0
  table(crit_2_features)
  dataOUT<-dataOUT[, crit_2_features]


  ## 5_C. imputation: knn  # input for knn function has 'genes' in rows.
  # visdat::vis_dat(dataOUT)
  # dataOUT %>% is.na %>% sum
  dataOUT<-impute.knn(dataOUT %>% t, param.impute.knn.k, colmax=0.99999, rowmax=0.9)$data %>% t %>% data.frame
  # visdat::vis_dat(dataOUT)





  ### Wrap up:

  ## Step_6: dataOUT

  ## 6_1: Push back covaraites: if exist, were seperated within Func_task_to_FeaturesMaster(), but were NOT return at the end
  Covariates_names<-task.int %>% getTaskFeatureNames() %>% .[str_sub(.,1, 9)=='Covariate']
  if(length(Covariates_names)>0) {
    #Covariates_names_edit<-gsub(paste0('Covariate', Assay.Analyte.sep), 'Covariate_', Covariates_names ) ## return original covaraites names
    dataOUT %<>% data.frame(.,  task.int %>% subsetTask(feature = Covariates_names) %>% getTaskData(target.extra = TRUE) %>% .$data)
  }


  ### 6_2: run preliminary cv.glmnet to find best lambda if too aggressive
  #%>% table
  preliminary_lasso_detect_aggressive_CV<-cv.glmnet(
    x = dataOUT %>% as.matrix,
    y = task.int %>%
      subsetTask(subset = (task.int$env$data %>% rownames) %in% (dataOUT %>% rownames))%>% ## not sure why I initially had task subset, if only target (y) is used. test, and remove if not required
      getTaskTargets(),
    alpha=1,
    nfolds = task.int %>% getTaskSize)
  # plot(preliminary_lasso_detect_aggressive_CV)
  alt_lambda<-cvglmnet_alternative_lambda.min_function(preliminary_lasso_detect_aggressive_CV, args$param.LASSO.n.features.arbitrary)
  #glmnet_sanity<-glmnet(dataOUT %>% as.matrix,task.int %>% getTaskTargets(), alpha=1)



  ## 6_3: Control
  control<-list(
    feature_names       = colnames(dataOUT), # this will NOT include target!
    alt_lambda          = alt_lambda,
    Features_Master_all = Features_and_DF_L$Features_Master_all
  ) # control must be a list

  ## 6_4: add outcome / target ## Note: dataOUT in predict() should INCLUDE the Targets variables!!
  dataOUT %<>% add_column(!!(task.int %>% getTaskTargetNames()) :=task.int %>% getTaskTargets() )


  # str(dataOUT)
  # task.int %>% getTaskFeatureNames() %>% tail()
  return(list(data = dataOUT, control = control)) # control must be a list
}



#' @export
#' @rdname Wrapper_Filter_2_train_predict

# PrePrcocess 3: Univ+clustering (Module/gene)

F_PreProc_3_UnivClust_Train_MaG<-function(data, target,
                                          args = list(
                                            param.Univ.filt.top.n.features   = param.Univ.filt.top.n.features,
                                            param.UnivClustRankTopN          = param.UnivClustRankTopN,
                                            param.cluster_method_KH          = param.cluster_method_KH,
                                            param.corrplot.n.clusters.k      = param.corrplot.n.clusters.k,
                                            param.corrplot.n.clusters.h      = param.corrplot.n.clusters.h,
                                            parame.gene.or.module            = parame.gene.or.module,
                                            param.LASSO.n.features.arbitrary = param.LASSO.n.features.arbitrary
                                          ) ) {
  # data = data_to_bake; target = target_to_bake; args = lrn_1$par.vals #args_UnivClust
  ## task=TaskAB_blocked; data=task %>% getTaskData(); target=task %>% getTaskTargetNames()

  ## Step_1: Setup: Decide which task type to create:
  # Classif :
  if( data %>% pull(!!target) %>% levels %>% length ==2 ) task.int<-makeClassifTask(data = data, target = target, positive='1')
  # Regr:
  if( data %>% pull(!!target) %>% levels %>% is.null)     task.int<-makeRegrTask(data = data, target = target)

  ##############################################################################################################
  ## Step_2: Run Univaraite (covariates will be adjusted for, but will NOT be returned)
  # task.int<-task; args = list(param.Univ.filt.top.n.features= param.Univ.filt.top.n.features, param.UnivClustRankTopN = param.UnivClustRankTopN, param.cluster_method_KH = param.cluster_method_KH, param.corrplot.n.clusters.k = param.corrplot.n.clusters.k, param.corrplot.n.clusters.h = param.corrplot.n.clusters.h, parame.gene.or.module = parame.gene.or.module)
  Features_and_DF_L<-Func_task_to_Passed_Features_DF_MaG(task.int, args) # Covariates are used to adjust univariate analysis, but not returned!
  ## Features_and_DF_L has two elements:
  ##  $Features_Master_all is used for feature ranking within clustering, as reference list
  ##  $DF_passed is used for the actual filtering

  Features_Master_all_P3<-Features_and_DF_L$Features_Master_all ## features master (for cluster purposes):
  # Features_Master_all_P3 %>% filter(assay_i_gene=='FACS.ZZZ.Pct_Senescent_CM_CD8')
  dataOUT_1<-Features_and_DF_L$DF_passed ## dataOUT: (non-imputed)
  ##############################################################################################################
  # Features_Master_all_P3$passed_univ_fw_abs %>% table
  # table(Features_Master_all_P3$assay_i, Features_Master_all_P3$passed_univ_fw_abs)

  ### Step_3: PRE-Impute/PRE-clustering: (crucial for clustering)

  ## 3_A. PRE-imput: remove SUBJECTS with 100% NAs
  dataTemp<-dataOUT_1 # task.int %>% subsetTask(features=Features_Master_Passed_Univ$gene) %>% getTaskData(target.extra = TRUE) %>% .$data

  crit_1_subjects<-dataTemp %>% t %>% data.frame %>% map(~!all(is.na(.x))) %>% unlist
  table(crit_1_subjects)
  dataTemp<-dataTemp[crit_1_subjects,]

  ## 3_B. PRE-clustering: remove FEATURES with sd=0 (can be done as early as here. should not affect knn.impute)
  crit_2_features<-dataTemp %>% map(sd, na.rm=TRUE) %>% unlist!=0
  table(crit_2_features)
  dataTemp<-dataTemp[, crit_2_features]

  ## 3_C. imputation: knn  # input for knn function has 'genes' in rows.
  # visdat::vis_dat(dataTemp)
  # dataTemp %>% is.na %>% sum
  dataTemp<-impute.knn(dataTemp %>% t, param.impute.knn.k, colmax=0.99999, rowmax=0.9)$data %>% t %>% data.frame
  # visdat::vis_dat(dataTemp)


  ### Step_4: Clusters, no filtering: (but for only the features that passed univ)

  ### 4_A. cluster raw:
  hclust_orig<-hclust(as.dist(1-cor(dataTemp, use='na.or.complete')))

  ### 4_B: cut tree
  # args$param.cluster_method_KH='method.k'; args$param.corrplot.n.clusters.k=100
  # args$param.cluster_method_KH='method.h'; args$param.corrplot.n.clusters.h=0.5
  if( length(hclust_orig$order) < args$param.corrplot.n.clusters.k) args$param.corrplot.n.clusters.k<-length(hclust_orig$order)
  if(args$param.cluster_method_KH=='method.k') cut_clusters<-cutree(hclust_orig, k = args$param.corrplot.n.clusters.k) #  used to be min(args$param.corrplot.n.clusters.k, nrow(dataTemp))
  if(args$param.cluster_method_KH=='method.h') cut_clusters<-cutree(hclust_orig, h = args$param.corrplot.n.clusters.h)
  # table(cut_clusters)
  # cut_clusters %>% length()

  ### 4_C: Update Feature_Master with cluster id
  # Features_Master_all_P3$gene[which( str_detect( Features_Master_all_P3$gene, 'RNAseq.Wholeblood..ZZZ. cell.adhesion') )]
  Features_Master_all_P3$gene<-Features_Master_all_P3$gene %>% make.names()
  # Features_Master_all_P3$gene[which( str_detect( Features_Master_all_P3$gene, 'RNAseq.Wholeblood..ZZZ..cell.adhesion') )]
  Features_Master_all_P3<-left_join(Features_Master_all_P3,
                                    cut_clusters %>% data.frame %>% rownames_to_column %>% setNames(c('name', 'cut_clusters_ordered')) ,
                                    by=c('assay_i_gene'='name') ) %>% ### actually a left join is enough. add NA for the features that did NOT passed univ
    group_by(cut_clusters_ordered) %>%
    mutate(rank.within.cluster.by.p = min_rank(p.value))
  # Features_Master_all_P3$cut_clusters_ordered %>% table

  # Features_Master_all_P3 %>% select(cut_clusters_ordered, rank.within.cluster.by.p)
  Features_Master_all_P3$rank.within.cluster.by.p[which(!Features_Master_all_P3$passed_univ_fw_abs)]<-NA
  # Features_Master_all_P3 %>% filter(cut_clusters_ordered==2) %>% arrange(rank.within.cluster.by.p) %>% select(-c(1:3))
  # Features_Master_all_P3 %>% filter(!is.na(cut_clusters_ordered))



  ## Step_5: Filter 2: keep only Cluster rank=1
  # Features_Master_all_P3$passed_UnivCLustRank1_fw_abs
  Features_Master_all_P3 %<>%
    mutate(passed_UnivCLustRank1_fw_abs=rank.within.cluster.by.p %>% map_lgl(~.x <= args$param.UnivClustRankTopN))

  Features_Master_Passed_ClusterRank<-Features_Master_all_P3 %>% filter(passed_UnivCLustRank1_fw_abs) %>% pull(assay_i_gene)
  # Features_Master_Passed_ClusterRank %>% select(c(8:10))
  # dataOUT_Non_Imputed<-task.int %>% subsetTask(features=Features_Master_Passed_ClusterRank$gene) %>% getTaskData()
  # names(dataTemp)
  dataOUT_3<-dataTemp %>% select(Features_Master_Passed_ClusterRank %>% make.names) #%>% add_column(!!(task.int %>% getTaskTargetNames()) :=task.int %>% getTaskTargets() )
  # names(dataOUT_3)





  ### Step_7: Wrap up:

  ## 7_A: Covariates: if exist, push back into dataOUT_3
  Covariates_names<-task.int %>% getTaskFeatureNames() %>% .[str_sub(.,1, 9)=='Covariate']
  if(length(Covariates_names)>0) {
    # Covariates_names_edit<-gsub(paste0('Covariate', Assay.Analyte.sep), 'Covariate_', Covariates_names ) ## return original covaraites names
    dataOUT_3 %<>% data.frame(.,  task.int %>% subsetTask(feature = Covariates_names) %>% getTaskData(target.extra = TRUE) %>% .$data)
  }


  ### 7_B: run preliminary cv.glmnet to find best lambda if too aggressive
  preliminary_lasso_detect_aggressive_CV<-cv.glmnet(
    x      = dataOUT_3 %>% as.matrix,
    y      = task.int %>% subsetTask(subset = (task.int$env$data %>% rownames) %in% (dataOUT_3 %>% rownames)) %>% getTaskTargets(),
    alpha  = 1,
    nfolds = task.int %>% getTaskSize)
  # plot(preliminary_lasso_detect_aggressive_CV)
  alt_lambda<-cvglmnet_alternative_lambda.min_function(preliminary_lasso_detect_aggressive_CV, args$param.LASSO.n.features.arbitrary)
  # glmnet_sanity<-glmnet(dataOUT_3 %>% as.matrix,task.int %>% getTaskTargets(), alpha=1)


  ## 7_C: If SUBJECT were removed (100% NA), also remove them from task (targets and covariates):
  Subjects_kept<-getTaskSubjectId(task.int) %in% rownames(dataOUT_3)
  if(getTaskSize(task.int) != nrow(dataOUT_3)) task.int %<>% subsetTask(subset = Subjects_kept)


  ## 7_D: Allocate all features+covariates to be passed to predict()
  control<-list(
    feature_names       = colnames(dataOUT_3), # this will NOT include target!
    alt_lambda          = alt_lambda,
    Features_Master_all = Features_Master_all_P3
  ) # control must be a list


  ## 7_E: Push back targets ## Note: dataOUT in predict() should INCLUDE the Targets variables!!
  dataOUT_3 %<>% add_column(!!(task.int %>% getTaskTargetNames()) :=task.int %>% getTaskTargets() )



  return(list(data = dataOUT_3, control = control))
}










#' @export
#' @rdname Wrapper_Filter_2_train_predict


# PreProc1_3: Predict, gene/Module, BOTH Univ only, AND Univ+Clust
F_PreProc13_BOTH_Predict_MaG<-function(data, target, args, control){
  #                 data = data_to_bake; target = target_to_bake; args = args_UnivOnly; control = baked_UnivOnly$control
  # args = args_UnivClust; control = baked_UnivClust$control
  ## task = task_i; data = task %>% getTaskData(target.extra = TRUE) %>% .$data; control = baked_UnivOnly$control

  # from Ext/CAV validation
  # data = data_to_bake; target = target_to_bake; control=Results.CAV.EXT.CV_1$features[[1]]; control = My_control

  # data=tasks_tib$data_to_bake[[1]]; target=tasks_tib$target_to_bake[[1]]; control=list(feature_names = Ref_feature_cohort$selected_features)

  # data include covariates, but does NOT include the targets (unknown for the predict!)


  ## 1: seperate covaraites: if exist, were seperated within Func_task_to_FeaturesMaster(), but were NOT return at the end
  if(args$parame.gene.or.module=='gene'){

    ## bypass: instead of applying same feature-selection as in Train (UnivOnly or UnivClust), take directly the names of the selected features from the pre-saved control slot.
    ## However, imputation is still required
    dataOUT_Predict<-data[, control$feature_names]# at gene level, do nothing. return only matched features, BUT exclude Y
    dataOUT_Predict<-impute.knn(dataOUT_Predict %>% t, param.impute.knn.k, colmax=0.99999, rowmax=0.9)$data %>% t %>% data.frame
  }

  if(args$parame.gene.or.module=='module'){
    ## ! add impute step as was done at gene level ! It is handled within nested function Func_Collapse_Data_by_module_names()

    ## step_1: remove covariates:
    dataOUT_Predict<-data
    Covariates_names<-colnames(dataOUT_Predict) %>% .[str_sub(.,1, 9)=='Covariate']
    if( length(Covariates_names)>0 ) {
      # dataOUT_Predict %<>% select(-one_of(!!Covariates_names)) # too slow for big datasets!!
      dataOUT_Predict %<>% .[, ! colnames(.) %in% Covariates_names]
      control$feature_names %<>% .[! . %in% Covariates_names]
    }
    # dim(dataOUT_Predict)

    ## Step_2: remove target (it also doesn't have to be pushed back at the end!)
    dataOUT_Predict %<>% .[, ! colnames(.) %in% target]

    ##############################################################################################################
    ## Step_3: collapse genes to modules, but ONLY to modules selected in train part that passed via control$ (covariates not included, but will be merged back later)
    ## make sure netierh covariates, nor target variables are at the passed data!!!
    # Train_feature_names<-modules_sets_combined_l %>% names %>% make.names %>% paste0('RNAseq.Wholeblood.ZZZ.')
    dataOUT_Predict<-Func_Collapse_Data_by_module_names(dataOUT_Predict, control$feature_names)
    ##############################################################################################################

    ## Step_4: features master (for cluster purposes):
    # Features_and_DF_L$Features_Master_all

    ## Step_5: dataOUT: (non-imputed!)
    # Data_Ncols<-Train_feature_names %>% length()
    # Data_Nrows<-nrow(data)
    # dataOUT_rand<-matrix(runif(Data_Ncols * Data_Nrows), Data_Nrows, Data_Ncols) %>% data.frame %>% set_colnames(Train_feature_names)
    # dataOUT_Predict<-data
    # dataOUT_Predict<-dataOUT_rand

    # dataOUT<-Features_and_DF_L$DF_passed


    ## 5_2: Push back covaraites: if exist, were seperated within Func_task_to_FeaturesMaster(), but were NOT return at the end
    if(length(Covariates_names)>0) {
      Covariates_names_edit<-gsub(paste0('Covariate', Assay.Analyte.sep), 'Covariate_', Covariates_names ) ## return original covaraites names
      dataOUT_Predict %<>% data.frame(.,  data[,Covariates_names_edit] )
    }
    # dim(dataOUT_Predict)
  } # end module if

  # task.int %>% getTaskFeatureNames() %>% tail()
  return(dataOUT_Predict) # return data.frame only
}
