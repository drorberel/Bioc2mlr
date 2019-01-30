#' Helper post_processing functions
#' @param all_kind all kind...
#'
#' @name post_processing_helpers
NULL

#' @return all kind...
#' @details TBA






#' @export
#' @rdname post_processing_helpers

# whole run mlr wrapper ----------------------------------------------------

whole_run_function<-function(task, lrn, lambda_opt, index){
  # task = task_j; lrn = lrn_1; lambda_opt = 'lambda.min'
  # task = task_j; lrn = bmr_tib$lrn_i[[1]]; lambda_opt = 'pushed_lambda'
  # lrn$par.vals
  # print(index)

  ## 1. Prepare data:
  data_to_bake  <-task %>% getTaskData()
  target_to_bake<-task %>% getTaskTargetNames()

  ## 2. Run pre-processing: Univariate + clustering feature selection
  baked_UnivClust<-F_PreProc_3_UnivClust_Train_MaG(data_to_bake, target_to_bake, args = lrn$par.vals) # ok to do with F_PreProc_3_UnivClust_Train_MaG() since no further ML is done

  ## 3. Extract selected features
  task_baked_UnivClust<-makeRegrTask(id = "baked_passed_univ", data = baked_UnivClust$data, target = "cpep_model_decayrate")
  #Model_UnivClust<-train(lrn, task) ## lrn will have defaults s, that is used to extract/predict specific coef, though all s values are used for calculation
  #alt_lambda<-Model_UnivClust$learner.model$control$alt_lambda[lambda_opt] %>% unlist
  alt_lambda  <-baked_UnivClust$control$alt_lambda[lambda_opt] %>% unlist
  if(!is.numeric(alt_lambda)) break

  # UnivOnlyCoef: will be extracted from saved diagnostic table in control.
  # Alternatively, could have been created independently using F_PreProc_1_UnivOnly_Train_MaG
  UnivOnlyCoef<-baked_UnivClust$control$Features_Master_all %>% filter(passed_univ_fw_abs) %>% pull(assay_i_gene) ## !! does not include coeffiecients, which are not part of the pre-processing feature selection. if considered, will be forced into final glmnet

  TargetName<-task %>% getTaskTargetNames
  Baked_x<-baked_UnivClust$data %>% select(-TargetName) %>% as.matrix
  Baked_y<-baked_UnivClust$data %>% pull(TargetName)

  ## 4.final feature selection via LASSO (glmnet)
  fit_baked<-glmnet(Baked_x, Baked_y, alpha=1)

  ## 5. selected model (coefficients)
  ## to do: allow user to select default (non-manually controlled) s from original learner!
  # alt_lambda<-lrn$next.learner$par.vals$s
  Coef_all<-coef(fit_baked, s = alt_lambda) %>% tidy # alternatively, could be taken from the above baked directly

  ## 6. model prediction (same training dataset, but imputed)
  response_direct_glmnet<-predict(fit_baked, Baked_x, s = alt_lambda) %>% data.frame %>% pull(X1) # done directly with glmnet, rather than mlr's since mlr's predict produce NA's for unknown reason

  ## 7. performance
  rmse_direct<-sqrt(mean((response_direct_glmnet - Baked_y) ^ 2))

  return(list(Coef_all = Coef_all, rmse_direct = rmse_direct, # initial
              fit_baked = fit_baked, alt_lambda_list = baked_UnivClust$control$alt_lambda, Baked_x = Baked_x, Baked_y = Baked_y,
              UnivOnlyCoef = UnivOnlyCoef)) # more generic
}
# out<-whole_run_function( task = task_j, lrn = lrn_1)






#' @export
#' @rdname post_processing_helpers


## If mlr's predict WOULD have NOT produce NA's:
whole_run_MLR_function<-function(task, lrn, adjust_alt_lambda){
  # task = task_j; lrn = lrn_1
  Model_UnivClust<-train(lrn, task) ## lrn will have defaults s, that is used to extract/predict specific coef, though all s values are used for calculation
  alt_lambda<-Model_UnivClust$learner.model$control$alt_lambda

  # coef:
  model_3_glmnet<-Model_UnivClust %>% getLearnerModel %>% getLearnerModel
  Coef_all<-coef(model_3_glmnet, s = alt_lambda) %>% tidy

  # prediction:
  Model_UnivClust_my_s<-Model_UnivClust # just to avoid over-ride
  if(adjust_alt_lambda) Model_UnivClust_my_s$learner$next.learner$par.vals$s<-alt_lambda

  Pred_UnivClust_my_s<-predict(Model_UnivClust_my_s, task_j) # self feed. overfitted
  rmse_mlr<-Pred_UnivClust_my_s %>% performance # with default 'measures' pre set parameters

  return(list(Coef_all = Coef_all, rmse_mlr = rmse_mlr, model = Model_UnivClust))
}

# out2<-whole_run_MLR_function( task = task_j, lrn = lrn_1)
