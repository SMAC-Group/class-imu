#' @export
comb.mat = function(n){
  c = rep(list(1:0), n)
  expand.grid(c)
}

#' @export
param_transform = function(model){

  np = model$plength

  # Initialise counter
  counter = 1

  for (j in 1:np){
    # is random walk?
    if (model$process.desc[j] == "RW"){
      model$theta[counter] = exp(model$theta[counter])
      counter = counter + 1
    }

    # is white noise?
    if (model$process.desc[j] == "WN"){
      model$theta[counter] = exp(model$theta[counter])
      counter = counter + 1
    }

    # is drift?
    if (model$process.desc[j] == "DR"){
      model$theta[counter] = exp(model$theta[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (model$process.desc[j] == "QN"){
      model$theta[counter] = exp(model$theta[counter])
      counter = counter + 1
    }

    # is AR1?
    if (model$process.desc[j] == "AR1"){
      model$theta[counter] = transform_phi(model$theta[counter])
      counter = counter + 1
    }

    # is SIGMA2?
    if (model$process.desc[j] == "SIGMA2"){
      model$theta[counter] = exp(model$theta[counter])
      counter = counter + 1
    }
  }
  model$theta
}

#' @export
inv_param_transform = function(model,starting_value){

  np = model$plength

  # Initialise counter
  counter = 1

  for (j in 1:np){
    # is random walk?
    if (model$process.desc[j] == "RW"){
      starting_value[counter] = log(starting_value[counter])
      counter = counter + 1
    }

    # is white noise?
    if (model$process.desc[j] == "WN"){
      starting_value[counter] = log(starting_value[counter])
      counter = counter + 1
    }

    # is drift?
    if (model$process.desc[j] == "DR"){
      model$theta[counter] = log(starting_value[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (model$process.desc[j] == "QN"){
      starting_value[counter] = log(starting_value[counter])
      counter = counter + 1
    }

    # is AR1?
    if (model$process.desc[j] == "AR1"){
      starting_value[counter] = inv_transform_phi(starting_value[counter])
      counter = counter + 1
    }

    # is SIGMA2?
    if (model$process.desc[j] == "SIGMA2"){
      starting_value[counter] = log(starting_value[counter])
      counter = counter + 1
    }
  }
  starting_value
}

#' @export
model_combination = function(model_max){

  # String description of model
  model_desc_max = model$desc

  # number of latent process in model max

  n_process_max = length(model_desc_max)


  # Build matrix of possible combination
  m = as.matrix(comb.mat(n_process_max))
  m = m[-nrow(m),]

  models_names = build_model_set(m,model_desc_max)

  n_models =length(models_names)

  all_model = list()


  for (i in 1:n_models){

    model_test = model

    model_test$desc = models_names[[i]]

    n_process = length(model_test$desc)

    model_test$starting = TRUE


    n_para = 0
    for (j in 1: n_process){
      if(model_test$desc[[j]] == "AR1"){
        n_para = n_para +  2
      }else{
        n_para = n_para + 1
      }
      model_test$plength =  n_para
    }

    process_desc_change = rep(NA,model_test$plength)
    theta_test = rep(NA,model_test$plength)
    obj_desc_list = list()

    counter = 1

    for (j in 1: n_process){
      if (model_test$desc[[j]] == "AR1"){
        process_desc_change[counter] = "AR1"
        process_desc_change[counter + 1] = "SIGMA2"
        theta_test[counter] = 0
        theta_test[counter + 1] = 1
        obj_desc_list[[j]] = c(1,1)

        counter = counter + 2
      }

      if (model_test$desc[[j]] == "WN"){
        process_desc_change[counter] = "WN"
        theta_test[counter] = 3
        obj_desc_list[[j]] = 1
        counter = counter + 1
      }

      if (model_test$desc[[j]] == "RW"){
        process_desc_change[counter] = "RW"
        theta_test[counter] = 4
        obj_desc_list[[j]] = 1
        counter = counter + 1
      }

      if (model_test$desc[[j]] == "QN"){
        process_desc_change[counter] = "QN"
        theta_test[counter] = 2
        obj_desc_list[[j]] = 1
        counter = counter + 1
      }

      if (model_test$desc[[j]] == "DR"){
        process_desc_change[counter] = "DR"
        theta_test[counter] = 5
        obj_desc_list[[j]] = 1
        counter = counter + 1
      }
      model_test$process.desc =  process_desc_change
      model_test$theta =  theta_test
      model_test$obj.desc =  obj_desc_list
    }
    all_model[[i]] = model_test
  }
  all_model
}

#' @export
model_selection = function(mimu, model, s_test = s_test, test_pval = FALSE){

  # model_max must be an object of type model

  # mimu must be an mimu obj

  # Number of replicates
  n_replicates = length(mimu)

  # Number of replicates fot compute the CV_WVIC

  s_valid = n_replicates - s_test

  # Matrix of possible combination of Time series to estimate
  pair = t(combn(n_replicates,s_test))

  # Number of possible combination
  n_replicates_permutation = (dim(pair)[1])

  # Create model object with all possible nested model in model_max
  model_est = model_combination(model_max = model)

  model_test = list()

  n_models = length(model_est)

  cv_wvic = rep(NA,n_models)

  obj_out_sample = matrix(NA,n_replicates_permutation, n_models)

  for (i in 1:n_models){
    desc = model_est[[i]]$desc

    np = model_est[[i]]$plength

    obj_desc = model_est[[i]]$obj.desc

    theta = model_est[[i]]$theta

    starting = model_est[[i]]$starting

    set.seed(2710)

    for (d in 1:n_replicates_permutation){

      mimu_est = list()
      class(mimu_est) = "mimu"
      mimu_test = mimu[-pair[d,sequence(s_test)]]

      for (s in 1:s_test){
        mimu_est[[s]] = mimu[[pair[d,s]]]
      }

      if(d == 1){
        para_gmwm = matrix(NA,np,s_test)
        N = rep(NA,s_test)
        for (k in 1:s_test){
          N[k] = length(mimu_est[[k]]$data)
          data = mimu_est[[k]]$data

          # Select G

          if(N[k] > 10000){
            G = 1e6
          }else{
            G = 20000
          }

          uni_gmwm = .Call('gmwm_gmwm_master_cpp', PACKAGE = 'gmwm', data, theta, desc, obj = obj_desc,
                           model.type = 'imu' , starting = model_est[[i]]$starting,
                           p = 0.05, compute_v = "fast", K = 1, H = 100, G = G,
                           robust=FALSE, eff = 1)
          para_gmwm[,k] = uni_gmwm[[1]]
        }
        starting_value = apply(para_gmwm, 1, mean)

        starting_value = inv_param_transform(model_est[[i]],starting_value)
      }

      out = optim(starting_value, mgmwm_obj_function, model = model_est[[i]], mimu = mimu_est)

      # transform e.g. variance = exp(out$par[1])

      model_test[[i]] = model_est[[i]]

      # Pass on the estimated paramters onto the model.
      model_test[[i]]$starting = FALSE
      model_test[[i]]$theta = out$par

      obj_out_sample[d,i] = mgmwm_obj_function(model_test[[i]]$theta, model_test[[i]], mimu_test)/s_valid

      model_test[[i]]$theta = param_transform(model_test[[i]])
    }
    cv_wvic[i] = mean(obj_out_sample[,i])
  }

  # Compute the average on all permutation of the out of sample WVIC
  mod_selected_cv = which.min(cv_wvic)

  ############ Just add if statement for test_pval

  if(test_pval == TRUE){
    wilcox_test_cv_wvic = rep(FALSE,n_models)

    model_selected_cv_size = model_test[[mod_selected_cv]]$plength

    for (i in 1:n_models){
      if(model_selected_cv_size > model_test[[i]]$plength){
        wilcox_test_cv_wvic[i] = wilcox.test(obj_out_sample[,i],obj_out_sample[,mod_selected_cv],paired = T,alternative = "greater")$p.val
      }
    }
    test_wilcox_result = wilcox_test_cv_wvic > .05

    if(sum(test_wilcox_result) > 0){
      index_select_wilcox = which(test_wilcox_result[1:n_models] == TRUE)

      if(length(index_select_wilcox) != 1){
        model_complexity = rep(NA,n_models)
        for(k in 1:n_models){
          model_complexity[k] = test_wilcox_result[k]*model_test[[k]]$plength
        }
        model_complexity[model_complexity == 0] = NA
        index_select_wilcox = which.min(model_complexity)
      }

    }else{
      index_select_wilcox = mod_selected_cv
    }
  }

  if(test_pval == FALSE){
    estimate = as.matrix(model_test[[mod_selected_cv]]$theta)
    rownames(estimate) = model_test[[mod_selected_cv]]$process.desc
    colnames(estimate) = "Estimates"

    #obj.value = model_test[[mod_selected_cv]]$value
    #names(obj.value) = "Value Objective Function"

    scales.num = rep(NA,length(mimu))

    for (i in 1:n_replicates){
      scales.num[i] = length(mimu[[i]]$scales)
    }

    max.scales.index = which.max(scales.num)
    max.scales = max(scales.num)
    scales.max = mimu[[max.scales.index]]$scales

    tau.max = 2^(1:max.scales)

    # WV implied by the parameter

    wv.implied = wv_theo(model_test[[mod_selected_cv]], tau.max)

    n_process = length(model_test[[mod_selected_cv]]$desc)

    #### Extact individual model for Theoretical decomposition

    model.desc.decomp.theo = list()

    # Initialise counter
    counter = 1

    for (i in 1:n_process){

      if (model_test[[mod_selected_cv]]$desc[i] == "RW"){
        model.desc.decomp.theo[[i]] = RW(gamma2 = model_test[[mod_selected_cv]]$theta[counter])
        counter = counter + 1
      }

      # is white noise?
      if (model_test[[mod_selected_cv]]$desc[i] == "WN"){
        model.desc.decomp.theo[[i]] =WN(sigma2 = model_test[[mod_selected_cv]]$theta[counter])
        counter = counter + 1
      }

      # is drift?
      if (model_test[[mod_selected_cv]]$desc[i] == "DR"){
        model.desc.decomp.theo[[i]] = DR(omega = model_test[[mod_selected_cv]]$theta[counter])
        counter = counter + 1
      }

      # is quantization noise?
      if (model_test[[mod_selected_cv]]$desc[i] == "QN"){
        model.desc.decomp.theo[[i]] = QN(q2 = model_test[[mod_selected_cv]]$theta[counter])
        counter = counter + 1
      }

      # is AR1?
      if (model_test[[mod_selected_cv]]$desc[i] == "AR1"){
        model.desc.decomp.theo[[i]] = AR1(phi = model_test[[mod_selected_cv]]$theta[counter], sigma2 = model_test[[mod_selected_cv]]$theta[counter + 1])
        counter = counter + 2
      }
    }

    # Compute individual theoretical wv
    decomp.theo = list()
    for (i in 1:n_process){
      model.decomp.theo = model.desc.decomp.theo[[i]]
      decomp.theo[[i]] =  wv_theo(model.decomp.theo, tau.max)
    }

    model.hat = model_test[[mod_selected_cv]]

  }else{

    estimate = as.matrix(model_test[[index_select_wilcox]]$theta)
    rownames(estimate) = model_test[[index_select_wilcox]]$process.desc
    colnames(estimate) = "Estimates"

    #obj.value = model_test[[mod_selected_cv]]$value
    #names(obj.value) = "Value Objective Function"

    scales.num = rep(NA,length(mimu))

    for (i in 1:n_replicates){
      scales.num[i] = length(mimu[[i]]$scales)
    }

    max.scales.index = which.max(scales.num)
    max.scales = max(scales.num)
    scales.max = mimu[[max.scales.index]]$scales

    tau.max = 2^(1:max.scales)

    # WV implied by the parameter

    wv.implied = wv_theo(model_test[[index_select_wilcox]], tau.max)

    n_process = length(model_test[[index_select_wilcox]]$desc)

    #### Extact individual model for Theoretical decomposition

    model.desc.decomp.theo = list()

    # Initialise counter
    counter = 1

    for (i in 1:n_process){

      if (model_test[[index_select_wilcox]]$desc[i] == "RW"){
        model.desc.decomp.theo[[i]] = RW(gamma2 = model_test[[index_select_wilcox]]$theta[counter])
        counter = counter + 1
      }

      # is white noise?
      if (model_test[[index_select_wilcox]]$desc[i] == "WN"){
        model.desc.decomp.theo[[i]] =WN(sigma2 = model_test[[index_select_wilcox]]$theta[counter])
        counter = counter + 1
      }

      # is drift?
      if (model_test[[index_select_wilcox]]$desc[i] == "DR"){
        model.desc.decomp.theo[[i]] = DR(omega = model_test[[index_select_wilcox]]$theta[counter])
        counter = counter + 1
      }

      # is quantization noise?
      if (model_test[[index_select_wilcox]]$desc[i] == "QN"){
        model.desc.decomp.theo[[i]] = QN(q2 = model_test[[index_select_wilcox]]$theta[counter])
        counter = counter + 1
      }

      # is AR1?
      if (model_test[[index_select_wilcox]]$desc[i] == "AR1"){
        model.desc.decomp.theo[[i]] = AR1(phi = model_test[[index_select_wilcox]]$theta[counter], sigma2 = model_test[[index_select_wilcox]]$theta[counter + 1])
        counter = counter + 2
      }
    }

    # Compute individual theoretical wv
    decomp.theo = list()
    for (i in 1:n_process){
      model.decomp.theo = model.desc.decomp.theo[[i]]
      decomp.theo[[i]] =  wv_theo(model.decomp.theo, tau.max)
    }

    model.hat = model_test[[index_select_wilcox]]

  }


  model_selected = structure(list(estimate = estimate,
                       decomp.theo = decomp.theo,
                       model = model,
                       cv_wvic = cv_wvic,
                       obj_out_sample = obj_out_sample,
                       mod_selected_cv = mod_selected_cv,
                       model_test = model_test,
                       model.hat = model.hat,
                       scales.max = scales.max,
                       wv.implied = wv.implied,
                       mimu = mimu), class = "mgmwm")
}
