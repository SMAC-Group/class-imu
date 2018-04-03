

#' @export
model_selection = function(model,mimu, s_est = NULL,
                           n_boot_ci_max = NULL, stationarity_test = FALSE, B_stationarity_test = NULL,
                           alpha_near_test = NULL, paired_test = FALSE,
                           alpha_paired_test = NULL, seed = 2710){

  # model_max must be an object of type model

  # mimu must be an mimu obj

  if(is.null(B_stationarity_test)){
    B = 100
  }else{
    B = B_stationarity_test
  }

  # Level of confidence for stationarity test
  if(is.null(alpha_near_test)){
    alpha_near_test = .05
  }else{
    alpha_near_test = alpha_near_test
  }

  # Level of confidence for Wilcoxon paired test
  if(is.null(alpha_paired_test)){
    alpha_paired_test = .05
  }else{
    alpha_paired_test = alpha_paired_test
  }

  # Number of replicates
  n_replicates = length(mimu)

  # Number of time series replicates for estimation
  if(is.null(s_est)){
    s_est = ceiling(n_replicates/2)
  }else{
    s_est = s_est
  }

  # Number of time series replicates for validation
  s_valid = n_replicates - s_est

  # Matrix of possible combination of Time series to estimate
  pair = t(combn(n_replicates,s_est))

  # Number of possible combination
  n_replicates_permutation = (dim(pair)[1])

  # Create model object with all possible nested model in model_max
  model_est = model_combination(model_max = model)

  model_test = list()

  n_models = length(model_est)

  cv_wvic = rep(NA,n_models)

  obj_out_sample = matrix(NA,n_replicates_permutation, n_models)

  ## Compute the maximum scales for signal of different length
  scales.num = rep(NA,length(mimu))

  for (i in 1:n_replicates){
    scales.num[i] = length(mimu[[i]]$scales)
  }

  max.scales.index = which.max(scales.num)
  max.scales = max(scales.num)
  scales.max = mimu[[max.scales.index]]$scales

  tau.max = 2^(1:max.scales)

  # Set seed for reproducibility
  set.seed(seed)

  pb <- progress_bar$new(
    format = "  Model :current of :total Models. Time remaining:  :eta",
    clear = FALSE, total = n_models, width = 100)

  for (i in 1:n_models){
    desc = model_est[[i]]$desc

    np = model_est[[i]]$plength

    obj_desc = model_est[[i]]$obj.desc

    theta = model_est[[i]]$theta

    starting = model_est[[i]]$starting


    N = rep(NA,n_replicates)
    param_starting = matrix(NA,n_replicates,np)

    for (k in 1:n_replicates){
      N[k] = length(mimu[[k]]$data)

      if(N[k] > 10000){
        G = 1e6
      }else{
        G = 20000
      }
      data = mimu[[k]]$data

      uni_gmwm = .Call('gmwm_gmwm_master_cpp', PACKAGE = 'gmwm', data, theta, desc, obj = obj_desc,
                       model.type = 'imu' , starting = starting,
                       p = 0.05, compute_v = "fast", K = 1, H = 100, G = G,
                       robust=FALSE, eff = 1)[[1]]

      param_starting[k,] = uni_gmwm
    }

    obj_value_starting_value = rep(NA, n_replicates)
    mgmwm_list = list()
    for (j in 1:n_replicates){
      starting_value = inv_param_transform(model_est[[i]], param_starting[j,])
      mgmwm_list[[j]] = optim(starting_value, mgmwm_obj_function, model = model_est[[i]], mimu = mimu)
      obj_value_starting_value[j] = mgmwm_list[[j]]$value
    }

    mgmwm_full = mgmwm_list[[which.min(obj_value_starting_value)]]$par

    for (d in 1:n_replicates_permutation){

      mimu_est = list()
      class(mimu_est) = "mimu"
      mimu_test = mimu[-pair[d,sequence(s_est)]]

      for (s in 1:s_est){
        mimu_est[[s]] = mimu[[pair[d,s]]]
      }

      out = optim(mgmwm_full, mgmwm_obj_function, model = model_est[[i]], mimu = mimu_est)

      # transform e.g. variance = exp(out$par[1])

      model_test[[i]] = model_est[[i]]

      # Pass on the estimated paramters onto the model.
      model_test[[i]]$starting = FALSE
      model_test[[i]]$theta = out$par

      obj_out_sample[d,i] = mgmwm_obj_function(model_test[[i]]$theta, model_test[[i]], mimu_test)/s_valid

      model_test[[i]]$theta = param_transform(model_test[[i]], model_test[[i]]$theta)
    }
    cv_wvic[i] = mean(obj_out_sample[,i])
    # Update progress bar
    pb$tick(tokens = list(what = "foo   "))
    Sys.sleep(1 / n_models)
  }

  # Compute the average on all permutation of the out of sample WVIC
  mod_selected_cv = which.min(cv_wvic)

  #Paired Wilcoxon test
  if(paired_test == TRUE){
    wilcox_test_cv_wvic = rep(FALSE,n_models)

    model_selected_cv_size = model_test[[mod_selected_cv]]$plength
    # Compute the Wilcoxon test
    for (i in 1:n_models){
      if(model_selected_cv_size >= model_test[[i]]$plength){
        wilcox_test_cv_wvic[i] = wilcox.test(obj_out_sample[,i],obj_out_sample[,mod_selected_cv],paired = T,alternative = "greater")$p.val
      }
    }
    test_wilcox_result = (wilcox_test_cv_wvic > alpha_paired_test)

    if(sum(test_wilcox_result) > 0){
      index_select_wilcox_list = which(test_wilcox_result[1:n_models] == TRUE)

      if(length(index_select_wilcox_list) != 1){
        model_complexity = rep(NA,n_models)
        for(k in 1:n_models){
          model_complexity[k] = test_wilcox_result[k]*model_test[[k]]$plength
        }
        model_complexity[model_complexity == 0] = NA
        index_select_wilcox = which.min(model_complexity)
      }else{
        paired_test = FALSE
      }
    }else{
      paired_test = FALSE
    }
  }


  # Output if no Wilcoxon test
  if(paired_test == TRUE){
    # Output if Wilcoxon test

    # Extract model selected through the cv-wvic
    model_hat_empty = model_test[[index_select_wilcox]]
    model_hat_empty$starting = TRUE

    theta_mgmwm = mgmwm(model_hat_empty,
                        mimu, stationarity_test = FALSE,
                        B_stationarity_test = B_stationarity_test,
                        alpha_near_test = alpha_near_test,
                        seed = seed)

    model_hat = theta_mgmwm$model
    model_hat$starting = FALSE
    model_hat$theta = theta_mgmwm$estimate


    ## Create the ouput for the selected model
    estimate = as.matrix(model_hat$theta)
    rownames(estimate) = model_hat$process.desc
    colnames(estimate) = "Estimates"

    #obj.value = model_test[[mod_selected_cv]]$value
    #names(obj.value) = "Value Objective Function"

    # List of models not rejected by the Wilcoxon test
    model_list_wilcox_test = list()
    if(length(index_select_wilcox_list) == 1){
      model_list_wilcox_test = model_hat
    }else{
      model_list_wilcox_test[[1]] = model_hat
      test_wilcox_result[index_select_wilcox] = FALSE
      index_select_wilcox_list = which(test_wilcox_result[1:n_models] == TRUE)
      for (j in 1:length(index_select_wilcox_list)){

        # First model of List is the selected one
        # Other models selected by test
        model_list_wilcox_test_empty = model_test[[index_select_wilcox_list[j]]]
        model_list_wilcox_test_empty$starting = TRUE

        theta_mgmwm_model_test = mgmwm(model_list_wilcox_test_empty,mimu, stationarity_test = FALSE,
                                       B_stationarity_test = B_stationarity_test,
                                       alpha_near_test = alpha_near_test,
                                       seed = seed)


        model_list_wilcox_test[[j+1]] = theta_mgmwm_model_test$model

        model_list_wilcox_test[[j+1]] $starting = FALSE
        model_list_wilcox_test[[j+1]] $theta = theta_mgmwm_model_test$estimate

      }
    }

    # Create output for full list of models not rejected by Wilcoxon test

    wv_implied = list()
    n_process = rep(NA,length(model_list_wilcox_test))
    wv_theo_latent_process = list()

    for (j in 1:length(model_list_wilcox_test)){
      # WV implied by the parameter
      wv_implied[[j]] = wv_theo(model_list_wilcox_test[[j]], tau.max)

      # Extact individual model for Theoretical decomposition

      n_process[j] = length(model_list_wilcox_test[[j]]$desc)


      model.desc.decomp.theo = list()
      # Initialise counter
      counter = 1

      for (i in 1:n_process[j]){

        if (model_list_wilcox_test[[j]]$desc[i] == "RW"){
          model.desc.decomp.theo[[i]] = RW(gamma2 = model_list_wilcox_test[[j]]$theta[counter])
          counter = counter + 1
        }

        # is white noise?
        if (model_list_wilcox_test[[j]]$desc[i] == "WN"){
          model.desc.decomp.theo[[i]] =WN(sigma2 = model_list_wilcox_test[[j]]$theta[counter])
          counter = counter + 1
        }

        # is drift?
        if (model_list_wilcox_test[[j]]$desc[i] == "DR"){
          model.desc.decomp.theo[[i]] = DR(omega = model_list_wilcox_test[[j]]$theta[counter])
          counter = counter + 1
        }

        # is quantization noise?
        if (model_list_wilcox_test[[j]]$desc[i] == "QN"){
          model.desc.decomp.theo[[i]] = QN(q2 = model_list_wilcox_test[[j]]$theta[counter])
          counter = counter + 1
        }

        # is AR1?
        if (model_list_wilcox_test[[j]]$desc[i] == "AR1"){
          model.desc.decomp.theo[[i]] = AR1(phi = model_list_wilcox_test[[j]]$theta[counter], sigma2 = model_list_wilcox_test[[j]]$theta[counter + 1])
          counter = counter + 2
        }
      }

      # Compute individual theoretical wv
      decomp.theo = list()
      for (i in 1:n_process[j]){
        decomp.theo[[i]] =  wv_theo(model.desc.decomp.theo[[i]], tau.max)
      }
      wv_theo_latent_process[[j]] = decomp.theo
    }

    model_selected = structure(list(estimate = estimate,
                                    wv_theo_latent_process = wv_theo_latent_process,
                                    model_hat = model_hat,
                                    model_test = model_test,
                                    mod_selected_cv = mod_selected_cv,
                                    model_list_wilcox_test = model_list_wilcox_test,
                                    scales.max = scales.max,
                                    wv_implied = wv_implied,
                                    mimu = mimu), class = "cvwvic")
  }else{
    # Extract model selected through the cv-wvic

    model_hat_empty = model_test[[mod_selected_cv]]
    model_hat_empty$starting = TRUE

    theta_mgmwm = mgmwm(model_hat_empty,mimu, stationarity_test = FALSE,
                        B_stationarity_test = B_stationarity_test,
                        alpha_near_test = alpha_near_test,
                        seed = seed)

    model_hat = theta_mgmwm$model
    model_hat$starting = FALSE
    model_hat$theta = theta_mgmwm$estimate


    ## Create the ouput for the selected model
    estimate = as.matrix(model_hat$theta)
    rownames(estimate) = model_hat$process.desc
    colnames(estimate) = "Estimates"

    #obj.value = model_test[[mod_selected_cv]]$value
    #names(obj.value) = "Value Objective Function"

    # WV implied by the parameter

    wv_implied = wv_theo(model_hat, tau.max)

    # Extact individual model for Theoretical decomposition

    model.desc.decomp.theo = list()

    n_process = length(model_hat$desc)

    # Initialise counter
    counter = 1

    for (i in 1:n_process){

      if (model_hat$desc[i] == "RW"){
        model.desc.decomp.theo[[i]] = RW(gamma2 = model_hat$theta[counter])
        counter = counter + 1
      }

      # is white noise?
      if (model_hat$desc[i] == "WN"){
        model.desc.decomp.theo[[i]] =WN(sigma2 = model_hat$theta[counter])
        counter = counter + 1
      }

      # is drift?
      if (model_hat$desc[i] == "DR"){
        model.desc.decomp.theo[[i]] = DR(omega = model_hat$theta[counter])
        counter = counter + 1
      }

      # is quantization noise?
      if (model_hat$desc[i] == "QN"){
        model.desc.decomp.theo[[i]] = QN(q2 = model_hat$theta[counter])
        counter = counter + 1
      }

      # is AR1?
      if (model_hat$desc[i] == "AR1"){
        model.desc.decomp.theo[[i]] = AR1(phi = model_hat$theta[counter], sigma2 = model_hat$theta[counter + 1])
        counter = counter + 2
      }
    }

    # Compute individual theoretical wv
    decomp.theo = list()
    for (i in 1:n_process){
      model.decomp.theo = model.desc.decomp.theo[[i]]
      decomp.theo[[i]] =  wv_theo(model.decomp.theo, tau.max)
    }

    # Cancel previous seed
    set.seed(as.numeric(format(Sys.time(),"%s"))/10)

    model_selected = structure(list(estimate = estimate,
                                    decomp.theo = decomp.theo,
                                    model_hat = model_hat,
                                    model_test = model_test,
                                    mod_selected_cv = mod_selected_cv,
                                    scales.max = scales.max,
                                    wv_implied = wv_implied,
                                    mimu = mimu), class = "mgmwm")
    }
  invisible(model_selected)
}


#' @export
plot.cvwvic = function(obj_list, process.decomp = FALSE,
                       add_legend_mgwmw = TRUE, legend_pos = NULL,
                       plot_type = "one", ylab_cvwvic = NULL){



  if(plot_type == "one"){

    if (is.null(ylab_cvwvic)){
      ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
    }else{
      ylab = ylab_cvwvic
    }

    plot(obj_list$mimu, add_legend = FALSE, ylab = ylab)
    U = length(obj_list$wv_theo_latent_process[[1]])
    col_wv = hcl(h = seq(100, 375, length = U + 1), l = 65, c = 200, alpha = 1)[1:U]

    if(process.decomp == TRUE){
      # Number of Latent proces

      # Plot lines of decomp theo
      for (i in 1:U){
        lines(t(obj_list$scales.max), obj_list$wv_theo_latent_process[[1]][[i]], col = col_wv[i])
      }
    }

    # Plot implied WV
    lines(t(obj_list$scales.max),obj_list$wv_implied[[1]], type = "l", lwd = 3, col = "#F47F24", pch = 1, cex = 1.5)
    lines(t(obj_list$scales.max),obj_list$wv_implied[[1]], type = "p", lwd = 2, col = "#F47F24", pch = 1, cex = 1.5)

    if(process.decomp == TRUE){
      legend_names = c("Implied WV", obj_list$model_hat$desc)
      col_legend = c("#F47F24",col_wv)
      p_cex_legend = c(1.5,rep(NA,U))
    }else{
      legend_names = c("Implied WV")
      col_legend = c("#F47F24")
      p_cex_legend = c(1.5)
    }


    if (is.null(legend_pos)){
      legend_pos = "bottomleft"
    }
    if (add_legend_mgwmw == TRUE){
      legend(legend_pos, legend_names, bty = "n", lwd = 1, pt.cex = 1.5, pch = p_cex_legend, col = col_legend)
    }
  }else if (plot_type == "all"){

    windows_col = ceiling(length(obj_list$wv_implied)/2)
    windows_row = ceiling(length(obj_list$wv_implied)/2)

    par(mfrow=c(windows_row,2), mar = c(4,4,4,2))

    for (i in 1:length(obj_list$wv_implied)){

      if (is.null(ylab_cvwvic)){
        ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
      }else{
        ylab = ylab_cvwvic
      }

      plot(obj_list$mimu, add_legend = FALSE, ylab = ylab)
      if(i ==1){
        title(main = "Model Selected")
      }
      U = length(obj_list$wv_theo_latent_process[[i]])
      col_wv = hcl(h = seq(100, 375, length = U + 1), l = 65, c = 200, alpha = 1)[1:U]

      if(process.decomp == TRUE){
        # Number of Latent proces

        # Plot lines of decomp theo
        for (j in 1:U){
          lines(t(obj_list$scales.max), obj_list$wv_theo_latent_process[[i]][[j]], col = col_wv[j])
        }
      }

      # Plot implied WV
      lines(t(obj_list$scales.max),obj_list$wv_implied[[i]], type = "l", lwd = 3, col = "#F47F24", pch = 1, cex = 1.5)
      lines(t(obj_list$scales.max),obj_list$wv_implied[[i]], type = "p", lwd = 2, col = "#F47F24", pch = 1, cex = 1.5)

      if(process.decomp == TRUE){
        legend_names = c("Implied WV", obj_list$model_list_wilcox_test[[i]]$desc)
        col_legend = c("#F47F24",col_wv)
        p_cex_legend = c(1.5,rep(NA,U))
      }else{
        legend_names = c("Implied WV")
        col_legend = c("#F47F24")
        p_cex_legend = c(1.5)
      }


      if (is.null(legend_pos)){
        legend_pos = "bottomleft"
      }
      if (add_legend_mgwmw == TRUE){
        legend(legend_pos, legend_names, bty = "n", lwd = 1, pt.cex = 1.5, pch = p_cex_legend, col = col_legend)
      }
    }
    par(mfrow=c(1,1), mar = c(4,4,4,4))
  }else if(plot_type == "merged"){
    process.decomp = FALSE

    if (is.null(ylab_cvwvic)){
      ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
    }else{
      ylab = ylab_cvwvic
    }

    plot(obj_list$mimu, add_legend = FALSE, transparency_wv = 0.4, transparency_ci = 0.05, ylab = ylab)

    U = length(obj_list$wv_implied)
    col_wv = hcl(h = seq(100, 375, length = U + 1), l = 65, c = 200, alpha = 1)[1:U]
    for (i in 1:length(obj_list$wv_implied)){


      # Plot implied WV
      lines(t(obj_list$scales.max),obj_list$wv_implied[[i]], type = "l", lwd = 3, col = col_wv[[i]], pch = 1, cex = 1.5)
      lines(t(obj_list$scales.max),obj_list$wv_implied[[i]], type = "p", lwd = 2, col = col_wv[[i]], pch = 1, cex = 1.5)

      legend_names = rep(NA,2*length(obj_list$wv_implied))
      col_legend = rep(NA,2*length(obj_list$wv_implied))
      if(i == 1){
        index = "selected"
      }else{
        index = i
      }
      # model_decomp_name = rep(NA,length(obj_list$model_list_wilcox_test[[i]]$desc))
      # for (j in 1:length(obj_list$model_list_wilcox_test[[i]]$desc)){
      #   model_decomp_name[j] = as.expression(bquote(paste0(obj_list$model_list_wilcox_test[[i]]$desc, collapse = '+'))))
      # }
      legend_names[(2*i)-1] = c(as.expression(bquote(paste(.("Implied WV "),  italic(M)[ .(index)]))))
      legend_names[2*i] = paste0(obj_list$model_list_wilcox_test[[i]]$desc, collapse = '+')
      col_legend[((2*i)-1):(2*i)] = c(col_wv[i], NA)
      p_cex_legend = rep(c(1.5,NA),length(obj_list$wv_implied))


      if (is.null(legend_pos)){
        legend_pos = "bottomleft"
      }
      if (add_legend_mgwmw == TRUE){
        legend(legend_pos, legend_names, bty = "n", lwd = 1, pt.cex = 1.5, pch = p_cex_legend, col = col_legend)
      }
    }
  }


}

