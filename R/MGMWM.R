

#' @export
mgmwm_obj_function = function(theta, model, mimu){

  M = length(model$process.desc)

  # Initialise counter
  counter = 1

  for (j in 1:M){
    # is random walk?
    if (model$process.desc[j] == "RW"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }

    # is white noise?
    if (model$process.desc[j] == "WN"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }

    # is drift?
    if (model$process.desc[j] == "DR"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (model$process.desc[j] == "QN"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }

    # is AR1?
    if (model$process.desc[j] == "AR1"){
      theta[counter] = transform_phi(theta[counter])
      counter = counter + 1
    }

    # is SIGMA2?
    if (model$process.desc[j] == "SIGMA2"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }
  }

  model$theta = theta

  # Step 1: compute theoretical WV
  tau = list()
  wv.theo = list()
  obj_len  = length(mimu)

  for (i in 1:obj_len) {
    tau[[i]] = 2^(1:length(mimu[[i]]$scales))
    # Compute theoretical wvar for each latent process
    wv.theo[[i]] = wv_theo(model, tau[[i]])
  }


  # Step 2: compute Omega
  Omega = list()
  for (i in 1:obj_len){
    Omega[[i]] = diag(1/(mimu[[i]]$ci_high - mimu[[i]]$ci_low)^2)
  }

  # Step 3: compute actual objective function
  out = 0
  for (i in 1:obj_len){
    dif_vect = wv.theo[[i]] - mimu[[i]]$variance
    out = out + t(dif_vect)%*%Omega[[i]]%*%dif_vect
  }
  out
}

#' @export
mgmwm = function(model, mimu, CI = FALSE, stationarity_test = FALSE, B_stationarity_test= 500,
                 alpha_ci = NULL, alpha_near_test = NULL, seed = 2710, n_boot_ci_max = 300){

  # Check if model is a ts object
  if(!is.mimu(mimu)){
    stop("`mimu` must be created from a `mimu` object. ")
  }

  # Check if mium is a valid object
  if(!is.mimu(mimu)){
    stop("`model` must be created from a `ts.model` object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }

  if(is.null(alpha_ci)){
    alpha_ci = .05
  }else{
    alpha_ci = alpha_ci
  }

  if(is.null(alpha_near_test)){
    alpha_near_test = .05
  }else{
    alpha_near_test = alpha_near_test
  }

  if(is.null(B_stationarity_test)){
    B = 100
  }else{
    B = B_stationarity_test
  }

  desc = model$desc

  n_process = length(model$desc)

  np = model$plength

  obj_desc = model$obj.desc

  theta = model$theta

  n_replicates = length(mimu)

  set.seed(seed)

  N = rep(NA,n_replicates)
  param_starting = matrix(NA,n_replicates,np)

  for (i in 1:n_replicates){
    N[i] = length(mimu[[i]]$data)

    if(N[i] > 10000){
      G = 1e6
    }else{
      G = 20000
    }
    data = mimu[[i]]$data

    uni_gmwm = .Call('gmwm_gmwm_master_cpp', PACKAGE = 'gmwm', data, theta, desc, obj = obj_desc,
                     model.type = 'imu' , starting = model$starting,
                     p = 0.05, compute_v = "fast", K = 1, H = 100, G = G,
                     robust=FALSE, eff = 1)[[1]]

    param_starting[i,] = uni_gmwm
  }

  obj_value_starting_value = rep(NA, n_replicates)
  mgmwm_list = list()
  for (i in 1:n_replicates){
    #starting_value = apply(param_starting, 2, median)

    starting_value = inv_param_transform(model, param_starting[i,])
    mgmwm_list[[i]] = optim(starting_value, mgmwm_obj_function, model = model, mimu = mimu)
    obj_value_starting_value[i] = mgmwm_list[[i]]$value
  }

  out = mgmwm_list[[which.min(obj_value_starting_value)]]


  # Create estimated model object
  model_hat = model

  # Pass on the estimated paramters onto the model.
  model_hat$starting = FALSE
  model_hat$theta = out$par

  model_hat$theta = param_transform(model_hat,model_hat$theta)


  # Create the near-stationnary test

  if(stationarity_test == TRUE){
    distrib_H0 = rep(NA,B)

    for (i in 1:B){
      sim.H0 = list()
      for (j in 1:n_replicates){
        sim.H0[[j]] = simts::gen_gts(N[j],model_hat)
      }
      simu.obj = make_wvar_mimu_obj(for_test = sim.H0, freq = 100, unit = "s", sensor.name = "MTiG - Gyro. X",
                                    exp.name = c("today", "yesterday", "a few days ago"))
      distrib_H0[i] = optim(out$par, mgmwm_obj_function, model = model, mimu = simu.obj)$value
    }

    # extract p_value from the test
    p_value = sum(distrib_H0 >= out$value)/B

    # decision rules from the test
    if(p_value >= alpha_near_test){
      test_res = " Data are stationary"
    }else{
      test_res = " Data are nearly-stationary"
    }
  }else{
    test_res = NA
    p_value = NA
  }


  # Extact the max number of scales.

  scales.num = rep(NA,n_replicates)

  for (i in 1:n_replicates){
    scales.num[i] = length(mimu[[i]]$scales)
  }

  max.scales.index = which.max(scales.num)
  max.scales = max(scales.num)
  scales.max = mimu[[max.scales.index]]$scales
  tau.max = 2^(1:max.scales)

  # WV implied by the parameter

  wv_implied = wv_theo(model_hat, tau.max)


  #### Extact individual model for Theoretical decomposition

  model.desc.decomp.theo = list()

  # Initialise counter
  counter = 1

  for (i in 1:n_process){

    if (desc[i] == "RW"){
      model.desc.decomp.theo[[i]] = RW(gamma2 = model_hat$theta[counter])
      counter = counter + 1
    }

    # is white noise?
    if (desc[i] == "WN"){
      model.desc.decomp.theo[[i]] =WN(sigma2 = model_hat$theta[counter])
      counter = counter + 1
    }

    # is drift?
    if (desc[i] == "DR"){
      model.desc.decomp.theo[[i]] = DR(omega = model_hat$theta[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (desc[i] == "QN"){
      model.desc.decomp.theo[[i]] = QN(q2 = model_hat$theta[counter])
      counter = counter + 1
    }

    # is AR1?
    if (desc[i] == "AR1"){
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

  estimate = as.matrix(model_hat$theta)
  rownames(estimate) = model_hat$process.desc
  colnames(estimate) = "Estimates"

  obj.value = out$value
  names(obj.value) = "Value Objective Function"

  # Cancel previous seed
  set.seed(as.numeric(format(Sys.time(),"%s"))/10)

  out = structure(list(estimate = estimate,
                       obj.value = obj.value,
                       decomp.theo = decomp.theo,
                       model = model,
                       model_hat = model_hat,
                       scales.max = scales.max,
                       p_value = p_value,
                       wv_implied = wv_implied,
                       test.result = test_res,
                       mimu = mimu), class = "mgmwm")
  invisible(out)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
######################################              GRAPH               ############################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


#' @export
plot.mgmwm = function(obj_list, process.decomp = FALSE,
                      add_legend_mgwmw = TRUE, legend_pos = NULL, ylab_mgmwm = NULL){

  mimu_obj_name = attr(obj_list[[7]], "exp.name")
  mimu_obj_name = paste("Empirical WV", mimu_obj_name)

  if (is.null(ylab_mgmwm)){
    ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
  }else{
    ylab = ylab_mgmwm
  }


  plot(obj_list$mimu, add_legend = FALSE,ylab = ylab)
  U = length(obj_list$decomp.theo)
  col_wv = hcl(h = seq(100, 375, length = U + 1), l = 65, c = 200, alpha = 1)[1:U]

  if(process.decomp == TRUE){
    # Number of Latent proces

    # Plot lines of decomp theo
    for (i in 1:U){
      lines(t(obj_list$scales.max), obj_list$decomp.theo[[i]], col = col_wv[i])
    }
  }
  # Plot implied WV
  lines(t(obj_list$scales.max),obj_list$wv_implied, type = "l", lwd = 3, col = "#F47F24", pch = 1, cex = 1.5)
  lines(t(obj_list$scales.max),obj_list$wv_implied, type = "p", lwd = 2, col = "#F47F24", pch = 1, cex = 1.5)

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
}

#' @export
ci_mgmwm = function(obj,alpha_ci = 0.025, n_boot_ci_max = NULL){

  if(is.null(n_boot_ci_max)){
    n_boot_ci_max = 300
  }else{
    n_boot_ci_max = n_boot_ci_max
  }

  #
  n_replicates = length(obj$mimu)

  np = obj$model_hat$plength

  # Set up starting value in model to false
  obj$model_hat$starting = TRUE

  I = iterpc(n_replicates, n_replicates, replace = TRUE)
  perm = getall(I)
  n_permutation = dim(perm)[1]

  distrib_param = matrix(NA,n_permutation,np)
  starting_value = inv_param_transform(obj$model_hat, obj$model_hat$theta)

  if(n_permutation < n_boot_ci_max){

    for (i in 1:n_permutation){
      sampled_imu_obj = list()
      for (j in 1:n_replicates){
        sampled_imu_obj[[j]] = obj$mimu[[perm[i,j]]]
        class(sampled_imu_obj) = "mimu"
      }
      distrib_param[i,] = optim(starting_value, mgmwm_obj_function, model = obj$model_hat, mimu = sampled_imu_obj)$par
      distrib_param[i,] = param_transform(obj$model_hat, distrib_param[i,])
    }
  }else{
    n_permutation = n_boot_ci_max
    for (i in 1:n_permutation){
      sampled_imu_obj = list()
      sampled_permutation = sample(1:n_replicates, n_replicates, replace = TRUE)
      for (j in 1:n_replicates){
        sampled_imu_obj[[j]] = mimu[[sampled_permutation[i]]]
        class(sampled_imu_obj) = "mimu"
      }
      distrib_param[i,] = optim(starting_value, mgmwm_obj_function, model = obj$model_hat, mimu = sampled_imu_obj)$par

      distrib_param[i,] = param_transform(obj$model_hat, distrib_param[i,])
    }
  }

  ci_low = rep(NA,np)
  ci_high = rep(NA,np)
  #Compute the empirical quantile
  for (k in 1:np){
    ci_low[k] = as.numeric(quantile(na.omit(distrib_param[,k]),(alpha_ci/2)))
    ci_high[k] = as.numeric(quantile(na.omit(distrib_param[,k]),(1-alpha_ci/2)))
  }

  out = structure(list(model_hat = obj$model_hat,
                       distrib_param = distrib_param,
                       ci_low =ci_low,
                       ci_high = ci_high), class = "ci_mgmwm")
  invisible(out)
}
