

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
mgmwm = function(model, mimu,CI = FALSE, stationarity_test = FALSE, B = 500,
                 fast = TRUE, alpha_ci = 0.05, alpha_near_test = 0.05){
  # Check if model is a ts object
  if(!is.mimu(mimu)){
    stop("`mimu` must be created from a `mimu` object. ")
  }

  # Check if mium is a valid object
  if(!is.mimu(mimu)){
    stop("`model` must be created from a `ts.model` object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }

  desc = model$desc

  n_process = length(model$desc)

  nr = length(mimu)

  np = model$plength

  obj_desc = model$obj.desc

  theta = model$theta

  set.seed(2710)

  N = rep(NA,nr)
  for (i in 1:nr){
    N[i] = length(mimu[[i]]$data)
  }

  if (fast == TRUE){
    index = which.max(N)
    N.fast = N[index]

    if(N.fast > 10000){
      G = 1e6
    }else{
      G = 20000
    }

    data = mimu[[index]]$data
    uni_gmwm = .Call('gmwm_gmwm_master_cpp', PACKAGE = 'gmwm', data, theta, desc, obj = obj_desc,
                     model.type = 'imu' , starting = model$starting,
                     p = 0.05, compute_v = "fast", K = 1, H = 100, G = G,
                     robust=FALSE, eff = 1)
    starting_value = uni_gmwm[[1]]
  }else{
    para_gmwm = matrix(NA,np,nr)
    for (i in 1:nr){
      data = mimu[[i]]$data

      # Select G

      if(N[i] > 10000){
        G = 1e6
      }else{
        G = 20000
      }

      uni_gmwm = .Call('gmwm_gmwm_master_cpp', PACKAGE = 'gmwm', data, theta, desc, obj = obj_desc,
                       model.type = 'imu' , starting = model$starting,
                       p = 0.05, compute_v = "fast", K = 1, H = 100, G = G,
                       robust=FALSE, eff = 1)
      para_gmwm[,i] = uni_gmwm[[1]]
    }

    starting_value = apply(para_gmwm, 1, mean)
  }

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
      starting_value[counter] = log(starting_value[counter])
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

  out = optim(starting_value, mgmwm_obj_function, model = model, mimu = mimu)


  # Compute the Confidence Intervals

  I = iterpc(nr, nr, replace = TRUE)
  perm = getall(I)
  n_permutation = dim(perm)[1]

  distrib_param = matrix(NA,n_permutation,length(out$par))

  if(CI == TRUE){
    if(n_permutation > 300){
      n_permutation = 300
    }else{
      n_permutation = n_permutation
    }

    for (i in 1:n_permutation){
      sampled_imu_obj = list()
      for (j in 1:nr){
        sampled_imu_obj[[j]] = mimu[[perm[i,j]]]
        class(sampled_imu_obj) = "mimu"
      }
      distrib_param[i,] = optim(out$par, mgmwm_obj_function, model = model, mimu = sampled_imu_obj)$par

      counter = 1
      for (j in 1:np){
        # is random walk?
        if (model$process.desc[j] == "RW"){
          distrib_param[,counter] = exp(distrib_param[,counter])
          counter = counter + 1
        }

        # is white noise?
        if (model$process.desc[j] == "WN"){
          distrib_param[,counter] = exp(distrib_param[,counter])
          counter = counter + 1
        }

        # is drift?
        if (model$process.desc[j] == "DR"){
          distrib_param[,counter] = exp(distrib_param[,counter])
          counter = counter + 1
        }

        # is quantization noise?
        if (model$process.desc[j] == "QN"){
          distrib_param[counter] = exp(distrib_param[,counter])
          counter = counter + 1
        }

        # is AR1?
        if (model$process.desc[j] == "AR1"){
          distrib_param[,counter] = transform_phi(distrib_param[,counter])
          counter = counter + 1
        }

        # is SIGMA2?
        if (model$process.desc[j] == "SIGMA2"){
          distrib_param[,counter] = exp(distrib_param[,counter])
          counter = counter + 1
        }
      }
    }

  }else{
    test_res = NA
    p_value = NA
  }


  # Create estimated model object
  model_hat = model

  # Pass on the estimated paramters onto the model.
  model_hat$starting = FALSE
  model_hat$theta = out$par



  # Initialise counter
  counter = 1

  for (j in 1:np){
    # is random walk?
    if (model$process.desc[j] == "RW"){
      model_hat$theta[counter] = exp(model_hat$theta[counter])
      counter = counter + 1
    }

    # is white noise?
    if (model$process.desc[j] == "WN"){
      model_hat$theta[counter] = exp(model_hat$theta[counter])
      counter = counter + 1
    }

    # is drift?
    if (model$process.desc[j] == "DR"){
      model_hat$theta[counter] = exp(model_hat$theta[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (model$process.desc[j] == "QN"){
      model_hat$theta[counter] = exp(model_hat$theta[counter])
      counter = counter + 1
    }

    # is AR1?
    if (model$process.desc[j] == "AR1"){
      model_hat$theta[counter] = transform_phi(model_hat$theta[counter])
      counter = counter + 1
    }

    # is SIGMA2?
    if (model$process.desc[j] == "SIGMA2"){
      model_hat$theta[counter] = exp(model_hat$theta[counter])
      counter = counter + 1
    }
  }


  # Create the near-stationnary test

  distrib_H0 = rep(NA,B)

  if(stationarity_test == TRUE){
    if(is.null(B)){
      B = 100
    }else{
      B = B
    }

    for (i in 1:B){
      sim.H0 = list()
      for (j in 1:nr){
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

  scales.num = rep(NA,length(mimu))

  for (i in 1:nr){
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
