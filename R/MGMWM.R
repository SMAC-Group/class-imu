
#' @export
mgmwm_obj_function = function(theta, model, mimu){


  # TRANSFORM e.g sigma = exp(variance)

  # Step 1: compute theoretical WV
  tau = list()
  wv.theo = list()
  obj_len  = length(mimu)
  for (i in 1:obj_len) {
    tau[[i]] = 2^(1:length(mimu[[i]]$scales))
    model$theta = theta
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
mgmwm = function(model, mimu, stationarity_test = FALSE, B = 500){
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

  para_gmwm = matrix(NA,np,nr)
  N = rep(NA,nr)
  for (i in 1:nr){
    N[i] = length(mimu[[i]]$data)
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

  # DO INVERSE TRANS. e.g sigma = log(variance)
  out = optim(starting_value, mgmwm_obj_function, model = model, mimu = mimu)

  # Create estimated model object
  model.hat = model

  # Pass on the estimated paramters onto the model.
  model.hat$starting = FALSE
  model.hat$theta = out$par

  # transform e.g. variance = exp(out$par[1])


  # Create the near-stationnary test

  distrib.H0 = rep(NA,B)

  if(stationarity_test == TRUE){
    if(is.null(B)){
      B = 100
    }else{
      B = B
    }
    for (i in 1:B){
      sim.H0 = list()
      for (j in 1:nr){
        sim.H0[[j]] = simts::gen_gts(N[j],model.hat)
      }
      simu.obj = make_wvar_mimu_obj(for_test = sim.H0, freq = 100, unit = "s", sensor.name = "MTiG - Gyro. X",
                                    exp.name = c("today", "yesterday", "a few days ago"))
      distrib.H0[i] = optim(starting_value, mgmwm_obj_function, model = model, mimu = simu.obj)$value
    }
  }

  p_value = sum(distrib.H0 >= out$value)/B
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

  wv.implied = wv_theo(model.hat, tau.max)


  #### Extact individual model for Theoretical decomposition

  model.desc.decomp.theo = list()

  # Initialise counter
  counter = 1

  for (i in 1:n_process){

    if (desc[i] == "RW"){
      model.desc.decomp.theo[[i]] = RW(gamma2 = model.hat$theta[counter])
      counter = counter + 1
    }

    # is white noise?
    if (desc[i] == "WN"){
      model.desc.decomp.theo[[i]] =WN(sigma2 = model.hat$theta[counter])
      counter = counter + 1
    }

    # is drift?
    if (desc[i] == "DR"){
      model.desc.decomp.theo[[i]] = DR(omega = model.hat$theta[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (desc[i] == "QN"){
      model.desc.decomp.theo[[i]] = QN(q2 = model.hat$theta[counter])
      counter = counter + 1
    }

    # is AR1?
    if (desc[i] == "AR1"){
      model.desc.decomp.theo[[i]] = AR1(phi = model.hat$theta[counter], sigma2 = model.hat$theta[counter + 1])
      counter = counter + 2
    }
  }

  # Compute individual theoretical wv
  decomp.theo = list()
  for (i in 1:n_process){
    model.decomp.theo = model.desc.decomp.theo[[i]]
    decomp.theo[[i]] =  wv_theo(model.decomp.theo, tau.max)
  }

  estimate = as.matrix(out$par)
  rownames(estimate) = model.hat$process.desc
  colnames(estimate) = "Estimates"

  obj.value = out$value
  names(obj.value) = "Value Objective Function"


  out = structure(list(estimate = estimate,
                       obj.value = obj.value,
                       decomp.theo = decomp.theo,
                       model = model,
                       model.hat = model.hat,
                       scales.max = scales.max,
                       mimu = mimu,
                       p_value = p_value,
                       wv.implied = wv.implied), class = "mgmwm")
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
plot.mgmwm = function(obj_list, process.decomp = FALSE){

  plot(mimu, add_legend = FALSE)

  if(process.decomp == TRUE){
    # Number of Latent proces
    U = length(obj_list$decomp.theo)

    hues = seq(100, 375, length = U + 1)
    col_wv = hcl(h = hues, l = 65, c = 200, alpha = 1)
    # Plot lines of decomp theo
    for (i in 1:U){
      lines(t(obj_list$scales.max), obj_list$decomp.theo[[i]], col = col_wv[i])
    }
  }
  # Plot implied WV
  lines(t(obj_list$scales.max),obj_list$wv.implied, type = "l", lwd = 3, col = "#F47F24", pch = 1, cex = 1.5)
  lines(t(obj_list$scales.max),obj_list$wv.implied, type = "p", lwd = 2, col = "#F47F24", pch = 1, cex = 1.5)
}
