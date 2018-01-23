
#' @export
mgmwm_obj_function = function(theta, model, mimu){

  # Step 1: compute theoretical WV
  tau = list()
  wv.theo = list()
  obj_len  = length(mimu)
  for (i in 1:obj_len) {
    tau[[i]] = 2^(1:length(mimu[[i]]$scales))
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
mgmwm = function(model, mimu){
  # Check if model is a ts object
  if(!is.mimu(mimu)){
    stop("`mimu` must be created from a `mimu` object. ")
  }

  # Check if mium is a valid object
  if(!is.mimu(mimu)){
    stop("`model` must be created from a `ts.model` object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }

  # Add starting value algo ->  call gmwm_gmwm_master_cpp function on each experiment in
  # the mimu object at hand, then we could average or perform some kind of random search
  # between the estimated values.
  desc = model$desc

  nr = length(mimu)

  obj = model$obj.desc

  np = model$plength

  starting = model$starting

  theta = model$theta


  para.gmwm = matrix(NA,np,nr)
  N = rep(NA,nr)
  for (i in 1:nr){
    N[i] = length(mimu[[i]]$data)
    data = mimu[[i]]$data

    out = .Call('gmwm_gmwm_master_cpp', PACKAGE = 'gmwm', data, theta, desc, obj,
                model.type = 'imu' , starting = model$starting,
                p = 0.05, compute_v = "bootstrap", K = 1, H = 100, G = 10000,
                robust=FALSE, eff = 1)
    para.gmwm[,i] = out[[1]]
  }
  starting.value = apply(para.gmwm, 1, mean)

  out2 = optim(starting.value, mgmwm_obj_function, model = model, mimu = mimu)

  estimate = out2$par
  rownames(estimate) = model$process.desc
  colnames(estimate) = "Estimates"

  class(obj) = "mgmwm"
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################



#' @export
plot.mgmwm = function(obj_list, split = FALSE, add_legend = TRUE, xlab = NULL,
                     ylab = NULL, col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL,
                     nb_ticks_y = NULL, legend_position = "bottomleft", ci_wv = NULL, point_cex = NULL,
                     point_pch = NULL, names = NULL){

  obj_name = attr(obj, "exp.name")
  obj_len  = length(obj_list)
  units = attr(obj, "unit")
  main = attr(obj, "sensor.name")

  # Check if passed objects are of the class wvar
  #is_wvar = sapply(obj_list, FUN = is, class2 = 'wvar')
  #if(!all(is_wvar == T)){
  #  stop("Supplied objects must be 'wvar' objects.")
  #}

  # Check length of time series argument
  if (obj_len == 0){
    stop('No object given!')
  }else if (obj_len == 1){
    # -> plot.wvar
    plot.wvar(..., nb_ticks_x = nb_ticks_x, nb_ticks_y = nb_ticks_y)
  }else{

    if (is.null(xlab)){
      if (is.null(units)){
        xlab = expression(paste("Scale ", tau, sep =""))
      }else{
        xlab = bquote(paste("Scale ", "(", .(units), ")", sep = " "))
      }
    }else{
      xlab = xlab
    }

    if (is.null(ylab)){
      ylab = bquote(paste("Wavelet Variance ", nu^2, sep = " "))
    }else{
      ylab = ylab
    }

    if (is.null(ci_wv)){
      ci_wv = rep(TRUE, obj_len)
    }else{
      ci_wv = rep(ci_wv, obj_len)
    }


    hues = seq(15, 375, length = obj_len + 1)
    # Line and CI colors
    if (is.null(col_wv)){
      col_wv = hcl(h = hues, l = 65, c = 200, alpha = 1)[seq_len(obj_len)]
    }else{
      if (length(col_wv) != obj_len){
        col_wv = hcl(h = hues, l = 65, c = 200, alpha = 1)[seq_len(obj_len)]
      }
    }

    if (is.null(col_ci)){
      col_ci = hcl(h = hues, l = 80, c = 100, alpha = 0.2)[seq_len(obj_len)]
    }else{
      if (length(col_ci) != obj_len){
        col_ci = hcl(h = hues, l = 80, c = 100, alpha = 0.2)[seq_len(obj_len)]
      }
    }

    # Range
    # Find x and y limits
    x_range = y_range = rep(NULL, 2)
    for (i in 1:obj_len){
      x_range = range(c(x_range, obj_list[[i]]$scales))
      y_range = range(c(y_range, obj_list[[i]]$ci_low, obj_list[[i]]$ci_high))
    }

    x_low = floor(log10(x_range[1]))
    x_high = ceiling(log10(x_range[2]))
    y_low = floor(log10(y_range[1]))
    y_high = ceiling(log10(y_range[2]))

    # Axes
    if (is.null(nb_ticks_x)){
      nb_ticks_x = 6
    }

    if (is.null(nb_ticks_y)){
      nb_ticks_y = 5
    }

    x_ticks = seq(x_low, x_high, by = 1)
    if (length(x_ticks) > nb_ticks_x){
      x_ticks = x_low + ceiling((x_high - x_low)/(nb_ticks_x + 1))*(0:nb_ticks_x)
    }
    x_labels = sapply(x_ticks, function(i) as.expression(bquote(10^ .(i))))
    x_at = 10^x_ticks
    x_actual_length = sum((x_at < x_range[2])*(x_at > x_range[1]))

    if (x_actual_length < (3 + as.numeric(split == FALSE))){
      x_low = floor(log2(x_range[1]))
      x_high = ceiling(log2(x_range[2]))
      x_ticks = seq(x_low, x_high, by = 1)
      if (length(x_ticks) > 8){
        x_ticks = seq(x_low, x_high, by = 2)
      }
      x_labels = sapply(x_ticks, function(i) as.expression(bquote(2^ .(i))))
      x_at = 2^x_ticks
    }
    y_ticks <- seq(y_low, y_high, by = 1)

    if (length(y_ticks) > nb_ticks_y){
      y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
    }
    y_labels = sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))
    y_at = 10^y_ticks

    # Legend position
    if (is.null(legend_position)){
      inter = rep(NA, obj_len)
      for (i in 1:obj_len){
        inter[i] = obj_list[[i]]$variance[1]
      }
      mean_wv_1 = mean(inter)
      if (which.min(abs(c(y_low, y_high) - log2(mean_wv_1))) == 1){
        legend_position = "topleft"
      }else{
        legend_position = "bottomleft"
      }
    }

    if (is.null(point_pch)){
      inter = rep(15:18, obj_len)
      point_pch = inter[1:obj_len]
    }else{
      if (length(point_pch) != obj_len){
        inter = rep(15:18, obj_len)
        point_pch = inter[1:obj_len]
      }
    }

    if (is.null(point_cex)){
      inter = rep(c(1.25,1.25,1.25,1.6), obj_len)
      point_cex = inter[1:obj_len]
    }else{
      if (length(point_pch) != obj_len){
        inter = rep(c(1.25,1.25,1.25,1.6), obj_len)
        point_cex = inter[1:obj_len]
      }
    }

    if (is.null(names)){
      names = obj_name
    }else{
      if (length(names) != obj_len){
        names = obj_name
      }
    }

    # Arguments passed into compare_wvar_split or compare_wvar_no_split
    graph_details = list(obj_list = obj_list, obj_len = obj_len, names = names, xlab = xlab,
                         ylab = ylab, col_wv = col_wv, add_legend = add_legend,
                         col_ci = col_ci, main = main, legend_position = legend_position,
                         ci_wv = ci_wv, point_cex = point_cex, point_pch = point_pch,
                         x_range = x_range, y_range = y_range, x_ticks = x_ticks,
                         x_labels = x_labels, y_labels = y_labels, x_at = x_at, y_at = y_at,
                         y_ticks = y_ticks, nb_ticks_x = nb_ticks_x, nb_ticks_y = nb_ticks_y)

    if (split == FALSE){
      # -> compare_wvar_no_split
      compare_wvar_no_split(graph_details)
    }else{
      # -> compare_wvar_split
      compare_wvar_split(graph_details)
    }
  }
}