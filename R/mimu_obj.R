#' @title Creates a Multivariate IMU object.
#'
#' @description This function creates a Multivariate IMU objects, combining multiple replicates of the
#' same inertial sensor.
#' @param X A \code{list} of vector \code{vector} value.
#' @param freq A \code{numeric} value.
#' @param units A \code{string}.
#' @param sensor.name A \code{string}.
#' @param freq A \code{string}
#' @author Stephane Guerrier
#' @export
#' @examples
#' n = 10^6
#' Xt = rnorm(n/4)
#' Yt = rnorm(n/2) + cumsum(rnorm(n/2, 0, 10^(-2)))
#' Zt = rnorm(n) + cumsum(rnorm(n, 0, 10^(-3)))
#' obj = make_wvar_mimu_obj(Xt, Yt, Zt, freq = 100, unit = "s",
#' sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))
make_wvar_mimu_obj = function(..., freq, unit = NULL, sensor.name = NULL, exp.name = NULL, for_test = NULL){

  if (is.null(for_test)){
    obj_list = list(...)
  }else{
    obj_list = for_test
  }

  obj_len  = length(obj_list)
  obj = list()
    for (i in 1:obj_len){
      inter = wvar(obj_list[[i]])
      inter$scales = inter$scales/freq
      obj[[i]] = c(data = list(obj_list[[i]]),inter)
    }

  if(is.null(sensor.name)){
    attr(obj, "") = sensor.name
  }else{
    attr(obj, "sensor.name") = sensor.name
  }

  if(is.null(freq)){
    attr(obj, "") = freq
  }else{
    attr(obj, "freq") = freq
  }

  if(is.null(unit)){
    attr(obj, "") = unit
  }else{
    attr(obj, "unit") = unit
  }

  if(is.null(exp.name)){
      attr(obj, "") = exp.name
  }else{
      attr(obj, "exp.name") = exp.name
  }

  class(obj) = "mimu"
  obj
}

#' @export
is.mimu = function(obj){class(obj) == "mimu"}


#' @export
plot.mimu = function(obj_list, split = FALSE, add_legend = TRUE, xlab = NULL,
                        ylab = NULL, col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL,
                        nb_ticks_y = NULL, legend_position = "bottomleft", ci_wv = NULL, point_cex = NULL,
                        point_pch = NULL, names = NULL, transparency_wv = NULL, transparency_ci = NULL){

  obj_name = attr(obj_list, "exp.name")
  obj_len  = length(obj_list)
  units = attr(obj_list, "unit")
  main = attr(obj_list, "sensor.name")

  # Check if passed objects are of the class mimu
  #is_mimu = class(obj_list)
  #if(!(is_mimu == T)){
  #  stop("Supplied object must be 'mimu' objects.")
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

    if (is.null(transparency_ci)){
      if(obj_len < 3){
        transparency_ci = 0.2
      }else{
        transparency_ci = 0.1
      }
    }else{
      transparency_ci = transparency_ci
    }

    if (is.null(transparency_wv)){
      if(obj_len < 3){
        transparency_wv = 1
      }else{
        transparency_wv = 0.8
      }
    }else{
      transparency_wv = transparency_wv
    }


    hues = seq(15, 375, length = obj_len + 1)
    # Line and CI colors
    if (is.null(col_wv)){
      col_wv = hcl(h = hues, l = 65, c = 200, alpha = transparency_wv)[seq_len(obj_len)]
    }else{
      if (length(col_wv) != obj_len){
        col_wv = hcl(h = hues, l = 65, c = 200, alpha = transparency_wv)[seq_len(obj_len)]
      }
    }

    if (is.null(col_ci)){
      col_ci = hcl(h = hues, l = 80, c = 100, alpha = transparency_ci)[seq_len(obj_len)]
    }else{
      if (length(col_ci) != obj_len){
        col_ci = hcl(h = hues, l = 80, c = 100, alpha = transparency_ci)[seq_len(obj_len)]
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

#' @export
wv_theo = function(model, tau){
  if (!class(model) == "ts.model"){
    # Error
  }

  # Number of scales
  J = length(tau)

  # Extract process description
  desc = model$desc

  # Number of latent processes
  M = length(desc)

  # Extract parameters
  theta = model$theta

  # Initialise counter
  counter = 1

  # Compute theoretical wvar for each latent process
  wv = matrix(NA, M, J)
  for (i in 1:M){
    # is random walk?
    if (desc[i] == "RW"){
      wv[i, ] = gmwm::rw_to_wv(gamma2 = theta[counter], tau = tau)
      counter = counter + 1
    }

    # is white noise?
    if (desc[i] == "WN"){
      wv[i, ] = gmwm::wn_to_wv(sigma2 = theta[counter], tau = tau)
      counter = counter + 1
    }

    # is drift?
    if (desc[i] == "DR"){
      wv[i, ] = gmwm::dr_to_wv(omega = theta[counter], tau = tau)
      counter = counter + 1
    }

    # is quantization noise?
    if (desc[i] == "QN"){
      wv[i, ] = gmwm::qn_to_wv(q2 = theta[counter], tau = tau)
      counter = counter + 1
    }

    # is AR1?
    if (desc[i] == "AR1"){
      wv[i, ] = gmwm::ar1_to_wv(phi = theta[counter], sigma2 = theta[counter + 1], tau = tau)
      counter = counter + 2
    }
  }
  # Summing them up and return
  apply(wv, 2, sum)
}
