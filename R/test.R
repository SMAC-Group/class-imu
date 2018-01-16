#' @title This funciton adds up two numbers
#'
#' @description What the functionkjh is doing
#' @param x A \code{numeric} value.
#' @param y A \code{numeric} value.
#' @return A \code{numeric} which adds \code{x} and \code{y}.
#' @author Stephane Guerrier
#' @export
#' @examples
#' add(3,4)
#' add(0,0)
add = function(x, y){
  x + y
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
