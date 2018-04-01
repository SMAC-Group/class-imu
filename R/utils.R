

#' @export
# ----- phi
transform_phi = function(phi_R){
  phi = 2*(inv_logit(phi_R)-1/2)
  phi
}

#' @export
# ----- phi
inv_transform_phi = function(phi){
  phi_R = log(2/(1-phi)-1)
  phi_R
}

#' @export
# ----- weights
# use inv_logit function to transform weights (to map R -> [0,1])
inv_logit = function(x){
  exp(x)/(1 + exp(x))
}

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
