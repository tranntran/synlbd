metrop = function(logA, init, nbatch = 10000, scale = 1e-4, X = X,
                  ub = Inf, lb = -Inf, lower_accept = 0, upper_accept = 1,
                  update_after = 10, adjust_scale_by = 2){
  dim = length(init)
  
  batch = matrix(rep(0, nbatch*dim),nrow = nbatch, ncol = dim)
  u = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
  accept = matrix(0, nrow = nbatch)
  
  ct = 0
  mean_acc = NA
  scale_ans = NA
  sigma = diag(rep(1, dim))*scale
  
  batch[1, ] = init
  for (i in 2:nbatch){
    propose_val = mvtnorm::rmvnorm(1, t(batch[i-1,]), sigma)
    check = TRUE
    check = check & min(X %*% t(propose_val)) >= lb
    check = check & max(X %*% t(propose_val)) <= ub
    
    if ((log(u[i]) < logA(c(propose_val)) - logA(batch[i-1, ])) & check){
      batch[i, ] = propose_val
      accept[i,] = 1
      ct = ct + 1
    } else {
      batch[i, ] = batch[i-1, ]
    }
    
    if (ct == update_after){
      sigma = cor(na.omit(batch)) + diag(dim)*0.00002 
      temp_acc = mean(accept[1:i])
      if (temp_acc < lower_accept){
        scale = scale/adjust_scale_by
      } else if (temp_acc > upper_accept) {
        scale = scale*adjust_scale_by
      }
      sigma = sigma*scale
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i]))
    scale_ans = cbind(scale_ans, scale)
  }
  mean_acc = mean_acc[-1]
  scale_ans = scale_ans[-1]
  
  out = list(batch = batch, accept_rate = mean(accept), scale = scale_ans, sigma = sigma)
  return(out)
}


KNG = function(init, ep, tau, sumX, X, Y, Cx, nbatch = 10000, scale = 1e-4, ub = Inf, 
               lb = -Inf, lower_accept = 0, upper_accept = 1, update_after = 10, 
               adjust_scale_by = 2){
  logA = function(beta){
    left = cbind(Y, X)%*% c(1, -beta)
    lessEq = (left <= 0)
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*Cx) - (1e-5)*(beta%*%beta)
    return(ans)
  }
  
  out = metrop(logA = logA, init = init, nbatch = nbatch, scale = scale, X = X, 
               ub = ub, lb = lb, lower_accept = lower_accept, upper_accept = upper_accept,
               update_after = update_after, adjust_scale_by = adjust_scale_by)
  
  return(out)
}


constrMetrop = function(logA, init, nbatch = 10000, scale = 1e-4, check_beta, check_data, 
                        ub = -Inf, lb = Inf, method = c("fixed", "varying"), 
                        type = c("upper", "lower"), lower_accept = 0, upper_accept = 1,
                        update_after = 10, adjust_scale_by = 2){
  method = match.arg(method)
  type = match.arg(type)
  
  dim = length(init)
  
  batch = matrix(rep(0, nbatch*dim),nrow = nbatch, ncol = dim)
  u = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
  accept = matrix(0, nrow = nbatch)
  
  ct = 0
  mean_acc = NA
  scale_ans = NA
  sigma = diag(rep(1, dim))*scale
  
  batch[1, ] = init
  for (i in 2:nbatch){
    check = TRUE
    
    if (method == "fixed") {
      propose_val = batch[i-1, ]
      propose_val[1] = propose_val[1] + rnorm(1, 0, sigma)
      check = check & ifelse(type == "lower", propose_val[1] < check_beta[1], propose_val[1] > check_beta[1])
      
    } else if (method == "varying") {
      if (dim == 1) {
        propose_val = batch[i-1, ] + rnorm(1, 0, sigma)
        
      } else {
        propose_val = t(mvtnorm::rmvnorm(1, t(batch[i-1,]), sigma))
        
      }
      check = check & ifelse(type == "lower",
                             all(check_data %*% propose_val < check_data %*% t(check_beta)),
                             all(check_data %*% propose_val > check_data %*% t(check_beta)))
    }
    
    check = check & max(as.matrix(check_data) %*% as.matrix(propose_val))<= ub
    check = check & min(as.matrix(check_data) %*% as.matrix(propose_val)) >= lb
    
    if ((log(u[i]) < logA(c(propose_val)) - logA(batch[i-1, ])) & check){
      batch[i, ] = propose_val
      accept[i,] = 1
      ct = ct + 1
    } else {
      batch[i, ] = batch[i-1, ]
    }
    
    if (ct == update_after){
      if (method == "fixed" | dim == 1) {
        sigma = 1
      } else {
        sigma = cor(na.omit(batch)) + diag(dim)*0.00002 
        
      }
      temp_acc = mean(accept[1:i])
      if (temp_acc < lower_accept){
        scale = scale/adjust_scale_by
      } else if (temp_acc > upper_accept) {
        #print(temp_acc)
        scale = scale*adjust_scale_by
      }
      sigma = sigma*scale
      #print(sigma)
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i]))
    scale_ans = cbind(scale_ans, scale)
  }
  mean_acc = mean_acc[-1]
  scale_ans = scale_ans[-1]
  
  out = list(batch = batch, accept_rate = mean(accept), scale = scale_ans, sigma = sigma)
  
  return(out)
}




constrKNG = function(init, ep, tau, sumX, X, Y, Cx, nbatch = 10000, scale = 1e-4,
                     check_beta, check_data = NULL, ub = Inf, lb = -Inf, 
                     method = c("fixed", "varying"), type = c("upper", "lower"), 
                     lower_accept = 0, upper_accept = 1, 
                     update_after = 10, adjust_scale_by = 2){
  if (is.null(check_data)){
    check_data = X
  } else {
    check_data = as.matrix(check_data)
    check_data = as.matrix(cbind(rep(1, nrow(check_data)), check_data))
  }
  
  
  logA = function(beta) {
    left = cbind(Y, X) %*% c(1,-beta)
    lessEq = (left <= 0)
    # 1e-5
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*Cx) - (1e-5)*(beta%*%beta)
    return(ans)
  }
  
  out = constrMetrop(logA = logA, init = init, nbatch = nbatch, scale = scale,
                     check_beta = check_beta, check_data = check_data,
                     ub = ub, lb = lb, method = method, type = type, 
                     lower_accept = lower_accept, upper_accept = upper_accept,
                     update_after = update_after, adjust_scale_by = adjust_scale_by)
  return(out)
}


stepwiseKNG = function(data, total_eps, median_eps = NULL, tau, scale = 1e-4, Cx,
                       nbatch = 10000, method = c("fixed", "varying"), 
                       ub = Inf, lb = -Inf, check_data = NULL, lower_accept = 0, 
                       upper_accept = 1, update_after = 10, adjust_scale_by = 2,
                       formula = NULL){
  method = match.arg(method)
  if(is.null(median_eps)){
    median_eps = ifelse(method == "fixed", 0.8, 0.9)
  }
  ep = total_eps*(1-median_eps)/(length(tau)-1)
  i = ncol(data)
  Y = data[,i]
  #R = max(abs(Y))
  #Y = Y/R
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  X[X > Cx] = Cx
  sumX = apply(X = X, 2, FUN = sum)
  m = ncol(X) - 1
  
  if (is.null(formula)) {
    vars = colnames(data)
    formula = paste(vars[i], " ~ .")
  }
  
  nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = 0.5)
  out = KNG(init = coef(nonpriv), ep = total_eps*median_eps, tau = 0.5, sumX = sumX, 
            X = X, Y = Y, Cx = Cx, nbatch = nbatch, scale = 1e-4, ub = ub, lb = lb, 
            upper_accept = upper_accept, lower_accept = lower_accept, 
            update_after = update_after, adjust_scale_by = adjust_scale_by)
  median_beta_kng = tail(out[[1]], 1)#*R
  accept_rate = out[[2]]
  scale_output = tail(out[[3]], 1)
  ans = t(median_beta_kng)
  
  tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
  tau_upper = sort(tau[tau > 0.5])
  
  if (length(tau_lower) > 0){
    check_beta = median_beta_kng
    for (i in 1:length(tau_lower)){
      new_tau = tau_lower[i]
      out = constrKNG(init = check_beta, ep = ep, tau = new_tau, sumX = sumX, X = X, Y = Y,
                      Cx = Cx, nbatch = nbatch, scale = scale, check_beta = check_beta, 
                      check_data = check_data, ub = ub, lb = lb, method = method, 
                      type = "lower", lower_accept = lower_accept, upper_accept = upper_accept, 
                      update_after = update_after, adjust_scale_by = adjust_scale_by)
      proposed_beta = tail(out[[1]], 1)#*R
      # message("Acceptance rate of ", new_tau, ": ", out[[2]])
      
      check_beta = proposed_beta
      ans = cbind(t(proposed_beta), ans)
      scale_output = cbind(tail(out[[3]], 1), scale_output)
      accept_rate = cbind(out[[2]], accept_rate)
      
    }
  }
  
  if (length(tau_upper) > 0){
    check_beta = median_beta_kng
    for (i in 1:length(tau_upper)){
      new_tau = tau_upper[i]
      out = constrKNG(init = check_beta, ep = ep, tau = new_tau, sumX = sumX, X = X, Y = Y, 
                      Cx = Cx, nbatch = nbatch, scale = scale, check_beta = check_beta, 
                      check_data = check_data, ub = ub , lb = lb, method = method, 
                      type = "upper", lower_accept = lower_accept, upper_accept = upper_accept, 
                      update_after = update_after, adjust_scale_by = adjust_scale_by)
      proposed_beta = tail(out[[1]], 1)#*R
      # message("Acceptance rate of ", new_tau, ": ", out[[2]])
      
      check_beta = proposed_beta
      ans = cbind(ans, t(proposed_beta))
      scale_output = cbind(scale_output, tail(out[[3]], 1))
      accept_rate = cbind(accept_rate, out[[2]])
    }
  }
  
  tau = sort(tau)
  colnames(ans) = tau
  colnames(scale_output) = tau
  rownames(scale_output) = "Scale"
  colnames(accept_rate) = tau
  
  return(list(ans, scale_output, accept_rate))
}


constrMetropSandwich = function(logA, init, nbatch = 10000, scale = 1e-4, lowerbeta, upperbeta, 
                                check_data, ub = Inf, lb = -Inf, method = c("fixed", "varying"),
                                lower_accept = 0, upper_accept = 1, update_after = 10, 
                                adjust_scale_by = 2){
  
  method = match.arg(method)
  dim = length(init)
  
  batch = matrix(rep(0, nbatch*dim),nrow = nbatch, ncol = dim)
  u = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
  accept = matrix(0, nrow = nbatch)
  
  ct = 0
  mean_acc = NA
  scale_ans = NA
  sigma = diag(rep(1, dim))*scale
  
  batch[1, ] = init
  for (i in 2:nbatch){
    check = TRUE
    
    if (method == "fixed") {
      propose_val = batch[i-1, ]
      propose_val[1] = propose_val[1] + rnorm(1, 0, sigma)
      check = check & (propose_val[1] < upperbeta[1]) & (propose_val[1] > lowerbeta[1])
      
    } else if (method == "varying") {
      if (dim == 1){
        propose_val = batch[i-1, ] + rnorm(1, 0, sigma)
      } else {
        propose_val = t(mvtnorm::rmvnorm(1, t(batch[i-1,]), sigma))
      }

      check = check & 
        all(check_data %*% propose_val < check_data %*% upperbeta) &
        all(check_data %*% propose_val > check_data %*% lowerbeta)
    }
    
    check = check & (max(as.matrix(check_data) %*% as.matrix(propose_val)) <= ub)
    check = check & (min(as.matrix(check_data) %*% as.matrix(propose_val)) >= lb)
    
    if ((log(u[i]) < logA(c(propose_val)) - logA(batch[i-1, ])) & check){
      batch[i, ] = propose_val
      accept[i,] = 1
      ct = ct + 1
    } else {
      batch[i, ] = batch[i-1, ]
    }
    
    if (ct == update_after){
      if (method == "fixed" | dim == 1) {
        sigma = 1
      } else {
        sigma = cor(na.omit(batch)) + diag(dim)*0.00002 
        
      }
      temp_acc = mean(accept[1:i])
      if (temp_acc < lower_accept){
        scale = scale/adjust_scale_by
      } else if (temp_acc > upper_accept) {
        scale = scale*adjust_scale_by
      }
      sigma = sigma*scale
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i]))
    scale_ans = cbind(scale_ans, scale)
  }
  mean_acc = mean_acc[-1]
  scale_ans = scale_ans[-1]
  
  out = list(batch = batch, accept_rate = mean(accept), scale = scale_ans, sigma = sigma)
  return(out)
}


constrKNGSandwich = function(init, ep, tau, sumX, X, Y, Cx, nbatch = 1000, scale = 1e-4,
                             lowerbeta, upperbeta, check_data = NULL, ub = Inf, 
                             lb = -Inf, method = c("fixed", "varying"),
                             lower_accept = 0, upper_accept = 1, 
                             update_after = 10, adjust_scale_by = 2){
  if (is.null(check_data)){
    check_data = X
  } else {
    check_data = as.matrix(check_data)
    check_data = as.matrix(cbind(rep(1, nrow(check_data)), check_data))
  }
  
  logA = function(beta) {
    left = cbind(Y, X) %*% c(1,-beta)
    lessEq = (left <= 0)
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*Cx) - (1e-5)*(beta%*%beta)
    return(ans)
  }
  
  out = constrMetropSandwich(logA = logA, init = init, nbatch = nbatch, scale = scale,
                             lowerbeta = lowerbeta, upperbeta = upperbeta, 
                             check_data = check_data, ub = ub, lb = lb, method = method,
                             lower_accept = lower_accept, upper_accept = upper_accept,
                             update_after = update_after, adjust_scale_by = adjust_scale_by)
  return(out)
}


sandwichKNG = function(data, total_eps, median_eps = NULL, main_tau_eps = NULL,
                       tau, main_tau, scale = 1e-4, sw_scale = 1e-4, Cx,
                       nbatch = 10000, method = c("fixed", "varying"), 
                       ub = Inf, lb = -Inf, check_data = NULL, 
                       lower_accept = 0, upper_accept = 1, 
                       update_after = 10, adjust_scale_by = 2, 
                       formula = NULL){
  
  main_tau_fac = as.factor(main_tau)
  sandwich_tau = tau[!tau %in% main_tau_fac]
  main_tau_order = tau[tau %in% main_tau_fac]
  
  if(is.null(main_tau_eps)){
    main_tau_eps = ifelse(method == "fixed", 0.8, 0.9)
  }
  
  data = as.matrix(data)
  i = ncol(data)
  Y = as.matrix(data[,i])
  #R = max(abs(Y))
  #Y = Y/R
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  X[X > Cx] = Cx
  sumX = apply(X = X, 2, FUN = sum)
  m = ncol(X) - 1
  
  if (is.null(formula)) {
    vars = colnames(data)
    formula = paste(vars[i], " ~ .")
  }
  
  out = stepwiseKNG(data = data, total_eps = total_eps*main_tau_eps, median_eps = median_eps, 
                    tau = main_tau_order, scale = scale, nbatch = nbatch, Cx = Cx,
                    method = method, ub = ub, lb = lb, check_data = check_data, 
                    lower_accept = lower_accept, upper_accept = upper_accept, 
                    update_after = update_after, adjust_scale_by = adjust_scale_by,
                    formula = formula)
  b = out[[1]]
  scale_output = out[[2]]
  accept_rate = out[[3]]
  main_tau = sort(main_tau)
  names(scale_output) = main_tau
  names(accept_rate) = main_tau
  colnames(b) = main_tau
  
  eps_sandwich = (1 - main_tau_eps) * total_eps / length(sandwich_tau)
  
  update_tau = main_tau
  sandwich_tau_lower = sort(sandwich_tau[sandwich_tau < 0.5], decreasing = TRUE)
  sandwich_tau_upper = sort(sandwich_tau[sandwich_tau > 0.5])
  if (length(sandwich_tau_lower) > 0){
    sandwich_scale_lower = rep(NA, length(sandwich_tau_lower))
    names(sandwich_scale_lower) = sandwich_tau_lower
    accept_rate_lower = rep(NA, length(sandwich_tau_lower))
    names(accept_rate_lower) = sandwich_tau_lower
    for (i in 1:length(sandwich_tau_lower)){
      curr_tau = sandwich_tau_lower[i]
      update_tau = sort(c(update_tau, curr_tau))
      idx = which(update_tau == curr_tau)
      lowertau = update_tau[idx-1]
      uppertau = update_tau[idx+1]
      
      lowerbeta = b[, which(colnames(b) == lowertau)]
      upperbeta = b[, which(colnames(b) == uppertau)]
      #main_tau_eps
      out = constrKNGSandwich(init = upperbeta, ep = eps_sandwich, tau = curr_tau, 
                              sumX = sumX, X = X, Y = Y, Cx = Cx, nbatch = nbatch, 
                              scale = sw_scale, lowerbeta = lowerbeta, upperbeta = upperbeta, 
                              check_data = check_data, ub = ub, lb = lb, method = method,
                              lower_accept = lower_accept, upper_accept = upper_accept, 
                              update_after = update_after, adjust_scale_by = adjust_scale_by)
      
      b_sandwich = tail(out[[1]], 1)
      sandwich_scale_lower[i] = tail(out[[3]], 1)
      accept_rate_lower[i] = out[[2]]
      b = cbind(b, t(b_sandwich))
      colnames(b)[ncol(b)] = curr_tau
      if(m == 0){
        b = t(as.matrix(b[, match(update_tau, colnames(b))]))
      } else {
        b = as.matrix(b[, match(update_tau, colnames(b))])
      }
    }
  }
  
  if (length(sandwich_tau_upper) > 0){
    sandwich_scale_upper = rep(NA, length(sandwich_tau_upper))
    names(sandwich_scale_upper) = sandwich_tau_upper
    accept_rate_upper = rep(NA, length(sandwich_tau_upper))
    names(accept_rate_upper) = sandwich_tau_upper
    for (i in 1:length(sandwich_tau_upper)){
      curr_tau = sandwich_tau_upper[i]
      update_tau = sort(c(update_tau, curr_tau))
      idx = which(update_tau == curr_tau)
      lowertau = update_tau[idx-1]
      uppertau = update_tau[idx+1]
      
      lowerbeta = b[, which(colnames(b) == lowertau)]
      upperbeta = b[, which(colnames(b) == uppertau)]
      #good sandwich startscale 1/R
      out = constrKNGSandwich(init = lowerbeta, ep = eps_sandwich, tau = curr_tau, 
                              sumX = sumX, X = X, Y = Y, Cx = Cx, nbatch = nbatch, 
                              scale = sw_scale, lowerbeta = lowerbeta, upperbeta = upperbeta, 
                              check_data = check_data, ub = ub, lb = lb, method = method,
                              lower_accept = lower_accept, upper_accept = upper_accept, 
                              update_after = update_after, adjust_scale_by = adjust_scale_by)
      b_sandwich = tail(out[[1]], 1)
      sandwich_scale_upper[i] = tail(out[[3]], 1)
      accept_rate_upper[i] = out[[2]]
      b = cbind(b, t(b_sandwich))
      colnames(b)[ncol(b)] = curr_tau
      if(m == 0){
        b = t(as.matrix(b[, match(update_tau, colnames(b))]))
      } else {
        b = as.matrix(b[, match(update_tau, colnames(b))])
      }
    }
  }
  
  #add case when there are no lower/upper quantiles
  scale_output = c(scale_output, sandwich_scale_upper, sandwich_scale_lower)
  scale_output = scale_output[order(names(scale_output))]
  accept_rate = c(accept_rate, accept_rate_lower, accept_rate_upper)
  accept_rate = accept_rate[order(names(accept_rate))]
  
  return(list(b, scale_output, accept_rate))
  
}
