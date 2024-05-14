# source code for stepwise and sandwich KNG from functions_final
# but with built in convergence check

metrop = function(logA, init, nbatch = 10000, nchains = 3, scale = 1e-4, X = X,
                  ub = Inf, lb = -Inf, lower_accept = 0, upper_accept = 1,
                  update_after = 10, adjust_scale_by = 2){
  dim = ncol(X)
  
  batch = rep(list(matrix(rep(NA, nbatch*dim),nrow = nbatch, ncol = dim)), nchains)
  u = list()
  for (chain in 1:nchains) {
    u[[chain]] = matrix(runif(nbatch, min = 0, max = 1))
  }
  accept = rep(list(matrix(0, nrow = nbatch)), nchains)
  
  ct = rep(list(0), nchains)
  mean_acc = rep(list(NULL), nchains)
  scale_ans = rep(list(NULL), nchains)
  sigma = rep(list(diag(rep(1, dim))*scale), nchains)
  
  for (j in 1:nchains) {
    batch[[j]][1, ] = init[[j]]
  }
  for (i in 2:nbatch){
    for (j in 1:nchains) {
      propose_val = mvtnorm::rmvnorm(1, t(batch[[j]][i-1,]), sigma[[j]])
      check = TRUE
      check = check & min(X %*% t(propose_val)) >= lb
      check = check & max(X %*% t(propose_val)) <= ub
      
      if ((log(u[[j]][i]) < logA(c(propose_val)) - logA(batch[[j]][i-1, ])) & check){
        batch[[j]][i, ] = propose_val
        accept[[j]][i,] = 1
        ct[[j]] = ct[[j]] + 1
      } else {
        batch[[j]][i, ] = batch[[j]][i-1, ]
      }
      
      if (ct[[j]] == update_after){
        sigma[[j]] = cor(na.omit(batch[[j]])) + diag(dim)*0.00002 
        temp_acc = mean(accept[[j]][1:i])
        if (temp_acc < lower_accept){
          scale = scale/adjust_scale_by
        } else if (temp_acc > upper_accept) {
          scale = scale*adjust_scale_by
        }
        sigma[[j]] = sigma[[j]]*scale
        ct[[j]] = 0
      }
      mean_acc[[j]] = cbind(mean_acc[[j]], mean(accept[[j]][1:i]))
      scale_ans[[j]] = cbind(scale_ans[[j]], scale)
    }
    
    if (i %% 1000 == 0) {
      conv_check = coda::gelman.diag(lapply(batch, function(x) coda::as.mcmc(x[1:i, ])))[[1]][, 1]
      print(conv_check)
      if (all(conv_check < 1.1)) {
        break
        break
        break
      }

    }
    
  }
  
  batch = lapply(batch, function(x) x[1:i, ])
  accept = lapply(accept, function(x) x[1:i, ])
  
  out = list(batch = batch, accept_rate = lapply(accept, mean), scale = scale_ans, sigma = sigma)
  return(out)
}

KNG = function(init, ep, tau, X, Y, Cx, nbatch = 10000, nchains = 3, scale = 1e-2, 
               h = 1, ub = Inf, lb = -Inf, lower_accept = 0, upper_accept = 1, 
               update_after = 10, adjust_scale_by = 2){
  logA = function(theta) {
    theta = as.matrix(theta)
    n = nrow(X)
    tmp = 0
    for (i in 1:n) {
      exp_val = exp(-(X[i,] %*% theta - Y[i])/h)
      tmp = tmp + (1/(1 + exp_val) - tau) %*% x[i, ]
    }
    delta = 2*max(1-tau, tau)*Cx
    k_bar = 1/4
    gamma = n*k_bar/h
    lambda = n*k_bar*Cx^2/h
    del = 2*lambda/ep
    ans = -ep / (2*delta) * norm(tmp, "M") - max(del - gamma, 0)*norm(as.matrix(theta), "F")^2/(2*n)
    return(ans)   
  }
  
  out = metrop(logA = logA, init = init, nbatch = nbatch, nchains = nchains, scale = scale, 
               X = X, ub = ub, lb = lb, lower_accept = lower_accept, upper_accept = upper_accept,
               update_after = update_after, adjust_scale_by = adjust_scale_by)
  
  return(out)
}


constrMetrop = function(logA, init, nbatch = 10000, nchains, scale = 1e-4, check_theta, 
                        check_data, ub = -Inf, lb = Inf, method = c("fixed", "varying"), 
                        type = c("upper", "lower"), lower_accept = 0, upper_accept = 1,
                        update_after = 10, adjust_scale_by = 2){
  method = match.arg(method)
  type = match.arg(type)
  
  dim = ncol(check_data)
  
  batch = rep(list(matrix(rep(NA, nbatch*dim),nrow = nbatch, ncol = dim)), nchains)
  u = list()
  for (chain in 1:nchains) {
    u[[chain]] = matrix(runif(nbatch, min = 0, max = 1))
  }
  accept = rep(list(matrix(0, nrow = nbatch)), nchains)
  
  ct = rep(list(0), nchains)
  scale_ans = rep(list(c(scale)), nchains)
  sigma = rep(list(diag(rep(1, dim))*scale), nchains)
  
  for (j in 1:nchains) {
    batch[[j]][1, ] = init[[j]]
  }
  
  for (i in 2:nbatch){
    for (j in 1:nchains) {
      scale = scale_ans[[j]][i-1]
      check = TRUE
      if (method == "fixed") {
        propose_val = batch[[j]][i-1, ]
        propose_val[1] = propose_val[1] + rnorm(1, 0, sigma[[j]])
        check = check & ifelse(type == "lower", 
                               propose_val[1] < check_theta[1], 
                               propose_val[1] > check_theta[1])
        
      } else if (method == "varying") {
        if (dim == 1) {
          propose_val = batch[[j]][i-1, ] + rnorm(1, 0, sigma[[j]])
          
        } else {
          propose_val = t(mvtnorm::rmvnorm(1, t(batch[[j]][i-1,]), sigma[[j]]))
          
        }
        check = check & ifelse(type == "lower",
                               all(check_data %*% propose_val < check_data %*% t(check_theta)),
                               all(check_data %*% propose_val > check_data %*% t(check_theta)))
      }
      
      check = check & max(as.matrix(check_data) %*% as.matrix(propose_val))<= ub
      check = check & min(as.matrix(check_data) %*% as.matrix(propose_val)) >= lb
      
      if ((log(u[[j]][i]) < logA(c(propose_val)) - logA(batch[[j]][i-1, ])) & check){
        batch[[j]][i, ] = propose_val
        accept[[j]][i,] = 1
      } else {
        batch[[j]][i, ] = batch[[j]][i-1, ]
      }
      ct[[j]] = ct[[j]] + 1
      
      if (ct[[j]] == update_after){
        if (method == "fixed" | dim == 1) {
          sigma[[j]] = 1
        } else {
          if (all(accept[[j]][1:i,] == 0)) {
            sigma[[j]] = sigma[[j]] + diag(dim)*scale
            
          } else {
            sigma[[j]] = cor(na.omit(batch[[j]])) + diag(dim)*0.00002 
          }
         
          
        }
        temp_acc = mean(accept[[j]][1:i])
        if (temp_acc < lower_accept){
          scale = scale/adjust_scale_by
        } else if (temp_acc > upper_accept) {
          #print(temp_acc)
          scale = scale*adjust_scale_by
        }
        sigma[[j]] = sigma[[j]]*scale
        #print(sigma)
        ct[[j]] = 0
      }
      scale_ans[[j]] = cbind(scale_ans[[j]], scale)
    }
    if (i %% 1000 == 0) {
      conv_check = coda::gelman.diag(lapply(batch, function(x) coda::as.mcmc(x[1:i, ])))#[[1]][, 1]
      cat('Run', i/1000, '\n')
      print(conv_check)
      if (all(conv_check[[2]] < 1.1)) {
        break
        break
        break
      }
      
    }
  }
  
  batch = lapply(batch, function(x) x[1:i, ])
  accept = lapply(accept, function(x) x[1:i, ])
  
  out = list(batch = batch, accept_rate = lapply(accept, mean), scale = scale_ans, sigma = sigma)  
  return(out)
}




constrKNG = function(init, ep, tau, X, Y, Cx, nbatch = 10000, scale = 1e-4, h = 1,
                     check_theta, check_data = NULL, ub = Inf, lb = -Inf, nchains = 3,
                     method = c("fixed", "varying"), type = c("upper", "lower"), 
                     lower_accept = 0, upper_accept = 1, 
                     update_after = 10, adjust_scale_by = 2){
  if (is.null(check_data)){
    check_data = X
  } else {
    check_data = as.matrix(check_data)
    check_data = as.matrix(cbind(rep(1, nrow(check_data)), check_data))
  }
  
  
  logA = function(theta) {
    theta = as.matrix(theta)
    n = nrow(X)
    tmp = 0
    for (i in 1:n) {
      exp_val = exp(-(X[i,] %*% theta - Y[i])/h)
      tmp = tmp + (1/(1 + exp_val) - tau) %*% x[i, ]
    }
    delta = max(1-tau, tau)*Cx
    k_bar = 1/4
    gamma = n*k_bar/h
    lambda = n*k_bar*Cx^2/h
    del = 2*lambda/ep
    ans = -ep / (2*delta) * norm(tmp, "M") - max(del - gamma, 0)*norm(as.matrix(theta), "F")^2/(2*n)
    return(ans)   
  }
  
  out = constrMetrop(logA = logA, init = init, nbatch = nbatch, scale = scale, 
                     check_theta = check_theta, check_data = check_data,
                     ub = ub, lb = lb, method = method, type = type, nchains = nchains,
                     lower_accept = lower_accept, upper_accept = upper_accept,
                     update_after = update_after, adjust_scale_by = adjust_scale_by)
  return(out)
}


stepwiseKNG = function(data, total_eps, median_eps = NULL, tau, scale = 1e-4, Cx,
                       nbatch = 10000, method = c("fixed", "varying"), nchains = 3,
                       ub = Inf, lb = -Inf, check_data = NULL, lower_accept = 0, 
                       upper_accept = 1, update_after = 10, adjust_scale_by = 2,
                       formula = NULL, h = 1){
  method = match.arg(method)
  if(is.null(median_eps)){
    median_eps = 0.4
  }
  ep = total_eps*(1-median_eps)/(length(tau)-1)
  i = ncol(data)
  Y = data[,i]
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  X[X > Cx] = Cx
  m = ncol(X) - 1
  
  if (is.null(formula)) {
    vars = colnames(data)
    formula = paste(vars[i], " ~ .")
  }
  
  nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = 0.5)
  out = KNG(init = list(coef(nonpriv), rep(0, m+1), rep(5, m+1)), ep = total_eps*median_eps, tau = 0.5, 
            X = X, Y = Y, Cx = Cx, nbatch = nbatch, scale = scale, ub = ub, lb = lb, h = h,
            upper_accept = 0.7, lower_accept = 0.2, nchains = nchains,
            update_after = 25, adjust_scale_by = 2)
  median_theta_kng = tail(out[[1]][[1]], 1)#*R
  accept_rate = out[[2]][[1]]
  scale_output = tail(t(out[[3]][[1]]), 1)
  scale = c(scale_output)
  ans = t(median_theta_kng)
  
  tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
  tau_upper = sort(tau[tau > 0.5])
  
  if (length(tau_lower) > 0){
    check_theta = median_theta_kng
    for (i in 1:length(tau_lower)){
      new_tau = tau_lower[i]
      print(new_tau)
      out = constrKNG(init = list(check_theta, rep(-2, m+1), rep(0, m+1)), ep = ep, 
                      tau = new_tau, X = X, Y = Y, h = h,
                      Cx = Cx, nbatch = nbatch, scale = scale, check_theta = check_theta, 
                      check_data = check_data, ub = ub, lb = lb, method = method, 
                      type = "lower", lower_accept = lower_accept, upper_accept = upper_accept, 
                      update_after = update_after, adjust_scale_by = adjust_scale_by,
                      nchains = nchains)
      proposed_theta = tail(out[[1]][[1]], 1)#*R
      # message("Acceptance rate of ", new_tau, ": ", out[[2]])
      
      #check_theta = proposed_theta
      ans = cbind(t(proposed_theta), ans)
      scale_output = cbind(tail(t(out[[3]][[1]]), 1), scale_output)
      accept_rate = cbind(out[[2]][[1]], accept_rate)
      
    }
  }
  
  if (length(tau_upper) > 0){
    check_theta = median_theta_kng
    for (i in 1:length(tau_upper)){
      new_tau = tau_upper[i]
      print(new_tau)
      out = constrKNG(init = list(check_theta, check_theta + 2, check_theta + 5), 
                      ep = ep, tau = new_tau, X = X, Y = Y, h = h,
                      Cx = Cx, nbatch = nbatch, scale = scale, check_theta = check_theta, 
                      check_data = check_data, ub = ub , lb = lb, method = method, 
                      type = "upper", lower_accept = lower_accept, upper_accept = upper_accept, 
                      update_after = update_after, adjust_scale_by = adjust_scale_by,
                      nchains = nchains)
      proposed_theta = tail(out[[1]][[1]], 1)#*R
      # message("Acceptance rate of ", new_tau, ": ", out[[2]])
      
      check_theta = proposed_theta
      ans = cbind(ans, t(proposed_theta))
      scale_output = cbind(scale_output, tail(t(out[[3]][[1]]), 1))
      accept_rate = cbind(accept_rate, out[[2]][[1]])
    }
  }
  
  tau = sort(tau)
  colnames(ans) = tau
  colnames(scale_output) = tau
  rownames(scale_output) = "Scale"
  colnames(accept_rate) = tau
  
  return(list(ans, scale_output, accept_rate))
}


constrMetropSandwich = function(logA, init, nbatch = 10000, nchains = 3, scale = 1e-4, 
                                lowertheta, uppertheta, check_data, ub = Inf, lb = -Inf, 
                                method = c("fixed", "varying"), lower_accept = 0, 
                                upper_accept = 1, update_after = 10, adjust_scale_by = 2){
  
  method = match.arg(method)
  dim = ncol(check_data)
  
  batch = rep(list(matrix(rep(NA, nbatch*dim),nrow = nbatch, ncol = dim)), nchains)
  u = list()
  for (chain in 1:nchains) {
    u[[chain]] = matrix(runif(nbatch, min = 0, max = 1))
  }
  accept = rep(list(matrix(0, nrow = nbatch)), nchains)
  
  ct = rep(list(0), nchains)
  scale_ans = rep(list(c(scale)), nchains)
  sigma = rep(list(diag(rep(1, dim))*scale), nchains)
  
  for (j in 1:nchains) {
    batch[[j]][1, ] = init[[j]]
  }
  
  for (i in 2:nbatch){
    for (j in 1:nchains) {
    scale = scale_ans[[j]][i-1]
    check = TRUE
    
    if (method == "fixed") {
      propose_val = batch[[j]][i-1, ]
      propose_val[1] = propose_val[1] + rnorm(1, 0, sigma[[j]])
      check = check & (propose_val[1] < uppertheta[1]) & (propose_val[1] > lowertheta[1])
      
    } else if (method == "varying") {
      if (dim == 1){
        propose_val = batch[[j]][i-1, ] + rnorm(1, 0, sigma[[j]])
      } else {
        propose_val = t(mvtnorm::rmvnorm(1, t(batch[[j]][i-1,]), sigma[[j]]))
      }
      
      check = check & 
        all(check_data %*% propose_val < check_data %*% uppertheta) &
        all(check_data %*% propose_val > check_data %*% lowertheta)
    }
    
    check = check & (max(as.matrix(check_data) %*% as.matrix(propose_val)) <= ub)
    check = check & (min(as.matrix(check_data) %*% as.matrix(propose_val)) >= lb)
    
    if ((log(u[[j]][i]) < logA(c(propose_val)) - logA(batch[[j]][i-1, ])) & check){
      batch[[j]][i, ] = propose_val
      accept[[j]][i,] = 1
      ct[[j]] = ct[[j]] + 1
    } else {
      batch[[j]][i, ] = batch[[j]][i-1, ]
    }
    
    if (ct[[j]] == update_after){
      if (method == "fixed" | dim == 1) {
        sigma[[j]] = 1
      } else {
        sigma[[j]] = cor(na.omit(batch[[j]])) + diag(dim)*0.00002 
        
      }
      temp_acc = mean(accept[[j]][1:i])
      if (temp_acc < lower_accept){
        scale = scale/adjust_scale_by
      } else if (temp_acc > upper_accept) {
        scale = scale*adjust_scale_by
      }
      sigma[[j]] = sigma[[j]]*scale
      ct[[j]] = 0
    }
    scale_ans[[j]] = cbind(scale_ans[[j]], scale)
    }
    
    if (i %% 1000 == 0) {
      conv_check = coda::gelman.diag(lapply(batch, function(x) coda::as.mcmc(x[1:i, ])))#[[1]][, 1]
      cat('Run', i/1000, '\n')
      print(conv_check)
      if (all(conv_check[[2]] < 1.1)) {
        break
        break
        break
      }
      
    }
  }
  
  batch = lapply(batch, function(x) x[1:i, ])
  accept = lapply(accept, function(x) x[1:i, ])
  
  out = list(batch = batch, accept_rate = lapply(accept, mean), scale = scale_ans, sigma = sigma)
  return(out)
}


constrKNGSandwich = function(init, ep, tau, X, Y, Cx, nbatch = 1000, scale = 1e-4,
                             lowertheta, uppertheta, check_data = NULL, ub = Inf, 
                             lb = -Inf, method = c("fixed", "varying"), h,
                             lower_accept = 0, upper_accept = 1, nchains = 3,
                             update_after = 10, adjust_scale_by = 2){
  if (is.null(check_data)){
    check_data = X
  } else {
    check_data = as.matrix(check_data)
    check_data = as.matrix(cbind(rep(1, nrow(check_data)), check_data))
  }
  
  logA = function(theta) {
    theta = as.matrix(theta)
    n = nrow(X)
    tmp = 0
    for (i in 1:n) {
      exp_val = exp(-(X[i,] %*% theta - Y[i])/h)
      tmp = tmp + (1/(1 + exp_val) - tau) %*% x[i, ]
    }
    delta = 2*max(1-tau, tau)*Cx
    k_bar = 1/4
    gamma = n*k_bar/h
    lambda = n*k_bar*Cx^2/h
    del = 2*lambda/ep
    ans = -ep / (2*delta) * norm(tmp, "M") - max(del - gamma, 0)*norm(as.matrix(theta), "F")^2/(2*n)
    return(ans)   
  }
  
  out = constrMetropSandwich(logA = logA, init = init, nbatch = nbatch, scale = scale,
                             lowertheta = lowertheta, uppertheta = uppertheta, nchains = nchains,
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
                       formula = NULL, nchains = 3, h = 1){
  
  main_tau_fac = as.factor(main_tau)
  sandwich_tau = tau[!tau %in% main_tau_fac]
  main_tau_order = tau[tau %in% main_tau_fac]
  
  if(is.null(main_tau_eps)){
    main_tau_eps = ifelse(method == "fixed", 0.7, 0.8)
  }
  
  data = as.matrix(data)
  i = ncol(data)
  Y = as.matrix(data[,i])
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  X[X > Cx] = Cx
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
                    formula = formula, h = h, nchains = nchains)
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
      print(curr_tau)
      update_tau = sort(c(update_tau, curr_tau))
      idx = which(update_tau == curr_tau)
      lowertau = update_tau[idx-1]
      uppertau = update_tau[idx+1]
      
      lowertheta = b[, which(colnames(b) == lowertau)]
      uppertheta = b[, which(colnames(b) == uppertau)]
      #main_tau_eps
      out = constrKNGSandwich(init = list(uppertheta, lowertheta, (uppertheta + lowertheta)/2), 
                              ep = eps_sandwich, tau = curr_tau, 
                              X = X, Y = Y, Cx = Cx, nbatch = nbatch, nchains = nchains, h = h, 
                              scale = sw_scale, lowertheta = lowertheta, uppertheta = uppertheta, 
                              check_data = check_data, ub = ub, lb = lb, method = method,
                              lower_accept = lower_accept, upper_accept = upper_accept, 
                              update_after = update_after, adjust_scale_by = adjust_scale_by)
      
      b_sandwich = tail(out[[1]][[1]], 1)
      sandwich_scale_lower[i] = tail(t(out[[3]][[1]]), 1)
      accept_rate_lower[i] = out[[2]][[1]]
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
      print(curr_tau)
      curr_tau = sandwich_tau_upper[i]
      update_tau = sort(c(update_tau, curr_tau))
      idx = which(update_tau == curr_tau)
      lowertau = update_tau[idx-1]
      uppertau = update_tau[idx+1]
      
      lowertheta = b[, which(colnames(b) == lowertau)]
      uppertheta = b[, which(colnames(b) == uppertau)]
      #good sandwich startscale 1/R
      out = constrKNGSandwich(init = list(lowertheta, uppertheta, (lowertheta+uppertheta)/2), 
                              ep = eps_sandwich, tau = curr_tau, 
                              X = X, Y = Y, Cx = Cx, nbatch = nbatch, nchains = nchains,
                              scale = sw_scale, lowertheta = lowertheta, uppertheta = uppertheta, 
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
