# This file contains the code for all methods used in comparing the utility of 
# synthetic data using simulation


### Note on parameters
# max_x = max(L_inf norm of x_i for i in 1...n)

### Loss functions

check_loss = function(theta, x, y, tau) {
  tmp = y - x%*%theta
  if (tmp <= 0) {
    ans = (tau-1)*tmp
  } else {
    ans = tau*tmp
  }
  return(ans)
}

check_loss_gradient = function(theta, x, y, tau) {
  ans = -tau*x 
  if (y - x%*%theta <= 0) {
    ans = ans + x
  }
  return(ans)
}

smooth_loss = function(theta, x, y, tau, alpha) {
  ans = tau*(y - x%*%theta) + alpha*log(1 + exp(-(y - x%*%theta)/alpha))
  return(ans)
}

smooth_loss_gradient = function(theta, x, y, tau, alpha) {
  ans = c(1/(1 + exp((y - x%*%theta)/alpha)) - tau) %*% x
  return(ans)
}

### Function for KNG with check loss
# Since it is not twice differentiable, the utility is not guaranteed.
kng_qr = function(theta, x, y, tau, epsilon = 1, max_x = 1) {
  n = nrow(x) 
  tmp = 0
  for (i in 1:n) {
    tmp = tmp - tau%*%x[i, ] + I(y[i]<= x[i, ]%*%theta)%*%x[i,]
  }
  delta = 2*(1-tau)*max_x
  ans = -epsilon / (2*delta) * norm(tmp, "M")
  return(ans)
}

estimate_kng_qr = function(data, total_eps, tau, max_x, nsteps = 100000, nchains = 4) {
  ep = total_eps/length(tau)
  n = nrow(data)
  i = ncol(data)
  y = matrix(data[, i], nrow = n)
  x = as.matrix(cbind(rep(1, nrow(data)), data))
  x = as.matrix(x[, -ncol(x)])
  x[x > max_x] = max_x
  ans = NULL
  start_points = matrix(0, ncol = ncol(x), nrow = nchains)
  for (t in tau) {
    tmp = fmcmc::MCMC(start_points, kng_qr, nsteps = nsteps, nchains = nchains, 
                      x = x, y = y, tau = t, epsilon = ep, max_x = max_x,
                      conv_checker = fmcmc::convergence_gelman(), 
                      kernel = fmcmc::kernel_ram())
    ans = rbind(ans, tail(as.matrix(tmp[[1]]), 1))
  }
  ans = t(ans)
  colnames(ans) = tau
  
  return(list(par = ans, tau = tau))
}

### Function for KNG with check loss and penalty terms
# The penalty term is added to regularize the loss function and makes it twice
# differentiable. The rate of penalty is currently set as 1/2, based on common 
# choices in the literature, but can be changed to another number.
# The penalty term helps KNG's performance become more stable, but incur higher bias.
kng_qr_reg = function(theta, x, y, tau, epsilon = 1, max_x = 1) {
  n = nrow(x) 
  tmp = 0
  for (i in 1:n) {
    tmp = tmp - tau%*%x[i, ] + I(y[i]<= x[i, ]%*%theta)%*%x[i,]
  }
  delta = 2*(1-tau)*max_x
  ans = -epsilon / (2*delta) * norm(tmp, "M") - 1/2*theta%*%theta
  return(ans)
}


estimate_kng_qr_reg = function(data, total_eps, tau, max_x, nsteps = 100000, 
                               nchains = 4) {
  ep = total_eps/length(tau)
  n = nrow(data)
  i = ncol(data)
  y = matrix(data[, i], nrow = n)
  x = as.matrix(cbind(rep(1, nrow(data)), data))
  x = as.matrix(x[, -ncol(x)])
  x[x > max_x] = max_x
  ans = NULL
  start_points = matrix(0, ncol = ncol(x), nrow = nchains)
  for (t in tau) {
    tmp = fmcmc::MCMC(start_points, kng_qr_reg, nsteps = nsteps, nchains = nchains, 
                      x = x, y = y, tau = t, epsilon = ep, max_x = max_x,
                      conv_checker = fmcmc::convergence_gelman(), 
                      kernel = fmcmc::kernel_ram())
    ans = rbind(ans, tail(as.matrix(tmp[[1]]), 1))
  }
  ans = t(ans)
  colnames(ans) = tau
  
  return(list(par = ans, tau = tau))
}


### Function for KNG with smooth loss (based on loss function of Song et al, 2012)
# Add a smoothing parameter alpha to make the loss function strong convex and 
# twice differentiable, resulting in the utility guarantee from KNG. As alpha 
# goes to 0, the smooth loss function will converge to the check loss. Higher 
# smoothing parameter will result in higher level of bias.
kng_qr_smooth = function(theta, x, y, tau, alpha, epsilon = 1, max_x = 1) {
  n = nrow(x)
  tmp = 0
  for (i in 1:n) {
    exp_val = exp(( y[i] - x[i, ]%*%theta)/alpha)
    tmp = tmp + (1 / (1 + exp_val) - tau) %*% x[i, ]
  }
  delta = max(abs(1-tau), abs(-tau))*max_x
  ans = -epsilon / (2*delta) * norm(tmp, "M")
  return(ans)
}

estimate_kng_qr_smooth = function(data, total_eps, tau, alpha, max_x, 
                                  nsteps = 100000, nchains = 4) {
  ep = total_eps/length(tau)
  n = nrow(data)
  i = ncol(data)
  y = matrix(data[, i], nrow = n)
  x = as.matrix(cbind(rep(1, nrow(data)), data))
  x = as.matrix(x[, -ncol(x)])
  x[x > max_x] = max_x
  ans = NULL
  start_points = matrix(0, ncol = ncol(x), nrow = nchains)
  for (t in tau) {
    tmp = fmcmc::MCMC(start_points, kng_qr_smooth, nsteps = nsteps, nchains = nchains, 
                      x = x, y = y, tau = t, alpha = alpha, epsilon = ep, 
                      max_x = max_x, conv_checker = fmcmc::convergence_gelman(), 
                      kernel = fmcmc::kernel_ram())
    ans = rbind(ans, tail(as.matrix(tmp[[1]]), 1))
  }
  ans = t(ans)
  colnames(ans) = tau
  
  return(list(par = ans, tau = tau))
}

### Function for OPM with smooth loss (based on loss function of Song et al, 2012)
# Traditionally, OPM is not used to generate private quantile regression estimates
# because of the twice differentiable requirement for loss function. The smooth 
# loss function helps mitigate this problem. However, as noted above, the higher
# the smoothing variable is, the higher the resulting bias is.
opm_qr = function(theta, x, y, tau, alpha, epsilon = 1, max_x = 1) {
  theta = matrix(theta, ncol = 1)
  n = nrow(x)
  p = ncol(x)
  xi = max(abs(1-tau), abs(tau))*sqrt(p)*max_x
  lambda = (p*max_x^2) / (4*alpha)
  
  delta = 2*lambda/epsilon
  
  # sample b
  norm = rgamma(n = 1, shape = p, rate = epsilon/(2*xi))
  vector = rnorm(n = p, m = 0, s = 1)
  b = vector*norm/sqrt(sum(vector^2))
  
  ans = 0
  for (i in 1:n){
    score_val = y[i] - x[i, ] %*% theta
    ans = ans + tau*score_val + alpha*log(1+ exp(-score_val/alpha))
  }
  ans = ans/n + delta/(2*n)*norm(theta, "F")^2 + b%*%theta/n
  
  return(ans)
}



estimate_opm_qr = function(data, total_eps, tau, alpha, max_x) {
  ep = total_eps/length(tau)
  n = nrow(data)
  i = ncol(data)
  y = matrix(data[, i], nrow = n)
  x = as.matrix(cbind(rep(1, nrow(data)), data))
  x = as.matrix(x[, -ncol(x)])
  x[x > max_x] = max_x
  ans = NULL
  start_points = rep(0, ncol(x))
  for (t in tau) {
    tmp = optim(start_points, opm_qr, x = x, y = y, tau = t, alpha = alpha, 
                epsilon = ep)$par
    ans = rbind(ans, tmp)
  }
  ans = t(ans)
  colnames(ans) = tau
  
  return(list(par = ans, tau = tau))
}

opm_qr1 = function(theta, x, y, tau, alpha, epsilon = 1, max_x = 1) {
  theta = matrix(theta, ncol = 1)
  n = nrow(x)
  p = ncol(x)
  xi = max(abs(1-tau), abs(tau))*sqrt(p)*max_x
  lambda = (p*max_x^2) / (4*alpha)
  
  delta = 2*lambda/epsilon
  
  # sample b
  norm = rgamma(n = 1, shape = p, rate = epsilon/(2*xi))
  vector = rnorm(n = p, m = 0, s = 1)
  b = vector*norm/sqrt(sum(vector^2))
  
  ans = 0
  for (i in 1:n){
    score_val = y[i] - x[i, ] %*% theta
    ans = ans + tau*score_val + alpha*log(1+ exp(-score_val/alpha))
  }
  # since loss function only involves sum, scale the resulting optimizing func
  # back up by n?
  ans = ans + delta/(2)*norm(theta, "F")^2 + b%*%theta
  
  return(ans)
}

estimate_opm_qr1 = function(data, total_eps, tau, alpha, max_x) {
  ep = total_eps/length(tau)
  n = nrow(data)
  i = ncol(data)
  y = matrix(data[, i], nrow = n)
  x = as.matrix(cbind(rep(1, nrow(data)), data))
  x = as.matrix(x[, -ncol(x)])
  x[x > max_x] = max_x
  ans = NULL
  start_points = rep(0, ncol(x))
  for (t in tau) {
    tmp = optim(start_points, opm_qr1, x = x, y = y, tau = t, alpha = alpha, 
                epsilon = ep)$par
    ans = rbind(ans, tmp)
  }
  ans = t(ans)
  colnames(ans) = tau
  
  return(list(par = ans, tau = tau))
}


opm_qr2 = function(theta, x, y, tau, alpha, epsilon = 1, max_x = 1) {
  theta = matrix(theta, ncol = 1)
  n = nrow(x)
  p = ncol(x)
  xi = max(abs(1-tau), abs(tau))*sqrt(p)*max_x
  lambda = (p*max_x^2) / (4*alpha)
  
  delta = 2*lambda/epsilon
  
  # sample b
  norm = rgamma(n = 1, shape = p, rate = epsilon/(2*xi))
  vector = rnorm(n = p, m = 0, s = 1)
  b = vector*norm/sqrt(sum(vector^2))
  
  ans = 0
  for (i in 1:n){
    score_val = y[i] - x[i, ] %*% theta
    ans = ans + tau*score_val + alpha*log(1+ exp(-score_val/alpha))
  }
  # since loss function only involves sum, scale the resulting optimizing func
  # back up by n?
  ans = ans + delta/(2*n)*norm(theta, "F")^2 + b%*%theta/n
  
  return(ans)
}


estimate_opm_qr2 = function(data, total_eps, tau, alpha, max_x) {
  ep = total_eps/length(tau)
  n = nrow(data)
  i = ncol(data)
  y = matrix(data[, i], nrow = n)
  x = as.matrix(cbind(rep(1, nrow(data)), data))
  x = as.matrix(x[, -ncol(x)])
  x[x > max_x] = max_x
  ans = NULL
  start_points = rep(0, ncol(x))
  for (t in tau) {
    tmp = optim(start_points, opm_qr2, x = x, y = y, tau = t, alpha = alpha, 
                epsilon = ep)$par
    ans = rbind(ans, tmp)
  }
  ans = t(ans)
  colnames(ans) = tau
  
  return(list(par = ans, tau = tau))
}


### Function for DP Adaptive Gradient Descent algorithm with smooth loss 
### (Lee & Kifer, 2018; Song et al, 2012)

# function for noisy max algorithm
# return index of the largest values in a list after adding random laplace noise
noisy_max = function(candidates, delta_f, epsilon) {
  candidates = candidates + rmutil::rlaplace(n = length(candidates), m = 0, s = delta_f/epsilon)
  return(which.max(candidates))
}

# function for finding the average gradient descent (Algorithm 2 in Lee & Kifer, 2018)
# used after the adaptive GD fails to find the descent direction and triggers 
# an increase in privacy budget required for the update
grad_avg = function(rho_old, rho_h, g, g_tilde, c_grad) {
  n = length(g)
  mu = rep(0, n)
  sgm = c_grad^2 / (2*(rho_h - rho_old)) * diag(n)
  g_tilde2 = g + MASS::mvrnorm(n = 1, mu = mu, Sigma = sgm)
  s_tilde = (rho_old*g_tilde + (rho_h - rho_old)*g_tilde2) / rho_h
  return(s_tilde)
}

# function for converting (epsilon, delta) privacy into rho-zCDP
# solve and check validity using equation 5 in Lee & Kifer, 2018
find_rho = function(epsilon_tot, delta_tot, tol = 1e-10) {
  rho = (-sqrt(log(1/delta_tot)) + sqrt(log(1/delta_tot) + epsilon_tot))^2
  check = abs(epsilon_tot - rho - 2*sqrt(rho*log(1/delta_tot))) <= tol
  
  while (!check) {
    epsilon_tot = epsilon_tot*1.01
    rho = -sqrt(log(1/delta_tot)) + sqrt(log(1/delta_tot) + epsilon_tot)
    check = (epsilon_tot - rho - 2*sqrt(rho*log(1/delta_tot))) >= 0
  }
  return(rho)
}

# define function f as shown in the algorithm
f = function(theta, x, y, tau, alpha) {
  ans = 0
  for (i in 1:nrow(x)) {
    ans = ans + smooth_loss(theta = theta, x = x[i, ], y = y[i], tau = tau, alpha = alpha)
  }
  return(ans)
}

# function for differentially private adaptive gradient descent 
dp_agd = function(x, y, tau, alpha, rho_nmax, rho_ng, epsilon_tot, delta_tot, gamma,
                  c_obj, c_grad, step_increase_rate, step_max_range, step_num, 
                  step_update_iter) {
  n = nrow(x)
  p = ncol(x)
  t = 1
  
  # initialize theta
  theta_t = matrix(0, p)
  
  # initialize list of possible step sizes, determined by step_max_range / step_num
  # this list will update every step_update_iter using 
  # new_step_max_range = (1 + step_increase_rate)*max(actual step size since last update)
  phi = seq(0, step_max_range, length.out = step_num)
  # initialize variable to save actual step size used
  actual_steps = NULL 
  
  # convert (epsilon, delta)-DP to rho-zCDP
  rho = find_rho(epsilon_tot = epsilon_tot, delta_tot = delta_tot)
  
  while (rho > 0) {
    print(t)
    i = 1
    g_t = 0
    for (j in 1:n) {
      tmp = loss_gradient(theta = as.matrix(theta_t[, t]), x = x[j, ], y = y[j], 
                          tau = tau, alpha = alpha)
      g_t = g_t + tmp / max(1, norm(tmp, "F") / c_grad)
    }
    g_t_tilde = g_t + MASS::mvrnorm(1, mu = rep(0, p), c_grad^2/(2*rho_ng)*diag(p))
    rho = rho - rho_ng
    g_t_tilde = g_t_tilde / norm(g_t_tilde, "F")
    
    while (i == 1) {
      omega = NULL
      for (a in phi) { # double check dim of g_t_tilde
        tmp = matrix(theta_t[, t] - a*g_t_tilde, ncol = 1)
        omega = c(omega, f(theta = tmp, x = x, y = y, tau = tau, alpha = alpha))
      }
      
      rho = rho - rho_nmax
      i = noisy_max(candidates = -omega, delta_f = c_obj, epsilon = sqrt(2*rho_nmax))
      print(i)
      if (i > 1) {
        if (rho > 0) {
          theta_t = cbind(theta_t, t(theta_t[, t] - phi[i]*g_t_tilde))
          actual_steps = c(actual_steps, phi[i])
        }
      } else {
        rho_old = rho_ng
        rho_ng = (1 + gamma)*rho_ng
        g_t_tilde = grad_avg(rho_old = rho_old, rho_h = rho_ng, g = g_t, 
                             g_tilde = g_t_tilde, c_grad = c_grad)
        rho = rho - (rho_ng - rho_old)
        if (rho <= 0) break
        print(paste("rho", rho))
      }
    }
    
    if (t %% step_update_iter == 0) {
      step_max_range = (1 + step_increase_rate)*max(actual_steps)
      phi = seq(0, step_max_range, length.out = step_num)
      actual_steps = NULL
    }
    
    t = t+1
    
  }
  
  return(list("par" = tail(t(theta_t), 1), "all" = theta_t))
}

estimate_dp_agd = function(data, total_eps, tau, alpha, max_x, delta = 1e-8) {
  ep = total_eps/length(tau)
  n = nrow(data)
  i = ncol(data)
  y = matrix(data[, i], nrow = n)
  x = as.matrix(cbind(rep(1, nrow(data)), data))
  x = as.matrix(x[, -ncol(x)])
  ans = NULL
  
  rho_nmax = find_rho(ep/120, delta)
  rho_ng = rho_nmax
  for (t in tau) {
    c_grad = max(abs(1-t), abs(t))*sqrt(ncol(x))*max_x
    c_obj = c_grad
    tmp = dp_agd(x = x, y = y, tau = t, alpha = alpha, rho_nmax = rho_nmax, 
                 rho_ng = rho_ng, epsilon_tot = ep, delta_tot = delta, gamma = 0.1, 
                 c_obj = c_obj, c_grad = c_grad, step_increase_rate = 0.1, 
                 step_max_range = 2, step_num = 20, step_update_iter = 10)$par
    ans = rbind(ans, tmp)
  }
  ans = t(ans)
  colnames(ans) = tau
  
  return(list(par = ans, tau = tau))
}


################################################################################

## Stepwise and Sandwich KNG code

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
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*Cx) - (1/2)*(beta%*%beta)
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



#what does blen do? should it be deleted?
#manipulate check_data here
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
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*Cx) - (1/2)*(beta%*%beta)
    return(ans)
  }
  
  out = constrMetrop(logA = logA, init = init, nbatch = nbatch, scale = scale,
                     check_beta = check_beta, check_data = check_data,
                     ub = ub, lb = lb, method = method, type = type, 
                     lower_accept = lower_accept, upper_accept = upper_accept,
                     update_after = update_after, adjust_scale_by = adjust_scale_by)
  return(out)
}

#instead of having lower_scale and upper_scale, scale should be input as a vector, with
#the same order as tau. If the same scale is used for all the quantiles, then it should be a
#vector of the same value
stepwiseKNG = function(data, total_eps, median_eps = NULL, tau, scale = 1e-2, Cx,
                       change_scale = NULL, change_quantile = NULL, 
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
  R = max(abs(Y))
  Y = Y/R
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  X[X > Cx] = Cx
  sumX = apply(X = X, 2, FUN = sum)
  m = ncol(X) - 1
  
  if (is.null(formula)) {
    vars = colnames(data)
    formula = paste(vars[i], " ~ .")
  }
  
  #set specific scale for the median
  nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = 0.5)
  out = KNG(init = coef(nonpriv)/R, ep = total_eps*median_eps, tau = 0.5, 
            sumX = sumX, X = X, Y = Y, Cx = Cx, nbatch = nbatch, scale = 1e-4/R, 
            ub = ub, lb = lb, upper_accept = upper_accept, lower_accept = lower_accept, 
            update_after = update_after, adjust_scale_by = adjust_scale_by)
  median_beta_kng = tail(out[[1]], 1)*R
  accept_rate = out[[2]]
  scale_output = tail(out[[3]], 1)
  ans = t(median_beta_kng)
  
  tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
  tau_upper = sort(tau[tau > 0.5])
  
  #names(scale) = tau
  
  if (length(tau_lower) > 0){
    check_beta = median_beta_kng
    for (i in 1:length(tau_lower)){
      new_tau = tau_lower[i]
      #tmp = b[, which(colnames(b) == lowertau)]
      # if (new_tau < 0.3){
      #   tmp = scale[1]
      # } else {
      #   tmp = scale[2]
      # }
      out = constrKNG(init = check_beta/R, ep = ep, tau = new_tau, sumX = sumX, X = X, 
                      Y = Y, Cx = Cx, nbatch = nbatch, scale = scale/R, 
                      check_beta = check_beta/R, check_data = check_data, ub = ub, 
                      lb = lb, method = method, type = "lower", lower_accept = lower_accept, 
                      upper_accept = upper_accept, update_after = update_after, 
                      adjust_scale_by = adjust_scale_by)
      proposed_beta = tail(out[[1]], 1)*R
      
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
      if (new_tau > change_quantile){
        tmp = change_scale
      } else {
        tmp = scale
      }
      out = constrKNG(init = check_beta/R, ep = ep, tau = new_tau, sumX = sumX, X = X, Y = Y, 
                      Cx = Cx, nbatch = nbatch, scale = tmp/R, check_beta = check_beta/R, 
                      check_data = check_data, ub = ub , lb = lb, method = method, 
                      type = "upper", lower_accept = lower_accept, upper_accept = upper_accept, 
                      update_after = update_after, adjust_scale_by = adjust_scale_by)
      proposed_beta = tail(out[[1]], 1)*R
      
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


constrMetropSandwich = function(logA, init, nbatch = 10000, scale = 1e-4, 
                                lowerbeta, upperbeta, 
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
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*Cx) - (1/2)*(beta%*%beta)
    return(ans)
  }
  
  out = constrMetropSandwich(logA = logA, init = init, nbatch = nbatch, scale = scale,
                             lowerbeta = lowerbeta, upperbeta = upperbeta, 
                             check_data = check_data, ub = ub, lb = lb, method = method,
                             lower_accept = lower_accept, upper_accept = upper_accept,
                             update_after = update_after, adjust_scale_by = adjust_scale_by)
  return(out)
}

#differnt lower and upper scale helps find scale easier (due to the long tail)
#scale needs to be a vector of the same order as tau
#quantile vector does not need to be in an order
sandwichKNG = function(data, total_eps, median_eps = NULL, main_tau_eps = NULL,
                       tau, main_tau, scale = c(1e-4, 1e-4), sw_scale = c(1e-4, 1e-4), 
                       change_scale, change_quantile, Cx,
                       sw_change_scale, sw_change_quantile,
                       nbatch = 10000, method = c("fixed", "varying"), ub = Inf, 
                       lb = -Inf, check_data = NULL, lower_accept = 0, upper_accept = 1, 
                       update_after = 10, adjust_scale_by = 2, formula = NULL){
  
  main_tau_fac = as.factor(main_tau)
  sandwich_tau = tau[!tau %in% main_tau_fac]
  main_tau_order = tau[tau %in% main_tau_fac]
  
  if(is.null(main_tau_eps)){
    main_tau_eps = ifelse(method == "fixed", 0.8, 0.9)
  }
  
  data = as.matrix(data)
  i = ncol(data)
  Y = as.matrix(data[,i])
  R = max(abs(Y))
  Y = Y/R
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
                    tau = main_tau_order, scale = scale, change_scale = change_scale, 
                    change_quantile = change_quantile, nbatch = nbatch, Cx = Cx,
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
      out = constrKNGSandwich(init = upperbeta/R, ep = eps_sandwich, tau = curr_tau, 
                              sumX = sumX, X = X, Y = Y, Cx = Cx, nbatch = nbatch, 
                              scale = sw_scale/R, lowerbeta = lowerbeta/R, upperbeta = upperbeta/R, 
                              check_data = check_data, ub = ub, lb = lb, method = method,
                              lower_accept = lower_accept, upper_accept = upper_accept, 
                              update_after = update_after, adjust_scale_by = adjust_scale_by)
      
      b_sandwich = tail(out[[1]], 1)*R
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
      
      if (curr_tau > sw_change_quantile){
        tmp = sw_change_scale
      } else {
        tmp = sw_scale
      }
      
      out = constrKNGSandwich(init = lowerbeta/R, ep = eps_sandwich, tau = curr_tau, 
                              sumX = sumX, X = X, Y = Y, nbatch = nbatch, scale = tmp/R, 
                              lowerbeta = lowerbeta/R, upperbeta = upperbeta/R, Cx = Cx,
                              check_data = check_data, ub = ub, lb = lb, method = method,
                              lower_accept = lower_accept, upper_accept = upper_accept, 
                              update_after = update_after, adjust_scale_by = adjust_scale_by)
      b_sandwich = tail(out[[1]], 1)*R
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

