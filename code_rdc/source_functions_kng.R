# Source functions for KNG, Stepwise KNG, and Sandwich KNG with convolution-type
# smooth quantile regression loss
# Note that the original logA function is now changed to kng_qr_logistic

metrop = function(logA, init, nbatch = 10000, scale = 1e-4, X = X, nchains = 3,
                  ub = Inf, lb = -Inf, lower_accept = 0, upper_accept = 1,
                  update_after = 10, adjust_scale_by = 2){
    dim = ncol(X)
    
    batch = rep(list(matrix(rep(NA, nbatch*dim),nrow = nbatch, ncol = dim)), nchains)
    u = list()
    for (i in 1:nchains) {
        u[[i]] = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
    }
    accept = rep(list(matrix(0, nrow = nbatch)), nchains)
    
    ct = rep(list(0), nchains)
    scale_ans = rep(list(scale), nchains)
    sigma = rep(list(diag(dim)*scale), nchains)
    
    for (j in 1:nchains) {
        batch[[j]][1, ] = init[[j]]
        
    }
    retry = TRUE
    while (retry) {
        retry = FALSE
        for (i in 2:nbatch){
            for (j in 1:nchains) {
                propose_val = mvtnorm::rmvnorm(1, t(batch[[j]][i-1,]), sigma[[j]])
                check = TRUE
                scale = scale_ans[[j]][i-1]
                check = check & min(X %*% t(propose_val)) >= lb
                check = check & max(X %*% t(propose_val)) <= ub
                
                if ((log(u[[j]][i]) < logA(c(propose_val)) - logA(batch[[j]][i-1, ])) & check){
                    batch[[j]][i, ] = propose_val
                    accept[[j]][i,] = 1
                    ct[[j]] = ct[[j]] + 1
                } else {
                    batch[[j]][i, ] = batch[[j]][i-1, ]
                }
                
                if(i %% update_after == 0) {
                temp_acc = mean(accept[[j]][1:i,])
                    if (temp_acc < lower_accept){
                        scale = (scale + runif(1, 0, scale))/adjust_scale_by
                    } else if (temp_acc > upper_accept) {
                        scale = scale*adjust_scale_by
                    }
                    #scale = ifelse(scale < 1e-4, 1e-4, scale)
                    sigma[[j]] = sigma[[j]]*scale
                }
                
                if (ct[[j]] == update_after){
                    sigma[[j]] = cor(na.omit(batch[[j]])) + diag(dim)*0.00002
                    if (temp_acc < lower_accept){
                        scale = scale/adjust_scale_by
                    } else if (temp_acc > upper_accept) {
                        scale = scale*adjust_scale_by
                    }
                    sigma[[j]] = sigma[[j]]*scale
                    ct[[j]] = 0
                }
                scale_ans[[j]] = c(scale_ans[[j]], scale)
    
            }
            
            if (i %% 1000 == 0) {
                tryCatch({
                    conv_check = coda::gelman.diag(coda::as.mcmc.list(lapply(batch, function(x) coda::as.mcmc(x[1:i,]))))[[2]]
                    cat("Run", i/1000, " ", c(conv_check), "\n")
                    if(conv_check <= 1.1) {
                        break
                        break
                        break
                        break
                    } else if (conv_check >= 15) {
                        return(NULL)
                        # break
                        # break
                        # break
                        # break
                    }}, error = function(cond) print('skip convergence check'))
                acc = lapply(accept, function(x) mean(x[1:i, ]))
                print(acc)
                print(lapply(batch, function(x) tail(x[1:i, ], 1)))
                for (nc in 1:nchains) {
                    cat(scale_ans[[nc]][i], '\n')
                    
                }
                
                if (any(acc < 1e-5)) {
                    idx = which(acc <1e-5)
                    for (id in idx) {
                        print(batch[[id]][1, ])
                        batch[[id]][1, ] = init[[min(which(acc >1e-5))]] + runif(dim, 0, 1)
                        print(batch[[id]][1, ])
                    }
                    retry = TRUE
                    break
                    break
                    break
                }
            }
    
    
        }
    }
    
    batch = lapply(batch, function(x) x[1:i, ])
    accept = lapply(accept, function(x) x[1:i, ])
    
    out = list(batch = batch, accept_rate = lapply(accept, mean), scale = scale_ans, sigma = sigma)
    return(out)
}


KNG = function(init, ep, tau, X, Y, Cx, nbatch = 10000, scale = 1e-4, h = 1, 
               ub = Inf, lb = -Inf, lower_accept = 0, upper_accept = 1, 
               update_after = 10, adjust_scale_by = 2, nchains = 3){
    logA = function(beta) {
        n = nrow(X)
        p = ncol(X)
        tmp = 0
        for (i in 1:n) {
            exp_val = exp(-(X[i,] %*% beta - Y[i])/h)
            tmp = tmp + (1/(1 + exp_val) - tau) %*% X[i, ]
        }
        delta = max(1-tau, tau)*Cx
        
        k_bar = 1/4
        lambda = n*k_bar*Cx^2/h*p
        del = 2*lambda#/epsilon
        gamma = n*k_bar/h # parameter for strong convexity
        
        ans = -ep / (2*delta) * norm(tmp, "M") - 1/2*norm(as.matrix(beta), "F")^2
        return(ans)   
    }
    
    out = metrop(logA = logA, init = init, nbatch = nbatch, scale = scale, X = X, nchains = nchains,
                 ub = ub, lb = lb, lower_accept = lower_accept, upper_accept = upper_accept,
                 update_after = update_after, adjust_scale_by = adjust_scale_by)
    
    return(out)
}


constrMetrop = function(logA, init, nbatch = 10000, scale = 1e-4, check_beta, check_data,
                        ub = -Inf, lb = Inf, method = c("fixed", "varying"), nchains = 3,
                        type = c("upper", "lower"), lower_accept = 0, upper_accept = 1,
                        update_after = 10, adjust_scale_by = 2){
    method = match.arg(method)
    type = match.arg(type)
    
    dim = ncol(check_data)
    
    batch = rep(list(matrix(rep(NA, nbatch*dim),nrow = nbatch, ncol = dim)), nchains)
    u = list()
    for (i in 1:nchains) {
        u[[i]] = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
    }
    accept = rep(list(matrix(0, nrow = nbatch)), nchains)
    
    ct = rep(list(0), nchains)
    scale_ans = rep(list(scale), nchains)
    sigma = rep(list(diag(rep(1, dim))*scale), nchains)
    
    for (j in 1:nchains) {
        batch[[j]][1, ] = init[[j]]
        
    }
    retry = TRUE
    while (retry) {
        retry = FALSE
        for (i in 2:nbatch){
            for (j in 1:nchains) {
                check = TRUE
                scale = scale_ans[[j]][i-1]
                if (method == "fixed") {
                    propose_val = batch[[j]][i-1, ]
                    propose_val[1] = propose_val[1] + rnorm(1, 0, sigma[[j]])
                    check = check & ifelse(type == "lower", propose_val[1] < check_beta[1], propose_val[1] > check_beta[1])
                    
                } else if (method == "varying") {
                    if (dim == 1) {
                        propose_val = batch[[j]][i-1, ] + rnorm(1, 0, sigma[[j]])
                        
                    } else {
                        propose_val = t(mvtnorm::rmvnorm(1, t(batch[[j]][i-1,]), sigma[[j]]))
                        
                    }
                    check = check & ifelse(type == "lower",
                                           all(check_data %*% propose_val < check_data %*% t(check_beta)),
                                           all(check_data %*% propose_val > check_data %*% t(check_beta)))
                }
                
                check = check & max(as.matrix(check_data) %*% as.matrix(propose_val))<= ub
                check = check & min(as.matrix(check_data) %*% as.matrix(propose_val)) >= lb
                
                if ((log(u[[j]][i]) < logA(c(propose_val)) - logA(batch[[j]][i-1, ])) & check){
                    batch[[j]][i, ] = propose_val
                    accept[[j]][i,] = 1
                    ct[[j]] = ct[[j]] + 1
                } else {
                    batch[[j]][i, ] = batch[[j]][i-1, ]
                }
                
                temp_acc = mean(accept[[j]][1:i,])
                if(temp_acc < 1e-2) {
                    if (temp_acc < lower_accept){
                        scale = (scale + runif(1, 0, scale))/adjust_scale_by
                    } else if (temp_acc > upper_accept) {
                        scale = scale*adjust_scale_by
                    }
                    #scale = ifelse(scale < 1e-4, 1e-4, scale)
                    sigma[[j]] = sigma[[j]]*scale
                }
                
                if (ct[[j]] == update_after){
                    if (method == "fixed" | dim == 1) {
                        sigma[[j]] = 1
                    } else {
                        sigma[[j]] = cor(na.omit(batch[[j]])) + diag(dim)*0.00002
                        
                    }
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
                scale_ans[[j]] = c(scale_ans[[j]], scale)
            }
            
            if (i %% 1000 == 0) {
                tryCatch({
                    conv_check = coda::gelman.diag(coda::as.mcmc.list(lapply(batch, function(x) coda::as.mcmc(x[1:i,]))))[[2]]
                    cat("Run", i/1000, " ", c(conv_check), "\n")
                    if(conv_check <= 1.1) {
                        break
                        break
                        break
                        break
                    } else if (conv_check >= 15) {
                        return(NULL)
                        # retry = TRUE
                        # break
                        # break
                        # break
                    }}, error = function(cond) print('skip convergence check'))
                acc = lapply(accept, function(x) mean(x[1:i, ]))
                print(acc)

                print(lapply(batch, function(x) tail(x[1:i, ], 1)))
                for (nc in 1:nchains) {
                    cat(scale_ans[[nc]][i], '\n')
                    
                }
                for (j in 1:nchains) {
                    if (acc[[j]] < lower_accept){
                        scale_ans[[j]][i] = (scale_ans[[j]][i])/adjust_scale_by
                    } else if (acc[[j]] > upper_accept) {
                        scale_ans[[j]][i] = scale_ans[[j]][i]*adjust_scale_by
                    }
                    #scale = ifelse(scale < 1e-4, 1e-4, scale)
                    sigma[[j]] = sigma[[j]]*scale_ans[[j]][i]
                    
                }
                for (nc in 1:nchains) {
                    cat(scale_ans[[nc]][i], '\n')
                    
                }
                
                if (any(acc < 1e-5)) {
                    idx = which(acc <1e-5)
                    for (id in idx) {
                        print(batch[[id]][1, ])
                        batch[[id]][1, ] = init[[min(which(acc >1e-5))]] + runif(dim, 0, 1)
                        print(batch[[id]][1, ])
                    }
                    retry = TRUE
                    break
                    break
                    break
                }
                
            }
        }
    }
    
    batch = lapply(batch, function(x) x[1:i, ])
    accept = lapply(accept, function(x) x[1:i, ])

    
    out = list(batch = batch, accept_rate = lapply(accept, mean), scale = scale_ans, sigma = sigma)
    
    return(out)
}



#what does blen do? should it be deleted?
#manipulate check_data here
constrKNG = function(init, ep, tau, X, Y, Cx, nbatch = 10000, scale = 1e-4, h = 1,
                     check_beta, check_data = NULL, ub = Inf, lb = -Inf,
                     method = c("fixed", "varying"), type = c("upper", "lower"),
                     lower_accept = 0, upper_accept = 1, nchains = 3,
                     update_after = 10, adjust_scale_by = 2){
    if (is.null(check_data)){
        check_data = X
    } else {
        check_data = as.matrix(check_data)
        check_data = as.matrix(cbind(rep(1, nrow(check_data)), check_data))
    }
    
    
    logA = function(beta) {
        n = nrow(X)
        p = ncol(X)
        tmp = 0
        for (i in 1:n) {
            exp_val = exp(-(X[i,] %*% beta - Y[i])/h)
            tmp = tmp + c(1/(1 + exp_val) - tau) %*% X[i, ]
        }
        delta = max(1-tau, tau)*Cx
        
        k_bar = 1/4
        lambda = n*k_bar*Cx^2/h*p
        del = 2*lambda#/epsilon
        gamma = n*k_bar/h # parameter for strong convexity
        
        ans = -ep / (2*delta) * norm(tmp, "M") - 1/2*norm(as.matrix(beta), "F")^2
        return(ans)   
    }
    
    out = constrMetrop(logA = logA, init = init, nbatch = nbatch, scale = scale,
                       check_beta = check_beta, check_data = check_data,
                       ub = ub, lb = lb, method = method, type = type, nchains = nchains,
                       lower_accept = lower_accept, upper_accept = upper_accept,
                       update_after = update_after, adjust_scale_by = adjust_scale_by)
    return(out)
}

#instead of having lower_scale and upper_scale, scale should be input as a vector, with
#the same order as tau. If the same scale is used for all the quantiles, then it should be a
#vector of the same value
stepwiseKNG = function(data, total_eps, median_eps = NULL, tau, scale = 1e-4, Cx,
                       nbatch = 10000, method = c("fixed", "varying"), h = 1,
                       ub = Inf, lb = -Inf, check_data = NULL, lower_accept = 0,
                       upper_accept = 1, update_after = 10, adjust_scale_by = 2,
                       formula = NULL){
    method = match.arg(method)
    if(is.null(median_eps)){
        median_eps = ifelse(method == "fixed", 0.3, 0.4)
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
    print('0.5')
    nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = 0.5)
    out = KNG(init = list(coef(nonpriv), rep(0, m+1), coef(nonpriv)/2), ep = total_eps*median_eps, 
              tau = 0.5, X = X, Cx = Cx, Y = Y, nbatch = nbatch, scale = scale, h = h, ub = ub, 
              lb = lb, upper_accept = 0.7, lower_accept = 0.1, update_after = 25, nchains = 3, 
              adjust_scale_by = 2)
    if (is.null(out)) {
        return(NULL)
    }
    median_beta_kng = tail(out[[1]][[1]], 1)#*R
    accept_rate = out[[2]][[1]]
    scale_output = tail(out[[3]][[1]], 1)
    scale = c(scale_output)
    ans = t(median_beta_kng)
    
    tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
    tau_upper = sort(tau[tau > 0.5])
    
    if (length(tau_lower) > 0){
        check_beta = median_beta_kng
        for (i in 1:length(tau_lower)){
            new_tau = tau_lower[i]
            print(new_tau)
            out = constrKNG(init = list(check_beta-1e-2, rep(0, m+1), check_beta/2), ep = ep, 
                            tau = new_tau, X = X, Y = Y, h = h, nbatch = nbatch, 
                            scale = scale, check_beta = check_beta, Cx = Cx, nchains = 3,
                            check_data = check_data, ub = ub, lb = lb, method = method,
                            type = "lower", lower_accept = lower_accept, upper_accept = upper_accept,
                            update_after = update_after, adjust_scale_by = adjust_scale_by)
            if (is.null(out)) {
                return(NULL)
            }
            proposed_beta = tail(out[[1]][[1]], 1)
            ans = cbind(t(proposed_beta), ans)
            scale_output = cbind(tail(out[[3]][[1]], 1), scale_output)
            accept_rate = cbind(out[[2]][[1]], accept_rate)
            
        }
    }
    
    if (length(tau_upper) > 0){
        check_beta = median_beta_kng
        for (i in 1:length(tau_upper)){
            new_tau = tau_upper[i]
            print(new_tau)
            out = constrKNG(init = list(check_beta+1e-2, rep(max(check_beta), m+1), check_beta*2), ep = ep, 
                            tau = new_tau, X = X, Y = Y, h = h, nbatch = nbatch, 
                            scale = scale, check_beta = check_beta, Cx = Cx, nchains = 3,
                            check_data = check_data, ub = ub , lb = lb, method = method,
                            type = "upper", lower_accept = lower_accept, upper_accept = upper_accept,
                            update_after = update_after, adjust_scale_by = adjust_scale_by)
            if (is.null(out)) {
                return(NULL)
            }
            proposed_beta = tail(out[[1]][[1]], 1)
            ans = cbind(ans, t(proposed_beta))
            scale_output = cbind(scale_output, tail(out[[3]][[1]], 1))
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


stepwiseKNG1 = function(data, total_eps, median_eps = NULL, tau, scale = 1e-4, Cx,
                       nbatch = 10000, method = c("fixed", "varying"), h = 1,
                       ub = Inf, lb = -Inf, check_data = NULL, lower_accept = 0,
                       upper_accept = 1, update_after = 10, adjust_scale_by = 2,
                       formula = NULL){
    method = match.arg(method)
    if(is.null(median_eps)){
        median_eps = ifelse(method == "fixed", 0.3, 0.4)
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
    print('0.5')
    nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = 0.5)
    out = KNG(init = list(coef(nonpriv)/2, coef(nonpriv)/2), ep = total_eps*median_eps, 
              tau = 0.5, X = X, Cx = Cx, Y = Y, nbatch = nbatch, scale = scale, h = h, ub = ub, 
              lb = lb, upper_accept = 0.7, lower_accept = 0.1, update_after = 25, nchains = 2, 
              adjust_scale_by = 2)
    if (is.null(out)) {
        return(NULL)
    }
    median_beta_kng = tail(out[[1]][[1]], 1)#*R
    accept_rate = out[[2]][[1]]
    scale_output = tail(out[[3]][[1]], 1)
    scale = c(scale_output)
    ans = t(median_beta_kng)
    
    tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
    tau_upper = sort(tau[tau > 0.5])
    
    if (length(tau_lower) > 0){
        check_beta = median_beta_kng
        for (i in 1:length(tau_lower)){
            new_tau = tau_lower[i]
            print(new_tau)
            out = constrKNG(init = list(check_beta-1e-2, check_beta/2), ep = ep, 
                            tau = new_tau, X = X, Y = Y, h = h, nbatch = nbatch, 
                            scale = scale, check_beta = check_beta, Cx = Cx, nchains = 2,
                            check_data = check_data, ub = ub, lb = lb, method = method,
                            type = "lower", lower_accept = lower_accept, upper_accept = upper_accept,
                            update_after = update_after, adjust_scale_by = adjust_scale_by)
            if (is.null(out)) {
                return(NULL)
            }
            proposed_beta = tail(out[[1]][[1]], 1)
            check_beta = proposed_beta
            ans = cbind(t(proposed_beta), ans)
            scale_output = cbind(tail(out[[3]][[1]], 1), scale_output)
            accept_rate = cbind(out[[2]][[1]], accept_rate)
            
        }
    }
    
    if (length(tau_upper) > 0){
        check_beta = median_beta_kng
        for (i in 1:length(tau_upper)){
            new_tau = tau_upper[i]
            print(new_tau)
            out = constrKNG(init = list(check_beta+1e-2, check_beta*1.25), ep = ep, 
                            tau = new_tau, X = X, Y = Y, h = h, nbatch = nbatch, 
                            scale = scale, check_beta = check_beta, Cx = Cx, nchains = 2,
                            check_data = check_data, ub = ub , lb = lb, method = method,
                            type = "upper", lower_accept = lower_accept, upper_accept = upper_accept,
                            update_after = update_after, adjust_scale_by = adjust_scale_by)
            if (is.null(out)) {
                return(NULL)
            }
            proposed_beta = tail(out[[1]][[1]], 1)
            check_beta = proposed_beta
            ans = cbind(ans, t(proposed_beta))
            scale_output = cbind(scale_output, tail(out[[3]][[1]], 1))
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

stepwiseKNG2 = function(data, total_eps, median_eps = NULL, tau, scale = 1e-4, Cx,
                        nbatch = 10000, method = c("fixed", "varying"), h = 1,
                        ub = Inf, lb = -Inf, check_data = NULL, lower_accept = 0,
                        upper_accept = 1, update_after = 10, adjust_scale_by = 2,
                        formula = NULL){
    method = match.arg(method)
    if(is.null(median_eps)){
        median_eps = ifelse(method == "fixed", 0.3, 0.4)
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
    print('0.5')
    nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = 0.5)
    out = KNG(init = list(coef(nonpriv)-1e-2, coef(nonpriv)/2), ep = total_eps*median_eps, 
              tau = 0.5, X = X, Cx = Cx, Y = Y, nbatch = nbatch, scale = scale, h = h, ub = ub, 
              lb = lb, upper_accept = 0.7, lower_accept = 0.1, update_after = 25, nchains = 2, 
              adjust_scale_by = 2)
    if (is.null(out)) {
        return(NULL)
    }
    median_beta_kng = tail(out[[1]][[1]], 1)#*R
    accept_rate = out[[2]][[1]]
    scale_output = tail(out[[3]][[1]], 1)
    scale = c(scale_output)
    ans = t(median_beta_kng)
    
    tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
    tau_upper = sort(tau[tau > 0.5])
    
    if (length(tau_lower) > 0){
        check_beta = median_beta_kng
        for (i in 1:length(tau_lower)){
            new_tau = tau_lower[i]
            print(new_tau)
            out = constrKNG(init = list(check_beta-1e-2, check_beta/2), ep = ep, 
                            tau = new_tau, X = X, Y = Y, h = h, nbatch = nbatch, 
                            scale = scale, check_beta = check_beta, Cx = Cx, nchains = 2,
                            check_data = check_data, ub = ub, lb = lb, method = method,
                            type = "lower", lower_accept = lower_accept, upper_accept = upper_accept,
                            update_after = update_after, adjust_scale_by = adjust_scale_by)
            if (is.null(out)) {
                return(NULL)
            }
            proposed_beta = tail(out[[1]][[1]], 1)
            ans = cbind(t(proposed_beta), ans)
            scale_output = cbind(tail(out[[3]][[1]], 1), scale_output)
            accept_rate = cbind(out[[2]][[1]], accept_rate)
            
        }
    }
    
    if (length(tau_upper) > 0){
        check_beta = median_beta_kng
        for (i in 1:length(tau_upper)){
            new_tau = tau_upper[i]
            print(new_tau)
            out = constrKNG(init = list(check_beta+1e-2, check_beta*1.25), ep = ep, 
                            tau = new_tau, X = X, Y = Y, h = h, nbatch = nbatch, 
                            scale = scale, check_beta = check_beta, Cx = Cx, nchains = 2,
                            check_data = check_data, ub = ub , lb = lb, method = method,
                            type = "upper", lower_accept = lower_accept, upper_accept = upper_accept,
                            update_after = update_after, adjust_scale_by = adjust_scale_by)
            if (is.null(out)) {
                return(NULL)
            }
            proposed_beta = tail(out[[1]][[1]], 1)
            ans = cbind(ans, t(proposed_beta))
            scale_output = cbind(scale_output, tail(out[[3]][[1]], 1))
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

constrMetropSandwich = function(logA, init, nbatch = 10000, scale = 1e-4, lowerbeta, upperbeta,
                                check_data, ub = Inf, lb = -Inf, method = c("fixed", "varying"),
                                lower_accept = 0, upper_accept = 1, update_after = 10,
                                adjust_scale_by = 2, nchains = 2){
    
    method = match.arg(method)
    dim = ncol(check_data)
    
    batch = rep(list(matrix(rep(NA, nbatch*dim),nrow = nbatch, ncol = dim)), nchains)
    u = list()
    for (i in 1:nchains) {
        u[[i]] = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
    }
    accept = rep(list(matrix(0, nrow = nbatch)), nchains)
    
    ct = rep(list(0), nchains)
    scale_ans = rep(list(scale), nchains)
    sigma = rep(list(diag(rep(1, dim))*scale), nchains)
    
    for (j in 1:nchains) {
        batch[[j]][1, ] = init[[j]]
        
    }
    retry = TRUE
    while (retry) {
        retry = FALSE
        for (i in 2:nbatch){
            for (j in 1:nchains) {
                check = TRUE
                scale = scale_ans[[j]][i-1]
                if (method == "fixed") {
                    propose_val = batch[[j]][i-1, ]
                    propose_val[1] = propose_val[1] + rnorm(1, 0, sigma[[j]])
                    check = check & (propose_val[1] < upperbeta[1]) & (propose_val[1] > lowerbeta[1])
                    
                } else if (method == "varying") {
                    if (dim == 1){
                        propose_val = batch[[j]][i-1, ] + rnorm(1, 0, sigma[[j]])
                    } else {
                        propose_val = t(mvtnorm::rmvnorm(1, t(batch[[j]][i-1,]), sigma[[j]]))
                    }
                    
                    check = check &
                        all(check_data %*% propose_val < check_data %*% upperbeta) &
                        all(check_data %*% propose_val > check_data %*% lowerbeta)
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
                
                temp_acc = mean(accept[[j]][1:i,])
                if(temp_acc < 1e-2) {
                    if (temp_acc < lower_accept){
                        scale = (scale + runif(1, 0, scale))/adjust_scale_by
                    } else if (temp_acc > upper_accept) {
                        scale = scale*adjust_scale_by
                    }
                    #scale = ifelse(scale < 1e-4, 1e-4, scale)
                    sigma[[j]] = sigma[[j]]*scale
                }
                
                if (ct[[j]] == update_after){
                    if (method == "fixed" | dim == 1) {
                        sigma[[j]] = 1
                    } else {
                        sigma[[j]] = cor(na.omit(batch[[j]])) + diag(dim)*0.00002
                        
                    }
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
                scale_ans[[j]] = c(scale_ans[[j]], scale)
            }
            
            if (i %% 1000 == 0) {
                tryCatch({
                    conv_check = coda::gelman.diag(coda::as.mcmc.list(lapply(batch, function(x) coda::as.mcmc(x[1:i,]))))[[2]]
                    cat("Run", i/1000, " ", c(conv_check), "\n")
                    if(conv_check <= 1.1) {
                        break
                        break
                        break
                        break
                    } else if (conv_check >= 15) {
                        return(NULL)
                        # retry = TRUE
                        # break
                        # break
                        # break
                    }}, error = function(cond) print('skip convergence check'))
                acc = lapply(accept, function(x) mean(x[1:i, ]))
                print(acc)
                
                print(lapply(batch, function(x) tail(x[1:i, ], 1)))
                for (nc in 1:nchains) {
                    cat(scale_ans[[nc]][i], '\n')
                    
                }
                for (j in 1:nchains) {
                    if (acc[[j]] < lower_accept){
                        scale_ans[[j]][i] = (scale_ans[[j]][i])/adjust_scale_by
                    } else if (acc[[j]] > upper_accept) {
                        scale_ans[[j]][i] = scale_ans[[j]][i]*adjust_scale_by
                    }
                    #scale = ifelse(scale < 1e-4, 1e-4, scale)
                    sigma[[j]] = sigma[[j]]*scale_ans[[j]][i]
                    
                }
                for (nc in 1:nchains) {
                    cat(scale_ans[[nc]][i], '\n')
                    
                }
                
                if (any(acc < 1e-5)) {
                    idx = which(acc <1e-5)
                    for (id in idx) {
                        print(batch[[id]][1, ])
                        batch[[id]][1, ] = init[[min(which(acc >1e-5))]] + runif(dim, 0, 1)
                        print(batch[[id]][1, ])
                    }
                    retry = TRUE
                    break
                    break
                    break
                }
                
        }
        
            }
    }
    batch = lapply(batch, function(x) x[1:i, ])
    accept = lapply(accept, function(x) x[1:i, ])
    
    out = list(batch = batch, accept_rate = lapply(accept, mean), scale = scale_ans, sigma = sigma)
    return(out)
}


constrKNGSandwich = function(init, ep, tau, X, Y, Cx, nbatch = 1000, scale = 1e-4, h=1,
                             lowerbeta, upperbeta, check_data = NULL, ub = Inf,
                             lb = -Inf, method = c("fixed", "varying"),
                             lower_accept = 0, upper_accept = 1, nchains = 2,
                             update_after = 10, adjust_scale_by = 2){
    if (is.null(check_data)){
        check_data = X
    } else {
        check_data = as.matrix(check_data)
        check_data = as.matrix(cbind(rep(1, nrow(check_data)), check_data))
    }
    
    logA = function(beta) {
        n = nrow(X)
        p = ncol(X)
        tmp = 0
        for (i in 1:n) {
            exp_val = exp(-(X[i,] %*% beta - Y[i])/h)
            tmp = tmp + (1/(1 + exp_val) - tau) %*% X[i, ]
        }
        delta = max(1-tau, tau)*Cx
        
        k_bar = 1/4
        lambda = n*k_bar*Cx^2/h*p
        del = 2*lambda#/epsilon
        gamma = n*k_bar/h # parameter for strong convexity
        
        ans = -ep / (2*delta) * norm(tmp, "M") - 1/2*norm(as.matrix(beta), "F")^2
        return(ans)   
    }
    
    out = constrMetropSandwich(logA = logA, init = init, nbatch = nbatch, scale = scale,
                               lowerbeta = lowerbeta, upperbeta = upperbeta, nchains = nchains,
                               check_data = check_data, ub = ub, lb = lb, method = method,
                               lower_accept = lower_accept, upper_accept = upper_accept,
                               update_after = update_after, adjust_scale_by = adjust_scale_by)
    return(out)
}

#differnt lower and upper scale helps find scale easier (due to the long tail)
#scale needs to be a vector of the same order as tau
#quantile vector does not need to be in an order
sandwichKNG = function(data, total_eps, median_eps = NULL, main_tau_eps = NULL,
                       tau, main_tau, Cx, scale = 1e-4, sw_scale = 1e-4, nbatch = 10000,
                       method = c("fixed", "varying"), ub = Inf, lb = -Inf, h=1,
                       check_data = NULL, lower_accept = 0, upper_accept = 1,
                       update_after = 10, adjust_scale_by = 2, formula = NULL){
    
    main_tau_fac = as.factor(main_tau)
    sandwich_tau = tau[!tau %in% main_tau_fac]
    main_tau_order = tau[tau %in% main_tau_fac]
    
    if(is.null(main_tau_eps)){
        main_tau_eps = ifelse(method == "fixed", 0.7, 0.7)
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
    
    out = stepwiseKNG1(data = data, total_eps = total_eps*main_tau_eps, median_eps = median_eps,
                      tau = main_tau_order, scale = scale, nbatch = nbatch, Cx = Cx, h=h,
                      method = method, ub = ub, lb = lb, check_data = check_data,
                      lower_accept = lower_accept, upper_accept = upper_accept,
                      update_after = update_after, adjust_scale_by = adjust_scale_by,
                      formula = formula)
    if (is.null(out)) {
        return(NULL)
    }
    b = out[[1]]
    scale_output = out[[2]]
    accept_rate = out[[3]]
    main_tau = sort(main_tau)
    names(scale_output) = main_tau
    names(accept_rate) = main_tau
    colnames(b) = main_tau
    sw_scale = mean(scale_output)/100
    
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
            update_tau = main_tau
            curr_tau = sandwich_tau_lower[i]
            print(curr_tau)
            update_tau = sort(c(update_tau, curr_tau))
            idx = which(update_tau == curr_tau)
            lowertau = update_tau[idx-1]
            uppertau = update_tau[idx+1]
            
            lowerbeta = b[, which(colnames(b) == lowertau)]
            upperbeta = b[, which(colnames(b) == uppertau)]
            out = NULL
            #main_tau_eps
            while(is.null(out)) {
                out = constrKNGSandwich(init = list((2*upperbeta + lowerbeta)/3, (upperbeta + lowerbeta)/2), 
                                        ep = eps_sandwich, tau = curr_tau, X = X, Y = Y, nbatch = nbatch, 
                                        scale = sw_scale, h = h, lowerbeta = lowerbeta, upperbeta = upperbeta, 
                                        Cx = Cx, check_data = check_data, ub = ub, lb = lb, method = method,
                                        lower_accept = lower_accept, upper_accept = upper_accept,
                                        update_after = update_after, adjust_scale_by = adjust_scale_by)
                
            }

            b_sandwich = tail(out[[1]][[1]], 1)
            sandwich_scale_lower[i] = tail(out[[3]][[1]], 1)
            accept_rate_lower[i] = out[[2]][[1]]
            b = cbind(b, t(b_sandwich))
            colnames(b)[ncol(b)] = curr_tau
            # if(m == 0){
            #     b = t(as.matrix(b[, match(update_tau, colnames(b))]))
            # } else {
            #     b = as.matrix(b[, match(update_tau, colnames(b))])
            # }
        }
    }
    
    if (length(sandwich_tau_upper) > 0){
        sandwich_scale_upper = rep(NA, length(sandwich_tau_upper))
        names(sandwich_scale_upper) = sandwich_tau_upper
        accept_rate_upper = rep(NA, length(sandwich_tau_upper))
        names(accept_rate_upper) = sandwich_tau_upper
        for (i in 1:length(sandwich_tau_upper)){
            update_tau = main_tau
            curr_tau = sandwich_tau_upper[i]
            print(curr_tau)
            update_tau = sort(c(update_tau, curr_tau))
            idx = which(update_tau == curr_tau)
            lowertau = update_tau[idx-1]
            uppertau = update_tau[idx+1]
            
            lowerbeta = b[, which(colnames(b) == lowertau)]
            upperbeta = b[, which(colnames(b) == uppertau)]
            out = NULL
            #good sandwich startscale 1/R
            while (is.null(out)) {
                out = constrKNGSandwich(init = list(upperbeta-1e-2, (upperbeta + lowerbeta)/2), 
                                        ep = eps_sandwich, tau = curr_tau, X = X, Y = Y, nbatch = nbatch, 
                                        scale = sw_scale, h = h, lowerbeta = lowerbeta, upperbeta = upperbeta, 
                                        Cx = Cx, check_data = check_data, ub = ub, lb = lb, method = method,
                                        lower_accept = lower_accept, upper_accept = upper_accept,
                                        update_after = update_after, adjust_scale_by = adjust_scale_by)
            }

            # if (is.null(out)) {
            #     return(NULL)
            # }
            b_sandwich = tail(out[[1]][[1]], 1)
            sandwich_scale_upper[i] = tail(out[[3]][[1]], 1)
            accept_rate_upper[i] = out[[2]][[1]]
            b = cbind(b, t(b_sandwich))
            colnames(b)[ncol(b)] = curr_tau
            # if(m == 0){
            #     b = t(as.matrix(b[, match(update_tau, colnames(b))]))
            # } else {
            #     b = as.matrix(b[, match(update_tau, colnames(b))])
            # }
        }
    }
    
    b = b[, sort(colnames(b))]
    
    #add case when there are no lower/upper quantiles
    scale_output = c(scale_output, sandwich_scale_upper, sandwich_scale_lower)
    scale_output = scale_output[order(names(scale_output))]
    accept_rate = c(accept_rate, accept_rate_lower, accept_rate_upper)
    accept_rate = accept_rate[order(names(accept_rate))]
    
    return(list(b, scale_output, accept_rate))
    
}
