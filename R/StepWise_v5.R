#source("splineBox1.9.R")

#install.packages("DescTools")
#install.packages("rlist")
#library(DescTools)
#library(rlist)

fit_model <- function(x, z, y, knots, order, max_iter = max_iter, epsilon = epsilon)
{
  # x : n x 1 vector
  # z : n x K matrix
  # y : n x 1 vector
  n = length(y)
  X = bsplines(x, knots = knots, order = order, derivative = 0)
  dimension = ncol(X)
  beta = rep(0, dimension)
  # z의 범주의 수
  num_cate = ncol(z)
  gamma = matrix(0, num_cate, dimension)
  BB = colSums(X^2)
  bb = matrix(0, num_cate, dimension)
  for(l in 1:num_cate)
    for(j in 1: dimension)
      bb[l, j] = (sum((z[,l] * X[, j])^2))
  active_gamma_index = which(bb != 0, arr.ind = T)
  residuals_old = y
  residuals = y
  RSS_old = sum(residuals^2) / (2 * n)
  for(iter in 1 : max_iter)
  {
    # update beta
    for(j in 1 : dimension)
    {
      partial_residuals = residuals + beta[j] * X[, j]
      #b = sum(X[, j]^2)
      c = sum(partial_residuals * X[, j]) / BB[j]
      beta[j] = c
      residuals = partial_residuals - beta[j] * X[, j]
    }
    # update gamma
    for(i in 1:dim(active_gamma_index)[1])
    {
      l = active_gamma_index[i, 1]
      j = active_gamma_index[i, 2]
      partial_residuals = residuals + z[,l] * X[, j] * gamma[l,j]
      #b = bb[l, j]
      c <- sum(partial_residuals * z[,l] * X[, j]) / bb[l, j]
      gamma[l,j] = c
      residuals = partial_residuals - z[,l] * gamma[l,j] * X[, j]
    }
    RSS = sum(residuals^2) / (2 * n)
    if(abs(RSS_old - RSS) <= epsilon)
      break
    else
      RSS_old = RSS
  }
  fitted_values = y - residuals
  # p = dimension + sum(gamma != 0)
  p = 2 * dimension
  # interior knots
  inter_knots = knots[(order + 1):(length(knots) - order)]
  BIC = n * log(sum(residuals^2) / n) + p * log(n)
  AIC = n * log(sum(residuals^2) / n) + p * 2
  # Results
  results = list()
  results$beta = beta
  results$gamma = gamma
  results$residuals = residuals
  results$fitted_values = fitted_values
  results$iteration = iter
  results$knots = inter_knots
  results$n = n
  results$BIC = BIC
  results$AIC = AIC
  return(results)
}


StepWise = function(x, z, y, max_knots = 10, order = 4,
                    criterion = "BIC", max_iter = 1000, epsilon = 1e-5)
{
  # cat("Start!!!", "\n")
  if (criterion == "BIC"){
    criter = 8
  } else if(criterion == "AIC"){
    criter = 9
  } else{
    print("Please select one of BIC and AIC as the criterion")
  }

  # 외부, 내부 매듭 분리
  n = length(x)
  all_knots = knots_quantile(x, dimension = max_knots + order,
                             order = order, type = "bs")
  # Left boundary knots
  left_knots <- all_knots[1:order]
  # Right boundary knots
  right_knots <- all_knots[(order + max_knots + 1):(order + max_knots + order)]
  # interior knots
  inter_knots <- all_knots[(order + 1):(order + max_knots)]

  # Results
  results = list()
  # NULL model (number of interior knots = 0)
  knots = c(left_knots, right_knots)
  null = fit_model(x, z, y, knots =  knots, order = order, max_iter = max_iter, epsilon = epsilon)
  pre_fit = null
  # cat("NULL : ", pre_fit$BIC, "\n")

  # naming
  selection_number = 1
  step = Forward(x, z, y, inter_knots, left_knots, right_knots,
                 pre_fit, order, criter, selection_number,
                 max_iter = max_iter, epsilon = epsilon)
  # cat("1 knots : ", step$knots, ", ", step$BIC, ", ", names(step$knots), "\n")

  if(unlist(pre_fit[criter]) < unlist(step[criter]))
  {
    results = pre_fit
    return(results)
  }
  else
    pre_fit = step

  iter = 3
  while(1)
  {
    selection_number = selection_number + 1
    sel_fit = Forward(x, z, y, inter_knots, left_knots, right_knots,
                      pre_fit, order, criter, selection_number,
                      max_iter = max_iter, epsilon = epsilon)

    # cat("sel_fit : ", sel_fit$knots,", ", unlist(sel_fit[criter]), ", ", names(sel_fit$knots), "\n")

    add_knot = sel_fit$knots[!(sel_fit$knots %in% pre_fit$knots)]
    existing_knot = pre_fit$knots
    selection_number = selection_number + 1

    rm_fit = Backward(x, z, y, inter_knots, left_knots, right_knots, add_knot,
                      existing_knot, order, criter, selection_number,
                      max_iter = max_iter, epsilon = epsilon)

    # cat("rm_fit : ", rm_fit$knots, ", ", unlist(rm_fit[criter]), ", ", names(rm_fit$knots), "\n")

    if(unlist(sel_fit[criter]) < unlist(rm_fit[criter]))
      step = sel_fit
    else
      step = rm_fit

    if(unlist(pre_fit[criter]) < unlist(step[criter]))
    {
      results = pre_fit
      break
    }
    else
    {
      pre_fit = step
      # cat("step_results : ", unlist(pre_fit$knots), ", ",
      #     unlist(pre_fit[criter]), ", ", names(pre_fit$knots), "\n")
      iter = iter + 1
    }
    # cat("iteration : ", iter, "\n")
  }
  return(results)
}

Forward <- function(x, z, y, inter_knots, left_knots, right_knots,
                    pre_fit, order, criter, selection_number,
                    max_iter = max_iter, epsilon = epsilon)
{
  if(pre_fit$knots[1] >= max(x))
  {
    extra_knots = inter_knots
    num_knots = NA
  }
  else
  {
    num_knots = pre_fit$knots
    extra_knots = inter_knots[!(inter_knots %in% num_knots)]
  }
  fit = list()
  BICs = rep(0, length(extra_knots))
  names(extra_knots) = rep(selection_number, length(extra_knots))
  for(i in 1 : length(extra_knots))
  {
    fit[[i]] = list()
    in_knots = sort(c(extra_knots[i], num_knots))
    in_knots = in_knots[!is.na(in_knots)]
    knots = c(left_knots, in_knots, right_knots)
    fit[[i]] = fit_model(x, z, y, knots =  knots, order = order,
                         max_iter = max_iter, epsilon = epsilon)
    BICs[i] = fit[[i]][criter]
  }
  sel_fit = fit[[which.min(BICs)]]
  return(sel_fit)
}


Backward = function(x, z, y, inter_knots, left_knots, right_knots, add_knot,
                    existing_knot, order, criter, selection_number,
                    max_iter = max_iter, epsilon = epsilon)
{

  index = as.numeric(names(existing_knot))
  deletion_knot = existing_knot[which.min(index)]

  # 오래된거 빼고 남은 knots
  tt = sort(c(existing_knot[-which.min(index)], add_knot))

  extra_knots = inter_knots[!(inter_knots %in% c(tt, deletion_knot))]

  names(extra_knots) = rep(selection_number, length(extra_knots))

  fit = list()
  BICs = rep(0, length(extra_knots))
  for(i in 1:length(extra_knots))
  {
    fit[[i]] = list()
    in_knots = sort(c(tt, extra_knots[i]))
    knots = c(left_knots, in_knots, right_knots)
    fit[[i]] = fit_model(x, z, y, knots =  knots, order = order,
                         max_iter = max_iter, epsilon = epsilon)
    BICs[i] = fit[[i]][criter]
  }
  rm_fit = fit[[which.min(BICs)]]
  return(rm_fit)
}


