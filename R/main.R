
#' @name fit_jpps
#' @description The JPPS curve (internal use). The function takes the observed data and
#' fits the model under our specified parameters and returns the residuals.
#' Assumes an ORDINIARY nonlinear least squares model.
#' @param x numeric vector; observations (e.g regressors)
#' @param y numeric vector; predicted value of y (y-hat; regressand)
#' @param par numeric vector; model parameter estimates
#' @return res; returns the residuals of the model
#' @export


fit_jpps <- function(x, y, par) {
  B = par["B"]
  D1 = par["D1"]
  D2 = par["D2"]
  D3 = par["D3"]
  C1 = par["C1"]
  C2 = par["C2"]
  C3 = par["C3"]
  yhat<-B*(1 - 1/(1 +((x/D1)^C1)+((x/D2)^C2)+((x/D3)^C3))) #JPPS equation
  res <- y-yhat #residuals
#  ret <- sum((res)^2) # optimize by minimizing residual sums of squares

  return(res)
}

#' @name fit_linear_regression
#' @description
#' @param a numeric vector; observations (e.g regressors)
#' @param b numeric vector; predicted value of y (y-hat; regressand)
#' @param par numeric vectors; model parameter estimates
#' @export

fit_linear_regression <-function(par, a, b){
  b0 <- par[1]
  b1 <- par[2]
  bhat <- b0 + b1 * a # linear model
  res <- b-bhat #residuals
  ret <- sum(abs(res)^2) #residual sums of squares
  return(ret)
}



#' @name optimize.jpps.ord
#' @description
#' @param par named numeric vector; JPPS parameter estimates
#' @param ind.obs numeric vector; independent observations that we are trying to fit (e.g. X)
#' @param dep.obs numeric vector; dependent observations that we are trying to fit (e.g. Y)
#' @export

optimize.jpps.ord <- function(params,
                              optim.method = "SANN",
                              ind.obs,
                              dep.obs){
  mod.new <- optim(par = params,

                   fn = function(x, y, par){
                    res <- fit_jpps(x, y, par)
                    ret <- sum((res)^2) # optimize by minimizing residual sums of squares
                    return(ret)
                    },

                   x = ind.obs, y = dep.obs,

                   method = optim.method,
                   control = list(maxit = 10000, trace = F))

  res <- fit_jpps(par = mod.new$par, x = ind.obs, y = dep.obs)

  # return
  ret <- list(
    params = mod.new$par, # save "converged" params
    residuals = res # save residuals associated with "converged" params
  )
  return(ret)

}

#' @name optimize.jpps.iteratively_reweight
#' @description
#' @param par named numeric vector; JPPS parameter estimates
#' @param optim.method string; optimizer to use in the optim function for global fit
#' @param optim.method string; optimizer to use in the optim function for weight calculations
#' @param ind.obs numeric vector; independent observations that we are trying to fit (e.g. X)
#' @param dep.obs numeric vector; dependent observations that we are trying to fit (e.g. Y)
#' @param tol numeric; tolerance to break iterative reweighting
#' @export


optimize.jpps.iteratively_reweight <- function(params,
                                               optim.method = "SANN",
                                               iter.method = "Nelder",
                                               tol,
                                               maxit,
                                               ind.obs,
                                               dep.obs){
  # TODO test that / catches

  # First fit init model
  mod.new <- optim(par = params,
                   fn = function(x, y, par){
                     res <- fit_jpps(x, y, par)
                     ret <- sum((res)^2) # optimize by minimizing residual sums of squares
                     return(ret)
                   },
                   x = ind.obs,
                   y = dep.obs,
                   method = optim.method,
                   control = list(maxit = 10000, trace = F))

  # perform iterative reweighting
  for (i in 1:maxit){
    # oldcof<-as.vector(mod.new$par)
    # B<-oldcof[1]
    # D1<-oldcof[2]
    # D2<-oldcof[3]
    # D3<-oldcof[4]
    # C1<-oldcof[5]
    # C2<-oldcof[6]
    # C3<-oldcof[7]
    # yhat<-B*(1 - 1/(1 + ((ind.obs/D1)^C1) + ((ind.obs/D2)^C2) + ((ind.obs/D3)^C3))) #JPPS eq.
    # res <- dep.obs - yhat #residuals

    # store old model coefficients
    names(mod.new$par) <- c("B", "D1", "D2", "D3", "C1", "C2", "C3")
    oldcof <- mod.new$par
    # get model residuals
    res <- fit_jpps(x = ind.obs, y = dep.obs, par = oldcof)

    # fit linear model to residuals for weights
    mod.res<-optim(par=c(.5,.5),
                   fn = fit_linear_regression,
                   a = ind.obs,
                   b = res,
                   control=list(abstol = tol),
                   method="Nelder") #this models any linear trend to the residuals; essentially the same as lm() function

    gam0 <- mod.res$par[1]
    gam1 <- mod.res$par[2]
    w.new <- abs(1/(gam0 + gam1 * dep.obs)) # estimate new weights

    #  mod.new<-optim(par=as.vector(mod.new$par),weight.func,x=dat$Age,y=dat$Avg.Size,w=w.new,control=list(abstol=tol),method="SANN") #run optimization again with new weights
    #
    mod.new <- optim(par= mod.new$par,
                     fn = function(par, x, y, w){
                       res <- fit_jpps(x, y, par)
                       ret <- sum(w*(res)^2) # now account for weights in sum of squares
                       return(ret)
                     },
                     x = ind.obs,
                     y = dep.obs,
                     w = w.new,
                     control=list(abstol=tol),
                     method = optim.method) #run optimization again with new weights


    newcof <- as.vector(mod.new$par)
    dif <- sum(sqrt((oldcof - newcof)^2))/length(newcof) #compare the difference; check tolerance with last pass of iterative weighting
    dif <- sqrt((oldcof[1] - newcof[1])^2)
    if(dif<=tol){break}

  }

  ret <- list(
    params = mod.new$par, # save "converged" params
    residuals = res # save residuals associated with "converged" params
  )

  return(ret)


}



































#----------------------------------
# PARKING LOT


#'
#' #' @description The JPPS curve (internal use). The function takes the observed data and
#' #' fits the model under our specified parameters and returns the residuals.
#' #' Assumes an WEIGHTED nonlinear least squares model.
#' #' @param x numeric vector; observations (e.g regressors)
#' #' @param y numeric vector; predicted value of y (y-hat; regressand)
#' #' @param par numeric vector; model parameter estimates
#' #' @param w numeric vector; weights for observations
#'
#'
#' optimize_jpps_fit.weighted <-function(x, y, par, w){
#'   B = par[1]
#'   D1 = par[2]
#'   D2 = par[3]
#'   D3 = par[4]
#'   C1 = par[5]
#'   C2 = par[6]
#'   C3 = par[7]
#'   yhat <- B*(1 - 1/(1 +((x/D1)^C1)+((x/D2)^C2)+((x/D3)^C3)))
#'   res <- y-yhat
#'   ret <- sum(w*(res)^2) # minimizing sums of squares with weighting matrix w
#'   return(ret)
#' }
#'
#'
#'
#'
#'
#'
#'
#'
