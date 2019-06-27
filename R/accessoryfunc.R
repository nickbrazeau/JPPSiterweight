
#' @name makeJPPSobj
#' @description OO structure for jpps object
#' @export

makeJPPS.firstorder_obj <- function(x){
  class(x) <- append(class(x), "jpps_firstorder")
  return(x)
}


#' @name makeJPPSobj
#' @description OO structure for jpps object
#' @export

makeJPPS.secondorder_obj <- function(x){
  class(x) <- append(class(x), "jpps_secondorder")
  return(x)
}

#' @name predictfirstorder overload
#' @description overload predict function
#' @export

predict.jpps_firstorder <- function(params, x){

  ret <- params["B"] *
         (
         1- 1/( 1 + (x/params["D1"])^params["C1"] +
         (x/params["D2"])^params["C2"] +
         (x/params["D3"])^params["D3"] )
         )

  ret <- data.frame(age = x, outcome = ret)

  return(ret)


}


#' @name predictsecondorder overload
#' @description overload predict function
#' @export

predict.jpps_secondorder <- function(params, x){

  ret <- params["B"] *
    (((x/params["D1"])^(params["C1"] - 1) * (params["C1"] * (1/params["D1"])) +
        (x/params["D2"])^(params["C2"] - 1) * (params["C2"] * (1/params["D2"])) +
        (x/params["D3"])^(params["C3"] - 1) * (params["C3"] * (1/params["D3"]))) /
       (1 + ((x/params["D1"])^params["C1"]) + ((x/params["D2"])^params["C2"]) +
          ((x/params["D3"])^params["C3"]))^2)


  ret <- data.frame(age = x, outcome = ret)

  return(ret)


}

