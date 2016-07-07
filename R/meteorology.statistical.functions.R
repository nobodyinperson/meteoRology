# Statistik-Funktionen

#' Root mean square error
#' @description Calculate the root mean square error of vectors \code{x}
#'  and \code{y}.
#' @param x,y numerical vectors of same length.
#' @param na.rm Remove NAs prior to calculation? Defaults to TRUE.
#' @details
#'  \code{sqrt( mean( ( x - y ) ^ 2 ) )}
#' @examples
#' rmse(1:10,21:30)
#' # [1] 20
#' set.seed(42) # We want to produce the same results!
#' rmse(rnorm(10),rnorm(10))
#' # [1] 2.108466
#' @export
rmse = function(x,y, na.rm = T) {
	if(length(x)!=length(y)) stop("'x' and 'y' need to be of same length!")
	sqrt( mean( ( x - y ) ^ 2 , na.rm = na.rm) )
}