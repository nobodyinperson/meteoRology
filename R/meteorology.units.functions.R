#### Units ####

#' Convert temperature in Celsius to Kelvin scale
#' @param T temperature in degrees Celsius
#' @details Unphysical input temperatures will be replaced by \code{NA}.
#' @return temperature converted to Kelvin scale
#' @seealso \code{\link{kel2cel}} for the opposite.
#' @export
cel2kel = function(T) {
	T[T<(-meteorology::celsius0)] <- NA
	T + meteorology::celsius0
}

#' Convert temperature in Kelvin to Celsius scale
#' @param T temperature in degrees Kelvin
#' @return temperature converted to Celsius scale
#' @details Negative input temperatures will be replaced by \code{NA}.
#' @seealso \code{\link{cel2kel}} for the opposite.
#' @export
kel2cel = function(T) {
	T[T<0] <- NA
	T - meteorology::celsius0
}

#' Convert degrees to radiant
#' @description Convert angles given in degrees to radiant.
#' @param x a numeric vector of angles in degrees
#' @return a numeric vector of respective angles in radiant.
#' @seealso \code{\link{rad2deg}} for the opposite.
#' @export
deg2rad = function(x) {
	2*pi/360*x
}
	
#' Convert radiant to degrees
#' @description Convert angles given in radiant to degrees.
#' @param x a numeric vector of angles in radiant
#' @return a numeric vector of respective angles in degrees. 
#' @seealso \code{\link{deg2rad}} for the opposite.
#' @export
rad2deg = function(x) {
	x*360/(2*pi)
}