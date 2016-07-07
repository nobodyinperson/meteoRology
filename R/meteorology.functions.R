#### Funktionen ####



#' Magnus formula for saturation water vapor pressure
#' 
#' @description Calculate the saturation water vapor pressure e_s(T) either over
#'  water or ice using the Magnus approximation.
#' @param T temperature, either in Kelvin.
#' @param over Calculate the saturation vapor pressure over "water" or over
#' "ice"?
#' @return The saturation water vapour pressure e_s(T) in Pascal [Pa].
#' @source \url{https://de.wikipedia.org/wiki/Sättigungsdampfdruck}
#' @export
e.magnus = function(T,over="water") {
	A = 6.112e2 # Pa
	if(over=="water") B = 17.62 # Einheitslos
	else              B = 22.46 # Einheitenlos
	if(over=="water") C = 243.12 # °C
	else              C = 272.62 # °C
	T <- meteorology::kel2cel(T) # In °Celsius umrechnen, weil so die Magnusformel definiert ist.
	A * exp ( ( B * T ) / ( C + T ) )
}

#' Meteorological wind direction
#' @param wind numerical vector of wind components c(u,v)
#' @param u lateral wind component
#' @param v longitudinal wind component
#' @return The meteorological wind direction of wind vector c(u,v), that is the angle from y-axis in mathematical negative direction to the inverted vector c(-u,-v).
#' @export
wind.direction = function(wind = c(u,v),u=wind[1],v=wind[2]) {
	# atan2 gives angle to x-axis
	# 90 - atan2 gives angle of vector to y-axis
	# 90 - atan2 + 180 "turns the vector around" and gives the angle of the inverted vector c(-u,-v) to the y-axis
	# the + 360 and modulo 360 is just to ensure the result to be positive
	( 90 - rad2deg( atan2( x = u, y = v ) ) + 180 + 360 ) %% 360
}
