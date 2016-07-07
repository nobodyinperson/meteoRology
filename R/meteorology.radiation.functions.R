#### Strahlungs-Funktionen ####

# Planck-Funktion
#' Planck's law
#' 
#' @description Calculate the spectral radiance of gray body per unit area of 
#' 	emitting surface, per unit solid angle and per spectral unit (frequency or 
#'  wavelength) using Planck's law.
#' 
#' @param epsilon emissivity, in the range c(0,1)
#' @param frequency frequency in [Hz]
#' @param wavelength wavelength in [m] (ignored if \code{frequency} is set)
#' @param temp temperature in [K]
#' 
#' @return Spectral radiance of body with emissivity e at absolute temperature T
#'  and frequency \code{frequency} or rather wavelength \code{wavelength}. The unit 
#'  is [W*m-2*sr-1*Hz] if \code{frequency} was given or [W*m-2*sr-1*m-1] if 
#'  wavelength \code{wavelength} was given.
#' @details If both \code{frequency} and wavelength \code{wavelength} are given, a 
#'  warning will be printed and the frequency taken.
#' @source \url{https://en.wikipedia.org/wiki/Planck\%27s_law}
#' @seealso \code{\link{brightness.temperature}} for backwards determination of 
#'  the brightness temperature by planck's law.
#' @export
planck=function(frequency=NULL,wavelength=NULL,T,epsilon=1) {
	if(is.numeric(frequency)) {
		if(is.numeric(wavelength)) warning("Both 'frequency' and 'wavelength' specified! Will use 'frequency'!")
		return( ( epsilon * 2 * meteorology::h * frequency ^ 3 ) / meteorology::c ^ 2 * 1 / ( exp( (meteorology::h * frequency ) / ( meteorology::k * T ) ) - 1 ) )
	}
	else {
		if(!is.numeric(wavelength)) stop("One of the arguments 'frequency' and 'wavelength' has to be specified!")
		return( epsilon * 2 * meteorology::h * meteorology::c ^ 2 / wavelength ^ 5 * 1 / ( exp( ( meteorology::h * meteorology::c ) / ( wavelength * meteorology::k * T ) ) - 1 ) )
	}
}

#' Create a Planck's law function
#' @description Create a Planck's law function and possibly set parameters.
#' @param frequency logical or numeric. Return a Planck's law function for frequency?
#' @param wavelength logical or numeric. Return a Planck's law function for wavelength? 
#'  Unused if both \code{frequency} and \code{wavelength} are given.
#' @param epsilon numeric. If possible, set the constant emissivity \code{epsilon} for the function. In the range c(0,1).
#' @param T numeric. If possible, set the constant Temperature \code{T} for the function. In the range c(0,1).
#' @return a Planck's law function
#' @details The function is compiled to bytecode with \code{cmpfun} from the \code{compiler} package
#'  if possible.\cr\cr
#'  If the package \code{functionutils} is available, all given numeric arguments are
#'  substituted into the returned function.
#' @examples
#' planck.function(frequency = T)
#' # function(frequency,T,epsilon) 
#' #    epsilon * ( 2 * meteorology::h * frequency ^ 3 ) / meteorology::c ^ 2 * 1 / ( exp( (meteorology::h * frequency ) / ( meteorology::k * T ) ) - 1 )
#' 
#' planck.function(frequency = T, epsilon = 0.5, T = 300)
#' # function (frequency) 
#' #	  0.5 * (2 * meteorology::h * frequency^3)/meteorology::c^2 * 1/(exp((meteorology::h * frequency)/(meteorology::k * 300)) - 1)
#' 
#' planck.function(frequency = 3e9)
#' # function (T, epsilon) 
#' #    epsilon * (2 * meteorology::h * 3e+09^3)/meteorology::c^2 * 1/(exp((meteorology::h * 3e+09)/(meteorology::k * T)) - 1)
#' @export
planck.function = function(frequency = NULL, wavelength = NULL, epsilon = NULL, T = NULL) {
	if(!is.null(frequency))
		fun <- function(frequency,T,epsilon) epsilon * ( 2 * meteorology::h * frequency ^ 3 ) / meteorology::c ^ 2 * 1 / ( exp( (meteorology::h * frequency ) / ( meteorology::k * T ) ) - 1 )
	else if(!is.null(wavelength))
		fun <- function(wavelength,T,epsilon) epsilon * 2 * meteorology::h * meteorology::c ^ 2 / wavelength ^ 5 * 1 / ( exp( ( meteorology::h * meteorology::c ) / ( wavelength * meteorology::k * T ) ) - 1 )
	else
		stop("Either frequency or wavelength must be specified.")

	# Substitute epsilon and T is possible
	if(requireNamespace("functionutils",quietly=T)) {
		if(is.numeric(epsilon)) 
			fun <- substitute.parameters(func = fun, params = list(epsilon = epsilon))
		if(is.numeric(T)) 
			fun <- substitute.parameters(func = fun, params = list(T = T))
		if(is.numeric(frequency))
			fun <- substitute.parameters(func = fun, params = list(frequency = frequency))
		if(is.numeric(wavelength))
			fun <- substitute.parameters(func = fun, params = list(wavelength = wavelength))
			
	}
	
	
	# Kompilieren, wenn möglich
	if(requireNamespace("compiler",quietly=T))
		fun <- compiler::cmpfun(fun)
	
	fun
}

#' Locate position of maximal radiation
#' @description Using Wien's law, determine the frequency or wavelength
#'  with the maximal radiation at temperature T.
#' @param T temperature in Kelvin.
#' @param frequency logical. Return the position as frequency? Defaults to FALSE.
#' @param wavelength logical. Return the position as wavelength? Defaults to TRUE.
#' @details By default, If both \code{frequency} and \code{wavelength} are \code{TRUE},
#'  the \code{wavelength} is returned.
#' @return The position of the maximal radiation at temperature T, either in Hertz [Hz] 
#'  for \code{frequency} or in metres [m] for \code{wavelength}.
#' @source \url{https://de.wikipedia.org/wiki/Wiensches_Verschiebungsgesetz}
#' @export
position.max.radiation = function(T,frequency=F,wavelength=T) {
	if(wavelength)
		return(2897.8e-6 / T)
	else if(frequency)
		return(5.878933e10 * T)
	else
		stop("Neither 'frequency', nor 'wavelength' is TRUE.")
}

#' Total body radiation at temperature T after Stefan-Boltzmann.
#' @description The total body radiation of a body with emissivity \code{epsilon} and temperature 
#'  \code{T} after Stefan-Boltzmann.
#' @param T temperature in Kelvin [K]
#' @param epsilon emissivity, in the range of c(0,1). Defaults to black body \code{epsilon=1}.
#' @return The total body radiation energy density of a body with emissivity 
#'  \code{epsilon} and temperature \code{T} after Stefan-Boltzmann in [W*m-2].
#' @seealso \code{\link{brightness.temperature.total}} for the opposite.
#' @export
total.radiation = function(T,epsilon=1) {
	epsilon * meteorology::sigma * T ^ 4
}

#' Brightness temperature
#' @description Calculate the brightness temperature from the total emitted radiation energy density.
#' @param I total radiation energy density in [W*m-2]
#' @param epsilon emissivity, in the range of c(0,1)
#' @return The brightness temperature in Kelvin.
#' @seealso \code{\link{total.radiation}} for the opposite.
#' @export
brightness.temperature.total = function(I,epsilon=1) {
	( I / ( epsilon * meteorology::sigma ) ) ^ ( 1/4 )
}


# Noch unfertig ####
# #' Radiation transfer through multiple homogenous layers
# #' @description NOT TESTED WELL !!! Calculate the resulting radiation energy density after 
# #'  passing through multiple homogenous layers with own absorptivity / thickness /
# #'  opacity / transmission properties. Except \code{I}, all arguments are vectors
# #'  of length of the amount of layers present, each element specifying a property of
# #'  the respective layer.
# #' @param I initial total radiation energy density in [W*m-2]
# #' @param B total radiation energy density of each layer in [W*m-2]. If not specified, \code{B} is calculated
# #'  with \code{\link{total.radiation}} from argument \code{T}.
# #' @param T layer temperature. Needed to determine \code{B} if not provided.
# #' @param t transmission. If not specified, \code{t} is calculated from argument
# #'  \code{tau}. See details.
# #' @param tau opacity. If not specified, \code{tau} is calculated from arguments
# #'  \code{alpha} and \code{s}. See details.
# #' @param s layer thickness in [m]. Needed together with \code{alpha} to determine 
# #'  \code{tau} if not provided.
# #' @param alpha layer absorbtivity in the range c(0,1). Needed together with \code{s} 
# #'  to determine \code{tau} if not provided. Also used to determine \code{B} with \code{T},
# #'  because absorptivity equals emissivity.
# #' @details
# #'  Radiation transport through homogenous layers is defined as follows:
# #'  \preformatted{ I_i = t * I_(i-1) + ( 1 - t ) * B_i  }
# #'  \describe{
# #'   \item{\code{I_i}}{resulting radiation energy density after passing the \code{i}-th layer}
# #'   \item{\code{t}}{transmission, defined as \code{exp( - tau ) = exp( - alpha * s )}}
# #'   \item{\code{B_i}}{own emitted radiation energy density of the \code{i}-th layer}
# #'   }
# #'  This function tries to gather all needed information from the given data. 
# #'  The radiation transfer is calculated recursively for all layers.
# #'  \cr\cr
# #'  All arguments except \code{I} are packed together in a \code{data.frame()}, 
# #'  so they must be either of same length or length 1!
# #' @return The resulting radiation energy density after passing through layers with 
# #'  specified absorptivity properties in [W*m-2].
# #' @source \url{https://en.wikipedia.org/wiki/Radiative_transfer#Local_thermodynamic_equilibrium}, 17.06.2015 14:21 CEST
# #' @export
# radiation.transport.homogenous = function(I,B=NA,t=NA,T=NA,tau=NA,s=NA,alpha=NA) {
# 	layers = data.frame(B=B,t=t,tau=tau,s=s,alpha=alpha,T=T) # Data frame from all arguments
# # 	print(which(apply(layers,2,function(x)all(is.na(x)))))
# 	layers <- layers[-as.vector(which(apply(layers,2,function(x)all(is.na(x)))))] # Skip columns with only NAs in it (not specified argument)
# 	
# 	if(!is.null(layers$B)) { # B is defined
# 		if(is.null(layers$t)) { # t is not defined
# 			if(is.null(layers$tau)) { # tau is not defined
# 				if(is.null(layers$s) | is.null(layers$alpha)) # The pair of s and alpha is not defined
# 					stop("Not enough information provided: Need either 't' or 'tau' or 's' and 'alpha'!")
# 				else
# 					layers$tau <- layers$s * layers$alpha
# 			}
# 			layers$t <- exp( - layers$tau ) # Calculate t from tau
# 		}
# 	} else { # B is not defined
# 		if(is.null(layers$T)) { # T is not defined
# 			stop("Not enough information provided: Need 'T'!")
# 		} else { # T is defined
# 			if(is.null(layers$alpha)) { # Alpha is not defined
# 				if(is.null(layers$s) | is.null(layers$tau)) # the pair of s and tau is not defined
# 					stop("Not enough information provided: Need either 'alpha' or 's' and 'tau' for 'B' calculation!")
# 				else
# 					layers$alpha = layers$tau / layers$s # Calculate alpha from s and tau
# 			} 
# 			layers$B <- total.radiation(T = layers$T, epsilon = layers$alpha)
# 		}
# 	}
# 	
# 	rth = function(I0,B,t) t*I0+(1-t)*B # basic function
# 	Ires <- I # Rekursionsanfang
# 	
# 	for(i in 1:length(row.names(layers))) # recursion!
# 		Ires = rth(I0 = I, B = layers$B[i], t = layers$t[i])
# 	
# 	return(Ires)
# }
