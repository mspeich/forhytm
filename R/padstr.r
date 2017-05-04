#' Auxiliary function for processing of date/time strings
#' @export
pad1 <- function(mon){
	monc <- as.character(mon)
	if(nchar(monc)==2){
		monout <- monc
	} else {
		monout <- paste0("0",monc)
	}
	return(monout)
}

#' Auxiliary function for processing of date/time strings
#'
#' @param mon One- or two-digit number indicating day or month number
#' @return A two-digit character string with a leading zero if the input was one digit
#' @export
padstr <- Vectorize(pad1)
