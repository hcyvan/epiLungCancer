
#' convert a string percentage into numeric
#' 
#' convert a string percentage into numeric
#' 
#' @param x a percentage string end with %
#' @return the numeric
#' @examples
#' percent2numeric('%0.01')
#' 
#' 0.1
#' @export
percent2numeric <- function(x){
  return(as.numeric(substr(x,0,nchar(x)-1))/100)
}
