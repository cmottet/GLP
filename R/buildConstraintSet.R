# Check that the parameters D, J1, and J2 are well defined
checkDJ1J2 <- function(D, J1, J2)
{
  if (!(is.finite(D) & (D >= 0)  & (floor(D) == D))){
    print("D must be a finite non-negative integer")
    return(FALSE)
  }

  if (!all(is.finite(J1) & (J1 >= 0) )){
    print("J1 must only contain finite non-negative real numbers")
    return(FALSE)
  }


  if (!all(is.finite(J2)  & (J2 > 0) & (floor(J2) == J2))) {
    print("J2 must only contain finite positive integers")
    return(FALSE)
  }

  if (max(J2) > D){
    print("D must be larger than the maximum value in the set J2")
    return(FALSE)
  }

  return(TRUE)
}

#' Constructs constraints functions of (26)
#'
#' Constructs the constraints functions \eqn{G_{j,1}} and  \eqn{G_{j,2}} defined in (26)
#' @param D Non-negative integer
#' @param J1 Vector of non-negative real values
#' @param J2 Vector of positive integers
#'
#' @return A list containing the constraints functions associated with J2 then J1
#' @export
#'
#' @examples
#' buildMomentDerivativeConstFunc(5, 1:2, 1:3)
#' buildMomentDerivativeConstFunc(1, 1:2, 1:3)
buildMomentDerivativeConstFunc <-  function(D,J1,J2)
{
  # Check D, J1, and J2
  if (!checkDJ1J2(D, J1, J2)) return(NA)

  output  <- list()
  J <- max(J2)
  J <- if (J == -Inf)  0 else J

  if (length(J2)>= 1)
    for (i in 1:length(J2))
      output[[i]] = eval(substitute (function(x) x^(J-j)/factorial(D-j),list(j=J2[i], D=D, J=J) ))

  K <- length(output)

  if (length(J1)>= 1)
    for (i in 1:length(J1))
      output[[K + i]] = eval(substitute (function(x) gamma(j+1)/gamma(j+D+1)*x^(j+J),list(j = J1[i], D=D, J=J)))


  return(output)
}


