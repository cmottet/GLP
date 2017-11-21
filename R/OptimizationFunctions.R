# Function to evaluate the constraints functions at point x
evalConst <- function(x,constFun){
  output <- matrix(0,nrow= length(constFun),ncol = length(x))

  for (i in 1:length(constFun)) output[i,] <- do.call(constFun[[i]],list(x = x))

  return(output)
}

# Function to evaluate the inner optimzation function at point x
innerOptimizationFunction <- function(x,lpdual,phase,constFun, objFun = NULL)
{
  newcol = evalConst(x,constFun)

  if (phase == 1) output = -lpdual%*%newcol
  if (phase == 2) output =  lpdual%*%newcol -  objFun(x)

  return(output)
}


#' Function to solve programs such as (24)
#'
#' Solves programs such as (24)  using the generalized linear programming approach described in Algorithm 1. We slightly generalize
#' this algorithm so that inequalities in (24) can be replaced by lower inequalities or equalities if desired.
#'
#' @param initBFS List containing a boolean \emph{feasible} indicating wether an initial feasible was found for (24), and if so, \emph{initBFS} should
#' also contain two vectors \emph{p} and \emph{x} representing a feasible distribution function.
#' @param objFun Function \eqn{H}
#' @param constFun List containing the functions \eqn{G_j}
#' @param constRHS Vector containing the bounds \eqn{gamma_j}
#' @param constDir Vector containing the direction of the constraints (\emph{e.g.} '=', '<=', ">=" )
#' @param constLambda Vector containing the the values of the parameters \eqn{\lambda_{j,M}}  (use default vector of 0's to solve classic generalized linear programs)
#' @param objLambda Scalar containing the value of \eqn{\lambda_M} (use default 0 to solve classic generalized linear programs)
#' @param gamma
#' @param C Upper bound of the support of feasible distribution functions (default 1e4)
#' @param IterMax Maximum number of iterations for the procedure (default = 100)
#' @param err Tolerance of the algorithm (default 1e-6)
#' @param rerr
#' @return A list containing
#' \item{p}{Vector containing the point masses of the optimal distribution function when the program is feasible and such distribution exists}
#' \item{x}{Vector containing the point supports of the optimal distribution function when the program is feasible and such distribution exists}
#' \item{s}{Optimal value of the variable \emph{s} when the program is feasible and an optimal solution  exists}
#' \item{lpdual}{Vector containing the optimal dual multipliers when the program is feasible and such an optimal solution exists}
#' \item{lB}
#' \item{uB}
#' \item{status}{Integer describing the final status of the procedure (0 => Reached a solution, 1 => the algorithm terminated by reaching the maximum number of iterations, 2 => the algorithm entered a cycle)}
#' \item{nIter}{Number of iterations reached when the procedure terminated}
#' \item{eps}{Opposite value of the inner optimization program when the procedure terminated}
#' \item{lastx}{Value minimizing the inner optimization program when the procedure terminated}
#'
#' @export
#' @importFrom lpSolve lp
#' @importFrom GenSA GenSA
#'

#' @examples
#' ####
#' #### Finding a the optimal probability measure (p,x) to the problem
#' ####
#' #### max P(X > d)
#' #### s.t. sum(p)  = 1
#' ####      sum(px) = 1
#' ####      sum(px^2) = 2
#' ####
#' #### where c is the 90th percentile of a standard exponential distribution
#' #### Note that the solution to this problem is known (see Theorem 3.3 of Bertsimas and Popescu)
#'
#' # Function and parameters for the integral of the objective
#' d <- qexp(0.9,rate)
#' objFun <- function(x) return(as.numeric(d <= x))
#'
#' # Function for the integrals of the constraints inequality
#' constFun = rep(list(
#' function(x) 1,
#' function(x) x,
#' function(x) x^2
#' ),2)
#'
#' # Direction of the inequality constraints
#' constDir <- rep(c("<=", ">="), each = 3)
#'
#' # Values on the RHS of each inequality
#' # here we choose the moment of order 0, 1, and 2 of an exponential distribution
#' rate <- 1
#' mu0 <- 1
#' mu1 <- 1/rate
#' mu2 <- 2/rate^2
#'
#' # Lambdas for the objective function and the constraints functions
#' constLambda <- rep(c(0,0,0),2)
#' objLambda <- 0
#'
#' # Get a basic feasible solution
#' initBFS  <-  phase1(constFun, constRHS,constDir, constLambda)
#'
#' # Check feasibility
#' with(initBFS, c(sum(p), sum(p*x), sum(p*x^2)))
#'
#' # Solve the optimization program of interest
#' output <- phase2(initBFS,
#'                  objFun,
#'                  constFun,
#'                  constRHS,
#'                  constDir,
#'                  constLambda,
#'                  objLambda,
#'                  C = 1000,
#'                  err = 5e-10)
#'
#' # Check that the output matches the analytical solution
#' CMsquare <- (mu2 - mu1^2)/mu1^2
#' delta <-  (d/mu1-1)
#'
#' data.frame(Algorithm = output$lB, Analytical = CMsquare/(CMsquare + delta^2))
phase2 <- function(initBFS,
                   objFun,
                   constFun,
                   constRHS,
                   constDir,
                   constLambda = rep(0,length(constRHS)),
                   objLambda = 0,
                   gamma = NULL,
                   C        = 1e4,
                   IterMax   = 100,
                   err       = 1e-6,
                   rerr      = 1e-4)

{

  ###
  ### Step 0 - Initialization of the optimization problem
  ###

  # Initialize output
  output <- list(p         = NA,
                 x         = NA,
                 s         = NA,
                 lastx     = NA,
                 lpdual    = rep(NA,length(constFun)),
                 status    = 1,
                 nIter     = 0,
                 lB        = Inf,
                 uB        = -Inf,
                 eps       = Inf)

  # If there is no initial feasible solution, the program in unbounded
  if (!initBFS$feasible) return(output)

  uB <-  Inf
  lB <- -Inf

  # Initialize the number
  N <- length(constRHS)

  # Initialize x
  x <- initBFS$x

  feasible <- initBFS$feasible

  # Initialize the objective function of the Master Program
  objectiveIn     <- rep(0,length(x)+1)
  objectiveIn[1]  <- if (objLambda == -Inf) 0 else objLambda # Same as setting s = 0 when lambda_M = -Inf
  objectiveIn[-1] <- sapply(x,objFun)

  # Initialize the constraints matrix of the Master Program
  constMat     <- matrix(0,length(objectiveIn), nrow = length(constRHS))
  constMat[,1] <- constLambda   # which we defined as the last one in f.con
  constMat[,-1] <- evalConst(x,constFun)

  for (k in 1:IterMax)
  {
    ###
    ### Master Program
    ###
    outOpt <- lp(direction = "max",
                 objectiveIn,
                 constMat,
                 constDir,
                 constRHS,
                 compute.sens = TRUE)

    if (outOpt$status !=0)
    {
      print(paste("The master problem does not converge. Error status " , outOpt$status))

      x <- x[-length(x)]
      status <- outOpt$status

      if (k ==1) xnew  = eps  = NA

      break
    }

    s      <- outOpt$solution[1]
    p      <- outOpt$solution[-1]
    lpdual <- outOpt$duals[1:N]
    lB     <- outOpt$objval

    ###
    ### Sub-program
    ###
    inOpt <- GenSA(par       = 0,
                  fn        = innerOptimizationFunction, # GenSA for GLOBAL max
                  lower     = 0,
                  upper     = C,
                  phase     = 2,
                  lpdual    = lpdual,
                  constFun = constFun,
                  objFun   = objFun)


    xnew <- inOpt$par
    eps  <- -inOpt$value

    if (is.null(gamma)) status <- as.integer(eps > err)
    else {
      uB <- min(gamma*eps +  lB, uB)   ;
      status <- as.integer( uB > lB + rerr*abs(lB) )
    }

    if ( !status ) break

    if (xnew %in% x)
    {
      print("Cycle") ;
      status  <- 2 ;
      break
    }

    # Add the column with xnew to the problem and iterate
    x  <- c(x,xnew)

    objectiveIn <- c(objectiveIn,objFun(xnew))
    constMat    <- cbind(constMat,evalConst(xnew,constFun))
  }

  # Record the optimal distribution function
  if (status == 0) {
    x <- x[p!=0]
    p <- p[p!=0]
  }

  #  Return the results of the phase 2 algorithm
  output <- list(p         = p,
                 x         = x,
                 s         = s,
                 lpdual    = lpdual,
                 feasible  = initBFS$feasible,
                 lB        = lB,
                 uB        = uB,
                 status    = status,
                 nIter     = k,
                 eps       = eps,
                 lastx     = xnew)

    return(output)



}


#' Function finding an initial feasible solution (provided it exists) to programs such as (24)
#'
#' Finds a feasible solution (provided it exists) to programs such as (24). It uses the generalized linear programming approach  described in Algorithm 2. We slightly generalize
#' this algorithm so that inequalities in (24) can be replaced by lower inequalities or equalities if desired.
#'
#'
#' @param x Non-negative scalar representing the first point support to enter the procedure (default 0)
#' @inheritParams phase2
#'
#' @return A list containing
#' \item{p}{Vector of point masses}
#' \item{x}{Vector of point supports}
#' \item{s}{Scalar}
#' \item{r}{Optimal value of the variable \emph{r} when the program is feasible and an optimal solution  exists}
#' \item{feasible}{Boolean indicating wether \emph{(p, x, s)} is a feasible solution to (24)}
#' \item{status}{Integer describing the final status of the procedure (0 => Reached a solution, 1 => the algorithm terminated by reaching the maximum number of iterations, 2 => the algorithm entered a cycle)}
#' \item{nIter}{Number of iterations reached when the procedure terminated}
#' \item{eps}{Opposite value of the inner optimization program when the procedure terminated}

#' @importFrom lpSolve lp
#' @importFrom GenSA GenSA
#' @export
#'
#' @examples
#'
#' ####
#' #### Finding a basic feasible solution (p, x, s) such that
#' ####
#' #### sum(p)  = 1
#' #### sum(px) = 1
#' #### sum(px^2) = 2 + s
#' ####
#' ####
#'
#' # Function for the integrals of the constraints inequality
#' constFun = rep(list(
#' function(x) 1,
#' function(x) x,
#' function(x) x^2 + s
#' ),2)
#'
#' # Direction of the inequality constraints
#' constDir <- rep(c("<=", ">="), each = 3)
#'
#' # Values on the RHS of each inequality
#' mu0 <- 1
#' mu1 <- 1
#' mu2 <- 2
#'
#' constRHS <- rep(c(mu0,mu1,mu2), 2)
#'
#' # Lambdas for the objective function and the constraints functions
#' constLambda <- rep(c(0,0,1),2)
#'
#' # Get a basic feasible solution
#' initBFS  <-  phase1(constFun, constRHS,constDir, constLambda)
#'
#' # Check feasibility
#' with(initBFS, c(sum(p), sum(p*x), sum(p*x^2) + s))
phase1 <- function(constFun,
                   constRHS,
                   constDir,
                   constLambda = rep(0,length(constRHS)),
                   x  = 0,
                   C        = 1e4,
                   IterMax   = 100)

{
  output <- list(p         = NA,
                 x         = NA,
                 s         = NA,
                 r         = NA,
                 feasible  = FALSE,
                 status    = 1,
                 nIter     = 0,
                 eps       = 0)

  ###
  ### Step 0 - Initialization of the optimization problem
  ###

  # Initialize the rhs direction  and constant vectors
  N <- length(constRHS)

  # Initialize the objective function of the Master Program
  objectiveIn     <- rep(0,length(x)+2)
  objectiveIn[1]  <- 1

  # Initialize the constraints matrix of the Master Program
  constMat <- matrix(0,length(objectiveIn),nrow = N)
  constMat[constDir == "<=",1] <- -1
  constMat[constDir == ">=",1] <- 1

  constMat[,2] <- constLambda

  constMat[,-c(1,2)] <- evalConst(x,constFun)

  for (k in 1:IterMax)
  {
    ###
    ### Master Program
    ###
    outOpt <- lp(direction = "min",
                 objectiveIn,
                 constMat,
                 constDir,
                 constRHS,
                 compute.sens = TRUE)


    if (outOpt$status!=0)
    {
      print(paste("The master problem does not converge. Error status " ,outOpt$status))
      break;
    }

    r      <- outOpt$solution[1]
    s      <- outOpt$solution[2]
    p      <- outOpt$solution[-c(1,2)]
    lpdual <- outOpt$duals[1:N]

    ###
    ### Sub-program
    ###
    inOpt <- GenSA(par       = 0,
                   fn        = innerOptimizationFunction,# GenSA for GLOBAL max
                   lower     = 0,
                   upper     = C,
                   lpdual    = lpdual,
                   phase     = 1,
                   constFun = constFun)

    xnew  <- inOpt$par
    eps   <- inOpt$value

    status <- as.integer(eps < 0 )

    if ( !status ) break
    if (xnew %in% x) {print("Cycle") ;  status = 2 ;  break}

    # Add the column with xnew to the problem and iterate
    #     print(c(xnew,eps))
    x <- c(x,xnew)

    objectiveIn <- c(objectiveIn,0)
    constMat    <- cbind(constMat,evalConst(xnew,constFun))
  }

  if (k == IterMax) x <- x[-length(x)]
  if (sum(p)> 0) {x <- x[p!=0]; p <- p[p!=0]  }

  output <- list(p         = p,
                 x         = x,
                 s         = s,
                 r         = r,
                 feasible  = (r == 0),
                 status    = status,
                 nIter     = k,
                 eps       = eps)
  return(output)
}


