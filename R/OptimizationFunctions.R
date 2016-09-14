
evalConst = function(x,constFun,paramConstFun){
  output = matrix(0,nrow= length(constFun),ncol = length(x))

  for (i in 1:length(constFun)) output[i,] = do.call(constFun[[i]],list(x = x,paramConstFun = paramConstFun))

  return(output)
}

innerOptimizationFunction = function(x,lpdual,phase,constFun,paramConsFun, objFun = NULL,paramObjFun = NULL)
{
  newcol = evalConst(x,constFun,paramConsFun)

  if (phase == 1) output = -lpdual%*%newcol
  if (phase == 2) output =  lpdual%*%newcol -  objFun(x,paramObjFun)

  return(output)
}


#' Title
#'
#' @param initBFS
#' @param objFun
#' @param constFun
#' @param constRHS
#' @param constDir
#' @param constLambda
#' @param objLambda
#' @param paramConsFun
#' @param paramObjFun
#' @param gamma
#' @param xf
#' @param IterMax
#' @param err
#' @param rerr
#' @param factor
#' @param scale
#' @param objFuncIndic
#'
#' @return
#' @export
#' @importFrom lpSolve lp
#' @importFrom GenSA GenSA
#'
#' @examples
#' ####
#' #### Finding a the optimal probability measure (p,x) to the problem
#' ####
#' #### max P(X > c)
#' #### s.t. sum(p)  = 1
#' ####      sum(px) = 1
#' ####      sum(px^2) = 2
#' ####
#' #### where c is the 90th percentile of a standard exponential distribution
#' #### Note that the solution to this problem is known (see Theorem 3.3 of Bertsimas and Popescu)
#'
#' # Function and parameters for the integral of the objective
#' objFun <- function(x,...) return(as.numeric(paramObjFun$c <= x))
#' paramObjFun <- list(c = qexp(0.9,rate))
#'
#' # Function for the integrals of the constraints inequality
#' constFun = rep(list(
#' function(x,...) 1,
#' function(x,...) x,
#' function(x,...) x^2
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
#' constLambda <- rep(c(0,0,1),2)
#' objLambda <- 0
#'
#' # Get a basic feasible solution
#' initBFS  <-  GLPPhase1(constFun, constRHS,constDir)
#'
#' # Check feasibility
#' with(initBFS, c(sum(p), sum(p*x), sum(p*x^2)))
#'
#' # Solve the optimization program of interest
#' output <- GLPPhase2(initBFS, objFun, constFun, constRHS,
#'                      constDir, constLambda, objLambda,
#'                      paramObjFun = paramObjFun, objFuncIndic = TRUE)
#'
#' # Check that the output matches the analytical solution
#' CMsquare <- (mu2 - mu1^2)/mu1^2
#' delta <-  (paramObjFun$c/mu1-1)
#'
#' data.frame(Algorithm = output$lB, Analytical = CMsquare/(CMsquare + delta^2))
GLPPhase2 = function(initBFS,
                     objFun,
                     constFun,
                     constRHS,
                     constDir,
                     constLambda,
                     objLambda,
                     paramConsFun = NULL,
                     paramObjFun = NULL,
                     gamma = NULL,
                     xf        = 1e4,
                     IterMax   = 100,
                     err       = 1e-6,
                     rerr      = 1e-4,
                     factor    = 1,
                     scale = 196,
                     objFuncIndic = FALSE)

{

  ###
  ### Step 0 - Initialization of the optimization problem
  ###

  # Initialize output
  output = list(p         = 0,
                x         = 0,
                r         = 0,
                lastx     = NA,
                lpdual   = rep(0,length(constFun)),
                feasible  = initBFS$feasible,
                status    = 1,
                nIter     = 0,
                dual_lB   = Inf,
                primal_uB = -Inf,
                eps       = Inf)

  if (!initBFS$feasible) return(output)

  uB =  Inf
  lB   = -Inf
  r         = 0

  N = length(constRHS)

  # Initialize the x
  x = initBFS$x
  feasible = initBFS$feasible

  # Initialize the objective function of the Master Program
  objectiveIn     = rep(0,length(x)+1)
  objectiveIn[1]  = objLambda
  objectiveIn[-1] = sapply(x,objFun,paramObjFun)

  # Initialize the constraints matrix of the Master Program
  constMat     = matrix(0,length(objectiveIn),nrow = length(constRHS))
  constMat[,1] = constLambda   # which we defined as the last one in f.con
  constMat[,-1] = evalConst(x,constFun,paramConsFun)

  for (k in 1:IterMax)
  {
    ###
    ### Master Program
    ###
    outOpt    = lp(direction = "max",
                   objectiveIn,
                   constMat,
                   constDir,
                   constRHS,
                   compute.sens = TRUE,
                   scale = scale)

    if (outOpt$status !=0)
    {
      print(paste("The master problem does not converge. Error status " ,outOpt$status))

      x = x[-length(x)] ; status = outOpt$status
      if (k ==1) xnew  = eps  =NA
      break
    }

    r      = outOpt$solution[1]
    p      = outOpt$solution[-1]
    lpdual = outOpt$duals[1:N]
    lB     = outOpt$objval
    ###
    ### Sub-program
    ###

    if (objFuncIndic & paramObjFun$c <= xf)
    {
      tmp = NULL

      tmp[[1]] = GenSA(par       = 0,
                       fn        = innerOptimizationFunction,# GenSA for GLOBAL max
                       lower     = 0,
                       upper     = paramObjFun$c,
                       phase     = 2,
                       lpdual    = lpdual,
                       constFun = constFun,
                       paramConsFun = paramConsFun,
                       objFun   = objFun,
                       paramObjFun   = paramObjFun)


      tmp[[2]] = GenSA(par       = 0,
                       fn        = innerOptimizationFunction,# GenSA for GLOBAL max
                       lower     = paramObjFun$c,
                       upper     = xf,
                       phase     = 2,
                       lpdual    = lpdual,
                       constFun = constFun,
                       paramConsFun = paramConsFun,
                       objFun   = objFun,
                       paramObjFun   = paramObjFun)


      inOpt = tmp[[ which.min(c(tmp[[1]]$value,tmp[[2]]$value)) ]]
    } else
      inOpt = GenSA(par       = 0,
                    fn        = innerOptimizationFunction,# GenSA for GLOBAL max
                    lower     = 0,
                    upper     = xf,
                    phase     = 2,
                    lpdual    = lpdual,
                    constFun = constFun,
                    paramConsFun = paramConsFun,
                    objFun   = objFun,
                    paramObjFun   = paramObjFun)


    xnew  = inOpt$par
    eps   = -inOpt$value
    #     print(c(xnew ,xf,lpdual))

    if (is.null(gamma)) status = eps > err
    else {uB  = min(gamma*eps +  lB,uB)   ; status = as.integer( uB > lB + rerr*abs(lB) )}
    #     print(c(lB,uB,gamma*eps ))
    if ( !status ) break
    if (xnew %in% x) {print("Cycle") ;  status = 2 ;  break}

    # Add the column with xnew to the problem and iterate
    x            = c(x,xnew)

    objectiveIn = c(objectiveIn,objFun(xnew,paramObjFun))
    constMat    = cbind(constMat,evalConst(xnew,constFun,paramConsFun))
  }

  if (status == 0) {
    x = x[p!=0]
    p = p[p!=0]
  }

  output = list(p         = p,
                x         = x,
                r         = r,
                lpdual    = lpdual,
                feasible  = feasible,
                lB        = min(factor*lB,factor*uB),
                uB        = max(factor*lB,factor*uB),
                status    = status,
                nIter     = k,
                eps       = eps,
                lastx     = xnew)
  #   print("------------------------")
  return(output)



}


#' Title
#'
#' @param constFun
#' @param constRHS
#' @param constDir
#' @param paramConsFun
#' @param x
#' @param xf
#' @param IterMax
#' @param scale
#'
#' @return a basic feasible solution
#' @importFrom lpSolve lp
#' @importFrom GenSA GenSA
#' @export
#'
#' @examples
#'
#' ####
#' #### Finding a basic feasible solution (p,x) such that
#' ####
#' #### sum(p)  = 1
#' #### sum(px) = 1
#' #### sum(px^2) = 2
#' ####
#' ####
#'
#' # Function objective to maximize
#' objFun = function(x,...)   # H for the optim with convexity constraint
#' {
#'   c = paramObjFun$c
#'   output = I(c <= x)
#'   return(output)
#' }
#'
#' # Function for the integrals of the constraints inequality
#' constFun = rep(list(
#' function(x,...) 1,
#' function(x,...) x,
#' function(x,...) x^2
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
#' objLambda <- 0
#'
#' # Get a basic feasible solution
#' initBFS  <-  GLPPhase1(constFun, constRHS,constDir)
#'
#' # Check feasibility
#' with(initBFS, c(sum(p), sum(p*x), sum(p*x^2)))
GLPPhase1 <- function(constFun,
                      constRHS,
                      constDir,
                      paramConsFun = NULL,
                      x  = NULL,
                      xf        = 1e4,
                      IterMax   = 100,
                      scale = 196)

{
  output <- list(p         = NA,
                x         = NA,
                r         = NA,
                feasible  = FALSE,
                bound     = Inf,
                status    = 1,
                nIter     = 0,
                eps       = 0,
                maxErrRel = 0)

  ###
  ### Step 0 - Initialization of the optimization problem
  ###

  # Initialize the rhs direction  and constant vectors
  N <- length(constRHS)

  # Initialize the x
  if (is.null(x))  x <- 0

  # Initialize the objective function of the Master Program
  objectiveIn     <- rep(0,length(x)+1)
  objectiveIn[1]  <- 1

  # Initialize the constraints matrix of the Master Program
  constMat <- matrix(0,length(objectiveIn),nrow = N)
  constMat[constDir == "<=",1] <- -1
  constMat[constDir == ">=",1] <- 1

  constMat[,-1] <- evalConst(x,constFun,paramConsFun)

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
                   compute.sens = TRUE,
                   scale = scale)


    if (outOpt$status!=0)
    {
      print(paste("The master problem does not converge. Error status " ,outOpt$status))
      break;
    }

    r      <- outOpt$solution[1]
    p      <- outOpt$solution[-1]
    lpdual <- outOpt$duals[1:N]
    bound  <- outOpt$objval

    ###
    ### Sub-program
    ###
    inOpt = GenSA(par       = 0,
                  fn        = innerOptimizationFunction,# GenSA for GLOBAL max
                  lower     = 0,
                  upper     = xf,
                  lpdual    = lpdual,
                  phase     = 1,
                  constFun = constFun,
                  paramConsFun = paramConsFun)

    xnew  <- inOpt$par
    eps   <- inOpt$value
    #     print(c(eps,r))
    status <- as.integer(eps < 0 )

    if ( !status ) break
    if (xnew %in% x) {print("Cycle") ;  status = 2 ;  break}

    # Add the column with xnew to the problem and iterate
    #     print(c(xnew,eps))
    x <- c(x,xnew)

    objectiveIn <- c(objectiveIn,0)
    constMat    <- cbind(constMat,evalConst(xnew,constFun,paramConsFun))
  }

  if (k == IterMax) x <- x[-length(x)]
  if (sum(p)> 0) {x <- x[p!=0]; p <- p[p!=0]  }

  output = list(p         = p,
                x         = x,
                r         = r,
                bound     = bound,
                feasible  = (r == 0),
                status    = status,
                nIter     = k,
                eps       = eps)
  return(output)
}


