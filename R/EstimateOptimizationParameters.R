# Clean and check that the parameters m and d
# are suited to be used in the function buildMomentAndDerivativesFunctions
cleanAndcheckdAndm <- function(d,m){

  if (!is.null(d)){
    if (any(round(d)!= d)) return("All derivatives order must a positive integer.")
    if (any(d > 3)) return("I am currently unable to estimate  density derivative or order higher than 2.")
    if (0 %in% d ) {
      d <- d[d!=0]
      if ( !(0 %in% m)) m <- c(m,0)
    }
    d <- sort(d, decreasing = TRUE)
  }

  if (!is.null(m)) m <- sort(m)

  output <- list(d = d, m= m)
  return(output)
}

# Function to obtain the bootstraped values of the functions contained in fboot
getBootstrapedValues <- function(sample, fboot, nboot = 1E3,mc.cores = 1)
{
  if (class(fboot) != "list") fboot <- list(fboot)

  FUN <- function(theta) bootstrap::bootstrap(x = sample, nboot, theta = theta)$thetastar
  bootstrapedValues <- parallel::mclapply(X =  fboot, FUN = FUN,mc.cores = mc.cores)
  output <- data.frame(t(plyr::ldply(bootstrapedValues)))

  return(output)
}

# Create the functions to estimate E[X^mI(x => a)] and f^(d)(a)
buildMomentAndDerivativesFunctions <- function(a, m = NULL, d = NULL, ...){
  # Initialize list of Functions
  nFunc <- length(d) + length(m)
  Func <- vector(nFunc,mode = "list")

  i <- 1

  # Build Functions for estimating the derivatives of of 1 <= order <= 3
  for (order in d){
    if (order == 1) Func[[i]] <- function(x) ks::kde(x = x, eval.points = a)$estimate
    if (order != 1) Func[[i]] <- eval(substitute(function(x) ks::kdde(x = x, eval.points = a, deriv.order = order-1)$estimate,list(order=order)))
    i <- i + 1
  }

  # Build Functions for moments
  for (order in m){
    if (order ==0) Func[[i]] <- function(x) 1 - sROC::kCDF(x = x, xgrid = a)$Fhat
    if (order !=0) Func[[i]] <- eval(substitute(function(x) mean(x^order*(x >= a), na.rm = TRUE),list(order=order)))
    i <- i + 1
  }

  output <- Func
  return(output)
}


#' Build confidence intervals for the quantities E[X^m I(X => a)] and  f^(d)(a)
#'
#' Using bootstrap, this function builds (1-alpha)-confidence intervals for the truncated moments E[X^m I(X => a)]
#' and derivatives f^(d)(a) for a given sample and threshold \emph{a}. The confidence intervals can be either hyperrectangles,
#' or ellipsoids.
#'
#' @param sample Vector containing the sample values of the random variable X
#' @param a Vector of values at which the derivatives and moments are estimated
#' @param m Vector of  real numbers
#' @param d Vector of positive integers
#' @param nboot The number of bootstrap samples desired
#' @param alpha Desired accuracy level for the confidence interval. Default value is 5\%
#' @param method A string either equal to \emph{hyperrectangle}, \emph{ellipsoid}, or \emph{both}
#'  defining the type of confidence intervals. Default value is \emph{hyperrectangle}
#' @param mc.cores Number of cores  used in the computation of the confidence intervals. Default value is 1.
#'  If \emph{mc.cores > 1}, parallel computing takes place with a total of \emph{mc.cores}
#' @param bootSample Logical value indicating wether the bootstraped sample should be returned by the function
#'
#' @return A list of with the same length as the vector \emph{a} of the argument list. Each cell of the list is also a list which contains
#' \item{hyperrectangle}{A data.frame containing the lower bounds and upper bounds of the derivatives and truncated moments. Only in
#'  the case when the argument \emph{method} = \emph{hyperrectangle} or \emph{both}}
#' \item{ellipsoid}{A list containing the vecor of means, the covariance matrix, and the radius of the ellipsoid describing the confidence interval. Only in
#' the case  when the argument \emph{method} = \emph{ellipsoid} or \emph{both}}
#' \item{a}{scalar value giving the threshold \emph{a} for the associated confidence interval}
#' \item{bootSample}{Matrix containing the bootstrapped samples used to build the confidence intervals.
#' Each row contains one bootstrap of the original sample}
#' @export
#' @importFrom dplyr '%>%'
#'
#' @examples
#' a <- c(0,1) ; m <- c(0, 1) ; d <- 1 ; nboot <- 1000
#' set.seed(100) ; sample <- rnorm(100, 0, 1)
#' hist(sample)
#' CI <- getCIMomentAndDerivatives(sample, a,m,d,nboot = nboot,mc.cores = 1, method = "both", bootSample = TRUE)
#'
#' library(ggplot2)
#' i <- 1
#' title <- paste0("95%-CI of the Estimated Parameters when a = ",a[i])
#' plot <- plotCI(CI[[i]]$bootSample, CI[[i]])
#' plot + geom_point(aes(x = dnorm(a[i]), y = 1-pnorm(a[i])), colour = "blue") +
#' geom_text(aes(x = 1.009*dnorm(a[i]), y = 1-pnorm(a[i]), label = "True Value"), colour = "blue")+
#' labs(x = "Derivative of order 1", y = "Moment of order 0") + ggtitle(title)
#'
#' i <- 2
#' title <- paste0("95%-CI of the Estimated Parameters when a = ",a[i])
#' plot <- plotCI(CI[[i]]$bootSample, CI[[i]])
#' plot + geom_point(aes(x = dnorm(a[i]), y = 1-pnorm(a[i])), colour = "blue") +
#' geom_text(aes(x = 1.009*dnorm(a[i]), y = 1-pnorm(a[i]), label = "True Value"), colour = "blue")+
#' labs(x = "Derivative of order 1", y = "Moment of order 0") + ggtitle(title)
#'
getCIMomentAndDerivatives = function(sample,a,m = NULL,d=NULL, nboot = 1E3, alpha = 0.05,
                                     method = c("hyperrectangle","ellipsoid", "both"),
                                     mc.cores = 1, bootSample = FALSE, bonferroni = FALSE)
{
  method <- match.arg(method)

  # Check conditions and clean up m and d
  tmp <- cleanAndcheckdAndm(d,m)
  d <- tmp$d
  m <- tmp$m

  BonferroniEllipse <- 1
  BonferroniRectangle  <- if (bonferroni) length(d) + length(m) else 1


  fboot <- lapply(a,buildMomentAndDerivativesFunctions,  m = m, d = d)
  bootstrapedMomentAndDerivatives <- lapply(fboot,getBootstrapedValues, sample = sample, nboot = nboot, mc.cores = mc.cores)

  output <- vector("list", length(a))
  for (i in 1:length(a)){

    # Assign name to the columns of the CI's
    namesCI <- NULL
    if (!is.null(d)) namesCI <- c(namesCI, paste0("d",d))
    if (!is.null(m)) namesCI <- c(namesCI, paste0("m",m))
    names(bootstrapedMomentAndDerivatives[[i]]) <- namesCI

    # Build the CI's based on the bootstrapped sample
    CIhyperrect <- apply(bootstrapedMomentAndDerivatives[[i]], 2, quantile, probs = c(alpha/(2*BonferroniRectangle),1-alpha/(2*BonferroniRectangle))) %>% data.frame
    CIellipsoid <-  list(mu = colMeans(bootstrapedMomentAndDerivatives[[i]], na.rm = TRUE),
                         Sigma = cov(bootstrapedMomentAndDerivatives[[i]], use = "complete.ob"),
                         radius = sqrt(qchisq(1-alpha/BonferroniEllipse, df = length(d) + length(m))))

    if (method == "hyperrectangle") output[[i]] <- list(hyperrectangle = CIhyperrect,a = a[i])
    if (method == "ellipsoid")   output[[i]] <- list(ellipsoid = CIellipsoid, a = a[i])
    if (method == "both")        output[[i]] <- list(hyperrectangle = CIhyperrect, ellipsoid = CIellipsoid, a = a[i])

    # If required, return the sample of derivatives and moments obtained by bootstrap
    if (bootSample) output[[i]] <- c(output[[i]], list(bootSample = bootstrapedMomentAndDerivatives[[i]]))
  }

  if (length(a) ==1) output <- output[[1]]
  return(output)
}

#' Plot hyperrectangular and ellipsoidal confidence intervals
#'
#' @param sample Vector of sample points
#' @param CI List containing either \emph{hyperrectangle}{, data.frame containing the lower bounds
#' and upper bounds of the derivatives and truncated moments.}, or \item{ellipsoid}{, a list containing the vector of means,
#' the covariance matrix, and the radius of the ellipsoid describing the confidence interval}
#'
#' @return an object of class ggplot2
#' @export
#' @import ggplot2
#' @importFrom gtools combinations
#'
#' @examples
#'set.seed(100)
#' n <- 300 ;  mu <- c(1,2,3) ; Sigma <- matrix(c(6,0,0,0,1,1,0,1,1),ncol =3) ; alpha = 0.05
#' sample <- MASS::mvrnorm(n,mu,Sigma) %>% as.data.frame
#' CIellipsoid <- list(mu = colMeans(sample, na.rm = TRUE), Sigma = cov(sample, use = "complete.ob"), radius = sqrt(qchisq(1-alpha/2, df = 2)))
#' CIrectangular <- apply(sample, 2, quantile, probs = c(alpha/2,1-alpha/2))
#' CI <- list(ellipsoid = CIellipsoid, hyperrectangle = CIrectangular)
#' plotCI(sample,CI)
#'
#' n <- 300 ;  mu <- c(1,2) ; Sigma <- matrix(c(6,1,1,1),ncol =2) ; alpha = 0.05
#' sample <- MASS::mvrnorm(n,mu,Sigma) %>% as.data.frame
#' CIellipsoid <- list(mu = colMeans(sample, na.rm = TRUE), Sigma = cov(sample, use = "complete.ob"), radius = sqrt(qchisq(1-alpha/2, df = 2)))
#' CIrectangular <- apply(sample, 2, quantile, probs = c(alpha/2,1-alpha/2))
#' CI <- list(ellipsoid = CIellipsoid, hyperrectangle = CIrectangular)
#' plotCI(sample,CI)

plotCI <- function(sample,CI)
{
  variable <- names(sample)
  nVar <- ncol(sample)
  allPairs <-  combinations(nVar,r = 2, repeats.allowed = F, v =  names(sample))
  dataPlot <- NULL
  for (i in 1:nrow(allPairs))
  {
    newData <- data.frame(xName = allPairs[i,1],
                          yName = allPairs[i,2],
                          x = sample[,allPairs[i,1]],
                          y = sample[,allPairs[i,2]])

    dataPlot <- rbind(dataPlot, newData)
  }
  plot <- ggplot(dataPlot) +
    geom_point(aes(x = x, y = y), colour = alpha("black", 1/3)) +
    facet_grid(yName ~ xName, scales = "free") +
    labs(x = "", y = "")

  #dataPlotHist <- dplyr::slice( dataPlot,which(dataPlot$xName == dataPlot$yName))
  #output <- plot + geom_histogram(data = dataPlotHist, aes(x = x))

  if ("hyperrectangle" %in% names(CI))
  {
    dataPlotRect <- NULL
    for (i in 1:nrow(allPairs))
    {
      CIhyperrect <- expand.grid(x = CI$hyperrectangle[,allPairs[i,1]] ,
                                 y = CI$hyperrectangle[,allPairs[i,2]])[c(1,3,4,2),]

      newData <- data.frame(xName = allPairs[i,1],
                            yName = allPairs[i,2],
                            CIhyperrect)

      dataPlotRect <- rbind(dataPlotRect, newData)
    }

    plot <- plot + geom_polygon(data = dataPlotRect, aes(x = x, y = y), colour = "red", fill = NA)
  }

  if ("ellipsoid" %in% names(CI))
  {
    dataPlotEllipse <- NULL
    segments <- 100

    for (i in 1:nrow(allPairs))
    {
      pairNames <- allPairs[i,]
      if (pairNames[1]!= pairNames[2]){
        Sigma <- CI$ellipsoid$Sigma[pairNames,pairNames ]
        chol_decomp <- chol(Sigma)
        radius <- CI$ellipsoid$radius
        center <- CI$ellipsoid$mu[pairNames]

        angles <- (0:segments) * 2 * pi/segments
        unit.circle <- cbind(cos(angles), sin(angles))
        CIellipse <- data.frame(t(center + radius * t(unit.circle %*% chol_decomp)))
        names(CIellipse) <- c("x", "y")

        newData <- data.frame(xName = allPairs[i,1],
                              yName = allPairs[i,2],
                              CIellipse)

        dataPlotEllipse <- rbind(dataPlotEllipse, newData)
      }
    }

    plot <- plot + geom_path(data = dataPlotEllipse, aes(x = x, y = y), colour = "green")
  }
  output <- plot
  return(output)
}
