#' Title
#'
#' @param d
#' @param m
#'
#' @return
#' @export
#'
#' @examples
#' d <- c(4,2,1)
#' m <- 0:3
#' cleanAndcheckdAndm(d,m)
#'
#' d <- c(1/2,2,1)
#' cleanAndcheckdAndm(d,m)
#'
#' d <- 0:3
#' m <- 1:3
#' cleanAndcheckdAndm(d,m)
#'
#'
cleanAndcheckdAndm <- function(d,m){

  if (any(round(d)!= d)) return("All derivatives order must a positive integer.")
  if (any(d > 3)) return("I am currently unable to estimate  density derivative or order higher than 2.")
  if (0 %in% d ) {
    d <- d[d!=0]
    if ( !(0 %in% m)) m <- c(m,0)
  }
  d <- sort(d, decreasing = TRUE)
  m <- sort(m)

  output <- list(d = d, m= m)
  return(output)
}

#' Title
#'
#' @param sample
#' @param fboot
#' @param nboot
#' @param ...
#'
#' @return
#' @export
#' @importFrom bootstrap bootstrap
#'
#' @examples
#' sample <- rnorm(1e4,0,1)
#' fboot <-  list(function(x) mean(x), function(x) sd(x))
#' bootstrapedValues <- getBootstrapedValues(sample,fboot, nboot = 1E4, mc.cores = 4)
#' head(bootstrapedValues)
#'
getBootstrapedValues <- function(sample, fboot, nboot = 1E3,mc.cores = 1,...)
{
  if (class(fboot) != "list") fboot <- as.list(fboot)

  FUN <- function(theta) bootstrap(x = sample, nboot, theta = theta,...)$thetastar
  bootstrapedValues <- parallel::mclapply(X =  fboot, FUN = FUN,mc.cores = mc.cores, ...)
  output <- data.frame(t(plyr::ldply(bootstrapedValues)))

  return(output)
}

#' Title
#'
#' @param a
#' @param m
#' @param d
#'
#' @return
#' @export
#' @importFrom ks kdde kcde
#'
#' @examples
#' a <- 0
#' m <- 0:2
#' d <- 1:3
#'
#' # Create the functions to estimate the derivatives and moments
#' MyFunc <- GLP:::buildMomentAndDerivativesFunctions(a,m,d)
#'
#' # Create a sample
#' sample <- rnorm(1E4, 0, 1 )
#' Truth <- c((a^2 - 1)*dnorm(a), -a*dnorm(a), dnorm(a), 1- pnorm(a), (2*pi)^(-1/2), 1/2)
#' data.frame(Truth = Truth , Estimates = sapply(MyFunc, function(f)do.call(f, list(sample = sample))))
buildMomentAndDerivativesFunctions <- function(a, m, d, ...){
  # Initialize list of Functions
  nFunc <- length(d) + length(m)
  Func <- vector(nFunc,mode = "list")

  i <- 1

  # Build Functions for estimating the derivatives of of 1 <= order <= 3
  for (order in d){
    Func[[i]] <- eval(substitute(function(sample) kdde(x = sample, eval.points = a, deriv.order = order-1, ...)$estimate,list(order=order)))
    i <- i + 1
  }

  # Build Functions for moments
  for (order in m){
    if (order ==0) Func[[i]] <- function(sample) 1-kcde(x = sample, eval.points = a)$estimate
    if (order !=0) Func[[i]] <- eval(substitute(function(sample) mean(sample^order*(sample >= a), na.rm = TRUE),list(order=order)))
    i <- i + 1
  }

  output <- Func
  return(output)
}

#' Title
#'
#' @param sample
#' @param a
#' @param m
#' @param d
#' @param nboot
#' @param alpha
#' @param method
#' @param mc.cores
#' @param ...
#'
#' @return
#' @export
#' @importFrom dplyr '%>%'
#'
#' @examples
#' a <- c(0,1) ; m <- c(0, 1) ; d <- 1 ; nboot <- 1000
#' set.seed(100) ; sample <- rnorm(100, 0, 1)
#' hist(sample)
#' CI <- getCIMomentAndDerivatives(sample, a,m,d,nboot = nboot,mc.cores = 4, method = "both", bootSample = TRUE)
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
                                     mc.cores = 1, bootSample = FALSE,...)
{

  method <- match.arg(method)

  # Check conditions and clean up m and d
  tmp <- cleanAndcheckdAndm(d,m)
  d <- tmp$d
  m <- tmp$m

  BonferroniEllipse <- length(a)
  BonferroniRectangle <- length(a)*(length(d) + length(m))
  fboot <- lapply(a,buildMomentAndDerivativesFunctions,  m = m, d = d)
  bootstrapedMomentAndDerivatives <- lapply(fboot,getBootstrapedValues, sample= sample, nboot = nboot, mc.cores = mc.cores, ...)

  output <- vector("list", length(a))
  for (i in 1:length(a)){
    names(bootstrapedMomentAndDerivatives[[i]]) <- c(paste0("d",d), paste0("m",m))
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

#' Title
#'
#' @param sample
#' @param CI
#'
#' @return
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



