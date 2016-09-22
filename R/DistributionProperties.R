####
#### Exponential distribution
####

# 1st order Derivative w.r.t. x
#' Title
#'
#' @param x
#' @param rate
#'
#' @return
#' @export
#'
#' @examples
ddexp   = function(x,rate)
{
  output = rep(0,length(x) )
  supp = x >=0
  xs   = x[supp]

  output[supp] = -rate*dexp(x,rate)
  return(output)
}

# 2nd order Derivative w.r.t. x
#' Title
#'
#' @param x
#' @param rate
#'
#' @return
#' @export
#'
#' @examples
d2exp   = function(x,rate)
{
  output = rep(0,length(x) )
  supp = x >=0
  xs   = x[supp]

  output[supp] = rate^2*dexp(x,rate)
  return(output)
}

# Convex point of the pdf
convpt.dlnorm = function(rate) 0

# UpperTuncated Moments
#' Title
#'
#' @param a
#' @param order
#' @param rate
#'
#' @return
#' @export
#'
#' @examples
UpperTruncMomExp = function(a,order,rate) #E[X^iI(x>a)]
{
  if (round(order)!= order)
  {
    print("The order must be an integer")
    break
  }

  result = 0

  for (i in 0:order)
    result = a^i*exp(-rate*a) + i/rate*result

  output = result

  return(output)
}
####
#### Log-Normal distribution
####

# 1st order Derivative w.r.t. x
ddlnorm   = function(x,meanlog,sdlog)
{
  output = rep(0,length(x) )
  supp = x >0
  xs   = x[supp]

  output[supp] = -1/x*dlnorm(xs,meanlog,sdlog)*(1 + (log(xs) - meanlog)/sdlog^2)
  return(output)
}

# 2nd order Derivative w.r.t. x
d2lnorm   = function(x,meanlog,sdlog)
{

  output = rep(0,length(x) )
  supp = x >0
  xs   = x[supp]

  y = (log(xs) - meanlog)/sdlog + sdlog
  output[supp] = -1/(xs^2*sdlog)*dlnorm(xs,meanlog,sdlog)*( (y^2 - 1)/sdlog + y)


  return(output)
}

# Convex point of the pdf
convpt.dlnorm = function(meanlog,sdlog) exp(meanlog + sdlog*(sqrt(sdlog^2 + 4) - 3*sdlog)/2 )

UpperTruncMomLogN = function(a,order,meanlog,sdlog) #E[X^iI(x>a)]
{
  if (round(order)!= order)
  {
    print("The order must be an integer")
    break
  }


  result = exp(meanlog*order +  1/2*(sdlog*order)^2)*(1-pnorm( (log(a) - meanlog)/sdlog,sdlog*order,1))
  output = result

  return(output)
}
####
#### GPD distribution
####

# 1st order Derivative w.r.t. x
ddGPD = function(x,xi,beta,u=0)
{

  output = rep(0,length(x) )
  supp = beta+xi*(x-u) >0
  xs   = x[supp]

  output[supp] =   -(xi+1)/(beta+xi*(xs-u))*dGPD(xs-u, xi,beta)

  return(output)
}

# UpperTuncated Moments
UpperTruncMomGPD = function(a,order,xi,beta,u=0) #E[X^iI(x>a)]
{
  if (round(order)!= order)
  {
    print("The order must be an integer")
    break
  }

  xf = u+ qGPD(1,xi,beta)

  if (order >= 1/xi) {
    result = Inf
  }else{


    result = NA

    if (order == 0 ) result = 1-pGPD(max(a,u)-u,xi,beta)
    if (order == 1)
    {
      term1 = 1-pGPD(a-u,xi,beta)
      term2 = max(a,u) - beta/(xi -1)*(1+ xi*(max(a,u) - u)/beta)
      result = term1*term2*(a<=xf)
    }
  }

  output = result

  return(output)
}


# CDF of log X when X ~ GPD
plogGPD = function(q,xi,beta,u=0)  pGPD(exp(q) -u,xi,beta)
dlogGPD = function(x,xi,beta,u=0)  dGPD(exp(x) -u,xi,beta)*exp(x)
ddlogGPD = function(x,xi,beta,u=0) exp(x)*(exp(x)*ddGPD(exp(x)-u,xi,beta) + dGPD(exp(x) -u,xi,beta) )

qlogGPD = function(p,xi,beta,u=0) log(qGPD(p,xi,beta) + u)


####
#### Pareto distribution
####


#pdf
dpareto  = function(x,scale=1,shape=1)
{
  pdf = rep(NA,length(x))

  pdf[x == 0] = 0
  pdf[x != 0] = shape*scale^shape/x^(shape+1)*(scale <= x)

  output = pdf
  return(output)
}

# CDF
ppareto  = function(x,scale=1,shape=1)
{
  p = rep(NA,length(x))

  p[x==0] = 0
  p[x!=0] = 1 - (scale/x)^shape*(scale <= x)

  output = p

  return(output)
}

# Quantile
qpareto = function(p,scale=1,shape=1)
{
  q = rep(NA,length(p))

  q[p==1] = Inf
  q[p==0] = scale
  q[0<p & p<1] = scale*(1-p[0<p & p<1])^(-1/shape)

  output = q

  return(output)
}
# 1st order Derivative
ddpareto = function(x,scale=1,shape=1)
{
  dpdf = rep(NA,length(x))

  dpdf[x == 0] = 0
  dpdf[x != 0] = -(shape+1)/x*dpareto(x,scale,shape)

  output = dpdf

  return(output)
}

# UpperTuncated Moments
UpperTruncMomPareto = function(a,order,scale=1,shape = 1) #E[X^iI(x>a)]
{
  uppMoments = rep(NA,length(order))

  # Case when the order is larger or equal than the shape
  uppMoments[(order >= shape)] = Inf

  # Case when the order is strictly smaller than the shape
  ind = which(order < shape)
  uppMoments[ind]  = shape*scale^shape/(shape-order[ind])*max(a,scale)^(order[ind]-shape)

  output = uppMoments
  return(output)
}


# # 1st order Derivative w.r.t. x
# ddpareto = function(x,shape,scale)
# {
#   output = -(shape+1)/(x+scale)*dpareto(x,shape,scale)
#   #dpareto is a function of the package actuar
#   return(output)
# }

####
#### Gamma distribution
####

# 1st order Derivative w.r.t. x
ddgamma  = function(x,shape,rate)
{

  output = rep(0,length(x) )
  supp = x>0
  xs   = x[supp]

  output[supp] =   ((shape-1)/xs - rate)*dgamma(xs,shape,rate)

  return(output)
}

# 2nd order Derivative w.r.t. x
d2dgamma = function(x,shape,rate)
{
  output = rep(0,length(x) )
  supp = x>0
  xs   = x[supp]

  output[supp] =  1/xs^2*( (shape-1 -rate*xs)^2 - shape + 1 )*dgamma(xs,shape,rate)
}

# Convex point of the pdf
convpt.dgamma =  function(shape,rate) ((shape - 1) + sqrt(shape-1))/rate

####
#### Log-Gamma distribution <=> log X ~ Gamma
####

# 1st order Derivative w.r.t. x
ddlgamma  = function(x,shapelog,ratelog)
{
  output = rep(0,length(x) )
  supp = x>1
  xs   = x[supp]

  output[supp] =  1/xs*(-1 -ratelog + (shapelog-1)/log(xs) )*dlgamma(xs,shapelog,ratelog)
  return(output)
}

# 2nd order Derivative w.r.t. x
d2dlgamma = function(x,shapelog,ratelog)
{
  r = ratelog  + 1
  s = shapelog - 1

  output = rep(0,length(x) )
  supp = x>1
  xs   = x[supp]

  output[supp] = 1/(xs*log(xs))^2*( log(xs)^2*r*(r+1) - log(xs)*s*(1 + 2*r) + s*(s-1))*dlgamma(xs,shapelog,ratelog)
  return(output)
}

# Convex point of the pdf
convpt.dlgamma = function(shapelog,ratelog)
{
  r = ratelog  + 1
  s = shapelog - 1

  c0 =  s*(s-1)
  c1 = -s*(1 + 2*r)
  c2 =  r*(r+1)

  sol_quad_eq = Re(polyroot(c(c0,c1,c2)))
  output = exp(max(sol_quad_eq))

  return(output)
}


####
#### Beta distribution
####
# 1st order Derivative w.r.t. x
ddbeta  = function(x,shape1,shape2)
{
  output = rep(0,length(x) )

  supp = 0<x & x <1
  xs   = x[supp]

  output[supp] = dbeta(xs,shape1,shape2)*((shape1 - 1)/xs - (shape2 - 1)/(1-xs))

  return(output)
}

d2dbeta  = function(x,shape1,shape2)
{
  output = rep(0,length(x) )

  supp = 0<x & x <1
  xs   = x[supp]

  output[supp] = dbeta(xs,shape1,shape2)*(((shape1 - 1)/xs + (shape2 - 1)/(1-xs))^2 -(shape2-1)/(1-xs)^2
                                          - (shape1-1)/xs^2)
  return(output)
}

####
#### "Hill Horror Plot" distribution
#### F(x) = 1 - 1/(x log x)
####

f_inv_horror_dist = function(d)
{
  # Define the interval in which to search for the root
  lower = sqrt(d)
  upper = max(exp(1),d)

  output = uniroot( function(x,d)x*log(x) -d, d = d, interval = c(lower,upper))$root
  return(output)
}



# Using the inversion CDF method, we simulate
# a sample
rhorror = function(n)
{
  U = runif(n)
  d = 1/(1-U)

  output = sapply(d,f_inv_horror_dist )
  return(output)
}

# Check that the simulation is accurate
# X = sort(simulate_horror_dist(1e5))
# Fn = ecdf(X)
# plot(Fn,xlim = c(f_inv(1),max(X)),log = "x" )
# lines(X,1- 1/(X*log(X)),col = "red"   )

qhorror   = function(p)
{
  q = rep(NA,length(p))
  q[p==1] = Inf

  ind = 0 <= p & p <1
  q[ind] =  sapply(1/(1-p[ind]),f_inv_horror_dist)

  output = q
  return(output)
}

phorror   = function(x)
{
  xmin = qhorror(0)
  ind_supp = which(x>= xmin)

  output = rep(0,length(x))
  z = x[ind_supp]
  output[ind_supp] = 1- 1/(z*log(z))

  return(output)
}

dhorror   = function(x){

  xmin = qhorror(0)
  ind_supp = which(x>= xmin)

  output = rep(0,length(x))
  z = x[ind_supp]
  output[ind_supp] = (log(z)+1)/(z*log(z))^2

  return(output)
}

ddhorror  = function(x) {
  xmin = qhorror(0)
  ind_supp = which(x>= xmin)
  output = rep(0,length(x))

  z = x[ind_supp]
  output[ind_supp] = (log(z) - 2*(log(z)+1)^2)/(z*log(z))^3
  return(output)

}

d2dhorror = function(x){
  xmin = qhorror(0)
  ind_supp = which(x>= xmin)
  output = rep(0,length(x))

  z = x[ind_supp]
  output[ind_supp] =(6*(log(z)+1)^3 -6*log(z) -7*log(z)^2 )/(z*log(z))^4
  return(output)

}

UpperTruncHorrorMoment = function(a,order,seed, eps = 1e-4)
{
  xmin = qhorror(0)
  A = pmax(a,xmin)

  if (order >= 1) output = Inf
  if (order < 1)
  {
    if (order == 0) output =  1 - phorror(A)
    else {
      rate = (1-order)*log(A)
      if (rate>=1)  output = as.numeric(A^(order - 1)/log(A) + order*gammainc(rate,0)[2])
      if (rate< 1){

        if (!missing(seed)) set.seed = seed
        X = rexp(1e7,rate)
        sig = sd(1/X*(X>=1))
        nsam = (2*qnorm(0.975)*sig/(eps*rate))^2

        if (!missing(seed)) set.seed = seed
        X = rexp(nsam,rate)
        output =  A^(order - 1)/log(A) + order*mean(1/X*(X>=1))/rate
      }
    }
  }
  return(output)
}

###
### Horror Hill distribution with log-transform
###

# Check that the simulation is accurate
# X = sort(simulate_horror_dist(1e5))
# Fn = ecdf(X)
# plot(Fn,xlim = c(f_inv(1),max(X)),log = "x" )
# lines(X,1- 1/(X*log(X)),col = "red"   )

qlhorror   = function(p)
{
  output = log(qhorror(p))
  return(output)
}

plhorror   = function(x)
{
  output = phorror(exp(x))
  return(output)
}

dlhorror   = function(x){
  output = exp(x)*dhorror(exp(x))
  return(output)
}

ddlhorror  = function(x) {
  output =  exp(2*x)*ddhorror(exp(x))+ exp(x)*dhorror(exp(x))
  return(output)

}

d2dlhorror = function(x){
  output = exp(3*x)*d2dhorror(exp(x)) + 3*exp(2*x)*ddhorror(exp(x)) + exp(x)*dhorror(exp(x))
  return(output)

}

UpperTrunclHorrorMoment = function(a,order,seed, eps = 1e-4)
{
  xmin = qlhorror(0)
  A = pmax(a,xmin)

  if (order <0){ print("The order must be non-negative") ; return(NA)} else
  {
    if (order == 0) output = 1 - plhorror(A)
    if (order == 1) {
      if (A > 1)  output = exp(-A) + as.numeric(gammainc(A,0)[2])
      if (A <= 1) output = exp(-A) + integrate(function(x)exp(-x)/x,lower = A,upper = Inf)$value
    }

    if ( !(order %in% c(0,1)) ) output = (1+ 1/(order-1))*as.numeric(gammainc(A,order)[2]) -  A^(order -1)*exp(-A)/(order-1)
  }
  return(output)
}
