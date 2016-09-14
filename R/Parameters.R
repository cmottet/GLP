getOptimParam = function(sample,a,m = NULL,d=NULL,mTruth = NULL,dTruth = NULL,seed = NULL,positive = TRUE)
{

  if (length(m)!=0 & length(d)!=0 &  0%in% m & 0%in% d) m = m[-(m==0)]
  if (any(d > 3)) {print("Cannot estimate  density derivative or order higher than 2") ; d = d[d <=2]}

  param = names =  NULL
  cover = matrix(TRUE,ncol = 1,nrow = length(a))

  # Set the seed so that all the parameters get
  # estimated on the same set of bootstrapepd samples
  if (is.null(seed)) seed = floor(runif(1,0,1e8))


  ## Adjust Bonferroni's correction based
  ## on the total number of CI's to estimate
  bonferroni = length(m) + length(d)

  ##
  ## Compute the derivatives estimates
  ##
  if (length(d) >=1)
  {
    names = c(names,paste(rep(paste("d",d,sep=""),each = 2),c("L","U"),sep="")  )

    for (order in d)
    {
      ##
      ## Compute the tail distribution estimation
      ##
      if (order == 0)
      {
        args = list(eval.points = a, positive = positive, func = kcde)
        tmp  = CI.bootstrap(sample,f1,args,bonferroni = bonferroni,seed = seed)
        newparam = matrix(cbind(1-tmp$CI$upper,1-tmp$CI$lower),nrow = length(a))
        param   = cbind(param, newparam)
      }

      ##
      ## Compute the moments estimates
      ##
      if (order >=1){
        args = list(eval.points = a, positive = positive,deriv.order = order-1, func = kdde)
        tmp = CI.bootstrap(sample,f2,args,bonferroni = bonferroni,seed = seed)
        newparam = t(apply((-1)^(order+1)*matrix(cbind(tmp$CI$lower,tmp$CI$upper),nrow = length(a)),1,sort))
        param = cbind(param,newparam)
      }

      if (!is.null(dTruth))
        cover = cover & apply(newparam,1,is.between,x = subset(dTruth,select = paste("d",order,sep ="") ))

    }
  }

  if( !(0 %in% d) ) {
    param   = cbind(param, matrix(rep(c(0,1),length(a)),ncol = 2,byrow = TRUE)) ; names = c(names,"d0L","d0U")
  }


  if (length(m) >=1)
  {
    names = c(names,paste(rep(paste("m",m,sep=""),each = 2),c("L","U"),sep="")  )

    for (order in m)
    {
      fboot =  eval(substitute(function(sample,args) matrix(mean(sample^i*(sample>=args))),list(i=order)))
      tmp = sapply(a,CI.bootstrap,data = sample,fboot = fboot,bonferroni = bonferroni,seed = seed)
      newparam = matrix(unlist(tmp[2,]),nrow = length(a),ncol = 2,byrow = TRUE)
      param    =  cbind(param,newparam)

      if (!is.null(mTruth))
        cover = cover & apply(newparam,1,is.between,x = subset(mTruth,select = paste("m",order,sep ="") ) )
    }
  }

  ##
  ## Structure output
  ##
  output = data.frame(matrix(param,nrow = length(a)))
  names(output) = names

  if (!is.null(mTruth) & !is.null(dTruth)) output$cover = cover

  return(output)
}


