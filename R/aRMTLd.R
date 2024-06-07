#' Title Comparing the difference of adjusted restricted mean time lost (aRMTLd) under competing risks
#'
#' @param time The follow-up time for right censored data.
#' @param status The status indicator, 0 = right censored, 1 = event of interest, 2 = competing event.
#' @param group The group indicator for comparison. The elements of this vector take either 1 or 0.
#' Normally, 0 = control group, 1 = treatment group.
#' @param covariates A matrix of covariates.
#' @param tau The value to specify the truncation time point for the RMTL calculation.
#'
#' @return The value of aRMTLd, the stand error and 95% confidence interval of aRMTLd。
#' @export
#'
#' @examples
#' library(aRMTL)
#' data(simdata)
#' time=simdata$time;status=simdata$status;group=simdata$group
#' covariates <- cbind(simdata$z1,simdata$z2,simdata$z3,simdata$z4)
#' aRMTLd(time,status,group,covariates,tau = 365)
aRMTLd <- function(time,status,group,covariates,tau){

  IPW_data <- data.frame(group,covariates)
  xnam <- paste0("X", 1:ncol(covariates))
  glm_mod <- glm(as.formula(paste("group ~ ", paste(xnam, collapse= "+"))) , data=IPW_data, family="binomial")
  pred = predict(glm_mod, type='response')
  weight = group/pred + (1-group)/(1-pred)

  data <- data.frame(time, status, group, weight, pred, covariates)
  data <- data[order(group),]
  Armtl=rep(0,2)
  Armtl_se=rep(0,2)
  var0=vector()
  var1=vector()

  for (gg in 1:2) {

    # Calculate cif
    dd=data[data$group==gg-1,]
    tj <- sort(unique(dd$time[dd$status!=0]))
    dj1 <- sapply(tj, function(x){sum(dd$weight[dd$time==x & dd$status==1])})
    dj2 <- sapply(tj, function(x){sum(dd$weight[dd$time==x & dd$status==2])})
    yj <- sapply(tj, function(x){sum(dd$weight[dd$time>=x])})
    hj1 = dj1/yj
    hj2 = dj2/yj
    Hj1 <- cumsum(hj1)
    Hj2 <- cumsum(hj2)
    Hj <- Hj1 + Hj2
    s  <- c(1,exp(-Hj))
    s <- s[-length(s)]
    cif1=cumsum(s*(hj1))
    cif2=cumsum(s*(hj2))

    rtime <- tj<=tau
    tj_r <- sort(c(tj[rtime],tau))
    time_diff <- diff(tj_r)
    s_r=s[rtime]
    cif1_r=cif1[rtime]
    cif2_r=cif2[rtime]
    rmtl=sum(time_diff*cif1_r)
    Armtl[gg]=rmtl

    # Calculate variance
    pre=data$pred
    X=as.matrix(cbind(1,data[,6:ncol(data)]))
    omega=matrix(0,ncol(X),ncol(X))
    n=nrow(data)
    for(k in 1:n)
    {
      omega=omega+pre[k]*(1-pre[k]) * X[k,] %o% X[k,]
    }
    omega=solve(omega/n)

    ##### Calculate hjk(t) #####
    # event 1
    yyj=matrix(NA,ncol(X),length(tj))
    ddj=matrix(NA,ncol(X),length(tj))
    h1=matrix(NA,ncol(X),length(tj))
    h2=matrix(NA,ncol(X),length(tj))
    wt=dd$weight
    X_group=as.matrix(cbind(1,dd[,6:ncol(data)]))
    a <- t( X_group * (wt-1) *(-1)^dd$group)

    for (i in 1:ncol(X)) {
      dd$a=a[i,]
      yyj[i,] <- sapply(tj, function(x){sum(dd$a[dd$time>=x])})
      ddj[i,] <- sapply(tj, function(x){sum(dd$a[dd$time==x & dd$status==1])})
      h1[i,] <- ddj[i,]
      h2[i,] <- dj1*(yyj[i,])
    }

    H1=matrix(NA,ncol(X),length(tj))
    H2=matrix(NA,ncol(X),length(tj))
    HH1=matrix(0,ncol(X),length(tj))
    for(i in 1:length(tj)){
      H1[,i]=h1[,i]/(yj[i])
      H2[,i]=h2[,i]/(yj[i])^2
    }
    HH1=t(apply(H1-H2,1,cumsum))
    HH1_r=HH1[,rtime]

    # event 2
    yyj=matrix(NA,ncol(X),length(tj))
    ddj=matrix(NA,ncol(X),length(tj))
    h1=matrix(NA,ncol(X),length(tj))
    h2=matrix(NA,ncol(X),length(tj))

    for (i in 1:ncol(X)) {
      dd$a=a[i,]
      yyj[i,] <- sapply(tj, function(x){sum(dd$a[dd$time>=x])})
      ddj[i,] <- sapply(tj, function(x){sum(dd$a[dd$time==x & dd$status==2])})
      h1[i,] <- ddj[i,]
      h2[i,] <- dj2*(yyj[i,])
    }

    H1=matrix(NA,ncol(X),length(tj))
    H2=matrix(NA,ncol(X),length(tj))
    HH2=matrix(0,ncol(X),length(tj))
    for(i in 1:length(tj)){
      H1[,i]=h1[,i]/(yj[i])
      H2[,i]=h2[,i]/(yj[i])^2
    }
    HH2=t(apply(H1-H2,1,cumsum))
    HH2_r=HH2[,rtime]

    ##### Calculate alpha1 #####

    X2 <- t(X * (data$group-pre))

    ALPHA1_1=matrix(NA,n,length(cif1_r))
    for (j in 1:n) {
      XX2=X2[,j]
      alpha1=rep(NA,length(cif1_r))
      X1=matrix(NA,ncol(X),length(cif1_r))
      for (i in 1:length(cif1_r)) {
        X1[,i]=HH1_r[,i]%*%omega
        alpha1[i]=X1[,i]%*%XX2
      }
      ALPHA1_1[j,]=alpha1
    }

    ALPHA1_2=matrix(NA,n,length(cif1_r))
    for (j in 1:n) {
      XX2=X2[,j]
      alpha1=rep(NA,length(cif1_r))
      X1=matrix(NA,ncol(X),length(cif1_r))
      for (i in 1:length(cif1_r)) {
        X1[,i]=HH2_r[,i]%*%omega
        alpha1[i]=X1[,i]%*%XX2
      }
      ALPHA1_2[j,]=alpha1
    }


    ##### Calculate dMjk #####
    # event 1
    Mi_1=matrix(0,nrow(dd),length(cif1_r))
    for (i in 1:nrow(dd)) {
      dNi1=ifelse(dd$time[i]==tj & dd$status[i]==1,dd$weight[i],0)
      dYi=ifelse(dd$time[i]>=tj,dd$weight[i],0)*hj1
      dMi=(dNi1-dYi)/(yj/n)
      Mi_1[i,]=cumsum(dMi)[rtime]
    }

    if(gg==1){Mi1=rbind(Mi_1,matrix(0,n-nrow(dd),length(cif1_r)))} else {Mi1=rbind(matrix(0,n-nrow(dd),length(cif1_r)),Mi_1)}

    # event 2
    Mi_2=matrix(0,nrow(dd),length(cif1_r))
    for (i in 1:nrow(dd)) {
      dNi2=ifelse(dd$time[i]==tj & dd$status[i]==2,dd$weight[i],0)
      dYi=ifelse(dd$time[i]>=tj,dd$weight[i],0)*hj2
      dMi=(dNi2-dYi)/(yj/n)
      Mi_2[i,]=cumsum(dMi)[rtime]
    }

    if(gg==1){Mi2=rbind(Mi_2,matrix(0,n-nrow(dd),length(cif1_r)))} else {Mi2=rbind(matrix(0,n-nrow(dd),length(cif1_r)),Mi_2)}

    fy1=ALPHA1_1+Mi1
    fy2=ALPHA1_2+Mi2

    ## Calculate variance
    var=rep(0,n)
    for (i in 1:n) {
      fy=c(0,fy1[i,]+fy2[i,])
      fy=fy[-length(fy)]
      var[i]=sum(cumsum( diff(c(0,fy1[i,])) * s_r)*time_diff) -   sum( cumsum(diff(c(0,cif1_r))* fy )*time_diff )
    }

    if(gg==1) {var0=var} else {var1=var}
    Armtl_se[gg]=sqrt(sum(var^2))/n

  }

  rmtl_diff=Armtl[2]-Armtl[1]
  rmtl_diff_se=sqrt(sum((var0-var1)^2))/n
  p=2*(1-pnorm(abs(rmtl_diff)/rmtl_diff_se))
  cil=rmtl_diff-(qnorm(1-0.05/2)*rmtl_diff_se)
  ciu=rmtl_diff+(qnorm(1-0.05/2)*rmtl_diff_se)

  result=data.frame(tau,Armtl[1],Armtl_se[1],Armtl[2],Armtl_se[2],rmtl_diff,rmtl_diff_se,p,cil,ciu)
  colnames(result) <- c("tau","aRMTL(Group = 0)","aRMTL_SE(Group = 0)",
                        "aRMTL(Group = 1)","aRMTL_SE(Group = 1)",
                        "aRMTLd(Group1 - Group0)","aRMTLd_SE","P_value",
                        "aRMTLd_CI_lower","aRMTLd_CI_upper")
  return(result)

}



