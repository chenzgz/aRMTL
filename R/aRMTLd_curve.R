#' Title The difference in adjusted restricted mean time lost (aRMTLd) curves under competing risks
#'
#' @param B The number of random variables generated of standard normal distribution. In example B=10 is chosen to calculate the time,
#' which can be 200,500 or 1000 for practical applications.
#' @param tL The beginning of aRMTLd curve time interval.
#' @param tU The end of aRMTLd curve time interval.
#' @param time The follow-up time for right censored data.
#' @param status The status indicator, 0 = right censored, 1 = event of interest, 2 = competing event.
#' @param group The group indicator for comparison. The elements of this vector take either 1 or 0.
#' Normally, 0 = control group, 1 = treatment group.
#' @param covariates A matrix of covariates.
#'
#' @return A six-column data box, where tt is listed as the limit time point,
#' diff is listed as the aRMTLd, ci_lower and ci_upper are the upper and lower bounds of the confidence interval,
#' and cb_lower and cb_upper are the lower bounds of the confidence band.
#' @export
#'
#' @examples
#' library(aRMTL)
#' data(simdata)
#' time=simdata$time;status=simdata$status;group=simdata$group
#' covariates <- cbind(simdata$z1,simdata$z2,simdata$z3,simdata$z4)
#' set.seed(2024)
#' B=10
#' tL <- ceiling(max(min(time[group==0&status==1]), min(time[group==1&status==1])))
#' tU <- 365
#' result <- aRMTLd_curve(B,tL,tU,time,status,group,covariates)
#' plot(result$tt,result$diff,type = "l",ylim = c(-105,35),xlim = c(0,365),xlab = "time(days)" , ylab = "aRMTLd",
#' col = "brown3",main = "Group1 - Group0: aRMTLd curve plot",lwd=3)
#' lines(result$tt,result$ci_lower,type = "l",col = "brown3",lwd=3,lty=2)
#' lines(result$tt,result$ci_upper,type = "l",col = "brown3",lwd=3,lty=2)
#' polygon(c(result$tt,rev(result$tt)),c(result$cb_lower,rev(result$cb_upper)),
#' col=rgb(1,0,0,0.15),border=FALSE)
#' legend("topleft", legend = "aRMTLd", col=c("brown3"),  lwd=3, cex=1, bty ="n", inset = c(0, 0))
#' abline(h=0,col = "gray30",lwd=3,lty=3)
aRMTLd_curve <- function(B,tL,tU,time,status,group,covariates){
  aRMTLd_sample <- function(time,status,group,covariates,tau){

    IPW_data <- data.frame(group,covariates)
    xnam <- paste0("X", 1:ncol(covariates))
    glm_mod <- glm(as.formula(paste("group ~ ", paste(xnam, collapse= "+"))) , data=IPW_data, family="binomial")
    pred = predict(glm_mod, type='response')
    weight = group/pred + (1-group)/(1-pred)

    data <- data.frame(time, status, group, weight, pred, covariates)
    data <- data[order(group),]
    Armtl=rep(0,2)
    rmtl_diff_se_sample=vector()
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
      #event 1
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

      #event 2
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

      ##### calculate alpha1 #####

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


      #####calculate dMjk #####
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

      sample=var*rnorm(n)
      if(gg==1) {var0=sample} else {var1=sample}

    }

    rmtl_diff=Armtl[2]-Armtl[1]
    rmtl_diff_se_sample=sum(var1-var0)/n

    result=data.frame(tau,rmtl_diff,rmtl_diff_se_sample)
    return(result)

  }
  tt <- c(tL,sort(unique(time[status==1 & time>=tL & time<=tU])),tU)
  sample <- matrix(0,B,length(tt))

  if (is.element("tcltk", installed.packages()[,1])==FALSE){
    install.packages("tcltk")
  }
  require(tcltk) # same as library statement

  pb <- tkProgressBar("Progress bar","completed %",0,100)
  star_time <- Sys.time()
  for (j in 1:B) {
    for (i in 1:length(tt)) {
      sample[j,i] <- aRMTLd_sample(time,status,group,covariates,tt[i])$rmtl_diff_se_sample
    }
    info <- sprintf("completed %d%%",round(j*100/B))
    setTkProgressBar(pb, j*100/B, sprintf("Progress bar (%s)",info),info)
  }
  end_time <- Sys.time()
  close(pb)
  run_time <- end_time - star_time

  se <- rep(0,length(tt))
  diff <- rep(0,length(tt))
  upper <- rep(0,length(tt))
  lower <- rep(0,length(tt))
  for (i in 1:length(tt)) {
    res <- aRMTLd(time,status,group,covariates,tt[i])
    diff[i] <- res[[6]]
    se[i] <- res[[7]]
    upper[i] <- res[[10]]
    lower[i] <- res[[9]]
  }
  qt <- apply(abs(t(sample)/se), 2, max)
  q95 <- quantile(qt,0.95)

  upper_band <- diff + se*q95
  lower_band <- diff - se*q95

  result <- data.frame(tt,diff,ci_lower=lower,ci_upper=upper,cb_lower=lower_band,cb_upper=upper_band)
  return(result)
}





