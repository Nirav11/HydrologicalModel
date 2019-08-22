constructDF <- function(x){
  
  DF <- iDF
  
  num  <- nrow(DF)
  
  DF$JD       <- (as.POSIXlt(as.Date(DF$Date, "%m/%d/%Y"))$yday)+1
  DF$dr       <- 1+0.033*cos(2*pi*DF$JD/365)
  DF$delta    <- 0.4093*sin(2*pi*DF$JD/365-1.405)
  DF$omegas   <- acos(-tan(phi)*tan(DF$delta)) 
  DF$Tsunset  <- DF$omegas/0.2618
  DF$N        <- 24/pi*DF$omegas
  DF$R0       <- 15.392*DF$dr*(DF$omegas*sin(phi)*sin(DF$delta)+cos(phi)*cos(DF$delta)*sin(DF$omegas))
  DF$PETt     <- 0.0023*DF$R0*(DF$Tavg + 17.8)*(sqrt((DF$Tmax)-(DF$Tmin))*30)
  DF$MtP      <- ifelse(DF$Tmax < Tb, 0, ifelse(e*(DF$Tavg+2.4) < 0, 0, e*(DF$Tavg+2.4)))
  
  DF$At       <- 0
  
  
  with(DF,
       At[2:num]  <- ifelse(Tmax[2:num]<Tb, ifelse((At[1:num-1] + Precip[2:num]-MtP[2:num])<0, 0, At[1:num-1] + Precip  - MtP[2:num]),
                            ifelse((At[1:num-1] - MtP[2:num])<0,0,At[1:num-1] - MtP[2:num])))
  DF$MtA       <- 0
  with(DF,
       
       MtA[2:num] <-ifelse((At[1:num-1]-At[2:num])<0,0,At[1:num-1]-At[2:num]))
  
  
  DF$PeffT     <-ifelse(DF$Tmax>=Tb,DF$Precip+DF$MtA,DF$MtA)
  
  
  for (i in 2:num) {
    DF$XUt[1]   <-XU0
    DF$Wat[i]   <-DF$PeffT[i]+DF$XUt[i-1]
    DF$EOt[i]   <-((DF$Wat[i]+x[2])/(2*x[1]))-sqrt(((DF$Wat[i]+x[2])/(2*x[1]))^2-(DF$Wat[i]*(x[2]/x[1])))
    
    DF$XUt[i]    <-(DF$EOt[i]^exp(-DF$PETt[i]/x[2]))
    
  }
  
  DF$Rt[2:num]   <-x[3]*(DF$Wat[2:num]-DF$EOt[2:num])
  DF$QUt[2:num]   <-(1-x[3])*(DF$Wat[2:num]-DF$EOt[2:num]) 
  
  for (i in 2:num) { 
    DF$XLt[1]     <-0
    DF$XLt[i]     <-(DF$Rt[i]+DF$XLt[i-1])/(1+x[4])
  }
  
  
  DF$QLt[2:num]         <- x[4]*DF$XLt[2:num]
  DF$Qt[2:num]          <- DF$QUt[2:num]+DF$QLt[2:num]
  DF$Residual[2:num]    <-DF$Qt[2:num]-DF$StrQobs[2:num]
  DF$SqResidual[2:num]  <- round((DF$Residual[2:num]^2),2)
  DF$SqD[2:num]         <- round(((DF$StrQobs[1:num-1]-mean(DF$StrQobs[1:num]))^2),2)
  
  return(DF)
  
}


NASH.Model <- function(y){
  
  DF <- constructDF(y)
  num  <- nrow(DF)
  
  SRavg <- mean(DF$SqResidual[2:num])
  
  SQDavg <- mean(DF$SqD[2:num])
  
  Nash <- 1-(SRavg/SQDavg)
  
  
  return(Nash)
}

RMSE.Model <- function(y){
  
  DF <- constructDF(y)
  num  <- nrow(DF)
  
  SRavg <- mean(DF$SqResidual[2:num])
  
  SQDavg <- mean(DF$SqD[2:num])
  
  RMSE <- sqrt(SRavg)
  
  return(RMSE)
}