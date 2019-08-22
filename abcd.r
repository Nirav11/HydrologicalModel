# x[1]   <- 00.8118
# x[2]   <- 500.0000
# x[3]   <- 0.6661
# x[4]   <- 00.0006
require(dplyr)
require(units)
require(lubridate)
require(hydromad)
require(ggplot2)
require(scales)
script.dir <- dirname(sys.frame(1)$ofile)
source(paste(script.dir,  "modelfx.R", sep="/"))

lat    <<- -4.2
phi    <<- lat*pi/180

Darea  <<- 2250


e   <<- 0
A0  <<- 0
Tb  <<- 2
XU0 <<- 300

iDF  <<- read.csv(paste(script.dir,"drainage2.csv", sep="/"))

o <- optim(c(0.8,500,0.66,0.0006), 
           NASH.Model, 
           control = list(fnscale = -1), 
           lower = c(0.1,0.1,0.1,0.1), 
           upper=c(.9,500,.9,.9), 
           method="L-BFGS-B")

D <- constructDF(o$par)
write.csv(D, paste(script.dir,"output.csv", sep="/"))

ggplot(D, mapping = aes(x = as.numeric(Date), y = XUt,XLt)) +
geom_line(aes(y=XUt),colour="red", size = .75)+
geom_line(aes(y=XLt/.25),colour="blue", size = .75)+
scale_y_continuous(name="Soil Moisture", sec.axis = sec_axis(~.*.25,name = "Groundwater Storage"))+labs(XUt = "Soil Moisture")#+theme(axis.text.x = element_text(angle = 90, hjust = 1))
 



ggplot(D, mapping = aes(x = as.numeric(Date), y = Tmax,Tmin,Tavg)) +
  geom_line(aes(y=Tmax),colour="red", size = .75)+
  geom_line(aes(y=Tmin),colour="blue", size = .75)+geom_line(aes(y=Tavg),colour="green", size = .75)+
  scale_y_continuous(name="Temp")

ggplot(D, mapping = aes(x = as.numeric(Date), y = Precip,PETt)) +
  geom_line(aes(y=Precip),colour="red", size = .75)+
  geom_line(aes(y=PETt/.5),colour="blue", size = .75)+
  scale_y_continuous(name="Precipitation(mm)", sec.axis = sec_axis(~.*.5,name = "PET(mm)"))

ggplot(D, mapping = aes(x = as.numeric(Date), y = StrQobs,Qt)) +
  geom_line(aes(y=StrQobs),colour="red", size = 1)+
  geom_line(aes(y=Qt/.5),colour="blue", size = .75)+
  scale_y_continuous(name="Flow(mm)")
