#Dead carbon fraction (reservoir effect) calculator for stalagmite datasets.
#Written by Franziska Lechleitner and Jens Fohlmeister

####Load packages
library(tidyverse)

####Load input data
#Atmospheric calibration curve (IntCal)
Intcal <- read_csv("input/Intcal13.csv") 

#Stalagmite 14C data
#FORMATTING NOTE FOR DATA: U-Th ages in yr BP (with BP referring to 1950CE), U-Th age error as 2 sigma.
stal_data <- read_csv("input/Stalagmite_input.csv") #Here list your stalagmite dataset


######Calculate atmospheric 14C values corresponding to stalagmite######
#Montecarlo simulation to find the best fitting atmospheric 14C age corresponding to stalagmite data.
#Because IntCal is at 5-10 year resolution, we first interpolate it to yearly resolution, so the matches with 
#the stalagmite data can be found easily.

#This bit of the code was taken from Jens Fohlmeister's Matlab code, credit goes to him.

anzahl = 10000; #Number of MC simulations (normally distributed)       
  
#Interpolate IntCal to yearly resolution
c14age_interp <- as.data.frame(approx(Intcal$Date, Intcal$c14age,       method = "linear", n = 50001), col.names = c("Date", "C14.age"))  
d14C_interp   <- as.data.frame(approx(Intcal$Date, Intcal$D14C,         method = "linear", n = 50001), col.names = c("Date", "D14C"))   
errd14C_interp<- as.data.frame(approx(Intcal$Date, Intcal$D14C.error,   method = "linear", n = 50001), col.names = c("Date", "D14C.error")) 
Intcal_interp <- merge(c14age_interp, d14C_interp,     by = "Date")
Intcal_interp <- merge(Intcal_interp, errd14C_interp,  by = "Date")
  
UThage <- as.array(stal_data$UThAge)
sigma <-  as.array(stal_data$UThAgeError/2) #Transform U-Th error to 1 sigma

#Produce 10000 realisations of ages around a stalagmite U-Th ages (normally distributed)
#This is used to estimate the uncertainty in DCF. By producing a large ensemble of ages 
#we can propagate the uncertainty from the U-Th dating to the DCF calculation.
x = as.data.frame(matrix(NA, nrow = anzahl, ncol = 1:dim(UThage))) 
for(n in 1:dim(UThage)) { 
  x[,n] <- floor(UThage[n] + sigma[n]  * rnorm(anzahl)); 
}

#Extract the IntCal values corresponding to the ensemble of 10000 possible ages and calculate the average and STDEV.
stats = as.data.frame(matrix(NA, nrow = 1:dim(UThage) , ncol = 5))
for(i in 1:dim(UThage)) {
x_new <- x[,i]
Intcal_new <- merge(x_new, Intcal_interp, by.x = 1, by.y = 'Date', sort = F) 
avg.Date    <- mean(Intcal_new$x)
avg.C14Cage <- mean(Intcal_new$C14.age)
avg.D14C    <- mean(Intcal_new$D14C)
STDEV.C14age <- sd(Intcal_new$C14.age)
STDEV.D14C   <- sd(Intcal_new$D14C)
stats [i,] <- c(avg.Date, avg.C14Cage, STDEV.C14age, avg.D14C, STDEV.D14C)
colnames(stats) <- c("avg.Date", "avg.C14Cage", "STDEV.C14age", "avg.D14C", "STDEV.D14C")
}


######Calculate reservoir age of the stalagmite#####
#By "reservoir age" we mean the difference between 14C age of the atmosphere and stalagmite.

res_age <- as.data.frame(stal_data$C14Age - stats$avg.C14Cage)
res_age_err <- sqrt(stal_data$C14AgeError^2 + stats$STDEV.C14age^2)
res_age_upper <- res_age + res_age_err
res_age_lower <- res_age - res_age_err
res.age <- data.frame(stal_data$SampleID, stal_data$UThAge, res_age, res_age_err, res_age_upper, res_age_lower)
colnames(res.age) <- c("Sample.ID", "UTh.Age", "reservoir.age", "reservoir.age.err", "reservoir.age.uci", "reservoir.age.lci")

#Plot reservoir age
ggplot()+
  geom_line(data = res.age,aes(y= reservoir.age.uci,x= UTh.Age),colour="grey")+
  geom_line(data = res.age,aes(y= reservoir.age.lci,x= UTh.Age),colour="grey")+
  geom_line(data=res.age,aes(y= reservoir.age,x= UTh.Age), colour="darkblue")+ 
  geom_point(data=res.age,aes(y= reservoir.age,x= UTh.Age), colour="darkblue")+ 
  xlab('U-Th age') +
  ylab('Reservoir age (yr)') 


#########Calculate DCF######################################################
#Calculate initial (decay corrected) stalagmite 14C
#14C decay constant
lambda <- 1/8267

#First convert the atmospheric 14C age to a14C and correct the stalagmite data for decay since deposition.
a14C_atm <- stats$avg.D14C/1000 + 1 
a14C_atm_sig <- stats$STDEV.D14C/1000

a14C_stal <- exp(-stal_data$C14Age/8033)
a14C_stal_sig <- a14C_stal * sqrt((stal_data$C14AgeError/stal_data$C14Age)^2) 
a14C_stal_ini <- a14C_stal * exp(UThage*lambda)
a14C_stal_ini_sig <- a14C_stal_ini * sqrt((a14C_stal_sig/a14C_stal)^2 + (sigma * lambda)^2)

dcf    <- (1- a14C_stal_ini/a14C_atm)*100;
dcferr <- (a14C_stal_ini/a14C_atm * sqrt((a14C_stal_ini_sig/a14C_stal_ini)^2 + (a14C_atm_sig/a14C_atm)^2)) *100

dcf_uc <- dcf + dcferr
dcf_lc <- dcf - dcferr
DCF <- data.frame(stal_data$SampleID, stal_data$UThAge, dcf, dcferr, dcf_uc, dcf_lc)
colnames(DCF) <- c("Sample.ID", "UTh.Age", "DCF", "DCF.err", "DCF.uci", "DCF.lci")

#Plot DCF
ggplot()+
  geom_line(data = DCF,aes(y= DCF.uci,x= UTh.Age),colour="grey")+
  geom_line(data = DCF,aes(y= DCF.lci,x= UTh.Age),colour="grey")+
  geom_line(data= DCF, aes(y= DCF,x= UTh.Age), colour="red")+ 
  geom_point(data=DCF, aes(y= DCF,x= UTh.Age), colour="red")+ 
  xlab('U-Th age') +
  ylab('DCF (%)') 

#Save output as csv file
write.csv(DCF, file = "DCF_output.csv")


