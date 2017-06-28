
# Author: Wells Bishop
# Class:  Physics 481 Final Project
# Title:  "Model Comparison and Analysis of Eugene Rainfall"

# Known bug: You can scroll through the different plots using arrows above plot window, however R seems to run
# out of system memory at a certain number of concurent plots. If unable to switch between plots, either clear
# all plots with the broom above the window, or go to session/restart R, this should solve it.


options(digits=12)
graphics.off()
library("stats")



normLikelihood = function( data, muList, sigmaList ) {
  # Returns a matrix of the likelihood of a normal distribution over the given parameter space.
  
  likelihood = matrix(0, nrow = length(muList), ncol = length(sigmaList))
  N = length(data)
  
  for (i in 1:length(muList)){
    for (j in 1:length(sigmaList)){
      
      Mu = muList[i]
      Sig = sigmaList[j]

      like = -0.5*N*log(2*pi) - 0.5*N*log(Sig^2) - 0.5*(1/Sig^2) * sum( (data - Mu)^2, nan.rm=TRUE )

      likelihood[j, i] = like

    }
  }
  likelihood = exp(likelihood-max(likelihood)) # Compute the log-likelihood and then raise it to the exponential when its returned.
  return( likelihood )
}
normLike = function(data, N, bounds){
  # Computes P(Data| m=Normal Distribution). 
  # Or the probability over all parameter space (for given bounds) that our data fits this model.

  mu.min = bounds[1] # Upper and lower bounds
  mu.max = bounds[2]
  sig.min = bounds[3]
  sig.max = bounds[4]
  
  dMu = (mu.max-mu.min) / N # Compute dMu and dSig
  dSig = (sig.max-sig.min) / N
  muList = seq(mu.min, mu.max, by=dMu)
  sigList = seq(sig.min, sig.max, by=dSig)
  
  like = normLikelihood(data, muList, sigList) # Run the above function over the given parameter space
  like = like * dMu * dSig # Multiply each point by dMu and dSig
  
  return(sum(like, na.rm=TRUE)) # Return the Riemman sum value over the parameter bounds
}






logisticLikelihood = function(data, muList, sList){
  # Returns a matrix of the likelihood at each paramter value for a logistic distribution.
  
  likelihood = matrix(0, nrow = length(muList), ncol = length(sList))
  N = length(data)
  
  for (i in 1:length(muList)){
    for (j in 1:length(sList)){
      
      Mu = muList[i]
      S = sList[j]
      
      like = sum(-1/S * (data - Mu)) - N*log(S) - 2* sum(log( 1 + exp( -(data - Mu)/S )), nan.rm=TRUE)

      likelihood[j, i] = like
    }
  }
  likelihood = exp(likelihood - max(likelihood))
  return( likelihood )
}
logisticLike = function(data, N, bounds){
  # Computes P(Data| m=Logistics).
  
  mu.min = bounds[1]
  mu.max = bounds[2]
  s.min = bounds[3]
  s.max = bounds[4]
  
  dMu = (mu.max-mu.min) / N
  dS = (s.max-s.min) / N
  muList = seq(mu.min, mu.max, by=dMu)
  sList = seq(s.min, s.max, by=dS)
  
  like = logisticLikelihood( data, muList, sList)
  like = like * dMu * dS

  return(sum(like, na.rm=TRUE))
}





cauchyLikelihood = function(data, muList, gammaList){
  # Returns a matrix of the likelihood at each paramter value for a Cauchy-Lorentz distribution
  
  likelihood = matrix(0, nrow = length(muList), ncol = length(gammaList))
  N = length(data)
  
  for (i in 1:length(muList)){
    for (j in 1:length(gammaList)){
      
      Mu = muList[i]
      Gamma = gammaList[j]
      
      like = -N*log(Gamma*pi) - sum(log(1 + ((data-Mu)/Gamma)^2), nan.rm=TRUE)
      likelihood[j, i] = like 
    }
  }
  likelihood = exp(likelihood - max(likelihood))
  return( likelihood )
}
cauchyLike = function(data, N, bounds){
  # Computes P(Data| m=Cauchy).
  
  mu.min = bounds[1] 
  mu.max = bounds[2]
  gamma.min = bounds[3]
  gamma.max = bounds[4]
  
  dMu = (mu.max-mu.min) / N
  dGamma = (gamma.max-gamma.min) / N
  muList = seq(mu.min, mu.max, by=dMu)
  gammaList = seq(gamma.min, gamma.max, by=dGamma)
  
  like = cauchyLikelihood(data, muList, gammaList)
  like = like * dMu * dGamma

  return(sum(like, na.rm=TRUE))
}




# ____________________________________________________________________________________________ #
# Eugene precipitation data from 1939-2014. 
# Organize data into useful form. Looking at the change in rainfall between January and July 
# of the same year.

data2 = read.csv("EUGpcpn.csv", header = FALSE, skipNul = TRUE)
data3 = numeric()
data.jan = numeric()
data.july = numeric()
data.change = numeric()

data.winter = numeric()
data.summer = numeric()
data.change2 = numeric()

for (k in seq(1, 912)){
  data3[k] = data2[k, 34]
  if (k %% 12 == 1){ # Pick out rainfall in January
    data.jan = c(data.jan, data3[k])
  }
  if (k %% 12 == 7){ # Pick out rainfall in July
    data.july = c(data.july, data3[k])
  }
}

m12 = seq(0,length(data3), by=12) # Another way to pick out months of data. m12=December, m1=January.
m1 = seq(1,length(data3), by=12)
m2 = seq(2,length(data3), by=12)

m6 =  seq(6,length(data3), by=12) # The three summer months
m7 =  seq(7,length(data3), by=12)
m8 =  seq(8,length(data3), by=12)

data.winter = (data3[m12]+data3[m1]+data3[m2])/3 # Taking the average of the summer and winter months
data.summer = (data3[m6]+data3[m7]+data3[m8])/3


for (i in 1:(length(data.jan))){ # Looks at the difference of rain between seasons
  delta = data.jan[i] - data.july[i]
  delta2 = data.winter[i] - data.summer[i]
  data.change = c(data.change, delta) # data.change is the data being analyzed
  data.change2 = c(data.change2, delta2)
}


data.change = data.change2 # Data.change2 corresponds to the updated and improved dataset of winter-summer months.



# ____________________________________________________________________________________________ #
# Model Comparison between Normal(Gaussian) , Logistic , and Student's t distributions.
# Change bounds to change which values the integral goes over.
# Change N to change how many sections the Riemman sum breaks each dimension up into. Greater than 2,000 is not recommended.

N = 100

show("Values for the Riemman sum of the three different models: ")
show("(Normal, Logistic, Cauchy)")


bounds.norm = c( -100, 500, 0.1, 500 ) # mu.min, mu.max, sig.min, sig.max
like.norm = normLike( data.change, N, bounds.norm )
show(like.norm)


bounds.logistic = c( -100, 500, 0.1, 500 ) # mu.min, mu.max, s.min, s.max
like.logistic = logisticLike( data.change, N, bounds.logistic )
show(like.logistic)


bounds.cauchy = c( -100, 500, 0.1, 500 ) # mu.min, mu.max, gamma.min, gamma.max
like.cauchy = cauchyLike( data.change, N, bounds.cauchy )
show(like.cauchy)





# ____________________________________________________________________________________________ #
# Plots #

# Histogram of original data and difference seasons:

#hist(data3*0.01, breaks = 40,col="gray", xlab = "Inches of Rain", main = "Monthly Rainfall (ALL)")
#hist(data.jan*0.01, breaks = 40,col="gray", xlab = "Inches of Rain", main = "Monthly Rainfall (JANUARY)")
hist(data.change*0.01, breaks = 35,col="gray", xlim=c(-4,4),xlab = "Inches of Rain", 
     main = "Difference in Rainfall Between an Average Winter and Summer Month")



#_____________________________________________________________________________________________________#
# ||| Normal distribution analysis:
# Specify range of parameters for Normal distribution contour plot over all parameters:
muList = seq(0, 100, length.out = 50)
sigmaList = seq(20, 200, length.out = 50)

like1 = normLikelihood( data.change, muList, sigmaList)

filled.contour(x=sigmaList, y=muList, z=like1/sum(like1) , 
               xlab="Sigma", ylab="Mu", main="Gaussian" ,
               color.palette = colorRampPalette(c("lightblue", "blue", "yellow", "red")))

# Create synthetic Normal data from our best parameter values from the contour plot above and compares to original:
data.norm.synth = rnorm(length(data.change) , mean=20, sd=135)

plot(sort(data.norm.synth), (1:length(data.norm.synth))/length(data.norm.synth),
     xlim = c(-400, 400), main="ECDF comparing Normal data")
lines(sort(data.change), (1:length(data.change))/length(data.change), col="red", lwd=3)
legend("topleft",legend=paste(c("Synthetic Data", "Real data")),pch=c("o","l"), col=c("black", "red"))



# _____________________________________________________________________________________________________#
# ||| Logistic distribution:
muList2 = seq(1, 100, length.out = 50)
sList = seq(20, 150, length.out = 50)

like2 = logisticLikelihood( data.change, muList2, sList)

filled.contour(x=sList, y=muList2, z=like2/sum(like2) , 
               xlab="S", ylab="Mu", main="Logistic",
               color.palette = colorRampPalette(c("lightblue", "blue", "yellow", "red")))

# Create synthetic data with parameters from contour plot:
data.logis.synth = rlogis(length(data.change) , location=28, scale=75) # FIx values!!!!!!!!

plot(sort(data.logis.synth), (1:length(data.logis.synth))/length(data.logis.synth),
     xlim = c(-400, 400), main="ECDF comparing Logistic data")
lines(sort(data.change), (1:length(data.change))/length(data.change), col="red", lwd=3)
legend("topleft",legend=paste(c("Synthetic Data", "Real data")),pch=c("o","l"), col=c("black", "red"))



#_________________________________________________________________________________________________________#
# ||| Cauchy distribution:
muList3 = seq(1, 100, length.out = 50)
gammaList = seq(20, 180, length.out = 50)

like3 = cauchyLikelihood( data.change, muList3, gammaList)

filled.contour(x=gammaList, y=muList3, z=like3/sum(like3) , 
               xlab="Gamma", ylab="Mu", main="Cauchy",
               color.palette = colorRampPalette(c("lightblue", "blue", "yellow", "red")))

# Create synthetic data with parameters from contour plot:
data.cauchy.synth = rcauchy(length(data.change) , location=40, scale=90)

plot(sort(data.cauchy.synth), (1:length(data.cauchy.synth))/length(data.cauchy.synth),
     xlim = c(-400, 400), main="ECDF comparing Cauchy data")
lines(sort(data.change), (1:length(data.change))/length(data.change), col="red", lwd=3)
legend("topleft",legend=paste(c("Synthetic Data", "Real data")),pch=c("o","l"), col=c("black", "red"))




#________________________________________________________________________________________________#
# Show that the Cauchy distribution with given parameters fits the original data histogram the best

k = seq(-400, 400, length.out = 800)

hist(data.change, breaks = 35,col="gray", xlim=c(-400,400),xlab = "Hundredths of an Inch of Rain", 
     main = "Difference in Rainfall Between an Average Winter and Summer Month")
lines(k,dcauchy(k, location = 40, scale=90)*2900, col="blue", lwd=4) # dcauchy is the Cauchy function in the "stats" package.



