genSimulationFilesPlain <- function(simulations=2,
                               samples=200,
                               predictors=1000,
                               infoVars=50,
                               SNRy=10,
                               SNRx=2,
                               kappa=5){
  #### input parameters
  # simulations:  the number of simulation sets to be generated
  # samples:      number of samples to be generated in each set
  # predictors:   number of predictors per sample
  # infovars:     number of variables to be assigned a non-zero coefficient in the
  #               calculation of the response variable (simple linear combination)
  # SNRy:         signal-to-noise ratio, used when adding random noise to the response variable
  # SNRx:         signal-to-noise ratio, used when adding random noise to the predictor variables
  # kappa:        coefficient in the logit exponential, it strongly regulates 
  #               class overlap in the two-class problem
  
  #set informative filename from list of parameters
  fileName <- paste(as.list(environment()), collapse="_")
    
  #initialise variables
  X <- array(0, c(simulations, samples, predictors))
  Xnoise <- array(0, c(simulations, samples, predictors))
  Y <- matrix(0, nrow=samples, ncol=simulations)
  Ynoise <- matrix(0, nrow=samples, ncol=simulations)
  yBin <- matrix(0, nrow=samples, ncol=simulations)
  yBinNoise <- matrix(0, nrow=samples, ncol=simulations)
  
  coeffs <- as.vector(rep(0, predictors))
  #fixed seed for reproducibility
  set.seed(1234)
  infoIndices <- sample(c(1:predictors), infoVars)
  #create a set of coefficients, the same for the whole set of simulations
  infoCoefs <- rnorm(infoVars, 1, 1)  
  coeffs[infoIndices] <- infoCoefs
  
  #debug#########################
  #View(coeffs)
  #varY<-0
  #############################################
  #generate dataset for each simulation
  for (simulation in 1:simulations) {
    #print(simulation)
    #new seeds for reproducibility
    set.seed(1234 + simulation*samples*4) #multiply by four because we sample four times below
                                          #CHECK THIS AGAIN
    #generate the samples
    for (predictor in 1:predictors) {
      #just sample from a normal distribution for each voxel
      curDim <- rnorm(samples, 0, 1)
      X[simulation, , predictor] <- curDim
      
      #also generate a noisy predictor
      noise <- rnorm(samples)
      #calculate the adj coefficient from variance of the predictor signal and the desired SNRx
      noiseCoeffX <- sqrt(var(curDim)/(SNRx * var(noise)))
      #generate the response with noise
      Xnoise[simulation, , predictor] <- X[simulation, , predictor] + noiseCoeffX*noise
    }
    #generate the response (from clean predictors!)
    Y[,simulation] <- X[simulation, ,] %*% coeffs
    
    #generate some noise for Y
    noise <- rnorm(samples)
    #calculate the adj coefficient from variance of the signal and the desired SNRy
    #varY<-varY+var(Y[,simulation])
    noiseCoeffY <- sqrt(var(Y[,simulation])/(SNRy * var(noise)))
    #generate the response with noise
    Ynoise[,simulation] <- Y[,simulation] + noiseCoeffY*noise
    #generate two-group response (kappa regulates the steepness and largely the overlap)
    yLogit <- 1/(1+exp(-kappa*Y[,simulation]))
    yBin[,simulation] <- rbinom(samples,1,yLogit)
    yLogitNoise <- 1/(1+exp(-kappa*Ynoise[,simulation]))
    yBinNoise[,simulation] <- rbinom(samples,1,yLogitNoise)
    
  }
  
  #save file containing continuous and binary responses, together with simulation parameters
  save(simulations, samples, predictors, infoVars, SNRy, SNRx, kappa, coeffs, X, Xnoise, Y, Ynoise, yBin, yBinNoise,
       file = paste(fileName, "Plain.rda", sep="_"))


#   #debug#########################
#   #just check behaviour of simple regression on first simulation set
# myModellm <- lm(Ynoise[,1] ~ X[1,,])
# plot(myModellm$fitted.values ~ Ynoise[,1])
# plot(myModellm$fitted.values ~ Y[,1])

#   plot(Y[,1] ~ Ynoise[,1])  
# 
#   boxplot(Y[,1] ~ yBin[,1])
#   boxplot(Ynoise[,1] ~ yBinNoise[,1])

}

#################################################################

genSimulationFiles2D <- function(simulations=2, #500
                                 samples=200, 
                                 predictors=400, #900
                                 SNRy=10,
                                 SNRx=2,
                                 kappa=5){
  
  
  #### input parameters
  # simulations:  the number of simulation sets to generate
  # samples:      number of samples to generate in each set
  # predictors:   number of predictors per sample, this should be a perfect cube
  # SNRy:          signal-to-noise ratio, used when adding random noise
  # SNRx:         signal-to-noise ratio, used when adding random noise to the predictor variables
  # kappa:        coefficient in the logit exponential, it strongly regulates 
  #               class overlap in the two-class problem
  
  #set informative filename from list of parameters
  fileName <- paste(as.list(environment()), collapse="_")
  
  quadrDim <- round(predictors^(1/2))
  
  #need to check for "quadratic-ness"
  if (quadrDim != (predictors^(1/2))) {
    stop("Volume is not quadratic!")
  }
  
  #initialise variables
  X <- array(0, c(simulations, samples, predictors))
  Xnoise <- array(0, c(simulations, samples, predictors))
  Y <- matrix(0, nrow=samples, ncol=simulations)
  Ynoise <- matrix(0, nrow=samples, ncol=simulations)
  yBin <- matrix(0, nrow=samples, ncol=simulations)
  yBinNoise <- matrix(0, nrow=samples, ncol=simulations)
  
  #simply use a gaussian spherical symmetry for the coefficient values
  #co-centred with the cube and SD=1/4 the quadratic predictor
  coeffs <- array(0, dim=c(quadrDim, quadrDim))
  quadrCentre <- c(quadrDim/2+.5, quadrDim/2+.5)
  coeffSD <- as.integer(quadrDim/4)
  for (coefX in 1:quadrDim){
    for (coefY in 1:quadrDim){
      distance <- as.numeric(dist(rbind(quadrCentre, c(coefX, coefY)), method = "euclidean"))
      coeffs[coefX, coefY] <- dnorm(distance, mean=0, sd=coeffSD)
    }
  }
  
  
  
  #   #debug#########################
  #   #varY<-0
  #   #############################################
  #generate dataset for each simulation
  for (simulation in 1:simulations) {
    #new seeds for reproducibility
    set.seed(1234 + simulation*samples*4) #multiply by four because we sample four times per iteration
    #generate the samples
    for (predictor in 1:predictors) {
      curDim <- rnorm(samples, 0, 1)
      X[simulation, , predictor] <- curDim
      
      #also generate a noisy predictor
      noise <- rnorm(samples)
      #calculate the adj coefficient from variance of the predictor signal and the desired SNRx
      noiseCoeffX <- sqrt(var(curDim)/(SNRx * var(noise)))
      #generate the response with noise
      Xnoise[simulation, , predictor] <- X[simulation, , predictor] + noiseCoeffX*noise
    }
    #generate the response
    Y[,simulation] <- X[simulation, ,] %*% as.vector(coeffs)
    #generate some noise
    noise <- rnorm(samples)
    #calculate the adj coefficient from variance of the signal and the desired SNRy
    noiseCoeffY <- sqrt(var(Y[,simulation])/(SNRy * var(noise)))
    #generate the response with noise
    Ynoise[,simulation] <- Y[,simulation] + noiseCoeffY*noise
    #generate two-group response (kappa regulates the steepness and largely the overlap)
    yLogit <- 1/(1+exp(-kappa*Y[,simulation]))
    yBin[,simulation] <- rbinom(samples,1,yLogit)
    yLogitNoise <- 1/(1+exp(-kappa*Ynoise[,simulation]))
    yBinNoise[,simulation] <- rbinom(samples,1,yLogitNoise)
    
  }
  
  #save file containing continuous and binary responses, together with simulation parameters
  save(simulations, samples, predictors, SNRy, kappa, coeffs, X, Xnoise, Y, Ynoise, yBin, yBinNoise,
       file = paste(fileName, "2D.rda", sep="_"))
}
  #debug#########################
  #just check behaviour of simple regression on first simulation set
  #myModel <- lm(Ynoise[,1] ~ X[1,,])
  # plot(myModel$fitted.values ~ Ynoise[,1])
  # plot(myModel$fitted.values ~ Y[,1])
  # plot(Y[,1] ~ Ynoise[,1])  
  
  # boxplot(Y[,1] ~ yBin[,1])
  # boxplot(Ynoise[,1] ~ yBinNoise[,1])
  
#######################################################################

genSimulationFiles2Dfixedpred <- function(simulations=2,
                                          samples=200,
                                          predictors=400,
                                          SNRy=10,
                                          SNRx=2,
                                          kappa=5){
  
  
  #### input parameters
  # simulations:  the number of simulation sets to generate
  # samples:      number of samples to generate in each set
  # predictors:   number of predictors per sample, this should be a perfect cube
  # SNRy:          signal-to-noise ratio, used when adding random noise
  # SNRx:         signal-to-noise ratio, used when adding random noise to the predictor variables
  # kappa:        coefficient in the logit exponential, it strongly regulates 
  #               class overlap in the two-class problem
  
  #set informative filename from list of parameters
  fileName <- paste(as.list(environment()), collapse="_")
  
  quadrDim <- round(predictors^(1/2))
  
  #need to check for "quadratic-ness"
  if (quadrDim != (predictors^(1/2))) {
    stop("Volume is not quadratic!")
  }
  
  #initialise variables
  X <- array(0, c(simulations, samples, predictors))
  Xnoise <- array(0, c(simulations, samples, predictors))
  Y <- matrix(0, nrow=samples, ncol=simulations)
  Ynoise <- matrix(0, nrow=samples, ncol=simulations)
  yBin <- matrix(0, nrow=samples, ncol=simulations)
  yBinNoise <- matrix(0, nrow=samples, ncol=simulations)
  
  #simply use a gaussian spherical symmetry for the coefficient values
  #co-centred with the cube and SD=1/4 the quadratic predictor
  coeffs <- array(0, dim=c(quadrDim, quadrDim))
  quadrCentre <- c(quadrDim/2+.5, quadrDim/2+.5)
  coeffSD <- as.integer(quadrDim/4)
  for (coefX in 1:quadrDim){
    for (coefY in 1:quadrDim){
      distance <- as.numeric(dist(rbind(quadrCentre, c(coefX, coefY)), method = "euclidean"))
      coeffs[coefX, coefY] <- dnorm(distance, mean=0, sd=coeffSD)
    }
  }
  
  for (predictor in 1:predictors) {
    #new seeds for reproducibility
    set.seed(1234) 
    curDim <- rnorm(samples, 0, 1)
    #predictors X are fixed over simulations
    #fix predictors in first simulation and give the same values to the predictors of the other simulations
    X[1, , predictor] <- curDim 
    #also generate a noisy predictor
    noise <- rnorm(samples)
    #calculate the adj coefficient from variance of the predictor signal and the desired SNRx
    noiseCoeffX <- sqrt(var(curDim)/(SNRx * var(noise)))
    #generate the response with noise
    Xnoise[1, , predictor] <- X[1, , predictor] + noiseCoeffX*noise
  }
  
  #   #debug#########################
  #   #varY<-0
  #   #############################################
  #generate dataset for each simulation
  for (simulation in 1:simulations) {
    #same predictors for all simulations 
    X[simulation, ,] <- X[1, ,]
    #generate the response
    Y[,simulation] <- X[simulation, ,] %*% as.vector(coeffs)
    #generate some noise
    noise <- rnorm(samples)
    #calculate the adj coefficient from variance of the signal and the desired SNRy
    noiseCoeffY <- sqrt(var(Y[,simulation])/(SNRy * var(noise)))
    #generate the response with noise
    Ynoise[,simulation] <- Y[,simulation] + noiseCoeffY*noise
    #generate two-group response (kappa regulates the steepness and largely the overlap)
    yLogit <- 1/(1+exp(-kappa*Y[,simulation]))
    yBin[,simulation] <- rbinom(samples,1,yLogit)
    yLogitNoise <- 1/(1+exp(-kappa*Ynoise[,simulation]))
    yBinNoise[,simulation] <- rbinom(samples,1,yLogitNoise)
    
  }
  
  #save file containing continuous and binary responses, together with simulation parameters
  save(simulations, samples, predictors, SNRy, kappa, coeffs, X, Xnoise, Y, Ynoise, yBin, yBinNoise,
       file = paste(fileName, "2Dfix.rda", sep="_"))
  
  #debug#########################
  #just check behaviour of simple regression on first simulation set
  #myModel <- lm(Ynoise[,1] ~ X[1,,])
  # plot(myModel$fitted.values ~ Ynoise[,1])
  # plot(myModel$fitted.values ~ Y[,1])
  # plot(Y[,1] ~ Ynoise[,1])  
  
  # boxplot(Y[,1] ~ yBin[,1])
  # boxplot(Ynoise[,1] ~ yBinNoise[,1])
  
}

#######################################################################

genSimulationFilesVols <- function(simulations=2, #500
                                   samples=200,
                                   predictors=8000,
                                   SNRx=2,
                                   SNRy=10,
                                   kappa=5){
 
  #### input parameters
  # simulations:  the number of simulation sets to generate
  # samples:      number of samples to generate in each set
  # predictors:   number of predictors per sample, this should be a perfect cube
  # SNRy:          signal-to-noise ratio, used when adding random noise
  # kappa:        coefficient in the logit exponential, it strongly regulates 
  #               class overlap in the two-class problem
  
  #set informative filename from list of parameters
  fileName <- paste(as.list(environment()), collapse="_")
  
  cubeDim <- round(predictors^(1/3))
  #   print(cubeDim)
  #   print(as.integer(cubeDim))
  #   #need to check for "cubic-ness"
  #   if (cubeDim != as.integer(cubeDim)) {
  #     stop("Volume is not a cube!")
  #   }
  
  #initialise variables
  X <- array(0, c(simulations, samples, predictors))
  Xnoise <- array(0, c(simulations, samples, predictors))
  Y <- matrix(0, nrow=samples, ncol=simulations)
  Ynoise <- matrix(0, nrow=samples, ncol=simulations)
  yBin <- matrix(0, nrow=samples, ncol=simulations)
  yBinNoise <- matrix(0, nrow=samples, ncol=simulations)
  
  #simply use a gaussian spherical symmetry for the coefficient values
  #co-centred with the cube and SD=1/8 the cube dimension
  coeffs <- array(0, dim=c(cubeDim, cubeDim, cubeDim))
  volCentre <- c(cubeDim/2+.5, cubeDim/2+.5, cubeDim/2+.5)
  coeffSD <- as.integer(cubeDim/8)
  for (coefX in 1:cubeDim){
    for (coefY in 1:cubeDim){
      for (coefZ in 1:cubeDim){
        distance <- as.numeric(dist(rbind(volCentre, c(coefX, coefY, coefZ)), method = "euclidean"))
        coeffs[coefX, coefY, coefZ] <- dnorm(distance, mean=0, sd=coeffSD)
      }
    }
  }
  
  ################ FOR LATER:
  # 1. add second smaller and less strong sphere
  # 2. also code the 2D case (might be of interest to plain image recognition people)
  
  #   #debug#########################
  #   #varY<-0
  #   #############################################
  #generate dataset for each simulation
  for (simulation in 1:simulations) {
    #new seeds for reproducibility
    set.seed(1234 + simulation*samples*4) #multiply by four because we sample four times per iteration
    #generate the samples
    for (predictor in 1:predictors) {
      curDim <- rnorm(samples, 0, 1)
      X[simulation, , predictor] <- curDim
      ################ FOR LATER:      
      #HERE I COULD ALSO USE CHOLESKY DECOMPOSITION OF A REALISTIC CORRELATION MATRIX
      #TO ADD A SIMULATION OF CORRELATION BETWEEN VOXELS
      #also generate a noisy predictor
      noise <- rnorm(samples)
      #calculate the adj coefficient from variance of the predictor signal and the desired SNRx
      noiseCoeffX <- sqrt(var(curDim)/(SNRx * var(noise)))
      Xnoise[1, , predictor] <- X[1, , predictor] + noiseCoeffX*noise
    }
    #generate the response
    Y[,simulation] <- X[simulation, ,] %*% coeffs
    #generate some noise
    noise <- rnorm(samples)
    #calculate the adj coefficient from variance of the signal and the desired SNRy
    #varY<-varY+var(Y[,simulation])
    noiseCoeffY <- sqrt(var(Y[,simulation])/(SNRy * var(noise)))
    #generate the response with noise
    Ynoise[,simulation] <- Y[,simulation] + noiseCoeffY*noise
    #generate two-group response (kappa regulates the steepness and largely the overlap)
    yLogit <- 1/(1+exp(-kappa*Y[,simulation]))
    yBin[,simulation] <- rbinom(samples,1,yLogit)
    yLogitNoise <- 1/(1+exp(-kappa*Ynoise[,simulation]))
    yBinNoise[,simulation] <- rbinom(samples,1,yLogitNoise)
    
  }
  
  #save file containing continuous and binary responses, together with simulation parameters
  save(simulations, samples, predictors, SNRx, SNRy, kappa, coeffs, X, Xnoise, Y, Ynoise, yBin, yBinNoise,
       file = paste(fileName, "Vols.rda", sep="_"))
  #    save(simulations, samples, predictors, infoVars, SNRy, kappa, coeffs,
  #         file = paste(fileName, "VolsTESTING.rda", sep="_"))
  
  #debug#########################
  #just check behaviour of simple regression on first simulation set
  #   myModel <- lm(Ynoise[,1] ~ X[1,,])
  #   plot(myModel$fitted.values ~ Ynoise[,1])
  #   plot(myModel$fitted.values ~ Y[,1])
  #   plot(Y[,1] ~ Ynoise[,1])  
  # 
  #   boxplot(Y[,1] ~ yBin[,1])
  #   boxplot(Ynoise[,1] ~ yBinNoise[,1])
  
}



#######################################################################
#analyse simulation including variance


analyseSimulation <- function(dataFile="3_200_1000_50_10_2_5_Plain.rda_OUT_10_100_0.4_100_0.1.rda") {
  
  library(pROC)
  library(Rmisc)
  
  data <- dataFile #save data file for later functions appearing
  
  load(dataFile)
  #   if (grepl(pattern = "_Plain", x = dataFile)){
  #     fileType <- fileTypePlain
  #   } else {
  #     fileType <- fileTypeVols
  #   }
  
  if (!exists("redSteps")){
    redSteps <- dim(predictionVector)[2]
  }
  
  
  
  
  #initalize variables for absolute errors
  AbsErrors <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save absolute errors after each simulation
  sumAbsErrors <- rep(0, redSteps) #vector to calculate sum of absolut errors
  AbsErrorsCI97.5 <- rep(0, redSteps)#vectors to save confidence intervals
  AbsErrorsCI2.5 <- rep(0, redSteps)
  AbsErrorsNoiseX <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save absolute errors after each simulation
  sumAbsErrorsNoiseX <- rep(0, redSteps)#vector to calculate sum of absolut errors
  AbsErrorsNoiseXCI97.5 <- rep(0, redSteps)#vectors to save confidence intervals
  AbsErrorsNoiseXCI2.5 <- rep(0, redSteps)
  AbsErrorsNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save absolute errors after each simulation
  sumAbsErrorsNoiseY <- rep(0, redSteps)#vector to calculate sum of absolut errors
  AbsErrorsNoiseYCI97.5 <- rep(0, redSteps)#vectors to save confidence intervals
  AbsErrorsNoiseYCI2.5 <- rep(0, redSteps)
  AbsErrorsNoiseXNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save absolute errors after each simulation
  sumAbsErrorsNoiseXNoiseY <- rep(0, redSteps)#vector to calculate sum of absolut errors
  AbsErrorsNoiseXNoiseYCI97.5 <- rep(0, redSteps) #vectors to save confidence intervals
  AbsErrorsNoiseXNoiseYCI2.5 <- rep(0, redSteps)
  
  #initialize variables for AUC
  AUCsb <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save AUCafter each simulation
  sumAUC <- rep(0, redSteps) #vector for sum of AUC and 
  CIAUC2.5 <- rep(0, redSteps)#vectors to save confidence intervals
  CIAUC97.5 <- rep(0, redSteps)
  AUCNoiseX <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save AUCafter each simulation
  sumAUCNoiseX <- rep(0, redSteps) #vector for sum of AUC and 
  CIAUCNoiseX2.5 <- rep(0, redSteps)#vectors to save confidence intervals
  CIAUCNoiseX97.5 <- rep(0, redSteps)
  AUCNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save AUCafter each simulation
  sumAUCNoiseY <- rep(0, redSteps) #vector for sum of AUC and 
  CIAUCNoiseY2.5 <- rep(0, redSteps)#vectors to save confidence intervals
  CIAUCNoiseY97.5 <- rep(0, redSteps)
  AUCNoiseXNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save AUCafter each simulation
  sumAUCNoiseXNoiseY <- rep(0, redSteps) #vector for sum of AUC and 
  CIAUCNoiseXNoiseY2.5<- rep(0, redSteps)#vectors to save confidence intervals
  CIAUCNoiseXNoiseY97.5<- rep(0, redSteps)
  
  
  #initialize variables for variance
  variance <- matrix(0, nrow = redSteps, ncol = simulations)#matrix save variance after each reduction step
  sumvariance <- rep(0, redSteps) #vector for sum of variance
  CIvariance2.5 <- rep(0, redSteps) #vector to save confidence intervals of variance
  CIvariance97.5 <- rep(0, redSteps)
  varianceNoiseX <- matrix(0, nrow = redSteps, ncol = simulations)#matrix save variance after each reduction step
  sumvarianceNoiseX <- rep(0, redSteps) #vector for sum of variance
  CIvariance2.5NoiseX <- rep(0, redSteps) #vector to save confidence intervals of variance
  CIvariance97.5NoiseX <- rep(0, redSteps)
  varianceNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix save variance after each reduction step
  sumvarianceNoiseY <- rep(0, redSteps) #vector for sum of variance
  CIvariance2.5NoiseY <- rep(0, redSteps) #vector to save confidence intervals of variance
  CIvariance97.5NoiseY <- rep(0, redSteps)
  varianceNoiseXNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix save variance after each reduction step
  sumvarianceNoiseXNoiseY <- rep(0, redSteps) #vector for sum of variance
  CIvariance2.5NoiseXNoiseY <- rep(0, redSteps) #vector to save confidence intervals of variance
  CIvariance97.5NoiseXNoiseY <- rep(0, redSteps)
  
  
  for (simulation in 1:simulations) {
    
    print(simulation)
    ########REGRESSION
    AbsErrors[,simulation] <- colSums(abs(t(predictionVector[simulation,,])-Y[,simulation]))/simulations
    AbsErrorsNoiseX[,simulation] <- colSums(abs(t(predictionVectorNoise[simulation,,])-Y[,simulation]))/simulations
    AbsErrorsNoiseY[,simulation] <- colSums(abs(t(predictionVector[simulation,,])-Ynoise[,simulation]))/simulations
    AbsErrorsNoiseXNoiseY[,simulation] <- colSums(abs(t(predictionVectorNoise[simulation,,])-Ynoise[,simulation]))/simulations
    
    #######CLASSIFICATION (VIA REGRESSION)
    for (reduction in 1:redSteps){
      #print(reduction)#DEBUG
      #sumAUC[reduction] <- sumAUC[reduction] + auc(yBin[,simulation], predictionVectorClass[simulation, reduction, ])
      #calculate AUC
      AUCsb[reduction, simulation] <- auc(yBin[,simulation], predictionVectorClass[simulation, reduction, ])
      AUCNoiseX[reduction, simulation] <- auc(yBin[,simulation], predictionVectorClassNoise[simulation, reduction, ])
      AUCNoiseY[reduction, simulation] <- auc(yBinNoise[,simulation], predictionVectorClass[simulation, reduction, ])
      AUCNoiseXNoiseY[reduction, simulation] <- auc(yBinNoise[,simulation], predictionVectorClassNoise[simulation, reduction, ])
      
      #calculate variance 
      variance[reduction, simulation] <- var(predictionVector[simulation , reduction ,])
      varianceNoiseX[reduction, simulation] <- var(predictionVectorNoise[simulation , reduction ,])
      varianceNoiseY[reduction, simulation] <- var(predictionVector[simulation , reduction ,] - Ynoise[,simulation])
      varianceNoiseXNoiseY[reduction, simulation] <- var(predictionVectorNoise[simulation , reduction ,] - Ynoise[,simulation] )
      
      
      #print(variance2.5) #DEBUG
    }
    
  } #end simulation loop
  
  
  
  for (reduction in 1:redSteps){
    #calculate confidence intervals for sum of absolut errors
    AbsErrorsCI <- CI(AbsErrors[reduction,]) #take all reductions steps over the simulations and calculate CI
    AbsErrorsCI97.5[reduction] <- AbsErrorsCI[1] #upper CI
    AbsErrorsCI2.5[reduction] <- AbsErrorsCI[3] #upper CI
    AbsErrorsCI <- CI(AbsErrorsNoiseX[reduction,]) #take all reductions steps over the simulations and calculate CI
    AbsErrorsNoiseXCI97.5[reduction] <- AbsErrorsCI[1] #upper CI
    AbsErrorsNoiseXCI2.5[reduction] <- AbsErrorsCI[3] #upper CI
    AbsErrorsCI <- CI(AbsErrorsNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    AbsErrorsNoiseYCI97.5[reduction] <- AbsErrorsCI[1] #upper CI
    AbsErrorsNoiseYCI2.5[reduction] <- AbsErrorsCI[3] #upper CI
    AbsErrorsCI <- CI(AbsErrorsNoiseXNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    AbsErrorsNoiseXNoiseYCI97.5[reduction] <- AbsErrorsCI[1] #upper CI
    AbsErrorsNoiseXNoiseYCI2.5[reduction] <- AbsErrorsCI[3] #upper CI
    
    #calculate confidence intervals for AUC
    CIAUC <- CI(AUCsb[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIAUC97.5[reduction] <- CIAUC[1] #upper CI
    CIAUC2.5[reduction] <- CIAUC[3] #lower CI
    CIAUC <- CI(AUCNoiseX[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIAUCNoiseX97.5[reduction] <- CIAUC[1] #upper CI
    CIAUCNoiseX2.5[reduction] <- CIAUC[3] #lower CI
    CIAUC <- CI(AUCNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIAUCNoiseY97.5[reduction] <- CIAUC[1] #upper CI
    CIAUCNoiseY2.5[reduction] <- CIAUC[3] #lower CI
    CIAUC <- CI(AUCNoiseX[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIAUCNoiseXNoiseY97.5[reduction] <- CIAUC[1] #upper CI
    CIAUCNoiseXNoiseY2.5[reduction] <- CIAUC[3] #lower CI
    
    #calculate confidence intervals for variance
    CIvariance <- CI(variance[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIvariance2.5[reduction] <-  CIvariance[3] #lower band
    CIvariance97.5[reduction] <- CIvariance[1] #upper band
    CIvariance <-  CI(varianceNoiseX[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIvariance2.5NoiseX[reduction] <-  CIvariance[3] #lower band
    CIvariance97.5NoiseX[reduction] <- CIvariance[1] #upper band
    CIvariance <-  CI(varianceNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIvariance2.5NoiseY[reduction] <-  CIvariance[3] #lower band
    CIvariance97.5NoiseY[reduction] <- CIvariance[1] #upper band
    CIvariance <-  CI(varianceNoiseXNoiseY[reduction,])  #take all reductions steps over the simulations and calculate CI
    CIvariance2.5NoiseXNoiseY[reduction] <-  CIvariance[3] #lower band
    CIvariance97.5NoiseXNoiseY[reduction] <- CIvariance[1] #upper band
  }
  
  sumAUC <- rowSums(AUCsb)/simulations
  sumAUCNoiseX <- rowSums(AUCNoiseX)/simulations
  sumAUCNoiseY <- rowSums(AUCNoiseY)/simulations
  sumAUCNoiseXNoiseY <- rowSums(AUCNoiseXNoiseY)/simulations
  
  sumAbsErrors <- rowSums(AbsErrors)/simulations
  sumAbsErrorsNoiseX <- rowSums(AbsErrorsNoiseX)/simulations
  sumAbsErrorsNoiseY <- rowSums(AbsErrorsNoiseY)/simulations
  sumAbsErrorsNoiseXNoiseY <- rowSums(AbsErrorsNoiseXNoiseY)/simulations
  
  sumvariance <- rowSums(variance)/simulations
  sumvarianceNoiseX <- rowSums(varianceNoiseX)/simulations
  sumvarianceNoiseY <- rowSums(varianceNoiseY)/simulations
  sumvarianceNoiseXNoiseY <- rowSums(varianceNoiseXNoiseY)/simulations
  #save these values
  #####convergenceIteration(data)#after how many redSteps do analyse plots converge, for plots 
  ####FOR LATER SAVE outcome of convergence####
  
  redStepsconv <- convergenceIteration(data) #after how many reduction steps does boosting iteration converge
  
  save(redStepsconv, AbsErrors, sumAbsErrors, sumAbsErrorsNoiseX, sumAbsErrorsNoiseY, sumAbsErrorsNoiseXNoiseY, 
       sumAUC, sumAUCNoiseX, sumAUCNoiseY, sumAUCNoiseXNoiseY, sumvariance, sumvarianceNoiseX, 
       file = paste(dataFile, "EVAL.rda"))
  
  nred <- redStepsconv[1,] #convergence of Iteration for AbsError
  nredNoiseX <- redStepsconv[2,]        
  nredNoiseY  <- redStepsconv[3,]       
  nredNoiseXNoiseY  <- redStepsconv[4,] 
  meanred <- sum(nred)/simulations #mean of reduction Steps for Absolute Error over simulations
  meanredNoiseX <- sum(nredNoiseX)/simulations
  meanredNoiseY <- sum(nredNoiseY)/simulations
  meanredNoiseXNoiseY <- sum(nredNoiseXNoiseY)/simulations
  
  nredAUC <- redStepsconv[5,]#convergence of Iteration for AUC
  nredAUCNoiseX <- redStepsconv[6,]
  nredAUCNoiseY <- redStepsconv[7,]
  nredAUCNoiseXNoiseY <- redStepsconv[8,]
  meanredAUC <- sum(nredAUC)/simulations #mean of reduction Steps for AUC over simulations
  meanredAUCNoiseX <- sum(nredAUCNoiseX)/simulations
  meanredAUCNoiseY <- sum(nredAUCNoiseY)/simulations
  meanredAUCNoiseXNoiseY <- sum(nredAUCNoiseXNoiseY)/simulations
  
  #Plots
  pdf(paste(dataFile, "plots.pdf")) #save plot
  plot(sumAbsErrors, ylim=range(sumAbsErrors, AbsErrorsCI97.5, AbsErrorsCI2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of absolute Errors", main = "Sum of absolute Errors")
  lines(AbsErrorsCI97.5, col="red") #add CIs to plot
  lines(AbsErrorsCI2.5, col="red")
  abline(h = sumAbsErrors[meanred], col = "green") #horizontal line that shows convergence
  plot(sumAbsErrorsNoiseX, ylim=range(sumAbsErrorsNoiseX, AbsErrorsNoiseXCI97.5, AbsErrorsNoiseXCI2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of absolute Errors", main = "Sum of absolute Errors for noisy X")
  lines(AbsErrorsNoiseXCI97.5, col="red") #add CIs to plot
  lines(AbsErrorsNoiseXCI2.5, col="red")
  abline(h = sumAbsErrorsNoiseX[meanredNoiseX], col = "green") #horizontal line that shows convergence
  plot(sumAbsErrorsNoiseY ,ylim=range(sumAbsErrorsNoiseY, AbsErrorsNoiseYCI97.5, AbsErrorsNoiseYCI2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of absolute Errors", main = "Sum of absolute Errors for noisy Y")
  lines(AbsErrorsNoiseYCI97.5, col="red") #add CIs to plot
  lines(AbsErrorsNoiseYCI2.5, col="red")
  abline(h = sumAbsErrorsNoiseY[meanredNoiseY], col = "green") #horizontal line that shows convergence
  plot(sumAbsErrorsNoiseXNoiseY ,ylim=range(sumAbsErrorsNoiseXNoiseY, AbsErrorsNoiseXNoiseYCI97.5, AbsErrorsNoiseXNoiseYCI2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of absolute Errors", main = "Sum of absolute Errors for noisy X and noisy Y")
  lines(AbsErrorsNoiseXNoiseYCI97.5, col="red") #add CIs to plot
  lines(AbsErrorsNoiseXNoiseYCI2.5, col="red")
  abline(h = sumAbsErrorsNoiseXNoiseY[meanredNoiseXNoiseY], col = "green") #horizontal line that shows convergence
  
  # #relative values (to the basic boosting) for first and last simulation
  # plot(AbsErrors[,simulations]/AbsErrors[1, simulations], ylim=range(AbsErrors[,simulations]/AbsErrors[1, simulations], AbsErrors[,1]/AbsErrors[1, 1]), type="l", main = "Relative value of absolute error to basic boosting (first and last simulation)", xlab = "Iteration", ylab = "Relative value")
  # lines(AbsErrors[,1]/AbsErrors[1, 1], type="l", col = "blue")
  # abline(h = 1, col = "red") #if line lays above the line the Abs Error of basic boosting is higher
  # plot(AbsErrorsNoiseX[,simulations]/AbsErrorsNoiseX[1, simulations], ylim=range(AbsErrorsNoiseX[,simulations]/AbsErrorsNoiseX[1, simulations], AbsErrorsNoiseX[,1]/AbsErrorsNoiseX[1, 1]), type="l", main = "Relative value of absolute error to basic boosting for noisy X (first and last simulation)", xlab = "Iteration", ylab = "Relative value")
  # lines(AbsErrorsNoiseX[,1]/AbsErrorsNoiseX[1, 1], type="l", col = "blue")
  # abline(h = 1, col = "red") #if line lays above the line the Abs Error of basic boosting is higher
  # plot(AbsErrorsNoiseY[,simulations]/AbsErrorsNoiseY[1, simulations], ylim=range(AbsErrorsNoiseY[,simulations]/AbsErrorsNoiseY[1, simulations], AbsErrorsNoiseY[,1]/AbsErrorsNoiseY[1, 1]), type="l", main = "Relative value of absolute error to basic boosting for noisy Y (first and last simulation)", xlab = "Iteration", ylab = "Relative value")
  # lines(AbsErrorsNoiseY[,1]/AbsErrorsNoiseY[1, 1], type="l", col = "blue")
  # abline(h = 1, col = "red") #if line lays above the line the Abs Error of basic boosting is higher
  # plot(AbsErrorsNoiseXNoiseY[,simulations]/AbsErrorsNoiseXNoiseY[1, simulations], ylim=range(AbsErrorsNoiseXNoiseY[,simulations]/AbsErrorsNoiseXNoiseY[1, simulations], AbsErrorsNoiseXNoiseY[,1]/AbsErrorsNoiseXNoiseY[1, 1]), type="l", main = "Relative value of absolute error to basic boosting for noisy X and noisy Y (first and last simulation)", xlab = "Iteration", ylab = "Relative value")
  # lines(AbsErrorsNoiseXNoiseY[,1]/AbsErrorsNoiseXNoiseY[1, 1], type="l", col = "blue")
  # abline(h = 1, col = "red") #if line lays above the line the Abs Error of basic boosting is higher
  # #classification AUC
  
 
  #relative values (to the basic boosting)  
  plot(sumAbsErrors/sumAbsErrors[1], main = "Relative value of absolute error to basic boosting ")
  plot(sumAbsErrorsNoiseX/sumAbsErrorsNoiseX[1], main = "Relative value of absolute error to basic boosting ")
  plot(sumAbsErrorsNoiseY/sumAbsErrorsNoiseY[1], main = "Relative value of absolute error to basic boosting ")
  plot(sumAbsErrorsNoiseXNoiseY/sumAbsErrorsNoiseXNoiseY[1], main = "Relative value of absolute error to basic boosting ")
  
  
  plot(sumAUC, ylim=range(sumAUC, CIAUC97.5, CIAUC2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of AUC", main = "Sum aof AUC")
  lines(CIAUC2.5, col="red") #add CIs to plot
  lines(CIAUC97.5, col="red")
  abline(h = sumAUC[meanredAUC], col = "green") #horizontal line that shows convergence
  plot(sumAUCNoiseX, ylim=range(sumAUCNoiseX, CIAUCNoiseX97.5, CIAUCNoiseX2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of AUC", main = "Sum aof AUC for noisy X")
  lines(CIAUCNoiseX2.5, col="red") #add CIs to plot
  lines(CIAUCNoiseX97.5, col="red")
  abline(h = sumAUCNoiseX[meanredAUCNoiseX], col = "green") #horizontal line that shows convergence
  plot(sumAUCNoiseY, ylim=range(sumAUCNoiseY, CIAUCNoiseY97.5, CIAUCNoiseY2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of AUC", main = "Sum aof AUC for noisy Y")
  lines(CIAUCNoiseY2.5, col="red") #add CIs to plot
  lines(CIAUCNoiseY97.5, col="red")
  abline(h = sumAUCNoiseY[meanredAUCNoiseY], col = "green") #horizontal line that shows convergence
  plot(sumAUCNoiseXNoiseY, ylim=range(sumAUCNoiseY, CIAUCNoiseXNoiseY97.5, CIAUCNoiseXNoiseY2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of AUC", main = "Sum aof AUC for noisy X and noisy Y")
  lines(CIAUCNoiseXNoiseY2.5, col="red") #add CIs to plot
  lines(CIAUCNoiseXNoiseY97.5, col="red")
  abline(h = sumAUCNoiseXNoiseY[meanredAUCNoiseXNoiseY], col = "green") #horizontal line that shows convergence
  #variance
  
  plot(sumvariance, ylim=range(sumvariance, CIvariance2.5, CIvariance97.5), col='black', type = "l", xlab = "Iteration", ylab = "Sum of Variance", main = "Sum of Variance")
  lines(CIvariance97.5, col="red") #add CIs to plot
  lines(CIvariance2.5, col="red")
  plot(sumvarianceNoiseX, ylim=range(sumvarianceNoiseX, CIvariance2.5NoiseX, CIvariance97.5NoiseX), col='black', type = "l", xlab = "Iteration", ylab = "Sum of Variance", main = "Sum of Variance for noisy X")
  lines(CIvariance2.5NoiseX, col="red") #add CIs to plot
  lines(CIvariance97.5NoiseX, col="red")
  plot(sumvarianceNoiseY, ylim=range(sumvarianceNoiseY, CIvariance2.5NoiseY, CIvariance97.5NoiseY), col='black', type = "l", xlab = "Iteration", ylab = "Sum of Variance", main = "Sum of Variance for noisy Y")
  lines(CIvariance2.5NoiseY, col="red") #add CIs to plot
  lines(CIvariance97.5NoiseY, col="red")
  plot(sumvarianceNoiseXNoiseY, ylim=range(sumvarianceNoiseXNoiseY, CIvariance2.5NoiseXNoiseY, CIvariance97.5NoiseXNoiseY), col='black', type = "l", xlab = "Iteration", ylab = "Sum of Variance" , main = "Sum of Variance for noisy Y and noisy X")
  lines(CIvariance2.5NoiseXNoiseY, col="red") #add CIs to plot
  lines(CIvariance97.5NoiseXNoiseY, col="red")
  dev.off() #end of saving plot
  #par(mfrow=c(1,1))
  
  
  
}

########################################################################


#######################################################################

convergenceIteration <- function(dataFile="3_200_1000_50_10_2_5_Plain.rda_OUT_10_1000_0.4_100_0.1.rda"){
  ###FOR LATER -> give back maximum value for no convergence
  ###FOR LATER -> same for noisy X/Y
  library(pracma)
  library(pROC)
  load(dataFile)
  
  
  grad <- 0 #variable for gradient of curve
  nred <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredNoiseX <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredNoiseY <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredNoiseXNoiseY <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredAUC <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredAUCNoiseX <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredAUCNoiseY <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredAUCNoiseXNoiseY <- rep(0, simulations) #variable to count how many reduction iterations were done
  
  AbsError <- matrix(0, nrow = redSteps, ncol = simulations) #AbsError to stop reduction iteration
  AbsErrorNoiseX <- matrix(0, nrow = redSteps, ncol = simulations) #AbsError to stop reduction iteration
  AbsErrorNoiseY <- matrix(0, nrow = redSteps, ncol = simulations) #AbsError to stop reduction iteration
  AbsErrorNoiseXNoiseY <- matrix(0, nrow = redSteps, ncol = simulations) #AbsError to stop reduction iteration
  sumAUC <- matrix(0, nrow=redSteps, ncol = simulations) #sumAUC to stop reduction iteration
  sumAUCNoiseX <- matrix(0, nrow=redSteps, ncol = simulations) #sumAUC to stop reduction iteration
  sumAUCNoiseY <- matrix(0, nrow=redSteps, ncol = simulations) #sumAUC to stop reduction iteration
  sumAUCNoiseXNoiseY <- matrix(0, nrow=redSteps, ncol = simulations) #sumAUC to stop reduction iteration
  
  
  for (simulation in 1:simulations) {
    print(simulation)
    stopred <- 0
    #calculate number of iteration when Absolute Error converges
    for (reduction in 1:redSteps){ #for sumabsErrors with no noisy X or Y
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      AbsError[reduction, simulation] <- sum(abs(t(predictionVector[simulation,reduction,])-Y[,simulation]))
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(AbsError[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.07){ 
          nred[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }
    }#end iteration reduction
    
    stopred <- 0
    for (reduction in 1:redSteps){ #for sumabsErrors with no noisy X or Y
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      AbsErrorNoiseX[reduction, simulation] <- sum(abs(t(predictionVectorNoise[simulation,reduction,])-Y[,simulation]))
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(AbsErrorNoiseX[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.07){ 
          nredNoiseX[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }
    }#end iteration reduction
    
    stopred <- 0
    for (reduction in 1:redSteps){ #for sumabsErrors with no noisy X or Y
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      AbsErrorNoiseY[reduction ,simulation] <- sum(abs(t(predictionVector[simulation, reduction,])-Ynoise[,simulation]))
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(AbsErrorNoiseY[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.07){ 
          nredNoiseY[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      } 
    }#end iteration reduction
    
    stopred <- 0
    for (reduction in 1:redSteps){ #for sumabsErrors with no noisy X or Y
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      AbsErrorNoiseXNoiseY[reduction ,simulation] <- sum(abs(t(predictionVectorNoise[simulation,,])-Ynoise[,simulation]))
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(AbsErrorNoiseXNoiseY[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.01){ 
          nredNoiseXNoiseY[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }
    }#end iteration reduction
    
    
    #calculate number of iteration when AUC converges
    stopred <- 0 
    for (reduction in 1:redSteps){
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      sumAUC[reduction, simulation] <- auc(yBin[,simulation], predictionVectorClass[simulation, reduction, ])
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(sumAUC[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.0001){
          nredAUC[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }  
    }#end reduction
  
  stopred <- 0 
  for (reduction in 1:redSteps){
    if(stopred == 1 ){  #if grad was small enough leave reduction iteration
      break
    }
    sumAUCNoiseX[reduction, simulation] <- auc(yBin[,simulation], predictionVectorClassNoise[simulation, reduction, ])
    if(as.integer(reduction/10)==(reduction/10)){
      #print(reduction)
      grad <- mean(gradient(sumAUCNoiseX[(reduction-9):reduction])) #mean of gradient of error curve
      #print(paste("grad:", grad)) #DEBUG
      if(abs(grad) < 0.0001){
        nredAUCNoiseX[simulation] <- reduction #save for checking how many reduction iterations were made
        # print("stop") #DEBUG
        stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
      }
    } 
  }#end reduction
  
  stopred <- 0 
  for (reduction in 1:redSteps){
    if(stopred == 1 ){  #if grad was small enough leave reduction iteration
      break
    }
    sumAUCNoiseY[reduction, simulation] <- auc(yBinNoise[,simulation], predictionVectorClass[simulation, reduction, ])
    if(as.integer(reduction/10)==(reduction/10)){
      #print(reduction)
      grad <- mean(gradient(sumAUCNoiseY[(reduction-9):reduction])) #mean of gradient of error curve
      #print(paste("grad:", grad)) #DEBUG
      if(abs(grad) < 0.0001){
        nredAUCNoiseY[simulation] <- reduction #save for checking how many reduction iterations were made
        # print("stop") #DEBUG
        stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
      }
    }
  }#end reduction
  
  stopred <- 0 
  for (reduction in 1:redSteps){
    if(stopred == 1 ){  #if grad was small enough leave reduction iteration
      break
    }
    sumAUCNoiseXNoiseY[reduction, simulation] <- auc(yBinNoise[,simulation], predictionVectorClassNoise[simulation, reduction, ])
    if(as.integer(reduction/10)==(reduction/10)){
      #print(reduction)
      grad <- mean(gradient(sumAUCNoiseXNoiseY[(reduction-9):reduction])) #mean of gradient of error curve
      #print(paste("grad:", grad)) #DEBUG
      if(abs(grad) < 0.0001){
        nredAUCNoiseXNoiseY[simulation] <- reduction #save for checking how many reduction iterations were made
        # print("stop") #DEBUG
        stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
      }
    }
  }#end reduction


  }#end simulation
  
  rbind(nred, nredNoiseX, nredNoiseY, nredNoiseXNoiseY, nredAUC, nredAUCNoiseX, nredAUCNoiseY, nredAUCNoiseXNoiseY) #save values as matrix to use them analyse function
  }#end function

#########################################
###FOR LATER: put function into plots
calculateCI <- function(input){
   }