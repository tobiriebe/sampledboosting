install.packages("devtools")
library(devtools)
install_github("mllg/batchtools")
library(batchtools)
#library(BatchExperiments)
library(ROCR)
library(mboost)
library(caret)
library(pROC)



setwd("/naslx/projects/ua341/di49suy")
#name ID an file direction
reg = makeRegistry(file.dir = "mytest")


################################################################################
####################################TREE########################################
################################################################################

sampledboosting.wrapper <- function(file, sampleRatio ){
  
  
  load(file)
  #this is not implementing inner folds (CV-like), just reduction over the whole training set of the current fold 
  #also weight is not implemented
  
  if (length(commandArgs(trailingOnly = TRUE))>0) {
    #####read the sample number to test upon from the call
    args <- commandArgs(trailingOnly = TRUE)
    
    nOuterFolds <- as.numeric(args[1])  #integer number of outer folds to use
    redSteps <- as.numeric(args[2])     #integer number of iterations to run for the sampled boosting
    sampleRatio <- as.numeric(args[3])  #ratio of voxels TO REMOVE
    fixedMstop <- as.numeric(args[4])   #mstop to be used when not implementing cv early stopping
    fixedNu <- args[5]                  #shrinkage coefficient for boosting
    file <- args[6]                 #file to load
    localRun <- FALSE
  } else {
    ###### local testing
    nOuterFolds <- 10 #number of folds
    redSteps <- 100 #number of steps to run for the reduction in voxels
    sampleRatio = sampleRatio #ratio of voxels to REMOVE
    fixedMstop <- 100
    fixedNu <- 0.1
    file <- file
    localRun <- TRUE
  }
  file <- load(file)
  ###########file types 
  fileTypePlain <- 1
  fileTypeVols <- 2
  
  #label for output
  labelOut <- ""
  addLabelOut <- ""
  
  #############load the data file with data and parameters#############
  #score <- read.table("CYP2D6ScoreTRAINING.txt")
  #load(file)
  #if (grepl(pattern = "_Plain", x = file)){
  # fileType <- fileTypePlain
  #} else {
  # fileType <- fileTypeVols
  #}
  n <- samples  #number of cases (samples)
  # simulations is already called 'simulations' in the file
  nVariables <- predictors 
  dataX <- X # X[simulation, sample, x]
  dataXnoise <- Xnoise # Xnoise[simulation, sample, x]
  #make a copy of the matrix, to keep for the final test (untouched predictors)
  originalX <- X
  originalXnoise <- Xnoise
  simulations <- simulations
  #### other variables that get loaded:
  #response variables:
  # Y
  # Ynoise
  # yBin
  # yBinNoise
  #...
  #####################################################################
  
  ########### output variables #####################
  # imgModelList <- as.array(0,c(simulations, redSteps, nVariables))  #this keeps the list of models that will be summed up to produce an image
  # imgModelListAdj <- as.array(0,c(simulations, redSteps, nVariables))  #same, with adjusted values (number of folds and reductions)
  predModelList <- array(0,c(simulations, nOuterFolds, redSteps, nVariables))  #this keeps the list of unaltered models to produce intermediate predictions
  #for monitoring
  predModelListNoise <- array(0,c(simulations, nOuterFolds, redSteps, nVariables))  
  predModelListClass <- array(0,c(simulations, nOuterFolds, redSteps, nVariables))  #this keeps the list of unaltered models to produce intermediate predictions
  #for monitoring
  predModelListClassNoise <- array(0,c(simulations, nOuterFolds, redSteps, nVariables))
  
  # predModelListAdj <- as.array(0,c(simulations, redSteps, nVariables))  #the same, but adjusted (number of folds)
  offsetFinal <- array(0,c(simulations, nOuterFolds, redSteps))
  offsetFinalNoise <- array(0,c(simulations, nOuterFolds, redSteps))
  offsetFinalClass <- array(0,c(simulations, nOuterFolds, redSteps))
  offsetFinalClassNoise <- array(0,c(simulations, nOuterFolds, redSteps))
  
  removedVoxels <- vector("list", redSteps)  #this keeps track of the voxels set to zero at each reduction
  removedVoxelsNoise <- vector("list", redSteps)
  
  #### variables used to plot trend for each fold as reductions increase for the classification case  
  # AUC0 <- matrix(0, nrow = redSteps)   #these are the AUC trends for each fold
  # modelWeight <- rep(0, redSteps)
  # monitorWeight <- rep(0, redSteps)
  predictionVector <- array(0, c(simulations, redSteps, n))
  predictionVectorNoise <- array(0, c(simulations, redSteps, n))
  predictionVectorClass <- array(0, c(simulations, redSteps, n))
  predictionVectorClassNoise <- array(0, c(simulations, redSteps, n))
  
  coefNumbers <- rep(0, redSteps) #the number of coefficient in the final linear model at each reduction (for monitoring)
  coefNumbersNoise <- rep(0, redSteps)
  
  set.seed(1234) #for reproducibility
  
  for (simulation in 1:simulations){
    
    print(paste("Simulation n.", simulation))   #DEBUG
    
    X <- as.matrix(dataX[simulation, , ])
    Xnoise <- as.matrix(dataXnoise[simulation, , ])
    #training is always performed on noisy response variables
    y <- Ynoise[, simulation]
    yClass <- yBinNoise[, simulation]
    
    #create Outer folds (list of indices, one list per fold, which specify the test sets)
    indexOuterList <- createFolds(Ynoise[, simulation], nOuterFolds) #sample(n)  
    
    #   #iteration to (progressively) eliminate selected voxels to produce images
    #   for (reduction in 1:redSteps){
    #make training/test sessions for this fold
    for (kkk in 1:nOuterFolds){
      
      #print(paste("Fold n.", kkk))   #DEBUG
      
      indvecOuter <- indexOuterList[[kkk]]  #indices for the test set of this Outer fold
      trainXOuter <- X[-indvecOuter,]
      trainXnoiseOuter <- Xnoise[-indvecOuter,]
      trainyOuter <- y[-indvecOuter]
      trainyOuterClass <- yClass[-indvecOuter]
      testXOuter <- X[indvecOuter,]
      testXnoiseOuter <- Xnoise[indvecOuter,]
      testyOuter <- y[indvecOuter]
      testyOuterClass <- yClass[indvecOuter]
      
      modelList <- matrix(0,nrow = nOuterFolds, ncol = nVariables) #this is the matrix with model coefficients
      modelListNoise <- matrix(0,nrow = nOuterFolds, ncol = nVariables) 
      modelListClass <- matrix(0,nrow = nOuterFolds, ncol = nVariables) #this is the matrix with model coefficients
      modelListClassNoise <- matrix(0,nrow = nOuterFolds, ncol = nVariables) 
      #colnames(modelList) <- colnames(X)
      #colnames(modelListClass) <- colnames(X)
      
      #start with all predictors
      trainX <- trainXOuter
      trainXnoise <- trainXnoiseOuter
      trainXClass <- trainXOuter
      trainXClassNoise <- trainXnoiseOuter
      
      #iteration to (progressively) eliminate selected voxels to produce images
      for (reduction in 1:redSteps){
        #     #make training sessions for this fold
        #     for (kkk in 1:nOuterFolds){
        
        #print(paste("Reduction n.", reduction)) #DEBUG
        
        #if implementing inner folds, here we need to reload training set from original
        #as the design matrix has been modified by the reduction inner cycles
        #    indvec <- testIndex 
        
        #       trainy <- trainyOuter
        #       trainyClass <- yClass
        #       testy <- y
        #       testyClass <- yClass
        
        #       predTest3 <- rep(0, n)   #variables to store prediction values
        #       predTest <- rep(0, length(testy))
        #       predTestTempReduced <- rep(0, length(testy))
        #       predTestAdj <- rep(0, n)   
        #       meanWeight <- 0 #mean weight of the models for the current reduction, to be used fo weighting the total model
        
        
        #make boosting model with FIXED mstop 
        model <- glmboost(trainX, trainyOuter,  center = FALSE, control = boost_control(mstop = fixedMstop, trace = FALSE, nu = fixedNu))
        modelNoisy <- glmboost(trainXnoise, trainyOuter, center = FALSE, control = boost_control(mstop = fixedMstop, trace = FALSE, nu = fixedNu))
        modelClass <- glmboost(trainXClass, trainyOuterClass, center = FALSE, control = boost_control(mstop = fixedMstop, trace = FALSE, nu = fixedNu))
        modelClassNoisy <- glmboost(trainXClassNoise, trainyOuterClass, center = FALSE, control = boost_control(mstop = fixedMstop, trace = FALSE, nu = fixedNu))
        
        #CI <- confint(model)
        #extract coefficients and offset of models
        cfs <- coef(model)
        cfsNoise <- coef(modelNoisy)
        #offset <- offset + attr(coef(model), "offset")/nOuterFolds
        offsetFinal[simulation, kkk, reduction] <- attr(coef(model), "offset")
        offsetFinalNoise[simulation, kkk, reduction] <- attr(coef(modelNoisy), "offset")
        modelList[kkk,as.numeric(substr(names(cfs), 2, 100))] <- cfs #* curWeight ############# I DON'T THINK KKK IS NEEDED HERE ####################
        modelListNoise[kkk,as.numeric(substr(names(cfsNoise), 2, 100))] <- cfsNoise
        predModelList[simulation, kkk, reduction, ] <- modelList[kkk,] ############# I DON'T THINK KKK IS NEEDED HERE ####################
        predModelListNoise[simulation, kkk, reduction, ] <- modelListNoise[kkk,]
        
        cfsClass <- coef(modelClass)
        cfsClassNoise <- coef(modelClassNoisy)
        #offsetClass <- offsetClass + attr(coef(modelClass), "offset")/nOuterFolds
        offsetFinalClass[simulation, kkk, reduction] <- attr(coef(modelClass), "offset")
        offsetFinalClassNoise[simulation, kkk, reduction] <- attr(coef(modelClassNoisy), "offset")
        modelListClass[kkk,as.numeric(substr(names(cfsClass), 2, 100))] <- cfsClass #* curWeight #history of coefficient values as selected at each reduction
        modelListClassNoise[kkk,as.numeric(substr(names(cfsClassNoise), 2, 100))] <- cfsClassNoise #* curWeight #history of coefficient values as selected at each reduction
        predModelListClass[simulation, kkk, reduction, ] <- modelListClass[kkk,] ############# I DON'T THINK KKK IS NEEDED HERE ##################
        predModelListClassNoise[simulation, kkk, reduction, ] <- modelListClassNoise[kkk,] ############# I DON'T THINK KKK IS NEEDED HERE ##################
        
        #build the final linear models for the current reduction
        lTempModel <- rep(0, nVariables)
        lTempOffset <- 0
        lTempModelNoise <- rep(0, nVariables)
        lTempOffsetNoise <- 0
        lTempModelClass <- rep(0, nVariables)
        lTempOffsetClass <- 0
        lTempModelClassNoise <- rep(0, nVariables)
        lTempOffsetClassNoise <- 0
        for (tempReduction in 1:reduction) {
          #basic model
          tempSingleModel <- predModelList[simulation, kkk, tempReduction, ]
          lTempModel <- lTempModel + tempSingleModel/reduction
          lTempOffset <- lTempOffset + offsetFinal[simulation, kkk, tempReduction]/reduction
          #noisy model
          tempSingleModelNoise <- predModelListNoise[simulation, kkk, tempReduction, ]
          lTempModelNoise <- lTempModelNoise + tempSingleModelNoise/reduction
          lTempOffsetNoise <- lTempOffsetNoise + offsetFinalNoise[simulation, kkk, tempReduction]/reduction
          #classification model
          tempSingleModelClass <- predModelListClass[simulation, kkk, tempReduction, ]
          lTempModelClass <- lTempModelClass + tempSingleModelClass/reduction
          lTempOffsetClass <- lTempOffsetClass + offsetFinalClass[simulation, kkk, tempReduction]/reduction
          #classification noisy model
          tempSingleModelClassNoise <- predModelListClassNoise[simulation, kkk, tempReduction, ]
          lTempModelClassNoise <- lTempModelClassNoise + tempSingleModelClassNoise/reduction
          lTempOffsetClassNoise <- lTempOffsetClassNoise + offsetFinalClassNoise[simulation, kkk, tempReduction]/reduction
        }
        
        # make the temporary predictions from the linear model based on coefficients, to later evaluate how final performance changes at each reduction step
        #lTempModel <- colSums(imgModelList)
        testXtemp <- matrix(originalX[simulation, indvecOuter, ], nrow = length(indvecOuter))  #use the original, unreduced design matrix!
        testXtempNoise <- matrix(originalXnoise[simulation, indvecOuter, ], nrow = length(indvecOuter))  
        predictionVector[simulation, reduction, indvecOuter] <- testXtemp %*% lTempModel + lTempOffset
        predictionVectorNoise[simulation, reduction, indvecOuter] <- testXtempNoise %*% lTempModelNoise + lTempOffsetNoise
        predictionVectorClass[simulation, reduction, indvecOuter] <- testXtemp %*% lTempModelClass + lTempOffsetClass
        predictionVectorClassNoise[simulation, reduction, indvecOuter] <- testXtempNoise %*% lTempModelClassNoise + lTempOffsetClassNoise
        
        #       #store the number of non-zero coefficients so far (for monitoring)
        #       coefNumbers[reduction] <- length(which(lTempModel!=0))
        
        #using some random coef which are not 0 in the current model, find predictors to be temporarily excluded
        voxelNum = round(length(which(colSums(modelList)!=0)) * sampleRatio)  #calculate number of predictors to remove
        voxelNumNoise = round(length(which(colSums(modelListNoise)!=0)) * sampleRatio)
        mostSelInd <- sample(which(colSums(modelList)!=0), voxelNum) #calculate random indices
        mostSelIndNoise <- sample(which(colSums(modelListNoise)!=0), voxelNumNoise) #calculate random indices
        voxelNumClass = round(length(which(colSums(modelListClass)!=0)) * sampleRatio)  #calculate number of predictors to remove
        mostSelIndClass <- sample(which(colSums(modelListClass)!=0), voxelNumClass) #calculate random indices
        voxelNumClassNoise = round(length(which(colSums(modelListClassNoise)!=0)) * sampleRatio)  #calculate number of predictors to remove
        mostSelIndClassNoise <- sample(which(colSums(modelListClassNoise)!=0), voxelNumClassNoise) #calculate random indices
        
        #print(paste("MostSel",mostSelInd))            ####### DEBUG ###########
        #print(paste("MostSelClass",mostSelIndClass))  ####### DEBUG ###########
        
        #start with all predictors
        trainX <- trainXOuter
        trainXnoise <- trainXnoiseOuter
        trainXClass <- trainXOuter
        trainXClassNoise <- trainXnoiseOuter
        #remove the selected ones
        for (indRow in 1:nrow(trainXOuter)) {
          #set the selected predictors to zero for the next iteration
          trainX[indRow,mostSelInd] <- 0
          trainXnoise[indRow,mostSelIndNoise] <- 0
          trainXClass[indRow,mostSelIndClass] <- 0
          trainXClassNoise[indRow,mostSelIndClassNoise] <- 0
        }
        
      } #end of reduction block
      
    } #end of fold block
    
  } #end of simulation block
  
  
  #make a name for the output file, using input file name and parameters passed
  #outLabel <- paste(file,"_OUT_", nOuterFolds, "_", redSteps, "_", sampleRatio, "_", fixedMstop, "_", fixedNu, ".rda", sep = "")
  #save all variables from input file, parameters and output
  list(y = y, Ynoise = Ynoise, yBin = yBin,yBinNoise = yBinNoise, 
       originalX = originalX, originalXnoise = originalXnoise, coeffs = coeffs, 
       predictors = predictors, kappa = kappa, samples = samples, 
       simulations = simulations, nOuterFolds = nOuterFolds, redSteps = redSteps, 
       sampleRatio = sampleRatio, fixedMstop = fixedMstop, fixedNu = fixedNu,
       offsetFinal = offsetFinal, predModelList = predModelList, offsetFinalClass = offsetFinalClass,
       predModelListClass = predModelListClass, predictionVector = predictionVector, 
       predictionVectorClass = predictionVectorClass,offsetFinalNoise = offsetFinalNoise, 
       predModelListNoise = predModelListNoise, predictionVectorNoise = predictionVectorNoise, 
       offsetFinalClassNoise = offsetFinalClassNoise, predModelListClassNoise = predModelListClassNoise,
       predictionVectorClassNoise = predictionVectorClassNoise)
  
  
} #end function

batchMap(fun = sampledboosting.wrapper, file = c( "1000_200_1000_50_10_2_5_Plain.rda", "1000_200_100_50_10_2_5_Plain.rda"), more.args = list(sampleRatio = c(0.1, 0.5)) )



################################################
################################################
################################################
################################################







# Submit the jobs to the batch system
submitJobs(reg, resources = list(walltime = 60L*60L*1L, memory = 1000L))

getStatus()


