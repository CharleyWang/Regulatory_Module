### Plot the output genes on a two dimentional plot
### the Plot is seperated by two lines '十'. seperation point is given by cutoff_forPirorClass
scatterDFplot <- function(df, cutOff_forPirorClass = 0) 
{
  #################plot
  plot(x=NA,y=NA,xlim = c(min(df$corDown),max(df$corDown)),
       ylim = c(min(df$corUP), max(df$corUP)), 
       #yaxt="n",xaxt = "n",bty = 'n',
       ylab = 'corUP', main = as.character(df$chr[1]),
       xlab = 'corDown',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
  abline(v = cutOff_forPirorClass, lty = 2)
  abline(h = cutOff_forPirorClass, lty = 2)
  
  points(df$corDown[which(df$state == 1)], df$corUP[which(df$state == 1)], col = 'navy', cex = 0.6)
  points(df$corDown[which(df$state == 4)], df$corUP[which(df$state == 4)], col = 'red', cex = 0.6)
  points(df$corDown[which(df$state == 2)], df$corUP[which(df$state == 2)], col = 'forestgreen', cex = 0.6)
  points(df$corDown[which(df$state == 3)], df$corUP[which(df$state == 3)], col = 'gold', cex = 0.6)
}

##Calculating median for each classes outputting from Markov model
##Order all median value and assign different class with specific numbaer from 1-4
##'1' represent downward correlation
##'2' represent no corelation with both side ## '3' represent upward rorrelation 
##'4' represent correlation with both side
##Require function 'categoryAssian'
categoryChangeDF<- function(df, cutOff_forPirorClass = 0) 
{
  medianUPcor1    <- median(df$corUP[which(df$state == 1)])
  medianDowncor1  <- median(df$corDown[which(df$state == 1)])
  medianUPcor2    <- median(df$corUP[which(df$state == 2)])
  medianDowncor2  <- median(df$corDown[which(df$state == 2)])
  medianUPcor3    <- median(df$corUP[which(df$state == 3)])
  medianDowncor3  <- median(df$corDown[which(df$state == 3)])
  medianUPcor4    <- median(df$corUP[which(df$state == 4)])
  medianDowncor4  <- median(df$corDown[which(df$state == 4)])
  newState <- rep(0, nrow(df))
  
  cat1 <- categoryAssian(medianDowncor1 , medianUPcor1, cutOff_forPirorClass)
  cat2 <- categoryAssian(medianDowncor2 , medianUPcor2, cutOff_forPirorClass)
  cat3 <- categoryAssian(medianDowncor3 , medianUPcor3, cutOff_forPirorClass)
  cat4 <- categoryAssian(medianDowncor4 , medianUPcor4, cutOff_forPirorClass)
  
  newState[which(df$state == 1)] <- cat1
  newState[which(df$state == 2)] <- cat2
  newState[which(df$state == 3)] <- cat3
  newState[which(df$state == 4)] <- cat4
  
  categoryOrder <- rep(4,4)
  categoryOrder[cat1] <- 1
  categoryOrder[cat2] <- 2
  categoryOrder[cat3] <- 3
  categoryOrder[cat4] <- 4
  #print(c(cat1, cat2, cat3, cat4))
  #print(categoryOrder)
  #print(head(df))
  df$state <- newState
  df[,9:12]<- data.frame(df$S1, df$S2, df$S3, df$S4)[,categoryOrder]
  #print(head(df))
  return(df)
}
## followed the function 'categoryChangeDF'
categoryAssian  <- function(med1, med2, cutOff_forPirorClass = 0) 
{
  if (med1 > cutOff_forPirorClass & med2 > cutOff_forPirorClass ) 
  {return (4)}
  if (med1 <= cutOff_forPirorClass & med2 > cutOff_forPirorClass) 
  {return (3)}
  if (med1 > cutOff_forPirorClass & med2 <= cutOff_forPirorClass) 
  {return (1)}
  if (med1 <= cutOff_forPirorClass & med2 <= cutOff_forPirorClass) 
  {return (2)}
}

########################################################################################
## Markov model building. Run 'timex' times, select the best model with highest loglik##
## Piror distribution is based on (piror) ##############################################                                            ##
#########################################################################################
markov <- function(x, timex = 10, cutOff_forPirorClass = median(c(x$corUP, x$corDown))) 
{
  print (paste( c('trainning on',  as.character(x$chr[1])) , sep = ' ') )
  print (paste('Piror correlation is ', cutOff_forPirorClass))
  result <- list()
  logLikhoodResult <- c()
  ## set piror
  piror <- rep(0, nrow(x))
  piror[which(x$corUP > cutOff_forPirorClass  & x$corDown > cutOff_forPirorClass)  ] <- 4
  piror[which(x$corUP <= cutOff_forPirorClass  & x$corDown <= cutOff_forPirorClass)] <- 2
  piror[which(x$corUP > cutOff_forPirorClass   & x$corDown <= cutOff_forPirorClass)] <- 3
  piror[which(x$corUP <= cutOff_forPirorClass  & x$corDown > cutOff_forPirorClass)]  <- 1
  pirorData <- data.frame(piror)
  ## for timex replicates , calculate log likelyhood and 
  for (i in 1:timex) 
  {
    mod <- mix(list(corUP ~ 1,corDown ~ 1), data = x,
               family = list(gaussian(), gaussian()),
               nstates = 4, 
               respstart = runif(16),
               #emcontrol=em.control(rand = F),
               prior =~ piror, initdata = pirorData
    )
    fm <- TryFit( fit(mod, verbose = FALSE, emc=em.control(rand = FALSE)) )
    
    ### if can't fit markov model, will return a vector from function 'TryFit '
    if (length(fm) > 1) 
    {
      result[[i]] <- c(0)
      logLikhoodResult[i] <- -1000000
    }else ## else put the final posterior to result
    {
      result[[i]] <- fm@posterior
      
      logLikhoodResult[i] <- logLik(fm)
    }
    
  }
  maxID       <- order(logLikhoodResult, decreasing = T)[1]
  
  df <- data.frame(x, result[[maxID]])
  ################# df category change
  df_final<- categoryChangeDF(df, cutOff_forPirorClass)
  
  ################plot df
  scatterDFplot(df_final, cutOff_forPirorClass)
  ###### return df with highest log likehood and max likehood
  return( list(df_final, max(logLikhoodResult)) )
}

#serving for function 'markov'
TryFit <- function(f){
  tryCatch(f,
           warning = function(e){ return( rep(NA,4) )}, 
           error = function(e){ return( rep(NA,4) )}) -> result
  return(result)
}
#########################################################################################