###################################################################################################
### Fit model for each chromosome based on function 'markov'######################################
### CorCutoff represents the candicate correlation cutoff   ######################################
### FittingTimes represents number of times running to find ######################################
### best candicat correlation                               ######################################
###################################################################################################

## require script 'markov.R' and 'markovPirorSelect.R'
markovFit <- function(df, fittingTimes = 5, 
                      corCutoff = seq(from = 0, to = 0.8, by = 0.05)) 
{
  chrBestCor <- list()
  fiRe <- data.frame(chr = 'chrZ',start = 1,
                     end = 1, strand = '+',
                     geneID = 'X', corUP = 1.1, corDown = 1.1,
                     state = 100, S1 = 1.1, S2 = 1.2, S3 = 1.3, S4 = 1.5)
  df$corUP   <- round(df$corUP, 3)
  df$corDown <- round(df$corDown, 3)
  logLikeSum <- list()
  
  aa <- sapply(unique(df$chr),
               function(x) 
               {
                 dfChr <- df[which(df$chr == x),]
                 ### fit corCutoff for pirior
                 bestCutOff <- markovPirorSelect(dfChr, fittingTimes = fittingTimes, 
                                                 corCutoff = corCutoff)
                 #bestCutOff <- corCutoff
                 print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
                 print(paste('trainning model based on selected correlaiton piror', bestCutOff))
                 print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
                 rrChr <- markov(dfChr, timex = fittingTimes, cutOff_forPirorClass = bestCutOff)
                 fiRe  <<-rbind(fiRe, rrChr[[1]])
                 logLike                         <- rrChr[[2]]
                 logLikeSum[[as.character(x)]]   <<- logLike
                 chrBestCor[[as.character(x)]]   <<- bestCutOff
               }
  )
  finalData <- fiRe[-1, ] 
  #print(chrBestCor)
  finalList <- list(finalData, 
                    logLikelihood = logLikeSum, 
                    corForTrain   = chrBestCor)
  
  print(paste("Done!"))
  return(finalList)
}