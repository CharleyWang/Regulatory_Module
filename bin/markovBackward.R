##################################################################
### Based on trainning model from markov, backward search model###
##################################################################
### Module can only make up by '1', '4', '3'. 
### The module should be end by '4' ||'3', if the next element is '1'
### Or should be end by '2'
backwardMarkov <- function( mrV, medianPosteirorV = 0.9 ) 
{
  finalModule <- data.frame(chr = 'chrZ',start = 1,end   = 1,numGene = 1, moduleID = 'pfxxxxx', pro = 0.5)
  finalModule <- finalModule[-1, ]
  colData   <- ncol(mrV)
  
  ### for each module candidate change its format
  moduleCheck <- function(moduleCandidate, cutoffGene = 2) 
  {
    numGene <- nrow(moduleCandidate)
    if (numGene >= cutoffGene) 
    {
      chr     <- moduleCandidate$chr[1]
      prob    <- moduleCandidate[,9:12]
      state   <- moduleCandidate$state
      geneID  <- paste(moduleCandidate$geneID, collapse = ',')
      #start   <- moduleCandidate$start[1]
      start   <- min(moduleCandidate$start)
      #end     <- moduleCandidate$end[nrow(moduleCandidate)]
      end     <- max(moduleCandidate$end)
      pro     <- mean( sapply(1:nrow(prob), function(u) prob[u, state[u]])  )
      data    <- data.frame(chr = chr,start = start,end   = end,numGene = numGene, 
                            moduleID = geneID, pro = pro)
      finalModule  <<- rbind(finalModule, data)
    }
  }
  
  aa <- sapply(unique(mrV$chr), 
               function(x) # for each chromosome
               {
                 dfChr <- mrV[which(mrV$chr == x), ]
                 pos   <- dfChr[, (colData - 3): colData]  
                 state <- dfChr[, colData  - 4]
                 
                 ## put the first gene into module candicate if its state != 2
                 ModuleTest  <- dfChr[1,]
                 if (state[1] == 2) { ModuleTest <- ModuleTest[-1, ] }
                 
                 ## from second genes to the end
                 for (i in 2:nrow(dfChr)) 
                 {
                   if (state[i] == 1 & (state[i-1] == 3|state[i-1] == 4)) 
                   {
                     if (nrow(ModuleTest) == 0) { ModuleTest <- rbind(ModuleTest, dfChr[i, ])
                     } else 
                     {
                       moduleCheck(ModuleTest);
                       ModuleTest <- ModuleTest[-c(1:nrow(ModuleTest)), ]
                       ModuleTest <- rbind(ModuleTest, dfChr[i, ])
                     }
                   } else if (state[i] == 1) {  ModuleTest <- rbind(ModuleTest, dfChr[i, ])}
                   
                   if (state[i] == 2 & nrow(ModuleTest) > 0) {                           
                     moduleCheck(ModuleTest);
                     ModuleTest <- ModuleTest[-c(1:nrow(ModuleTest)), ] }
                   
                   if (state[i] == 3 | state[i] == 4) { ModuleTest <- rbind(ModuleTest, dfChr[i, ]) }
                   
                   if (i == nrow(dfChr) & (state[i] == 3 | state[i] == 4)) 
                   {moduleCheck(ModuleTest);} 
                 }
                 return(0)
               }
  )
  return(finalModule)
}
