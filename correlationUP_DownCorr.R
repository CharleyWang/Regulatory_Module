######################################################################
### calculated the correlation with the elements located on upstream##
### or downstream                                                   ##
######################################################################

#!!! requirement
## colI in df should be chromosome
## colII to colIII should be start and end position of elements
## colIV to colV   should be information such as ID or strand of gene
## colVI to the end should be sequencing signal for different stages or tissues


## numGeneCheck indicates the number of genes on upstream or downstream are
## considered here for correlation calculation (average correlation is calculated)

correlationUP_DownCorr <- function(df, numGeneCheck = 1) 
{
  #rank the genes based on chromosome and start sites
  orderID <-  order(df$chr, df$start) 
  dfOrder <-  df[orderID, ]
  
  #Initial a data frame, final result will put inside this data frame
  corUP_Down <- data.frame(
    chr    = 'chrZ', start  = 1, end = 1, strand = '+',
    geneID = 'PF', corUP  =  0, corDown=  0)
  
  ## iterating chromosome, calculate correlation with up or downward 
  aa <- sapply(unique(dfOrder$chr), 
               function(x)
               {
                 dfChr <- dfOrder[which(dfOrder$chr == x), ]
                 startRow <- numGeneCheck + 1
                 endRow   <- nrow(dfChr) - numGeneCheck
                 
                 ###for beginning $numGeneCheck genes, only calculate downward correlation
                 ###corrlation with upstream is defined as 0 here
                 for (j in 1:numGeneCheck) 
                 {
                   downMatrix<- dfChr[j: (j+numGeneCheck)       , 5:ncol(dfChr)]
                   rownames(downMatrix) <- downMatrix[,1]
                   downMatrix     <- downMatrix[,-1]
                   corDown        <- cor(t(downMatrix))
                   corDOWNv       <- corDown[1, ]
                   corDOWNval <- (sum(corDOWNv[!is.na(corDOWNv)]) - 1)/numGeneCheck
                   corDF <- data.frame(
                     chr    = dfChr$chr[j],    start  = dfChr$start[j],  end    = dfChr$end[j],
                     strand = dfChr$strand[j], geneID = dfChr$geneID[j], corUP  = 0,
                     corDown= corDOWNval
                   )
                   corUP_Down <<- rbind(corUP_Down, corDF)
                 }
                 ############for the genes in the middle part, calculate both up and down stream
                 for (i in startRow : endRow) 
                 {
                   upMatrix  <- dfChr[(i - numGeneCheck ) : i, 5:ncol(dfChr)]
                   downMatrix<- dfChr[i: (i+numGeneCheck)    , 5:ncol(dfChr)]
                   rownames(upMatrix)   <- upMatrix[,1]
                   rownames(downMatrix) <- downMatrix[,1]
                   upMatrix       <- upMatrix[,-1]
                   downMatrix     <- downMatrix[,-1]
                   corUP          <- cor(t(upMatrix))
                   corDown        <- cor(t(downMatrix))

                   corUPv         <- corUP[numGeneCheck + 1,]
                   corDOWNv       <- corDown[1, ]
                   # minus 1(-1) to remove the correlation with itself
                   corUPval   <- (sum(corUPv[!is.na(corUPv)]) - 1)/numGeneCheck
                   corDOWNval <- (sum(corDOWNv[!is.na(corDOWNv)]) - 1)/numGeneCheck
                   
                   corDF <- data.frame(
                     chr    = dfChr$chr[i],
                     start  = dfChr$start[i],
                     end    = dfChr$end[i],
                     strand = dfChr$strand[i],
                     geneID = dfChr$geneID[i],
                     corUP  = corUPval,
                     corDown=corDOWNval
                   )
                   corUP_Down <<- rbind(corUP_Down, corDF)
                 }
                 ##### for genes in the end, only calcuated the correlation with upstream
                 for (m in (nrow(dfChr) - numGeneCheck + 1) : nrow(dfChr)) 
                 {
                   upMatrix  <- dfChr[(m - numGeneCheck ) : m, 5:ncol(dfChr)]
                   rownames(upMatrix)   <- upMatrix[,1]
                   upMatrix       <- upMatrix[,-1]
                   corUP          <- cor(t(upMatrix))
                   #print(corUP)
                   
                   corUPv         <- corUP[numGeneCheck + 1,]
                   corUPval   <- (sum(corUPv[!is.na(corUPv)]) - 1)/numGeneCheck
                   corDF <- data.frame(
                     chr    = dfChr$chr[m],
                     start  = dfChr$start[m],
                     end    = dfChr$end[m],
                     strand = dfChr$strand[m],
                     geneID = dfChr$geneID[m],
                     corUP  = corUPval,
                     corDown= 0
                   )
                   corUP_Down <<- rbind(corUP_Down, corDF)
                 }
               }
  )
  corUP_Down <- corUP_Down[-1,]
  corUP_Down$chr <- factor(corUP_Down$chr)
  return(corUP_Down)
}