#Sequencing Data input
#load RNA-seq data
dfRNA <- read.table('gene.exp.signal',
                         header = T, sep = '\t')
# load ATAC-Seq data
dfATAC    <- read.table('ATAC.signal.all.intergenic',
                        header = F, sep = '\t') 
# re-type column name
dfATAC   <-  data.frame(chr = dfATAC[,1], start = dfATAC[,2], end = dfATAC[,3],
                        strand = rep('ATACp', nrow(dfATAC)), geneID = dfATAC[,4], 
                        t5 = dfATAC[,5], t10 = dfATAC[,6], t15 = dfATAC[,7],
                        t20 = dfATAC[,8], t25 = dfATAC[,9], t30 = dfATAC[,10], 
                        t35 =dfATAC[,11], t40 = dfATAC[,12])
# combine RNA-Seq and ATAC data into on matrix
dfINPUTdata <- rbind(dfATAC, dfRNA)
# order the row based on chromosome and start site
df          <- dfINPUTdata[ order(dfINPUTdata$chr, dfINPUTdata$start),   ]

# load Hidden Markov package
# original paper link: 'https://cran.r-project.org/web/packages/depmixS4/vignettes/depmixS4.pdf'
library("depmixS4")


###load function 'correlationUP_DownCorr'
setwd('~/Documents/work/GeneExpModule_ATAC_GENE/final_Package/')
source('correlationUP_DownCorr.R')
corUP.Down <- correlationUP_DownCorr(df, numGeneCheck = 1)

#########################################################################################
##############################Training markov model######################################
source('markov.R')
source('markovPirorSelect.R')
source('markovFit.R')
markovModel <- markovFit(corUP.Down, 20, corCutoff = seq(from = 0, to = 0.8, by = 0.05) )
##In first list from the out put of 'markovFit', using the number '1-4' represents true state
##'1' represent downward correlation
##'2' represent no corelation with both side ## '3' represent upward rorrelation 
##'4' represent correlation with both side

#########################################################################################

### based on status trained from markov model
### extract true gene modules. all modules should be seperated by '3||4, 1' or '2'
source('~/Documents/work/GeneExpModule_ATAC_GENE//markov/markovBackward.R')
module      <- backwardMarkov(markovModel[[1]])
write.table(module, file = '~/Documents/work/GeneExpModule_ATAC_GENE/moduleRegion//module.bed', 
            row.names = F, col.names = T, quote = F,sep = '\t')
