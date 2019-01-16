#########################################################################################
##For piror correlation cutoff, select best cutoff based on loglik#######################
#########################################################################################

markovPirorSelect <- function(df, fittingTimes = 20, 
                              corCutoff = seq(from = 0, to = 0.8, by = 0.05)) 
{
  print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
  print('$ Select best correlation cutoff for piror distribution $')
  print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
  
  likhod <- c()
  for (i in corCutoff)
  {
    print('--------------------------------------------------------------------------')
    print(paste('fit model for corCutoff = ', i, sep = ' ')) 
    
    moT <- markov(df, timex = fittingTimes, cutOff_forPirorClass = i)
    likhod <- c(likhod, moT[[2]])
    
    print(paste('log likelyhood is', moT[[2]], sep = ' '))
  }
  ff <- data.frame(corCutoff, likhod)
  bestCutOff <- ff$corCutoff[order(ff$likhod, decreasing = T)[1]]
  
  print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
  print(paste('The best piror fitting correlatin cutoff is', bestCutOff, sep = ' ') )
  print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
  return(bestCutOff)
}
###################################################################################################