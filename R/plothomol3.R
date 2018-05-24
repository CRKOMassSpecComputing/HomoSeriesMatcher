######################################################################################################################################
# Modified from Martin Loos' "nontarget" package found on GitHub at https://github.com/blosloos/nontarget/blob/master/R/plothomol.R
#
# Modified by Christopher Ruybal and Karl Oetjen
# Modifed on 5/1/2016
####################################################################


plothomol3 <-
  
  function (homol, xlim = FALSE, ylim = FALSE, plotlegend = TRUE) 
  {
    if (xlim[1] != FALSE) {
      if (length(xlim) > 2) {
        stop("xlim not correct!")
      }
    }
    if (ylim[1] != FALSE) {
      if (length(xlim) > 2) {
        stop("xlim not correct!")
      }
    }
    if (length(homol[[5]]) < 1) {
      stop("no homologue series found!")
    }
    sc <- close.screen()
    if (sc[1] != FALSE) {
      for (m in 1:length(sc)) {
        close.screen(sc[m])
      }
    }
    plot.new()
    if (xlim[1] != FALSE & ylim[1] == FALSE) {
      plot.window(xlim = xlim, ylim = c(min(homol[[1]][, 1]), 
                                        max(homol[[1]][, 1])))
    }
    if (ylim[1] != FALSE & xlim[1] == FALSE) {
      if (plotlegend == TRUE) {
        plot.window(xlim = c(min(homol[[1]][, 3]), max(homol[[1]][, 
                                                                  3]) * 1.2), ylim = ylim)
      }
      else {
        plot.window(xlim = c(min(homol[[1]][, 3]), max(homol[[1]][, 
                                                                  3])), ylim = ylim)
      }
    }
    if (xlim[1] == FALSE & ylim[1] == FALSE) {
      if (plotlegend == TRUE) {
        plot.window(xlim = c(min(homol[[1]][, 3]), max(homol[[1]][, 
                                                                  3]) * 1.2), ylim = c(min(homol[[1]][, 1]), max(homol[[1]][, 
                                                                                                                            1])))
      }
      else {
        plot.window(xlim = c(min(homol[[1]][, 3]), max(homol[[1]][, 
                                                                  3])), ylim = c(min(homol[[1]][, 1]), max(homol[[1]][, 
                                                                                                                      1])))
      }
    }
    if (xlim[1] != FALSE & ylim[1] != FALSE) {
      plot.window(xlim = xlim, ylim = ylim)
    }
    
    box()
    axis(1)
    axis(2)
    title(xlab = "Retention time", ylab = "m/z")
    points(homol[[1]][, 3], homol[[1]][, 1], cex = 0.3, pch = 19, 
           col = "lightgrey")
    
    dat1 <- data.frame(mz = homol[[1]][,1],
                       RT = homol[[1]][,3],
                       PeakID = homol[[1]][,4],
                       GroupID = homol[[1]][,5])
    
    
  for (rrr in 1:length(dat1$mz)){
    if ((dat1$GroupID[[rrr]]==0) & (!is.na(namelist$Name[rrr]))) {
      textxy(dat1$RT[rrr],dat1$mz[rrr],dat1$PeakID[rrr])
     }
    }
    
    
    this <- round(homol[[3]][, 3], digits = 2)
    this_rt <- round(homol[[3]][, 4], digits = 2)
    that <- levels(as.factor(this))
    
    colo <- rainbow(length(that))
    
    for (i in 1:length(homol[[5]])) {
      for (j in 2:length(homol[[5]][[i]])) {
        for (k in 1:length(homol[[5]][[i]][j - 1])) {
          for (m in 1:length(homol[[5]][[i]][j])) {
            
            lines(c(homol[[1]][homol[[5]][[i]][[j - 1]][k], 
                               3], homol[[1]][homol[[5]][[i]][[j]][m], 3]), 
                  c(homol[[1]][homol[[5]][[i]][[j - 1]][k], 
                               1], homol[[1]][homol[[5]][[i]][[j]][m], 
                                              1]), col = colo[that == this[i]], lwd = 1.8)
            
            points(homol[[1]][homol[[5]][[i]][[j - 1]][k], 
                              3], homol[[1]][homol[[5]][[i]][[j - 1]][k], 
                                             1], col = colo[that == this[i]], pch = 19, 
                   cex = 0.5)
            
            
           
            
            if (is.na(namelist$Name[homol[[1]][homol[[5]][[i]][[j - 1]][k],4]])) {
            
            textxy(homol[[1]][homol[[5]][[i]][[j - 1]][k], 
                              3], homol[[1]][homol[[5]][[i]][[j - 1]][k],1], 
                   homol[[1]][homol[[5]][[i]][[j - 1]][k],4])
            } else{
            textxy(homol[[1]][homol[[5]][[i]][[j - 1]][k], 
                                3], homol[[1]][homol[[5]][[i]][[j - 1]][k],1], 
                     homol[[1]][homol[[5]][[i]][[j - 1]][k],4],col=colo[that == this[i]])
            }  
              
              
            
            points(homol[[1]][homol[[5]][[i]][[j]][m], 
                              3], homol[[1]][homol[[5]][[i]][[j]][m], 1], 
                   col = colo[that == this[i]], pch = 19, cex = 0.5)
            
            if (is.na(namelist$Name[homol[[1]][homol[[5]][[i]][[j]][m], 4]])) {
             textxy(homol[[1]][homol[[5]][[i]][[j]][m], 
                                3], homol[[1]][homol[[5]][[i]][[j]][m], 1],
                     homol[[1]][homol[[5]][[i]][[j]][m], 4])
          
            
            }else {
              textxy(homol[[1]][homol[[5]][[i]][[j]][m], 
                              3], homol[[1]][homol[[5]][[i]][[j]][m], 1],
                   homol[[1]][homol[[5]][[i]][[j]][m], 4],col=colo[that == this[i]])
            }
            
        }
      }
    }
  }  
    

  ####Legend for M/Z only############
    
    if (plotlegend == TRUE) {
      plot.window(xlim = c(0, 1), ylim = c(min(as.numeric(that)), 
                                           max(as.numeric(that))))
      lines(c(0.94, 0.94), c(min(as.numeric(that)), max(as.numeric(that))), 
            col = "lightgrey", lwd = 6)
      
      it <- 2
      
      for (i in 1:length(that)) {
        points(0.94, as.numeric(that[i]), pch = 19, col = colo[i])
        text(0.94, as.numeric(that[i]), labels = that[i], 
             col = colo[i], cex = 0.65, pos = it)
        text(0.94, max(as.numeric(that)), labels = expression(paste(Delta,"m/z")) , 
             cex = 0.65, pos = 3)
      }
    }
      
    ######### For RT Shift####  
    
      if (plotlegend == TRUE) {
        
        plot.window(xlim = c(0, 1), ylim = c(min(as.numeric(this_rt)), 
                                             max(as.numeric(this_rt))))
        
        lines(c(0.98, 0.98), c(min(as.numeric(this_rt)), max(as.numeric(this_rt))), 
              col = "lightgrey", lwd = 6)
        
        it_rt <- 4
        
        for (i in 1:length(this_rt)) {
          points(0.98, as.numeric(this_rt[i]), pch = 19, col = colo[i])
          text(0.98, as.numeric(this_rt[i]), labels = this_rt[i], 
               col = colo[i], cex = 0.65, pos = it_rt)
          text(0.98, max(as.numeric(this_rt)),labels = expression(paste(Delta,"RT")) 
               ,cex = 0.65, pos = 3)
        
        }
      }
}


    
  
