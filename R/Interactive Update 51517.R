# Interactive plot for homol series ######################################################
#
# Christopher Ruybal and Karl Oetjen
# Created 5/17/2016
# Modified 5/17/2016
##########################################################################################
# Uses outputs from Martin Loos' "nontarget" and "nontargetData" packages found on GitHub

plot_interact <-
  
  function (homol) 
   {
    
    
    
    dat1 <- data.frame(mz = homol[[1]][,1],
                       RT = homol[[1]][,3],
                       PeakID = homol[[1]][,4],
                       GroupID = homol[[1]][,5])
    
    I <- which(dat1$GroupID == 0)
    
    dat2 <- data.frame(mz = homol[[1]][I,1],
                       RT = homol[[1]][I,3], col = 0)  
    
    
    list2 <- matrix(data=NA,nrow=strtoi(max(summary(homol[[5]])[,1])),ncol=length(homol[[5]]))
    
    
    for (ll in 1:length(homol[[5]])) {
      for (mm in 1:length(homol[[5]][[ll]])) {
        list2[mm,ll] <- homol[[5]][[ll]][[mm]][1]
      }
    }
    
    
    this <- round(homol[[3]][, 3], digits = 2)
    that <- levels(as.factor(this))
    
    colo <- rainbow(length(that))
    
    
    
    df <- NULL
    for (p in 1:length(homol[[5]])) {
      temp_df <- data.frame(x=dat1$RT[list2[,p]], y=dat1$mz[list2[,p]], col = p, MZ = this[p],PeakID = paste("PeadkID: " ,dat1$PeakID[list2[,p]]))
      df <- rbind(df,temp_df) 
     }
    
    
    names(df)[1] = "RT"
    names(df)[2] = "mz"
    
    
    #plot_ly(data = df, x = RT, y = mz, mode = "markers+lines",
    #        color = factor(col,labels = this)) 
      
    #plot_ly(data = dat2, x = RT, y = mz, mode = "markers",  marker=list(color="grey" , size=5 , opacity=0.5))
    
    f <- list(
      family = "sans-serif",
      size = 12,
      color = "#000"
    )
    l <- list(
      font = f,
      bgcolor = NA,
      bordercolor = NA,
      
      borderwidth = 2
      
    )
   
    yy <- list(title = "m/z")
    
    
    idx = !is.na(df$RT)
    
    df2 <- df[idx,]
    
    
    # Get unique MZ
    mzU <- unique(df2$MZ)
    
    cols <- grDevices::rainbow(length(unique(df$MZ)))
    
    p <- plot_ly()
    
    # for each unique MZ
    for (ii in 1:length(mzU)){
      
      
      
      p <- add_trace(p,data = filter(df2, df2$MZ == mzU[ii]),
                     type = "scatter",
                     x = ~RT,
                     y = ~mz,
                     mode = "markers",
                     color = ~factor(MZ,as.character(mzU), mzU),
                     colors = grDevices::rainbow(length(unique(df$MZ))),
                     legendgroup = as.character(mzU[ii]),
                     name = as.character(mzU[ii]),
                     hoverinfo = "none",
                     showlegend = TRUE) 
      
      
      p <- add_trace(p,data = filter(df2, MZ == mzU[ii]),
                     type = "scatter",
                     x = ~RT,
                     y = ~mz,
                     text = df2$PeakID[which(df2$MZ== mzU[ii])],
                     split = ~col,
                     color = ~factor(MZ,as.character(mzU), mzU),
                     colors = grDevices::rainbow(length(unique(df$MZ))),
                     mode = "lines",
                     legendgroup = as.character(mzU[ii]),
                     name = as.character(mzU[ii]),
                     hoverinfo = "x+y+text",
                     showlegend = FALSE)
      
    }
    
    p <- add_trace(p,data = dat2, 
                   type = "scatter",
                   x = ~RT, 
                   y = ~mz, 
                   mode = "markers",  
                   marker=list(color="grey" , size=5 , opacity=0.5),
                   name="No Series",
                   hoverinfo = "x+y",
                   showlegend=TRUE)
    
    p
    # Set the font and size for axis labels 
    f1 <- list(
      family = "Arial, sans-serif",
      size = 18,
      color = "black"
    )
    
    
    # Set the font and size for the tick text 
    f2 <- list(
      family = "Old Standard TT, serif",
      size = 14,
      color = "black"
    )
    
    # Y axis setup 
    axy <- list(
      title = "m/z",
      titlefont = f1,
      showticklabels = TRUE,
      tickangle = 0,
      tickfont = f2,
      autotick = TRUE,
      ticks = "inside",
      tick0 = NA,
      dtick = NA,
      ticklen = 5,
      tickwidth = 2,
      showline = TRUE,
      mirror = "ticks",
      gridcolor = toRGB("white"),
      gridwidth = 0.5,
      linecolor = toRGB("black"),
      linewidth = 1
    )
    
    # X axis setup 
    axx <- list(
      title = "RT (minutes)",
      titlefont = f1,
      showticklabels = TRUE,
      tickangle = 0,
      tickfont = f2,
      autotick = TRUE,
      ticks = "inside",
      tick0 = NA,
      dtick = NA,
      ticklen = 5,
      tickwidth = 2,
      showline = TRUE,
      mirror = "ticks",
      gridcolor = toRGB("white"),
      gridwidth = 0.5,
      linecolor = toRGB("black"),
      linewidth = 1
    )
    
    layout(p,xaxis = axx, yaxis = axy) 
  }