# Chad Hammerquist
# A set of functions to plot in R and create a pplotj document.
# These functions have most of the capability of R plot and of PublishPlot
# Using these functions plots the data in R and creates a .ppltj file
# the common functions are available: plot, lines, points, text as pp.plot, etc

# still needs work, not all options are possible yet

pp.plot=function (x, y = NULL, File="R_PPlotJ.ppltj",type = "p", xlim = NULL, ylim = NULL,
                log = "", main = NULL, xlab = NULL, ylab = NULL,lty=1,
                col="black",lwd=1,cex=1,pch=1) {
  # make the plot
  plot(x=x, y = y, type = type, xlim = xlim, ylim = ylim,
                  log = log, main = main, xlab = xlab, ylab = ylab,
                  col=col,lwd=lwd,cex=cex,lty=lty)

  # figure out labels
  xlabel = if (!missing(x)){
    deparse(substitute(x))
  }

  ylabel = if (!missing(y)){
    deparse(substitute(y))
  }

  #  get x-y coordinates and labels
  xy = xy.coords(x, y, xlabel, ylabel, log)
  if (is.null(xlab)){
    xlab = xy$xlab
  }
  if (is.null(ylab)){
    ylab = xy$ylab
  }

  # get limits
  if(is.null(xlim)){
    xlim = range(xy$x[is.finite(xy$x)])
    auto_x = "true" # should the axis automatically adjust
  }else{
    auto_x="false"
  }

  if(is.null(ylim)){
    ylim = range(xy$y[is.finite(xy$y)])
    auto_y = "true" # should the axis automatically adjust
  }else{
    auto_y="false"
  }

 # what about log?
  logX = "false"
  logY = "false"
  if(log=="x"||log=="xy"){
    logX = "true"
  }
  if(log=="y"||log=="xy"){
    logY="true"
  }

  # Create header and other settingss
  header = rbind("#PublishPlotJ Document",
              "#version	1",
              "#beginClass::JNPlotComments",
              " Created in R using package:ChadsTools",
              "#beginClass::JNPlotFrame",
              "#setStyle	true,0.0,0.0,true",
              "#setBGColor	(1.0 1.0 1.0 1.0)",
              "#setTicks	0.2,1.5,false",
              "#setPlotSize	725,510")


  LabelX = rbind("#beginClass::JNPlotLabelX",
                 paste('#setText',xlab,sep="\t"),
                 "#setFont	Arial,0,18.0",
                 "#setObjColor	(0.0 0.0 0.0 1.0)")

  LabelY = rbind("#beginClass::JNPlotLabelY",
                 paste('#setText',ylab,sep="\t"),
                 "#setFont	Arial,0,18.0",
                 "#setObjColor	(0.0 0.0 0.0 1.0)")



  AxisLabelX = rbind("#beginClass::JNXAxisLabels",
                     paste("#setRange",paste(xlim[1],xlim[2],0.1111*(xlim[2]-xlim[1]),0,sep=","),sep="\t"),
                     paste("#setOptions",paste("",auto_x,",",logX,",false",sep=""),sep="\t"),
                     "#setTickSides	true,false,true,true",
                     "#setFont	Arial,0,18.0",
                     "#setObjColor	(0.0 0.0 0.0 1.0)")

  AxisLabelY = rbind("#beginClass::JNYAxisLabels",
                     paste("#setRange",paste(ylim[1],ylim[2],0.1111*(ylim[2]-ylim[1]),0,sep=","),sep="\t"),
                     paste("#setOptions",paste("",auto_y,",",logY,",false",sep=""),sep="\t"),
                     "#setTickSides	true,false,true,true",
                     "#setFont	Arial,0,18.0",
                     "#setObjColor	(0.0 0.0 0.0 1.0)")

  # Parse color, line type,and etc
  color = as.numeric(col2rgb(col,TRUE)/255)
  line_width = 0.2*lwd
  line_type="none"
  if(type=="l"||type=="b"||type=="o"){
    if((lty == 1)||(lty == "solid")){
      line_type="solid"
    }else if((lty == 2)||(lty == "dashed")){
      line_type = "dashed"
    }else if((lty == 3)||(lty == "dotted")){
      line_type = "dotted"
    }else if((lty == 4)||(lty == "dotdash")){
      line_type = "dot-dash"
    }else{
      line_type="solid"
      warning("line type not available, reverting to solid line")
    }
  }

  # points
  point_type="none"
  if(type=="p"||type=="b"||type=="o"){
    if((pch == 0)){
      point_type="square"
    }else if((pch == 1)){
      point_type="circle"
    }else if((pch == 2)){
      point_type="triangle"
    }else if((pch == 5)){
      point_type="diamond"
    }else if((pch == 6)){
      point_type="inverted triangle"
    }else{
      point_type="circle"
      warning("point type not available, reverting to circle")
    }
  }
  # create array header
  ArrayHeader = rbind("#beginClass::JNPlotArray",
                      paste('#setColor\t"(',paste(color,collapse=" "),')"',sep=""),
                      paste("#setLineType",line_type,sep="\t"),
                      paste("#setLineWidth",line_width,sep="\t"),
                      paste("#setSymbolType",point_type,sep="\t"),
                      paste("#setSymbolSize",1.25*cex,sep="\t"),
                      paste('#setSymbolLineColor\t"(',paste(color,collapse=" "),')"',sep=""),
                      '#setName	"pp.plot 1"')

  # Full header
  Full_header = rbind(header,LabelX,LabelY,AxisLabelX,AxisLabelY,ArrayHeader)

  # file connection
  file_connection = file(File,"wb",encoding="UTF-8")

  # Write header
  cat(Full_header,file=file_connection,sep="\n")


  # write data
  write.table(cbind(xy$x,xy$y),sep="\t",file=file_connection,append=TRUE,col.names = FALSE,row.names = FALSE,
              eol="\n",fileEncoding="UTF-8")

  #variable
  assign("THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC",File,envir = .GlobalEnv)

  # close
  close(file_connection)

}


# the Line function
pp.lines = function (x, y = NULL,lty=1,col="black",lwd=1,cex=1){
  # make sure the file exist
  if(!exists("THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC")){
    stop("pplot must be called before pp.line")
  #}else if(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC==""){
    #stop("pplot must be called before pp.line")
  }else if(!file.exists(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC)){
    stop("Can't find .ppltj file. pp.plot must be called first")
  }
  #plot data
  lines(x,y,lty=lty,col=col,lwd=lwd,cex=cex)
  # format data
  xy = xy.coords(x, y)
  # Parse color, line type,and etc
  color = as.numeric(col2rgb(col,TRUE)/255)
  line_width = 0.2*lwd
  line_type="none"
  if((lty == 1)||(lty == "solid")){
    line_type="solid"
  }else if((lty == 2)||(lty == "dashed")){
    line_type = "dashed"
  }else if((lty == 3)||(lty == "dotted")){
    line_type = "dotted"
  }else if((lty == 4)||(lty == "dotdash")){
    line_type = "dot-dash"
  }else if((lty ==0)||(lty == "blank")){
    line_type = "filled"
  }else{
    line_type="solid"
    warning("line type not available, reverting to solid line")
  }
  # make array header
  ArrayHeader = rbind("","#beginClass::JNPlotArray",
                      paste('#setColor\t"(',paste(color,collapse=" "),')"',sep=""),
                      paste("#setLineType",line_type,sep="\t"),
                      paste("#setLineWidth",line_width,sep="\t"),
                      '#setName	"pp.line"')


  # file connection
  file_connection = file(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC,"ab",encoding="UTF-8")

  # Cat
  cat(ArrayHeader,file=file_connection,sep="\n",append=TRUE)

  # write data
  write.table(cbind(xy$x,xy$y),file=file_connection,
              append=TRUE,col.names = FALSE,row.names = FALSE,sep="\t",
              eol="\n",fileEncoding="UTF-8")
  # close
  close(file_connection)

}

# the point function
pp.points = function (x, y = NULL,pch=1,col="black",cex=1,type="p"){
  # make sure the file exist
  if(!exists("THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC")){
    stop("pplot must be called before pp.text")
    #}else if(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC==""){
    #stop("pplot must be called before pp.line")
  }else if(!file.exists(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC)){
    stop("Can't find .ppltj file. pplot must be called first")
  }

  #plot data
  points(x,y,pch=pch,col=col,cex=cex,type=type)

  # format data
  xy = xy.coords(x, y)

  # Parse color, line type,and etc
  color = as.numeric(col2rgb(col,TRUE)/255)
  if((pch == 0)){
    point_type="square"
  }else if((pch == 1)){
    point_type="circle"
  }else if((pch == 2)){
    point_type="triangle"
  }else if((pch == 5)){
    point_type="diamond"
  }else if((pch == 6)){
    point_type="inverted triangle"
  }else{
    point_type="circle"
    warning("point type not available, reverting to circle")
  }

  # check plot type
  line_type="none"
  if(type=="p"||type=="b"||type=="o"){
    if(type=="b"||type=="o"){
      line_type = "solid"
    }
  }else{
    warning("plot type is not available, reverting to points only")
  }
  # make array header
  ArrayHeader = rbind("","#beginClass::JNPlotArray",
                      paste('#setColor\t"(',paste(color,collapse=" "),')"',sep=""),
                      paste("#setLineType",line_type,sep="\t"),
                      paste("#setSymbolType",point_type,sep="\t"),
                      paste("#setSymbolSize",1.25*cex,sep="\t"),
                      paste('#setSymbolLineColor\t"(',paste(color,collapse=" "),')"',sep=""),
                      '#setName	"point 1"')

  # file connection
  file_connection = file(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC,"ab",encoding="UTF-8")

  # Cat
  cat(ArrayHeader,file=file_connection,sep="\n",append=TRUE)

  # write data
  write.table(cbind(xy$x,xy$y),file=file_connection,
              append=TRUE,col.names = FALSE,row.names = FALSE,sep="\t",
              eol="\n",fileEncoding="UTF-8")
  # close
  close(file_connection)

}

# the text function (not finished)
pp.text = function (x, y = NULL,labels=seq_along(x),col="black",cex=1){
  # make sure the file exist
  if(!exists("THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC")){
    stop("pp.plot must be called before pp.text")
    #}else if(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC==""){
    #stop("pplot must be called before pp.line")
  }else if(!file.exists(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC)){
    stop("Can't find .ppltj file. pp.plot() must be called first")
  }
  #plot data
  text(x,y,labels=labels,col=col,cex=cex)
  # format data
  xy = xy.coords(x,y)
  # Parse color, line type,and etc
  color = as.numeric(col2rgb(col,TRUE)/255)

  # make array header
  ArrayHeader = rbind("","#beginClass::JNPlotText",
                      paste('#setPoint\t',x[1],",",y[1],collapse="",sep=""),
                      paste('#setText\t',labels[1],collapse="",sep=""),
                      paste('#setFont\tArial,0,',cex*18,collapse="",sep=""),
                      paste('#setObjColor\t"(',paste(color,collapse=" "),')"',sep=""))

  # file connection
  file_connection = file(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC,"ab",encoding="UTF-8")

  # Cat
  cat(ArrayHeader,file=file_connection,sep="\n",append=TRUE)

  # write data
  write.table(cbind(xy$x,xy$y),file=file_connection,
              append=TRUE,col.names = FALSE,row.names = FALSE,sep="\t",
              eol="\n",fileEncoding="UTF-8")
  # close
  close(file_connection)

}

# the text function (not finished)
pp.arrows = function (x0, y0,x=x0,y=y0,col="black",cex=1){
  # make sure the file exist
  if(!exists("THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC")){
    stop("pp.plot must be called before pp.arrow")
  }else if(!file.exists(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC)){
    stop("Can't find .ppltj file. pp.plot() must be called first")
  }
  #plot data
  arrows(x0,y0,x,y,col=col,lwd=cex)

  # Parse color, line type,and etc
  color = as.numeric(col2rgb(col,TRUE)/255)

  # make array header
  ArrayHeader = rbind("","#beginClass::JNPlotArrow",
    "#setShape	0",
    paste('#setColor\t"(',paste(color,collapse=" "),')"',sep=""),
    paste('#setLineWidth\t',cex*0.1,collapse="",sep=""),
    paste('#setSymbolLineWidth\t',cex*0.1,collapse="",sep=""),
    paste('#setSymbolSize\t',cex*1.25,collapse="",sep=""),
    paste('#setSymbolLineColor\t"(',paste(color,collapse=" "),')"',sep=""),
    paste('#setSymbolFillColor\t"(',paste(color,collapse=" "),')"',sep=""),
    "#setName\t R arrow")

  # file connection
  file_connection = file(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC,"ab",encoding="UTF-8")

  # Cat
  cat(ArrayHeader,file=file_connection,sep="\n",append=TRUE)

  # write data
  write.table(cbind(rbind(x0,x),rbind(y0,y)),file=file_connection,
              append=TRUE,col.names = FALSE,row.names = FALSE,sep="\t",
              eol="\n",fileEncoding="UTF-8")
  # close
  close(file_connection)

}

# Close out connection
pp.close=function(){
  if(exists("THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC")){
    rm("THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC")
  }
}

# the Line function
pp.fanplot = function (x,y,ylim,percent=NULL,ind=NULL){
  # make sure the file exist
  if(!exists("THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC")){
    stop("pplot must be called before pp.fanplot")
    #}else if(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC==""){
    #stop("pplot must be called before pp.line")
  }else if(!file.exists(THE_FILE_NAME_FOR_PPLOT_AND_PPLINES_ETC)){
    stop("Can't find .ppltj file. pp.plot must be called first")
  }

  # get size of y data
  ysize = dim(y)
  if(is.null(ysize)){
    stop(" y needs to be a matrix")
  }

  # make ends work
  bot = rep(ylim[1],ysize[2])
  X = c(x[1],x,x[length(x)])
  Y = rbind(bot,y,bot)

  # set up colors
  ncols = floor(ysize[2]/2)
  if(2*ncols == ysize[2]){
    col = rev(grey.colors(ncols))
    col = c(col,rev(col))
  }else{
    col = rev(grey.colors(ncols+1))
    col = c(col,rev(col))
  }

  # loop through and plot filled curves
  for(k in ysize[2]:2){
    pp.lines(X,Y[,k],lty="blank",col=col[k])
  }
  pp.lines(X,Y[,1],lty="blank",col="white")

   # plot percentile lines
  if(!is.null(ind)){
    for(k in ind){
      pp.lines(x,y[,k],lwd=0.5,col=grey(0.15))
    }
  }
  # plot percentile text
  xt  = max(x)+.05*abs(diff(range(x)))
  if(!is.null(percent)){
    for(k in ind){
      pp.text(xt,y[ysize[1],k],percent[k],cex=0.75)
    }
  }

}
