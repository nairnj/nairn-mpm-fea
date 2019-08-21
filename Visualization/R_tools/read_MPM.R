# Chad Hammerquist
# This contains two functions needed to read files from
# OSParticulas MPM software developed by John Nairn.

# Contains 2 functions , read binary and read global

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note: New and better way to read binary archive files:
#: command = paste(ExtractMPM -T,mat_option,data_option,file,collapse="",sep="")
#: z = pipe(command)
#: out = read.table(z)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This is a function to read in .global files
# Doesn't have error checking (yet)
ReadMPMGlobal=function(filename=NULL,Col=NULL){
  # Read in data, allow for the use of clipboard
  if((tolower(filename)=="clipboard")||is.null(filename)){
    temp=read.delim("clipboard",header=TRUE,stringsAsFactors=FALSE) # read in data
  }
  else{
    temp=read.delim(filename,header=TRUE,stringsAsFactors=FALSE) # read in data
  }
  if(!is.null(Col)){
    column_names=c("Time",temp[1,Col]) # extract column names
    temp=temp[,c(1,Col)] # remove all but selected columns
  }
  else{
    column_names=c("Time",temp[1,]) # Chage the name of the first col
    column_names=column_names[-2]
  }
  temp=temp[-1,] # remove first line
  colnames(temp)=column_names # name columns
  # Convert to numeric
  nd=dim(temp)
  for(k in 1:nd[2]){
    temp[,k]=as.numeric(temp[,k])
  }
  #
  return(temp)
}



# Read MPM binary archive files
# Not finished, kind of rough
# cant handle temperature
# and no cracks
ReadMPM = function(filename,DataOption=""){
  # Open Connection
  out = list()
  file_size = file.info(filename)$size
  con = file(description = filename, open = "rb")

  # Read Header
  out$ver = readChar(con,4,TRUE) # version of file
  nchar1 = as.integer(paste(c("0x",readBin(con,raw(),1)),collapse="")) #number of characters in archive format
  out$AFormat = readChar(con,nchar1,TRUE) #Archive format
  nchar2 = as.integer(paste(c("0x",readBin(con,raw(),1)),collapse="")) #number of characters in crack format
  out$CFormat = readChar(con,nchar2,TRUE) #Crack Archive format
  out$Dim = as.integer(readChar(con,1,TRUE))  # Dimension of MPM 3D or 2D
  lf = 64-nchar1-nchar2-7
  if(lf>0){
    leftover = readBin(con,"raw",n=lf) #Advance file
  }

  # Get Particle Size in bytes
  StoredFlag = unlist(strsplit(out$AFormat,split = ""))[-1]
  FormatExpanded = StoredFlag=="Y"
  historyoptions = c("Y",1,2,3,4,5,6,7,8,9,":",";","<","=",">","?")
  historysize = c(1,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4)
  if(out$Dim == 2){
    PSize = sum((FormatExpanded)*c(64,16,32,32,32,0,8,8,8,0,16,8,0,24,8,4,8,16))
    #CFormatExpanded = unlist(strsplit(out$CFormat,split = ""))[-1] == "Y"
    #PSize = PSize+sum((CFormatExpanded)*c(88,16,16,20))

  }else{
    PSize = sum((FormatExpanded)*c(88,24,48,48,48,0,8,8,8,0,0,8,0,32,8,4,24,24))
  }
  num_of_hist =sum(historysize*(historyoptions==StoredFlag[13]))
  PSize = PSize+8*num_of_hist  # history variable
  NPart = floor((file_size-64)/PSize)
  out$NPart = NPart

  ## Pre-allocate
  size1 = out$Dim
  size2 = size1*2
  size2d = size1*2*8
  size1d = size1*8
  sizecd = (size1+1)*8
  sizerd = ((size1==2)+(size1==3)*3)*8
  sizehv = 8*num_of_hist
  # See what you data you want
  Property = unlist(tolower(strsplit(DataOption,split=" ")))
  compare = c("mass","velocity","stress","strain","plasticstrain",
              "workenergy","plasticenergy","totalshear","strainenergy","history","damagenormal")
  select_prop = !is.na(pmatch(compare,Property))
  select_prop = select_prop&c(1,FormatExpanded[c(2:5,7,9,11:13,18)])

  # Default property
  out$material = matrix(0,NPart,1)
  out$position = matrix(0,NPart,size2)


  # Mass?
  if(select_prop[1]){
    out$mass = matrix(0,NPart,1)
  }
  # Velocity?
  if(select_prop[2]){
    out$velocity = matrix(0,NPart,size1)
  }
  # Stress?
  if(select_prop[3]){
    out$stress = matrix(0,NPart,size2)
  }
  # Strain?
  if(select_prop[4]){
    out$strain = matrix(0,NPart,size2)
  }
  # Plastic Strain?
  if(select_prop[5]){
    out$pstrain = matrix(0,NPart,size2)
  }

  # work energy
  if(select_prop[6]){
    out$work= matrix(0,NPart,1)
  }
  # plastic energy
  if(select_prop[7]){
    out$penergy= matrix(0,NPart,1)
  }

  #total  shear strain
  if(select_prop[8]){
    out$totalshearstrain= matrix(0,NPart,2)
  }

  # strain energy
  if(select_prop[9]){
    out$senergy= matrix(0,NPart,1)
  }

  # history
  if(select_prop[10]){
    out$history= matrix(0,NPart,num_of_hist)
  }


  # damagenormal
  if(select_prop[11]){
    out$dnormal= matrix(0,NPart,size1)
  }
  ##  Loop through and read

  for(k in 1:NPart){

    temp = readBin(con,"raw",n=4,size=1) # Discard this part
    # Mass?
    if(select_prop[1]){
      out$mass[k] = readBin(con,"double",n=1,size=8)
    }else{
      temp = readBin(con,"raw",n=8,size=1) # Discard this part
    }
    # Material number always read
    out$material[k] = readBin(con,"int",n=1,size=2)
    # Discard this part
    if(size1 == 2){
      temp1  = readBin(con,"raw",n=18,size=1)
    }else{
      temp1  = readBin(con,"raw",n=26,size=1)
    }

    # Position
    out$position[k,] = readBin(con,"double",n=size2,size=8)

    # Velocity?
    if(select_prop[2]){
      out$velocity[k,]  = readBin(con,"double",n=size1,size=8)
    }else if(FormatExpanded[2]){
      temp2 = readBin(con,"raw",n=size1d)
    }
    # Stress?
    if(select_prop[3]){
      out$stress[k,] = readBin(con,"double",n=size2,size=8)
    }else if(FormatExpanded[3]){
      temp3 = readBin(con,"raw",n=size2d)
    }
    # Strain?
    if(select_prop[4]){
      out$strain[k,] = readBin(con,"double",n=size2,size=8)
    }else if(FormatExpanded[4]){
      temp4 = readBin(con,"raw",n=size2d)
    }
    # Plastic Strain?
    if(select_prop[5]){
      out$pstrain[k,] = readBin(con,"double",n=size2,size=8)
    }else if(FormatExpanded[5]){
      temp5 = readBin(con,"raw",n=size2d)
    }
    # work energy
    if(select_prop[6]){
      out$work[k] = readBin(con,"double",n=1,size=8)
    }else if(FormatExpanded[7]){
      temp6 = readBin(con,"raw",n=8,size=1)
    }
    if(FormatExpanded[8]){
      temp7 = readBin(con,"raw",n=8,size=1) # discard temperature
    }
    # plastic energy
    if(select_prop[7]){
      out$penergy[k] = readBin(con,"double",n=1,size=8)
    }else if(FormatExpanded[9]){
      temp8 = readBin(con,"raw",n=8,size=1)
    }

    #total  shear strain
    if(select_prop[8]){
      out$totatlshearstrain[k,] = readBin(con,"double",n=2,size=8)
    }else if(FormatExpanded[11]){
      temp9 = readBin(con,"raw",n=16,size=1)
    }
    # strain energy
    if(select_prop[9]){
      out$senergy[k] = readBin(con,"double",n=1,size=8)
    }else if(FormatExpanded[12]){
      temp10 = readBin(con,"raw",n=8,size=1)
    }

    # history
    if(select_prop[10]){
      out$history[k,] = readBin(con,"double",n=num_of_hist,size=8)
    }else if(any(historyoptions==StoredFlag[13])){
      temp11 = readBin(con,"raw",n=sizehv,size=1)
    }
    # Particle concentration
    if(FormatExpanded[14]){
      temp12 = readBin(con,"raw",n=sizecd,size=1)
    }
    # Particle heat energy
    if(FormatExpanded[15]){
      temp13 = readBin(con,"raw",n=8,size=1)
    }
    # Particle element crossings
    if(FormatExpanded[16]){
      temp14 = readBin(con,"raw",n=4,size=1)
    }
    # Rotational Strain
    if(FormatExpanded[17]){
      temp15 = readBin(con,"raw",n=sizerd,size=1)
    }
    # Damage normal
    if(select_prop[11]){
      out$dnormal[k,] = readBin(con,"double",n=size1,size=8)
    }else if(FormatExpanded[18]){
      temp16 = readBin(con,"raw",n= size1d,size=1)
    }
  }

  # Close connection and output
  close(con)
  return(out)
}


