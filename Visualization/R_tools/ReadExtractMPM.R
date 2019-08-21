# Chad Hammerquist
# Function to import read MPM
ReadExtractMPM=function(file,extract,mat=0){
  # make prefix
  prefix = "ExtractMPM -T"
  data_option = paste("-q",extract," ",collapse ="")
  if(mat[1]==0){
    mat_option = ""
  }else if(mat[1]<1){
    mat_option = paste(" -m",abs(mat)," ",collapse="")
  }else{
    mat_option = paste(" -M",abs(mat)," ",collapse="")
  }
  # combine into one name
  command = paste(prefix,mat_option,data_option,file,collapse="",sep="")

  # Pipe
  z = pipe(command)
  out = read.table(z)
  return(out)

}
