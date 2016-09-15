# Chad Hammerquist
# This file combines to movies into single frames with an optional spacer

# Get packages
require(tiff)
require(png)
require(abind)

# parameters

# path to parent folder
parent = "/Volumes/Numerical\ Modeling/Modeling_Files/Nairn_Research/MPM/Water\ Column"

# path from parent to folder with frames
frames1 = "MPM+PS-FLIP/Frames/"
frames2 = "MPM+PS-XPIC2/Frames/"

# path from frames back to parent (same for both)
backup = "../.."

# output folder (must exist, will overwrite)
output = "movie/stitch_"

# file to place between the images
divide = "scale.png"

# combine horizontal (2) or vertical (1)
axis = 2

# Which pixels to include

# horizontal direction
hlim = 363:735     # images (same for both)
dlim = 1:60        # divider

# vertical direction with top=0 and bottom=max
vlim = 245:695

# Set working directory
setwd(parent)
 
# Go into first frames and get list of files
setwd(frames1)
files1=list.files()
setwd(backup)
 
# Go into secopnd frames and get list of files
setwd(frames2)
files2=list.files()
setwd(backup)

 # Read space picture
spacer = readPNG(divide)

 # Loop through and read file
for(k in 1:length(files1)){
  # read
  part1 = readTIFF(paste(frames1,files1[k],sep=""))
  part2 = readTIFF(paste(frames2,files2[k],sep=""))
  
  # concatenate RGB arrays of the pictures along the second margin
  combine = abind(part1[vlim,hlim,1:3],spacer[vlim,dlim,1:3],part2[vlim,hlim,1:3],along = axis)
  
  # write out
  writePNG(combine,target = paste(output,k,".png",sep=""))
  
}