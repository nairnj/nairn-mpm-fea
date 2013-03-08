This version of NairnMPM includes:

Rotation of a cylinder(3D) or circle (2D)
Dislocation Based material model
work based frictional heating



To upload modified documents to the branch:

Move into the nairn-mpm-fea-Fagan folder (using cygwin on windows). Then type

svn commit -m "short note about the recent changes"


To include additional files, be sure to add them before submitting by:

svn add "/path to file/file to add"

To make the executable, use the makefile by typing:

make SYSTEM=cygwin

for a windows machine, or type

make SYSTEM=gpu 

for the gpu cluster on CSIRO