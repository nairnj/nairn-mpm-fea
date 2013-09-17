This version of NairnMPM includes:

Rotation of a cylinder(3D) or circle (2D)
Dislocation Based material model
work based frictional heating
hard coded heat input


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


To download the code to a new location:

svn checkout https://nairn-mpm-fea.googlecode.com/svn/branches/Fagan nairn-mpm-fea-Fagan --username timcfagan@gmail.com.

To extract the ExtractMPM file go to NairnMPM/tools then type:

g++ -o ExtractMPM ExtractMPM.cpp

or if on cluster

g++ -o ExtractMPM ExtractMPM_Cluster.cpp





To extract particle data using the new ExtractMPMUltimate6, type:

./ExtractMPMUltimate6 -V -q hist -o <file name to create> -s <file extension>


The -q gives the data to show, which must be stored to be converted. 
The options can be temp, mat, hist etc look in documentation for others. 
The variable I've added to this Extract version is "hist" for the stored history variables. View these vtk output files in paraview.