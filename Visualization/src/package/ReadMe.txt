--------------------------------
NairnFEAMPM Package for Windows
--------------------------------

This folder contains a complete and compiled nairn-fea-mpm computational mechanics package. It can run material point method (MPM) calculations using NairnMPM or finite element analysis (FEA) using NairnFEA on a Windows PC. The calculations are set up and visualized using the main NairnFEAMPMViz.exe Windows application. The included PublishPlotJ.exe file can be used to reopen any 2D plot saved with the "Save Plot File..." command in NairnFEAMPMViz.exe for later editing. Editing in PublishPlotJ.exe also provides additional options for customizing your ploted results.

--------------------------------
Quick Start
--------------------------------

This downloaded and extracted NairnFEAMPMViz folder can be moved to any location on your computer. A good choice is to move it to the "Program Files" or "Program Files (x86)" folder.

To start the application, double click on the NairnFEAMPMViz.exe file. A convenient option for future runs is to right click on NairnFEAMPMViz.exe, choose menu command "Create short cut," and then move that short cut to your desktop. Double clicking the new short cut will start the application.

Once the application is running, you can learn to use it by opening the help window using the "Help" menu command. To get some immediate results, you can choose one of the examples in the File->Examples submenu. After picking an example, a command file will open. Save that new file to a folder and then choose the "Run MPM/FEA Analysis" menu command to run an example calculation. The examples or MPM calculations except those with "FEA" in their name.

For more detailed documentation, consult the OSUPDocs web site at:

https://osupdocs.forestry.oregonstate.edu/index.php/Main_Page

--------------------------------
Other Information
--------------------------------

This application is written in Java but provided here wrapped into an exe file using the Launch4j tool. Although it runs as a Windows application, it lacks some features of native Windows application. The main limitation is that the only way to start the application is to double click NairnFEAMPMViz.exe (or double click a short created as explained above). You cannot double click files created by the application. You have to open files using the File menu after starting the application.

The key files types for working with NairnFEAMPMViz.exe by extension are:

1. fmcmd files - these are command files for setting up calculations or for writing scripts.
2. mpm files - these are output summaries of MPM calculations. The output will also incLude a folder of more results.
3. fea files - these are output of FEA calculations.

