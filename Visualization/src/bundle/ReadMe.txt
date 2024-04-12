The NairnFEAMPMViz Bundle Folder
--------------------------------
This "bundle" folder is used to store resources such that development can test the $(bundle) options while running Java app within Eclipse. Basically the src folder becomes a jar file folder and therefore this bundle folders acts like it is the same folder at the runnable jar file.

Windows
-------

The items needed in this folder when developing in Windows are:

1. NairnFEA.exe and NairnMPM.exe compiled code engines
2. NairnFEA.dtd and NairnMPM.dtd XML format definition files
3. xerces-c_###.dll
4. vcruntime###.dll
5. vcomp###.dll
6. msvcp###.dll
7. ExtractMPM.exe to test features that use this tool

where "###" indicates a version number such as "3_1" for xerces or "140" for dlls from Visual Studio 14. Note that the required (and redistributable) dlls are in the folder

   Microsoft Visual Studio ##\VC\redist\x64

In Windows, this bundle folder is only intended for testing during developing. See the CreateExePackage.txt file for the best way to create a standalone package for installation on other Windows PCs.

MacOS
-----

The items needed when developing on a Mac are:

1. NairnFEA and NairnMPM compiled code engines
2. NairnFEA.dtd and NairnMPM.dtd XML format definition files
3. ExtractMPM to test features that use this tool

Note that these files no not include the xerces library or c++ libraries. I was not able to get MacOS to use xerces in this bundle folder and when using clang for OpenMP development, the engines need special c++ libraries. As a result MacOS development needs prior installs of xerces and clang-mp.

A downloadable Mac version of NairnFEAMPMViz can be provided with the following folder:

     NairnFEAMPMViz
     ... NairnFEAMPMViz.jar (exported runnable jar file)
     ... ReadMe.txt (if desired to brief instructions)
     ... bundle
     ....... (all bundle files listed above)

Like above, someone that downloads such a folder can only run it if their Mac has prior installs of xerces and clang-mp. A better option for Mac users is to download NairnFEAMPM and follow its instructions for setting up the system.


