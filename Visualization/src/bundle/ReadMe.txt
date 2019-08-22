The NairnFEAMPMViz Bundle Folder
--------------------------------
This "bundle" folder is used to store resources that can be distributed with the NairnFEAMPMViz.jar file to create a stand-alone package that can run NairnFEA and NairnMPM simulations on other computers without needing to install any developer tools or to compile any code.

For distribution, create a new folder and move the NairnFEAMPMViz.jar file into that folder. Next, populate this bundle folder with the following, platform-dependent files:

Windows
-------
1. NairnFEA.exe and NairnMPM.exe compiled code engines
2. NairnFEA.dtd and NairnMPM.dtd XML format definition files
3. xerces-c_###.dll
4. vcruntime###.dll
5. vcomp###.dll
6. msvcp###.dll

where "###" indicates a version number such as "3_1" for xerces or "140" for dlls from Visual Studio 14. Note that the required (and redistributable) dlls are in the folder

   Microsoft Visual Studio ##\VC\redist\x64

MacOS
-----
1. NairnFEA and NairnMPM compiled code engines
2. NairnFEA.dtd and NairnMPM.dtd XML format definition files
3. libxerces-c-###.dylib to the lib folder
4. A symbolic link libxerces-c.dylib to the library in #3 to the lib folder

where "###' refers to the compiled xerces version (such as "3.1")

After this bundle folder is fully populated, copy it the entire folder to the same folder with the NairnFEAMPMViz.jar file. That final folder can be copied to other computers and run by launching the NairnFEAMPMViz.jar file. Make sure all code preferences are set to $(bundle) to direct NairnFEAMPMViz to use the code engines in the bundle folder.

Eclipse Development
-------------------
When developing and running NairnFEAMPMViz directly in Eclipse, the bin folder is a stand in for the jar file. Because this bundle folder is in the same folder as the bin folder, it will act just like a bundle folder that is distributed in the same folder as the NairnFEAMPMViz.jar file. If you populate the bundle folder in Visualization\src folder, it can be used to test code development that interacts with the bundle files.