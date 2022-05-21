Installation Steps:
1. Copy the BCL folder from the disk image to a convenient location, for example, /Applications
2. Place the license (obtainable from http://www.meilerlab.org/servers/bcl-academic-license) for the bcl applications in 
   the BCL folder (Applications/BioChemicalLibrary-3.1.0)
3. (Optional) If you want to be able to run the bcl from any directory, add it to your path. To do this, open a terminal 
   (/Applications/Utilities/Terminal), and then type the following command
   echo "export PATH=/Applications/BioChemicalLibrary-3.1.0:${PATH}" >> ~/.bashrc
   
The remaining sections assume that you are either in the bcl directory from a terminal, or that you've added the bcl to
your path (by following step 3, above)

Usage:
  bcl <application_group>:<application> [Other parameters] [Flags] [@filenames]
  see bcl -help for more information about the application groups and application
  Note: Be patient! The first time it is loaded, the bcl may take several seconds or even a minute to start up. 
     
Help/Usage:
  To learn more about an application, including example command lines, and learn what you must cite if you use it in 
  your research, type: 
    bcl ApplicationName -readme
  To learn about the options (flags/parameters) for e.g. Jufo, type 
    bcl Jufo -help
  To see a list of all applications and a brief description of what they do, type: 
    bcl all:Help
  "help" for any flag or parameter to obtain more help for just that flag or parameter, e.g.
    bcl Jufo help 
    OR
    bcl Jufo -output help
    
Submitting bug reports:
  Bug reports and licensing questions may be forwarded to bcl-academic-support@meilerlab.org.  When submitting a bug 
  report, please include the following information:
  1. BCL version number (if not sure, type bcl.exe -help | grep "BCL v")
  2. Platform 
  3. Mac OS version (Click the Apple icon in the upper left hand corner of your menu bar and select About this Mac)
  5. The command you tried to run
  6. What you expected to happen when running that command (if it was not an error)
  7. What actually happened (error messages, crash, no output, etc.)
  8. If any log messages were written, please send us the log file, e.g. 
     bcl SomeCommand ... -logger File log.txt 
     then attach log.txt to the output.  If it is larger than 1MB, please compress it first using gzip, bzip2, or tar

Known Bugs / Workarounds:
   If you are getting any error messages about opencl and want to remove them; or the bcl is crashing when you run
   multiple bcl instances on the same machine, then try these steps
    1. Try running the command with -opencl Disable set (e.g. bcl.exe Fold ... -opencl Disable) 
    2. As a last resort, try renaming the opencl_kernels folder (so that the kernels are not found), e.g. 
       mv opencl_kernels opencl_kernels.back
       If the problem vanishes but you need opencl support, submit a bug report to bcl-academic-support@meilerlab.org
       with the make and model of gpu that your machine has.  If you're unsure of this information, follow these steps
       A. Click the Apple icon in the upper left hand corner of your menu bar and select About this Mac
       B. Click More info
       C. Click System Report button (if available)
       D. Click Hardware, then Graphics/Displays tab in the resulting window
          Send us the make and model name (e.g. NVIDIA GeForce GT 330M) of the gpu card
       Also send us the output of 
       bcl ApplicationName -opencl help 
       (where ApplicationName is replaced by the name of the application you were trying to run)
