Installation Steps:
  1. Run the installer and agree to licensing terms
     -> For convenience, we recommend enabling options to add the bcl to your path 
  2. Place the license (bcl.license obtainable from http://www.meilerlab.org/servers/bcl-academic-license) for
     the bcl applications in the BCL folder (by default C:\Programs (x86)\bcl_3.1\)
     -> You may now run the bcl from within the installation folder or by typing the full path to the bcl from any other 
     directory, using Windows Powershell or any other terminal emulator.
   
The remaining sections assume that you are either in the bcl directory from a terminal, or that you've added the bcl to
your path (by allowing the installer to do so, in step 1)

Usage:
  bcl.exe <application_group>:<application> [Other parameters] [Flags] [@filenames]
  see bcl.exe -help for more information about the application groups and application
  Note: Be patient! The first time it is loaded, the bcl may take several seconds or even a minute to start up. 
     
Help:
  To learn more about an application, including example command lines, and learn what you must cite if you use it in 
  your research, type: 
    bcl.exe ApplicationName -readme
  To learn about the options (flags/parameters) for e.g. Jufo, type 
    bcl.exe Jufo -help
  To see a list of all applications and a brief description of what they do, type: 
    bcl.exe all:Help
  "help" for any flag or parameter to obtain more help for just that flag or parameter, e.g.
    bcl.exe Jufo help 
    OR
    bcl.exe Jufo -output help
  
Submitting bug reports:
  Bug reports and licensing questions may be forwarded to bcl-academic-support@meilerlab.org.  When submitting a bug 
  report, please include the following information:
  1. BCL version number (if not sure, type bcl.exe -help | grep "BCL v")
  2. Platform 
  3. Windows version (if unsure, follow the directions at http://windows.microsoft.com/en-us/windows/which-operating-system)
  4. The command you tried to run
  5. What you expected to happen when running that command (if it was not an error)
  6. What actually happened (error messages, crash, no output, etc.)
  7. If any log messages were written, please send us the log file, e.g. 
     bcl.exe SomeCommand ... -logger File log.txt 
     then attach log.txt to the output.  If it is larger than 1MB, please compress it first into a zip file
