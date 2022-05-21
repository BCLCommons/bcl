Installation Steps:
  1. Run the installer (bcl-3.1.0-Linux-x86_64.sh) and agree to licensing terms
  2. Place the license (bcl.license obtainable from http://www.meilerlab.org/servers/bcl-academic-license) for
     the bcl applications in the BCL folder 
     -> You may now run the bcl from within the installation folder 
  The next steps are optional, but enable use of the BCL from directories other than the installation directory:
  3. Add the bcl folder to your ld library path, e.g. 
     (tcsh/csh) setenv LD_LIBRARY_PATH /path/to/bcl-3.1.0-Linux-x86_64:${LD_LIBRARY_PATH}
     (bash/sh)  export LD_LIBRARY_PATH=/path/to/bcl-3.1.0-Linux-x86_64:${LD_LIBRARY_PATH}
     replacing /path/to with the path to where the bcl-3.1.0 folder was installed
  4. Changing your LD_LIBRARY_PATH manually only allows the current shell to run the bcl; to avoid needing to repeat 
     step 3 with every new shell you want to use the bcl with, append the setenv command from step 3 to your ~/.cshrc file 
     if you're using csh or tcsh as your default shell.  For example:
       echo "setenv LD_LIBRARY_PATH /path/to/bcl-3.1.0-Linux-x86_64:${LD_LIBRARY_PATH}" >> ~/.cshrc
     If you're using bash/sh, append the export command to your ~/.bashrc file instead, e.g.
       echo "export LD_LIBRARY_PATH=/path/to/bcl-3.1.0-Linux-x86_64:${LD_LIBRARY_PATH}" >> ~/.bashrc
     -> You may now run the bcl from any directory by typing the full path to the bcl.exe
  If typing the full path to the bcl gets old, you may also add it to your path variable, for example,
  5. (in tcsh/csh shell) echo "setenv PATH /path/to/bcl-3.1.0-Linux-x86_64:${PATH}" >> ~/.cshrc
     (in bash/sh shell) echo "export PATH=/path/to/bcl-3.1.0-Linux-x86_64:${PATH}" >> ~/.bashrc
     Then type source ~/.cshrc (or .bashrc) and 
     -> you can just type bcl.exe from any directory to run the bcl
     
The remaining sections assume that you are either in the bcl directory from a terminal, or that you've added the bcl to
your path (by following steps 1-5, above)

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
  3. Linux distribution (command line: cat /etc/*-release)
  4. Linux kernel version (command line: uname -a)
  5. The command you tried to run
  6. What you expected to happen when running that command (if it was not an error)
  7. What actually happened (error messages, crash, no output, etc.)
  8. If any log messages were written, please send us the log file, e.g. 
     bcl.exe SomeCommand ... -logger File log.txt 
     then attach log.txt to the output.  If it is larger than 1MB, please compress it first using gzip, bzip2, or tar

Known Bugs / Workarounds:
   If you are getting any error messages about opencl and want to remove them; or the bcl is crashing when you run
   multiple bcl instances on the same machine, then try these steps
    1. Try running the command with -opencl Disable set (e.g. bcl.exe Fold ... -opencl Disable) 
    2. As a last resort, try renaming the opencl_kernels folder (so that the kernels are not found), e.g. 
       mv opencl_kernels opencl_kernels.back
       If the problem vanishes but you need opencl support, submit a bug report to bcl-academic-support@meilerlab.org
       with the output of the following commands:
       ls /etc/OpenCL/vendors/
       cat /etc/OpenCL/vendors/*.icd
       bcl.exe ApplicationName -opencl Help 
       (where ApplicationName is replaced by the name of the application you were trying to run)
