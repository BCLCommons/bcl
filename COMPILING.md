Document Notes
-------------------------------------------------------------------------------------

This document summarizes important points about developing and building the BCL.  For more detailed information
see additional files in the documentation/ directory.  For higher-level code documentation for the BCL, 
some documents have been written and placed in the include/*/ directories to give an overview of the classes in
each namespace.

Introduction to building, compiling, and developing in the BCL
-------------------------------------------------------------------------------------

Pre-requisites:  
GCC/G++ 4.6.3 or greater installed and in your PATH (clang++ 4.2 or above on apple; mingw compiler on windows)  
CMake 2.8.2 or higher  
Python 2.4 or higher (2.7 recommended)  

External Packages and Attribution:
-------------------------------------------------------------------------------------

External packages included in source distribution, that are not developed by the MeilerLab:  
- zlib     http://zlib.net/ - compression library  
- bzip2    http://www.bzip.org/ - compression library  

Linux64 additional packages:  
- Freeocl (http://code.google.com/p/freeocl/) used for OpenCL GPU acceleration on Linux  
	 - Installation steps: download freeocl package (0.3.6 or later), compile  
	 - install headers into extern/noarch/freeocl/0.3.6/include  
	 - mkdir -p extern/noarch/freeocl/0.3.6/etc/OpenCL/vendors/  
	 - echo "libFreeOCL.so" > extern/noarch/freeocl/0.3.6/etc/OpenCL/vendors/freeocl.icd  

Win32 additional optional packages:
1. ATI/AMD; allows use of OpenCL GPU acceleration functionality on Windows
	 - Installation steps: download AMD SDK, install headers into extern/noarch/ati/2.5/include
	 - Compile and install OpenCL.lib into extern/windows/win32/ati/2.5/lib/OpenCL.lib
2. pthreads (http://sourceware.org/pthreads-win32/); allows use of POSIX threads on Windows
	 - Copy headers (pthread.h, sched.h, and semaphore.h) into extern/noarch/pthreads/2.8.0/include/
	 - Compile/install pthreadGC2-w32.dll into extern/win32/pthreads/2.8.0/bin/mingw

Internally we also use libmysql and libmysqlpp, but they are not included the source release due to the need for
additional, site-specific configuration and lack of necessity for all but the most demanding cheminformatics applications.
If you have a strong interest in getting such functionality into the BCL, please let us know; we have some interest in
migrating to SQL-lite but knowing that others would be interested in this too might provide the impetus needed to get it
moving forward.

Development Environment - Eclipse
-------------------------------------------------------------------------------------
Eclipse (http://www.eclipse.org/) is the primary IDE used for BCL development and is the easiest way for new developers
to begin building and using the BCL. This section assumes you have installed Eclipse with CDT (C++ dev package), if not,
follow the steps at http://www.eclipse.org/ first.

Importing the BCL into an eclipse workspace (external users):
1. Open eclipse
2. Select File -> Import -> General -> Existing Projects Into Workspace
3. Choose the root directory containing the downloaded bcl source (same directory as this ReadMe is located in)
4. Click Finish

Importing the BCL into an eclipse workspace from GitHub:
1. Request access to a github fork of the BCL project or create a fork from https://github.com/BCLCommons/bcl
2. Open eclipse
3. Select File -> Import -> Git -> Projects from Git and click next
4. Select Clone URL and click next
5. Type the URL of your fork under URL
6. Type your username and password and click next
7. On the branch selection page click next
8. Select the destination directory and click next
9. Allow the checkout to process and select import existing Eclipse project and click next
10. Click Finish

Building the BCL for the first time (Linux64):
1. Select Project -> Build Configurations -> Set Active -> linux64_release
2. Build the project (Project -> Build All, or Ctrl + B)

Building the BCL for the first time (Win32):
1. Edit compiler names in cmake/toolchain/i686-w64-mingw32.cmake for your configuration
	 - usually this just means removing the ${TARGET_PLATFORM}- prefix for all executables in that file
2. Select Project -> Build Configurations -> Set Active -> win32_release
3. Build the project (Project -> Build All, or Ctrl + B)

Building the BCL for the first time (Apple64):
1. Edit compiler names in cmake/toolchain/x86_64-apple-darwin10.cmake for your configuration
	 - usually this just means removing the ${TARGET_PLATFORM}- prefix for all executables in that file, and changing
		 gcc to clang, g++ to clang++
2. Select Project -> Build Configurations -> Set Active -> apple64_release
3. Build the project (Project -> Build All, or Ctrl + B)


Where are my executables?
-------------------------------------------------------------------------------------
After building using eclipse, bcl executables are located in:  
build/{arch}_(debug|release)/bin/bcl-(example|apps|example-apps)-static.exe ex. build/linux64_release/bin/bcl-apps-static.exe  
bcl-apps-static.exe - contains all applications, without unit tests  
bcl-example-static.exe - unit tests for internal classes only, run with bcl-example-static.exe Examples  
bcl-example-apps-static.exe - applications and unit tests for applications  


Build Configuration Syntax
-------------------------------------------------------------------------------------
{arch}(_apps?)(_debug|release) ex. linux64_apps_release  
{arch} is the system architecture; currently win32, apple64, and linux64 are supported  
apps if present; only build bcl-apps-static.exe; otherwise, build applications and examples (unit tests)  
debug - builds debuggable (with gdb or similar) versions; this is usually much slower and requires a lot more memory; prefer release for all but the hardest-to-find segfaults  
release - builds executables without debug symbols, at optimization level O2  


Building shared libraries
-------------------------------------------------------------------------------------
On clusters or for use of the BCL as an external library in other applications, it is useful to build the BCL as a
shared library.  In eclipse, this can be done by modifying your build configuration as follows:
1. Select Project -> Properties -> Build Configurations
2. Click C++ Build from the window
3. Select the configuration that you'd like to build a shared library for next to the Configuration option
4. Click the behavior tab
5. Next to build, there is a field with the text "static".  Change "static" to "shared"
6. Build the BCL; the result libbcl.so is in build/linux64_release/lib/libbcl.so
7. The applications in build/linux64_release/bin/ now lack the -static.exe because they depend on libbcl.so


Development Environment - Non-Eclipse
-------------------------------------------------------------------------------------
Any IDE that can run arbitrary commands as part of the build process can be used to build the BCL, because ultimately
the building is left up to the BCL's internal build scripts, cmake, and the selected toolchain.  Here are the sequence
of build steps for linux64_release and what each command does.
1. Open a shell and `cd` into the bcl directory
2. Set environment variables: TOOLCHAIN_ALIAS=linux64, TOOLCHAIN=x86_64-unknown-linux-gnu, BUILD_TYPE_ALIAS=release,
	 BCL_DIRECTORY=`pwd`
3. make directory build/${TOOLCHAIN_ALIAS}_${BUILD_TYPE_ALIAS} _if it does not already exist. On linux:_  
	 `mkdir -p build/${TOOLCHAIN_ALIAS}_${BUILD_TYPE_ALIAS}`  
4. `./cmake/scripts/CheckCmakeLists.py ./ -o`  
	 -- this command updates the CMakeLists files with any new source files that have been added or removed.  It does not
			automatically add new folders under example/ to example/CMakeLists or from source/ to source/CMakeLists.txt; nor
			does it create new namespace CMakeLists.  In short, if you create a new namespace or folder, you'll need to
			manually add it to the associated CMakeLists.txt files in the next-level up directories.  
5. `./scripts/code/CreateNamespaceForwardHeaders.py ./ -o`  
	 -- this command updates the .fwd.hh file in each namespace folder (include/*/) _with any new classes that have been
			declared.  It may also create a bcl_(namespace).depends.fwd.hh for cross-namespace dependencies.  These all serve
			to minimize the number and scope of interdepencies and number of files that need to be recompiled when classes
			are added or changed._  
6. `./scripts/build/check_pump_make.csh -k static`  
	 -- this command runs cmake (if needed) and make with the arguments given.  If distcc is available, it runs the commands
			under pump with as many processes as your distcc server has declared available.

Notes:
- In step 6, you may change static to shared if shared objects are desired.
- Win32 -- change TOOLCHAIN_ALIAS to win32, TOOLCHAIN to i686-w64-mingw32
- Mac OS X -- change TOOLCHAIN_ALIAS to apple64, TOOLCHAIN to x86_64-apple-darwin-10

Disclaimer for Windows users:
The above toolchain has only been tested on Linux and (to a lesser extent) Mac OS X.  Minor tweaks may be necessary
in Cygwin; and substantial re-writing of the build system may be necessary for use on non-Cygwin based Windows.


Building using Command line
-------------------------------------------------------------------------------------
There are scripts that can be run to assist in building the BCL in the following directory:  
	scripts/build/
For example, to build _on linux:_  
	`cd bcl/`  
	`./scripts/build/build_cmdline.linux.sh`
