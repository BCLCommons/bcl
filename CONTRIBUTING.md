-------------------------------------------------------------------------------------

How to Contribute to BCL Development
-------------------------------------------------------------------------------------

We welcome anyone who wants to contribute to the project, whether by adding features, fixing bugs, or improving
the documentation. We appreciate that folks who use the code may want to give back to its development.

First, begin by opening an issue on the BCL Commons Github page that describes the change you want to make. 
This gives everyone a chance to discuss it before you put in a lot of work. For bug fixes, we will confirm 
that the observed behavior is actually a bug (as opposed to an intended behavior or user-specific error) and
that the proposed fix is correct. For new features, we will check whether the desired functionality can be accomplished 
with existing code (i.e. to try and prevent redundancy). If not, and if the proposed feature has the potential to benefit
the BCL Commons community, then we are happy to help discuss possible designs for it.

Once everyone is in agreement, the next step is to create a pull request with the code changes. For larger features,
feel free to create the pull request even before the implementaton is finished so as to get early feedback on the
code. When doing this, tag the pull request with "not ready for review". This allows everyone to track progress on the pull 
request and provide feedback. 

The core developers will review the pull request and may suggest changes. Push the changes to the branch that
is being pulled from, and they will automatically be added to the pull request. Core developers will run the 
full set of BCL example checks to validate the implementation. You may also be asked to write a new example for 
the feature you are adding. We are happy to help as needed. Once the tests are passing and everyone is 
satisfied with the code, the pull request will be merged.

-------------------------------------------------------------------------------------

Formatting and development in the BCL
-------------------------------------------------------------------------------------
BCL formatting and naming guidelines are detailed in documentation/Formatting.pdf

BCL Development is an active an ongoing process.  Tutorials for using some of the core BCL classes and functionality
are available on the BCL Commons GitHub page or upon request (send an email to bcl-academic-support@meilerlab.org) to 
parties with interest in developing additions or enhancements to the BCL.

-------------------------------------------------------------------------------------

Helper / External code-related applications
-------------------------------------------------------------------------------------

Applications that can be called using eclipse external tools configurations (Run -> ExternalTools) or on the command
line.

1. BCL Code formatting
	Eclipse name: FixBCLGuidelinesViolations
	Command: Report issues only - ./scripts/code/BclCleaner.sh ./
					 Report and fix unambiguous issues - ./scripts/code/BclCleaner.sh ./ -o
	How to use in eclipse: Click the bcl folder (or any file thereof) in
												 eclipse, then click run -> external tools -> FixBCLGuidelinesViolations
	What it does - fixes operator spacing according to BCL guidelines. Sorts headers, fixes header guards, warns about
								 some issues that it cannot unambiguously fix.  Also checks doxy comment links between classes and unit
								 tests (examples), but these require manualy fixing.
	What it does not do - fix excessively long lines, check variable naming conventions, fix all indentation issues, find
												seg faults, make coffee, etc.
	Similar command: ReportBCLGuidelinesViolations; runs the same script, but does not make any changes, just outputs the
									 issues that around found to the screen with the recommended fixes
2. Remove Unnecessary headers from CPPs
	Eclipse name: RemoveUnnecessaryHeadersFromCPPs
	Command: scripts/deheader/run_deheader {selected files or folders}
	How to use in eclipse:
		Select any set .cpp files (usually just those you've been working on), then click
		Run -> External Tools -> RemoveUnnecessaryHeadersFromCPPs
	What it does - exhaustively tries removing each include header from each source file, recompiles, checks for errors.
								 Discards includes that are not required for compilation, which may ultimately speed up compilation and
								 reduce interdependencies
3. RemoveCMakeCaches
	Eclipse name: RemoveCMakeCaches
	Command: /bin/rm -f ./build/*/CMakeCache.txt ./build/*/Makefile
	How to use in eclipse: run -> external tools -> RemoveCMakeCaches
	What it does - removes cmake caches; this forces cmake to rerun.  This is useful if, say, you update your environment
								 to have a new gcc or distcc, or add a new library.  CMake will only check for such changes if you blow
								 away its cache files; which this command makes it convenient.  Rarely, cmake may also miss the fact
								 that it needs to recompile a file because another file changed (typically with static instantiations
								 that are never used).
4. CreateCMakeFilesForAllConfigurations
	Eclipse name: CreateCMakeFilesForAllConfigurations
	Command: ./cmake/scripts/create_cmake_lists_arch_if_necessary.csh
	Usage in eclipse: run -> external tools -> CreateCMakeFilesForAllConfigurations
	What it does - used when testing changes to the core cmake list build setup primarily; to ensure that cmake still runs
								 properly against all the primary configurations used in the lab.
