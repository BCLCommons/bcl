This folder contains scripts that are one of three things:
1. Building or maintaining the code (Folders: build, code, debug, deheader, nightly_build_helper)
   * build - scripts related to building the bcl on the command line or in eclipse
   * code - analysis and fixing of the bcl codebase itself
2. Common workflows in which the BCL is the integral component
3. Analysis scripts used to generate data, such as histograms, that are used or are primarily useful for the BCL

Rules:
* All scripts belong in a subfolder of scripts, not directly inside scripts
* General utility scripts that are useful outside the BCL and not in any way tied to the BCL do NOT belong here
