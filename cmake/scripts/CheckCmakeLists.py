#!/usr/bin/python
'''
Created on Feb 19, 2011

This script checks CmakeLists for the bcl project, optionally updating CMakeLists.txt with new source files and removing
non-existent files.  

@author: mendenjl
'''

import os
import sys
import os.path
from BclCmakeFile import *

# print usage info 
def usage():
  print "\nusage: CheckCmakeLists.py bcl-path [-o]\n"
  print "-o if given, update the CMakeLists"

def main():

  # do not update by default
  should_update = False

  # get the arguments given by the user
  args = sys.argv[1::1]
  if (len(args) != 1 and len(args) != 2) or args[0] == "h" or args[0] == "help":
    usage()
    sys.exit(1)

  # if there was a second argument and it was -o, update the cmake lists
  if len(args) == 2 and args[1] == '-o':
    should_update = True

  # save the original working directory
  original_directory = os.getcwd()

  # move to the bcl path
  os.chdir(args[0])

  # get a dictionary (map) from directories a list of associated source files
  source_files_and_directories = {}
  # get a dictionary (map) from directories a list of cmake files
  cmake_list_paths = {}
  source_directories = [ "example", "source", "apps"]
  for directory in source_directories:
    for root, dirs, files in os.walk(directory):
      # add source files (.cpp) and generated source files (.cpp.in)
      source_files_and_directories[ root] = [ str(file).strip('.in') for file in files if os.path.isfile(root + os.sep + file) and (file.endswith(".cpp") or file.endswith(".cpp.in")) ]
      if os.path.isfile(root + os.sep + "CMakeLists.txt"):
        cmake_list_paths[ root] = root + os.sep + "CMakeLists.txt"
      if '.svn' in dirs:
        dirs.remove('.svn')  # don't visit svn directories

  for source_directory, directory_sources in source_files_and_directories.items():
    # does the cmake list file even exist yet?
    if source_directory not in cmake_list_paths.keys():
      if len(directory_sources) > 0:
        print source_directory + " needs a CMakeLists.txt file containing:\n" + '\n'.join(directory_sources) + '\n'
    else:
      # open the cmake lists file
      CheckCmakeLists(source_directory, source_files_and_directories, should_update)

  # return to the original directory
  os.chdir(original_directory)

if __name__ == '__main__':
  main()
