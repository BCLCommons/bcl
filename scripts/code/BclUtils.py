'''
Created on Apr 14, 2010
@brief Utilities for bcl code analysis and editing
@author: mendenjl
'''

import string
import sys
import os
from CodeFileUtils import stripComments, getFilesFromDirectory

def BCLNameFromFilename(line):
  split_line = line.split('_')
  if len(split_line) < 2:
    return line

  new_name = split_line[0] + "::"
  split_line = split_line[1::1]

  for index in xrange(len(split_line)):
    if str(split_line[index][0]).isdigit():
      split_line[index] = split_line[index].upper()
    else:
      split_line[index] = split_line[index].capitalize()

  new_name += ''.join(split_line)
  return new_name

def getDependencies(knownHeaders):
  knownFilesSet = set(knownHeaders)
  knownHeaders = sorted(knownFilesSet)
  pos = 0
  while pos < len(knownHeaders):
    headerFile = knownHeaders[ pos]
    #print "scanning " + headerFile
    if os.path.exists(headerFile) and os.path.isfile(headerFile):
      ifile = open(headerFile, 'r')
      lines = ifile.readlines()
      lines = stripComments(lines)
      ifile.close()
      for line in lines:
        if line.startswith("#include \""):
          line = line[10:].strip('" \t.' + os.sep)
          #line = line.replace('.fwd.hh', '.h').strip()
          outer_directory = ""
          directory = ""
          filename = ""
          #print "found include " + line
          if line.find(os.sep) > 0:
            directories_and_file = line.split(os.sep)
            filename = directories_and_file[ -1]
            directory = directories_and_file[ -2] + os.sep
            if len(directories_and_file) > 2:
              if(line.startswith("example_")):
                outer_directory = "example" + os.sep
              outer_directory += (os.sep).join(directories_and_file[0:-2:1])
            elif line.endswith('.cpp'):
              outer_directory = "source"
            else:
              outer_directory = "include"
          elif line.startswith("bcl."):
            filename = line
            if filename.endswith(".cpp"):
              outer_directory = "source"
            else:
              outer_directory = "include"
          elif line.startswith("bcl_app"):
            filename = line
            if line != "bcl_app_examples.cpp":
              directory = "apps" + os.sep
            else:
              directory = "example" + os.sep
          elif line == 'apps.cpp':
            directory = "apps" + os.sep
            filename = line
          elif line.startswith("bcl_"):
            # directory could be include or source
            directory = line.split("_")[ 1].split('.')[0] + os.sep
            filename = line
            if filename.endswith(".cpp"):
              outer_directory = "source"
            else:
              outer_directory = "include"
          elif line.startswith("example_"):
            directory = "example" + os.sep
            filename = line
          elif line.startswith("example."):
            directory = "example" + os.sep
            filename = line
          else:
            continue
          new_filename = outer_directory + os.sep + directory + filename
          new_filename = new_filename.strip(os.sep)
          if new_filename not in knownFilesSet:
            #print "will scan " + new_filename
            knownHeaders.append(new_filename)
            knownFilesSet.add(new_filename)

          if filename.endswith(".h"):
            forward_header = new_filename.replace('.h', '.fwd.hh')
            if forward_header not in knownFilesSet:
              knownHeaders.append(forward_header)
              knownFilesSet.add(forward_header)
          elif filename.endswith(".fwd.hh"):
            hdr = new_filename.replace('.fwd.hh', '.h')
            if hdr not in knownFilesSet:
              #print "will scan " + new_filename
              knownHeaders.append(hdr)
              knownFilesSet.add(hdr)

          # add the cpp as well
          if outer_directory == "include":
            outer_directory = "source"
          elif outer_directory == "source":
            outer_directory = "include"
          if filename.endswith(".h"):
            filename = filename.replace(".h", ".cpp")
          elif filename.endswith(".fwd.hh"):
            filename = filename.replace('.fwd.hh', '.cpp')
          elif filename.endswith(".cpp"):
            filename = filename.replace(".cpp", ".h")
          new_filename = outer_directory + os.sep + directory + filename
          new_filename = new_filename.replace('//', '/').strip(os.sep)
          if new_filename not in knownFilesSet:
            knownHeaders.append(new_filename)
            knownFilesSet.add(new_filename)
            #print "will scan " + new_filename
    elif os.path.exists(headerFile) and os.path.isdir(headerFile):
      #print "scanning directory " + headerFile
      for suffix in ['cpp', '.h', '.fwd.hh']:
        more_headers = getFilesFromDirectory(headerFile, suffix)
        for file in more_headers:
          file_str = str(file)
          if file_str not in knownFilesSet:
            #print "will scan " + file_str
            knownHeaders.append(file_str)
            knownFilesSet.add(file_str)
    pos += 1

  #print "known headers pre os.path.exists: " + '\n'.join(knownHeaders) + "\n\n"
  knownHeaders = [ hdr for hdr in knownHeaders if os.path.exists(hdr) and os.path.isfile(hdr)]
  return knownHeaders
