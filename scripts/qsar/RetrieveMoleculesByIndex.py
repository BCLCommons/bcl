#!/usr/bin/python
'''
Created on Nov 17, 2010

@author: mendenjl
'''

import os
import sys
import os.path
import bz2, gzip
from symbol import exec_stmt

def usage():
  print "SplitFileAfter.py input-file output_file ids"
  sys.exit(1)

search_type = 0

def main():
  arguments = sys.argv[1::1]

  if len(arguments) < 3:
    usage()

  ids = []
  try:
    ids = [ int(x) for x in arguments[2:]];
  except:
    usage()

  id_set = set()
  for id in ids:
    id_set.add(id)

  input_filename = arguments[0]
  output_file = arguments[1]

  file = None
  if input_filename.endswith('.bz2'):
    file = bz2.BZ2File(input_filename, 'r')
  elif input_filename.endswith('gz'):
    file = gzip.GzipFile(input_filename, 'r')
  else:
    file = open(input_filename, 'r')
  ofile = open(output_file, 'w')
  line_number = 0
  molecule_number = 0
  are_writing = molecule_number in id_set
  while 1:
    line = file.readline()
    if line is None or len(line) == 0:
      break
    if are_writing:
      ofile.write(line)
    if line == "$$$$\n":
      molecule_number += 1
      are_writing = molecule_number in id_set
    line_number += 1
  file.close()
  ofile.close()

if __name__ == '__main__':
    main()
