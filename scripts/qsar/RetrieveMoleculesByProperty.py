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
  print "RetrieveMoleculesByProperty.py input_file output_file property_name"
  sys.exit(1)

def main():
  arguments = sys.argv[1::1]

  if len(arguments) < 3:
    usage()

  property = '> <' + arguments[2] + '>'

  input_filename = arguments[0]
  output_file = arguments[1]

  file = None
  if input_filename.endswith('.bz2'):
    file = bz2.BZ2File(input_filename, 'r')
  elif input_filename.endswith('gz'):
    file = gzip.GzipFile(input_filename, 'r')
  else:
    file = open(input_filename, 'r')
  
  current_mol = ""
  property_val = ""
  vals_seen_before = set()
  while 1:
    line = file.readline()
    if line is None or len(line) == 0:
      break
    current_mol += line
    if line.startswith(property):
      property_val = file.readline()
      current_mol += property_val
      property_val = property_val.strip().replace(' ','_').replace('\t','_')
    elif line == "$$$$\n":
      open_mode = 'a'
      if property_val not in vals_seen_before:
        open_mode = 'w'
        vals_seen_before.add(property_val)
      ofile = open(output_file + property_val + '.sdf', open_mode)
      ofile.write(current_mol)
      ofile.close()
      current_mol = ""
  file.close()

if __name__ == '__main__':
    main()
