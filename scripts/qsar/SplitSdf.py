#!/usr/bin/python
'''
Created on Nov 17, 2010

@author: mendenjl
'''

import os, io
import sys
import os.path
import gzip, bz2

def usage():
  print "SplitSdf.py input-file output-file-base [optional: number of molecules per file] {Compression:gz|bz2}"
  sys.exit(1)

def DetermineNumberDigitsInFilename(number_line_numbers):
  i = 1
  while number_line_numbers > 0:
    i += 1
    number_line_numbers /= 10
  return i

def MakeFilenameDigitStr(number_digits, index):
  if index == 0:
    return ''.join([ "0" for i in xrange(number_digits - 1)])
  number_0s = number_digits - DetermineNumberDigitsInFilename(index)
  return ''.join([ "0" for i in xrange(number_0s)]) + str(index)

def main():
  arguments = sys.argv[1::1]

  if len(arguments) > 4:
    usage()

  input_filename = arguments[0]
  search_string = '$$$$\n'
  n = 1
  if len(arguments) >= 3:
    try:
      n = int(arguments[2])
    except:
      usage()
  compression=0
  if len(arguments) == 4:
    if arguments[3] == 'gz':
      compression=1
    elif arguments[3] == 'bz2': 
      compression=2
    else:
      print "unrecognized compression type: " + arguments[3]
      sys.exit(1)
  output_file_base = arguments[1]
  output_file_suffix = 'sdf' + ('.' + arguments[3] if len(arguments) == 4 else '')

  filea = open(input_filename, 'r')
  filea = None
  if input_filename.endswith('.bz2'):
    filea = bz2.BZ2File(input_filename, 'r')
  elif input_filename.endswith('gz'):
    filea = io.BufferedReader(gzip.open(input_filename, 'r'))
  else:
    filea = open(input_filename, 'r')
  line_numbers_to_split_at = []
  line_number = 0
  occurrence_count = 0
  sys.stdout.write('Counting molecules: ')
  sys.stdout.flush()
  while 1:
    line = filea.readline()
    if line is None or len(line) == 0:
      if len(line_numbers_to_split_at) > 0 and line_numbers_to_split_at[-1] != line_number - 1:
        line_numbers_to_split_at.append(line_number - 1)
      break
    if line == search_string:
      occurrence_count += 1
      if (occurrence_count % 100)==0:
        sys.stdout.write('\rCounting molecules: '+ str(occurrence_count+len(line_numbers_to_split_at) * n)+' '*20)
        sys.stdout.flush()
      if occurrence_count == n:
        line_numbers_to_split_at.append(line_number)
        occurrence_count = 0
    line_number += 1
  filea.close()
  print "Found " + str(len(line_numbers_to_split_at)) + " molecules to split at in " + input_filename

  n_digits = DetermineNumberDigitsInFilename(len(line_numbers_to_split_at))
  filea = None
  if input_filename.endswith('.bz2'):
    filea = bz2.BZ2File(input_filename, 'r')
  elif input_filename.endswith('gz'):
    filea = io.BufferedReader(gzip.open(input_filename, 'r'))
  else:
    filea = open(input_filename, 'r')
  line_number = 0
  file_number = 0
  for splitting_line_number in line_numbers_to_split_at:
    fname=output_file_base + '.' + MakeFilenameDigitStr(n_digits, file_number) + '.' + output_file_suffix
    new_file = open(fname, 'w') if compression == 0 else ( gzip.GzipFile(fname, 'w') if compression == 1 else bz2.BZ2File(fname, 'w'))
    while line_number <= splitting_line_number:
      line_number += 1
      new_file.write(filea.readline())
    new_file.close()
    print fname
    file_number += 1
  filea.close()

if __name__ == '__main__':
    main()
