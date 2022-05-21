#!/usr/bin/python
'''
Created on Nov 17, 2010

@author: mendenjl
'''

import os
import sys
import os.path
from symbol import exec_stmt

search_styles = [ "exact", "starts", "ends", "any"]

def usage():
  print "SplitFileAfter.py input-file search-string N output-file-base output-file-suffix [" + "/".join(search_styles) + "]"
  sys.exit(1)

search_type = 0

def matches(search_string, line):
  if search_type == 0:
    if len(line) == len(search_string) and search_string == line:
        return 1
  elif search_type == 1:
    if line.startswith(search_string):
      return 1
  elif search_type == 2:
    if line.endswith(search_string):
      return 1
  elif search_type == 3:
    if line.find(search_string) >= 0:
      return 1
  return 0

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

  if len(arguments) != 6:
    usage()

  try:
    int(arguments[ 2]);
  except:
    usage()

  input_filename = arguments[0]
  search_string = arguments[1]
  n = int(arguments[2])
  output_file_base = arguments[3]
  output_file_suffix = arguments[4].strip('.')
  search_type = findIndex(search_styles, arguments[ 5])
  if search_type == -1:
    usage()

  print "search type " + arguments[5] + " " + str(search_type)
  print "search string " + search_string
  print "n " + str(n)

  if search_type == 0 or search_type == 2:
    search_string += '\n'
  filea = open(input_filename, 'r')
  line_numbers_to_split_at = []
  line_number = 0
  occurrence_count = 0
  while 1:
    line = filea.readline()
    if line is None or len(line) == 0:
      if len(line_numbers_to_split_at) > 0 and line_numbers_to_split_at[-1] != line_number - 1:
        line_numbers_to_split_at.append(line_number - 1)
      break
    if matches(search_string, line) == 1:
      occurrence_count += 1
      if occurrence_count == n:
        line_numbers_to_split_at.append(line_number)
        occurrence_count = 0
    line_number += 1
  filea.close()
  print "Found " + search_string + " " + str(len(line_numbers_to_split_at)) + " times in " + input_filename

  n_digits = DetermineNumberDigitsInFilename(len(line_numbers_to_split_at))
  print '\n' + '\n'.join([ MakeFilenameDigitStr(n_digits, i) for i in xrange(len(line_numbers_to_split_at) + 1)])

  filea = open(input_filename, 'r')
  line_number = 0
  file_number = 0
  for splitting_line_number in line_numbers_to_split_at:
    new_file = open(output_file_base + '.' + MakeFilenameDigitStr(n_digits, file_number) + '.' + output_file_suffix, 'w')
    while line_number <= splitting_line_number:
      line_number += 1
      new_file.write(filea.readline())
    new_file.close()
    file_number += 1
  filea.close()

if __name__ == '__main__':
    main()
