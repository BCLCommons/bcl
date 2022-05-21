#!/usr/bin/python

import os
import sys
import os.path
import bz2, gzip
from symbol import exec_stmt

def usage():
  print "SimplifyDecisionTree.py [-c] input-file"
  print "-c print coefficients only; do not pretty print the equation"
  sys.exit(1)

def main():
  arguments = sys.argv[1::1]

  pretty_print=True
  if len(arguments) == 0 or len(arguments) > 2:
    usage()
  if len(arguments) == 2:
    pretty_print = False
    if arguments[0] == '-c':
      arguments[0] = arguments[1]

  input_filename = arguments[0]
  input_f = None
  if input_filename.endswith('.bz2'):
    input_f = bz2.BZ2File(input_filename, 'r')
  elif input_filename.endswith('gz'):
    input_f = gzip.GzipFile(input_filename, 'r')
  else:
    input_f = open(input_filename, 'r')
  lines = [x.strip() for x in input_f.readlines()]
  input_f.close()
  
  matrix_lines = [x for x in range(len(lines)) if lines[x].startswith('bcl::linal::Matrix<float>')]
  bias_vector_lines = [x for x in range(len(lines)) if lines[x].find('Vector<float>') >= 0]
  if len(matrix_lines) != 1 or len(bias_vector_lines) == 0:
    print "Not a valid logistic regression model file! Had ",len(matrix_lines)," matrices and ",len(bias_vector_lines)," bias lines"
    sys.exit(1)
  matrix_size_line_num = matrix_lines[0]+1
  bias_vector_line_number = bias_vector_lines[-1]+2
  n_outputs = int(lines[matrix_size_line_num].split()[0])
  n_inputs = int(lines[matrix_size_line_num].split()[1])
  coefficients_matrix = []
  for output_n in range(n_outputs):
    coefficients_matrix.append([float(x) for x in lines[matrix_size_line_num+1+output_n].split()])
  biases = [float(x) for x in lines[bias_vector_line_number].split()]
  all_range_lines = [[float(y) for y in x.translate(None,',[]').split(' ') if len(y)] for x in lines if x.startswith('[')]
  input_ranges_from = all_range_lines[0:n_inputs]
  input_range_to = all_range_lines[n_inputs]
  output_ranges_from = all_range_lines[n_inputs+1:n_inputs+n_outputs+1]
  output_ranges_to = all_range_lines[-1]
  #print "Biases: ",biases
  #print "input_ranges_from: ",input_ranges_from
  #print "input_range_to: ",input_range_to
  output_range_to_width = 1.0
  input_range_to_width = input_range_to[1]-input_range_to[0]
  input_range_to_min = input_range_to[0]
  # convert input range max to expansion factor
  for x in range(len(input_ranges_from)):
    input_ranges_from[x][1] = input_range_to_width / ( input_ranges_from[x][1]-input_ranges_from[x][0]) if input_ranges_from[x][1]-input_ranges_from[x][0] > 1e-38 else 0.0
  #print "output_ranges_from: ",output_ranges_from
  #print "output_ranges_to: ",output_ranges_to
  for output in range(n_outputs):
    st = "Equation for output" + str(output + 1) + ": " if pretty_print else ""
    output_range_expansion = (output_ranges_from[output][1] - output_ranges_from[output][0])/output_range_to_width
    offset_accumulator = biases[output] - output_ranges_to[0] 
    for inp in range(n_inputs):
      offset_accumulator += coefficients_matrix[output][inp] * ( input_range_to_min - input_ranges_from[inp][1] * input_ranges_from[inp][0])
      st += str( coefficients_matrix[output][inp] * input_ranges_from[inp][1] * output_range_expansion) + ( " * X" + str(inp) if pretty_print else " ")
      if pretty_print: 
        st += ' + '
    st += str( offset_accumulator * output_range_expansion + output_ranges_from[output][0])
    if pretty_print:
      st = st.replace( '+ -', '- ')
    print st
  
if __name__ == '__main__':
    main()
