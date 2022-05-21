#!/usr/bin/python

import os
import sys
import os.path

def getAllArguments(originalString, startPosa):
  startPos = startPosa
  argList = []
  while originalString[startPos] != '(':
    startPos += 1
  parenDepth = 1
  startPos += 1
  argStartPos = startPos
  while 1:
    if originalString[startPos] == '(':
      parenDepth += 1
    elif originalString[startPos] == ')':
      parenDepth -= 1
      if parenDepth == 0:
        if len(argList) > 0 or startPos > argStartPos:
          argList.append(originalString[argStartPos:startPos])
        break
    elif originalString[startPos] == ',' and parenDepth == 1:
      argList.append(originalString[argStartPos:startPos])
      argStartPos = startPos + 1
    startPos += 1

  for arg in argList:
    arg = arg.strip()
  return argList

def usage():
  print " Given a bcl-decision tree model, prints a if() else style representation of the tree"
  print " usage: SimplifyDecisionTree.py model-file [threshold] [descriptors-file (one per line)]"
  sys.exit(1)

def spaceEnd(originalString, i):
  j = i
  while j < len(originalString) and originalString[j] == ' ':
    j += 1
  return j

class decisionTreeNode:

  def __init__(self):
    self.decision_indx = 0
    self.decision_val = 0
    self.count = 0
    self.decision = ""
    self.branch_a = None
    self.branch_b = None

def main():
  arguments = sys.argv[1::1]

  if len(arguments) == 0 or len(arguments) > 3:
    usage()

  parent_node = decisionTreeNode()
  cutoff = 0.5
  feature_list = []
  if len(arguments) >= 2:
    cutoff = float(arguments[1])
  if len(arguments) >= 3:
    features_file = open(arguments[2], 'r')
    features_str = ""
    open_count = 0
    close_count = 0
    for line in features_file.readlines():
      open_count += line.count('(')
      close_count += line.count(')')
      features_str += line.strip()
      if open_count > 0 and open_count == close_count:
        break
    features_file.close()
    #feature_list = [x for x in getAllArguments(features_str, 0) if not x.startswith('Define(')]
    feature_list = getAllArguments(features_str, 0)

  fname = arguments[0]
  input_f = open(fname, 'r')
  skip_line = 3
  on_decision = True
  prev_cnt = 0
  prev_split_indx = 0
  prev_split_val = 0
  indent = 0
  used_indices = set()
  descriptor_strings = ["No", "Yes"]
  for line in input_f.readlines():
    if skip_line:
      skip_line -= 1
      continue
    if on_decision:
      on_decision = False
      skip_line = 2
      indent = spaceEnd(line, 0) - 4
      split_line = line.rstrip().split(' ')
      prev_cnt = split_line[-1]
      prev_split_val = split_line[-2]
      prev_split_indx = split_line[-3]
      if len(feature_list) and prev_split_indx != "nan":
        prev_split_indx = feature_list[int(prev_split_indx)]
    else:
      on_decision = True
      skip_line = 3
      if prev_split_indx == "nan":
        values = line.strip()
        values = ' '.join([ descriptor_strings[float(x) >= cutoff] for x in values.split()])
        print indent * ' ' + "Leaf\tWeight\t" + prev_cnt + '\t' + values + '\t' + line.strip()
      else:
        used_indices.add(prev_split_indx)
        print indent * ' ' + "Decision\tWeight\t" + prev_cnt + "\tFeature\t" + prev_split_indx + "\t<\t" + prev_split_val

  print "All used indices: " + '\t'.join([str(x) for x in sorted(used_indices)])
if __name__ == '__main__':
    main()
