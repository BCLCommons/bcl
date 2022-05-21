'''
Created on Apr 14, 2010
@brief Classes used for holding statistics and information on example checks for each file
@author: mendenjl
'''

import string
from CodeFileUtils import *
from BclUtils import *

class BCLDeveloperStatistic:

  name = ""
  number_of_examples_complete = 0
  number_of_examples_incomplete = 0
  number_of_classes = 0 # actually just the classes that require an example
  number_of_functions = 0
  number_of_example_checks = 0
  number_of_functions_distributed = 0.0
  number_of_examples_complete_distributed = 0.0
  number_of_examples_incomplete_distributed = 0.0
  number_of_classes_distributed = 0.0
  number_of_example_checks_distributed = 0.0

  def ConsiderFile(self, file):
    if len(file.classes) == 0 and file.example_unnecessary == 0:
      if self.name in file.authors:
        self.number_of_classes += 1
        self.number_of_classes_distributed += 1.0 / float(len(file.authors))
        self.number_of_functions += file.number_functions
        self.number_of_functions_distributed += float(file.number_functions) / float(len(file.authors))
    elif len(file.classes):
      for class_i_wrote in file.classes.values():
        if self.name in class_i_wrote.authors and class_i_wrote.example_unnecessary == 0:
          self.number_of_classes += 1
          self.number_of_classes_distributed += 1.0 / float(len(class_i_wrote.authors))
          self.number_of_functions += class_i_wrote.number_functions
          self.number_of_functions_distributed += float(class_i_wrote.number_functions) / float(len(class_i_wrote.authors))

  def ConsiderExample(self, example_i_wrote):
    if example_i_wrote.number_of_checks > 0:
      if example_i_wrote.status == "complete":
        self.number_of_examples_complete += len(example_i_wrote.included_files)
        self.number_of_examples_complete_distributed += float(len(example_i_wrote.included_files)) / float(len(example_i_wrote.authors))
      else:
        self.number_of_examples_incomplete += len(example_i_wrote.included_files)
        self.number_of_examples_incomplete_distributed += float(len(example_i_wrote.included_files)) / float(len(example_i_wrote.authors))
      self.number_of_example_checks += example_i_wrote.number_of_checks
      self.number_of_example_checks_distributed += float(example_i_wrote.number_of_checks) / float(len(example_i_wrote.authors))

  def __init__(self):
    self.number_of_examples_complete = 0
    self.number_of_examples_incomplete = 0
    self.number_of_classes = 0
    self.number_of_example_checks = 0
    self.number_of_examples_complete_distributed = 0.0
    self.number_of_examples_incomplete_distributed = 0.0
    self.number_of_classes_distributed = 0.0
    self.number_of_example_checks_distributed = 0.0

class BCLNamespaceStatistic:

  name = ""
  number_of_examples_complete = 0
  number_of_examples_incomplete = 0
  number_of_classes = 0 # actually just the classes that require an example
  number_of_functions = 0
  number_of_example_checks = 0

  def ConsiderFile(self, file):
    if file.namespace != self.name:
      return
    if len(file.classes) == 0 and file.example_unnecessary == 0:
      self.number_of_classes += 1
      self.number_of_functions += file.number_functions
    elif len(file.classes):
      for class_i_wrote in file.classes.values():
        if class_i_wrote.example_unnecessary == 0:
          self.number_of_classes += 1
          self.number_of_functions += class_i_wrote.number_functions

  def ConsiderExample(self, example):
    if example.number_of_checks > 0:
      if example.status == "complete":
        self.number_of_examples_complete += len(example.included_files)
      else:
        self.number_of_examples_incomplete += len(example.included_files)
      self.number_of_example_checks += example.number_of_checks

  def __init__(self):
    self.number_of_examples_complete = 0
    self.number_of_examples_incomplete = 0
    self.number_of_classes = 0
    self.number_of_example_checks = 0
    self.number_of_functions = 0

class BCLExampleFile:

  filename = ""
  short_filename = ""
  included_files = set([])
  namespace = ""
  status = ""
  date = ""
  authors = []
  reviewer_info = []
  number_of_checks = 0
  number_of_functions = 0
  valid_statuses = [ "empty", "incomplete", "complete"]

  def ExtractFromFile(self, filenm):
    self.filename = filenm
    file = open(self.filename, 'r')
    # load all lines from the file, then strip out multiline comments 
    ourLines = stripMultilineComments(file.readlines())
    file.close()
    self.short_filename = self.filename[self.filename.rfind('/') + 1:]
    if len(self.short_filename.split('_')) < 2 or self.short_filename.split('_')[0] != "example":
      self.filename = ""
      self.short_filename = ""
      return
    self.namespace = self.short_filename.split('_')[1].split('.')[0]

    #print self.filename
#    print self.short_filename
#    print self.namespace

    self.included_files = set([])
    for i in xrange(len(ourLines)):
      if ourLines[i].startswith("// include the header of the class"):
        i += 1
        while i < len(ourLines) and ourLines[i].startswith('#include \"') or (ourLines[i].startswith('//#include \"') and ourLines[i].endswith('.cpp\"')):
          include = ourLines[i].strip('/')[10:-1:1].split('/')[-1]
          include_namespace = ""
          split_filename = include.split('.')[0].split('_')
          if len(split_filename) > 1:
            include_namespace = split_filename[1] + os.sep
          #elif len(split_filename) == 1:
          #  include_namespace = split_filename[0] + os.sep
          self.included_files.add(include_namespace + include)
          i += 1
        break;

#    print "included: " + '\t'.join(self.included_files)

    lines, strings = extractStrings(stripComments(ourLines), "")
    is_example = 0
    for line in lines:
      if line.startswith("BCL_Example") and len(line) > len("BCL_Example") and "RICA_".find(line[len("BCL_Example")]) >= 0:
        self.number_of_checks += 1
      elif line.endswith("public ExampleInterface"):
        is_example = 1
#    if not is_example:
#      print "not an example"
#    else:
#      print "# checks: " + str(self.number_of_checks)
    if is_example == 0:
      self.__init__()
      return

    #print "our lines: " + '\n'.join(ourLines)

    for block in getDoxyBlocksContainingTag(ourLines, "example"):
      #print "block keys\t" + '\t'.join(block.keys())
      #print "block values\n" + '\n'.join(block.values())
      if 'author' in block and len(block['author'].strip(', ')):
        self.authors.extend(block['author'].strip(', ').split(', '))
      if 'date' in block:
        if len(self.date):
          self.date += ' ' + block['date']
        else:
          self.date += block['date']
      if 'remarks' in block:
        split_remarks = block['remarks'].split('\n')
        for remark in split_remarks:
          if remark.strip().startswith('status '):
            if len(self.status):
              self.status += ' ' + remark.strip()[7:].strip()
            else:
              self.status = remark.strip()[7:].strip().lower()
          elif remark.strip().startswith('reviewed '):
            self.reviewer_info.append(remark.strip())

  def GetErrorString(self, class_dict):
    # errors in example comments:
    #   No "// include the header of the class which this example is for" or #include following this comment
    #   No @remark status
    #   Remark claims anything other than empty but there are no example checks
    #   no authors
    #   no date
    errors = []
    if len(self.included_files) == 0:
      errors.append("target class is needed |")
    if self.number_of_checks != 0 and len(self.authors) == 0:
      errors.append("anonymous example |")

    # handle incorrect status
    if len(self.status) == 0:
      errors.append("no status |")
    elif Contains(BCLExampleFile.valid_statuses, self.status) == 0:
      errors.append("invalid status | " + self.status)
    elif self.number_of_checks == 0 and self.status == 'complete':
      errors.append("status should be empty |")
    elif self.number_of_checks != 0 and self.status == "empty":
      errors.append("status should be complete or incomplete |")

    for file in self.included_files:
      if file not in class_dict:
        errors.append("example claims to test a file that could not be found | " + file)
      elif len(class_dict[file].example_files) == 0 and not file.endswith('.cpp'):
        errors.append("file " + file + " forgot to link to example |")
      elif len(class_dict[file].example_files) == 1:
        for ex_file in class_dict[file].example_files:
          if ex_file != self.short_filename:
            errors.append("example claims to test this file but file does claims it is tested in " + ex_file + " | " + file)

    return errors

  def __init__(self):
    self.filename = ""
    self.short_filename = ""
    self.included_files = set([])
    self.namespace = ""
    self.status = ""
    self.date = ""
    self.authors = []
    self.reviewer_info = []
    self.number_of_checks = 0
