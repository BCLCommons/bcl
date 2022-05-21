'''
Created on Apr 14, 2010
@brief Classes used for holding information about bcl classes related to example checks
@author: mendenjl
'''

import string
from CodeFileUtils import *
from BclUtils import *
from BclExampleFile import *
from copy import deepcopy

def ExtractExampleFromBlock(block):
  if 'see' in block:
    block_str = block['see'].strip()
    for line in block_str.split('\n'):
      if line.startswith('@link ') and line.endswith('.cpp @endlink'):
        return line[len('@link '):-len('@endlink')].strip()
  if 'remarks' in block:
    split_remarks = block['remarks'].split('\n')
    for remark in split_remarks:
      if remark.strip() == 'example unnecessary':
        return 'unnecessary'
  return ""

class BCLClassComment:

  class_name = ""
  authors = set([])
  date = []
  example_unnecessary = 0
  example_file = ""
  had_comment = 0
  number_functions = 0
  functions = []

  def ExtractFromCommentBlock(self, block):
    if 'class' in block:
      self.class_name = block['class']
    if 'author' in block and len(block['author'].strip(', ')):
      self.authors |= set(block['author'].strip(', ').split(', '))
    if 'date' in block:
      self.date += '\n' + block['date']
    if 'remarks' in block:
      split_remarks = block['remarks'].split('\n')
      for remark in split_remarks:
        if remark.strip().lower() == 'example unnecessary':
          self.example_unnecessary = 1
    if len(block) > 0:
      self.had_comment = 1
    ex = ExtractExampleFromBlock(block)
    if ex == 'unnecessary':
      self.example_unnecessary = 1
    else:
      self.example_file = ex

  def __init__(self):
    self.class_name = ""
    self.authors = set([])
    self.date = []
    self.example_unnecessary = 0
    self.example_file = ""
    self.had_comment = 0
    self.number_functions = 0
    self.functions = []

class BCLClassFile:

  filename = ""
  short_filename = ""
  namespace = ""
  classes = {}
  examples = []
  authors = set([])
  example_files = set([])
  example_unnecessary = 0
  inapplicable_class_blocks = []
  uncommented_classes = []
  number_functions = 0
  functions = []

  def ExtractFromFile(self, filenm, example_dictionary):
    self.filename = filenm
    file = open(self.filename, 'r')
    # load all lines from the file, then strip out multiline comments 
    ourLines = stripMultilineComments(file.readlines())
    file.close()
    self.short_filename = self.filename[self.filename.rfind('/') + 1:]
    if len(self.short_filename.split('_')) < 2:
      self.namespace = self.short_filename.split('.')[0]
    else:
      self.namespace = self.short_filename.split('_')[1].split('.')[0]
    for block in getDoxyBlocksContainingTag(ourLines, "file"):
      file_example_str = ExtractExampleFromBlock(block)
      if len(file_example_str) > 0:
        if file_example_str == 'unnecessary':
          self.example_unnecessary = 1
        else:
          self.example_files.add(file_example_str)
      if 'author' in block and len(block['author'].strip(', ')):
        self.authors |= set(block['author'].strip(', ').split(', '))

    class_names = set(getClassesAndStructs(ourLines))
    class_blocks = getDoxyBlocksContainingTag(ourLines, "class")

    for block in class_blocks:
      class_name = block['class']
      if class_name in class_names:
        new_class_comment = BCLClassComment()
        new_class_comment.ExtractFromCommentBlock(block)
        if len(new_class_comment.class_name):
          self.classes[ class_name] = new_class_comment
          if len(new_class_comment.example_file):
            self.example_files.add(new_class_comment.example_file)
        for author in new_class_comment.authors:
          self.authors.add(author)
      else:
        self.inapplicable_class_blocks.append(block)

    for name in class_names:
      if name not in self.classes.keys():
        self.uncommented_classes.append(name)

    for ex in self.example_files:
      if ex in example_dictionary:
        self.examples.append(example_dictionary[ex])

  def NeedsAnExample(self):
    if self.example_unnecessary == 1:
      return 0

    if len(self.examples) > 0:
      number_checks = 0
      for example in self.examples:
        number_checks += example.number_of_checks
      return number_checks == 0

    if len(self.classes) == 0:
      return self.number_functions > 0

    for class_type in self.classes.values():
      if class_type.example_unnecessary == 0 and class_type.number_functions > 0:
        return 1

    return 0

  def TotalNumberFunctions(self):
    number_funcs = self.number_functions
    for classa in self.classes.values():
      number_funcs += classa.number_functions
    return number_funcs

  def DistributeFunctions(self):
    if len(self.examples) == 1:
      for example in self.examples:
        example.number_of_functions += self.TotalNumberFunctions()

  def GetErrorString(self):
    # errors in example comments:
    #   No "// include the header of the class which this example is for" or #include following this comment
    #   No @remark status
    #   Remark claims anything other than empty but there are no example checks
    #   no author
    #   no date
    errors = []
    if len(self.example_files) > 1:
      errors.append("different example filenames were linked |")

    for example in self.examples:
      if self.namespace == 'bcl':
        if self.short_filename not in example.included_files:
          errors.append("linked example did not include this class | " + example.short_filename)
      elif self.namespace + os.sep + self.short_filename not in example.included_files:
        errors.append("linked example did not include this class | " + example.short_filename)
    for block in self.inapplicable_class_blocks:
      errors.append("incorrect class name in class comment | " + block['class'])
    for class_name in self.uncommented_classes:
      errors.append("uncommented class | " + class_name)
    for class_type in self.classes.items():
      if len(class_type[1].authors) == 0:
        errors.append("class needs author | " + class_type[0])
    if len(self.authors) == 0 and len(self.classes) == 0 and self.number_functions > 0:
      errors.append("file needs author |")
    return errors

  def ConsiderTag(self, tag):
    # for class function, add the corresponding # of functions to the class
    if tag.holder_type == 0:
      class_name = tag.holder_name[ tag.holder_name.rfind(':') + 1:]
      if class_name in self.classes.keys():
        self.classes[class_name].number_functions += 1
        self.classes[class_name].functions.append(class_name + " | " + tag.function_name + tag.signature + " | " + ', '.join(self.classes[class_name].authors) + " |")
      # if the class could not be found, it is probably because it was never commented
      #else:
        #print "error, could not find class: " + class_name
    elif tag.function_name not in self.classes.keys():
      # the if statment here eliminates instantiations of s_Instance, which is found in header files and identified as a
      # function if the class is templated
      new_function_name = "| " + tag.function_name + tag.signature + " | " + ', '.join(self.authors) + " |"
      if new_function_name not in self.functions:
        self.number_functions += 1
        self.functions.append(new_function_name)

  def GetAllFunctions(self):
    functions = deepcopy(self.functions)
    for classes in self.classes.values():
      functions.extend(classes.functions)
    return functions

  def __init__(self):
    self.filename = ""
    self.short_filename = ""
    self.namespace = ""
    self.classes = {}
    self.examples = []
    self.authors = set([])
    self.example_files = set([])
    self.example_unnecessary = 0
    self.inapplicable_class_blocks = []
    self.uncommented_classes = []
    self.number_functions = 0
    self.functions = []
