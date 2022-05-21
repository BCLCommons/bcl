import string
from CodeFileUtils import *

access_types = ["none", "public", "protected", "private"]
implementation_types = ["default", "virtual", "pure virtual"]
function_holder_types = ["class", "namespace", "none"]

class Ctag:

  holder_type = 0
  holder_name = ""
  file_name = ""
  function_name = ""
  signature = ""
  function_type = ""
  access_type = 0
  implementation_type = 0

  def __init__(self, line):
    self.holder_type = 0
    self.holder_name = ""
    self.file_name = ""
    self.function_name = ""
    self.signature = ""
    self.function_type = ""
    self.access_type = 0
    self.implementation_type = 0

    pseudo_signature_start_pos = line.find('\t/^')
    pseudo_signature_end_pos = line.rfind(';"') + 2
    line = line[0:pseudo_signature_start_pos] + line[pseudo_signature_end_pos:]

    parts = line.strip(' \r').split('\t')
    parts = [ part for part in parts]
    if len(parts) == 0 or parts[0].startswith('!_TAG'):
      return

    if len(parts) < 4 or len(parts) > 8:
      print "Error in ctag parser, given line was: " + line + " which has # fields: " + str(len(parts))
      i = 0
      for part in parts:
        print str(i) + ": " + part
        i += 1
      return

    self.function_name = parts[0]

    file = parts[1]
    file_parts = file.split('/')[-2:]
    if file_parts[0] == 'include' or file_parts[0] == 'source':
      file_parts = [file_parts[1]]
    self.file_name = '/'.join(file_parts)

    self.function_type = parts[2]

    holder, self.signature = partition(parts[-1], ':')[::2]
    if holder != 'signature':
      print "Error in ctag parser, expected signature, but received: " + holder + " in line " + line
      return

    if len(parts) == 4:
      self.holder_name = ""
      self.holder_type = 2
    else:
      holder, self.holder_name = partition(parts[3], ':')[::2]
      if holder == "struct":
        holder = "class"
      self.holder_type = findIndex(function_holder_types, holder)
      for part in parts[4:-1:1]:
        holder, value = partition(part, ':')[::2]
        if holder == "implementation":
          self.implementation_type = findIndex(implementation_types, value)
        elif holder == "access":
          self.access_type = findIndex(access_types, value)
        elif holder != "file":
          print "Error in ctag parser, unknown tag:" + holder + " with value " + value + " in line " + line
          return

