#!/usr/bin/python2
'''
Created on Dec 5, 2012

Functions:
Creates .fwd.hh files for all namesapces in the bcl
For forward headers with typedefs that require external namespaces, 
also creates a .depends.fwd.hh file with the required classes

@author: mendenjl
'''

import os
import sys
import os.path
import datetime
from CodeFileUtils import *
import cStringIO
from curses.ascii import isspace, ispunct
from macpath import getmtime

copyright_block = []

def prune_terminal_comment(line):
  cpp_comment_pos = line.find('//')
  c_comment_pos = line.find('/*')
  if c_comment_pos < 0 and cpp_comment_pos < 0:
    return line
  start_comment_pos = 0
  if c_comment_pos >= 0:
    if cpp_comment_pos >= 0:
      start_comment_pos = min(cpp_comment_pos, c_comment_pos)
    else:
      start_comment_pos = cpp_comment_pos
  else:
    start_comment_pos = cpp_comment_pos
  return line[:start_comment_pos].rstrip()

def get_indent(line, start_pos = 0):
  i = start_pos
  while i < len(line) and line[i] == ' ':
    i += 1
  return i

def align_trailing_words(lines, n_words, indent):
  if n_words <= 0:
    return lines
  nlines = []
  longest_word_len = [ 0 ] * (n_words + 1)
  lines_split = [ line.rstrip().rsplit(' ', n_words) for line in lines]
  already_indented = []
  for i in xrange(len(lines_split)):
    lines_split[i][0] = lines_split[i][0].rstrip(' ')
    if lines_split[i][0][-1] == '\n':
      lines_split[i][0] += ' ' * indent
      already_indented.append(True)
    else:
      already_indented.append(False)
  for line_number in xrange(len(lines)):
    for i in xrange(len(lines_split[line_number])):
      actual_line_len = len(lines_split[line_number][i])
      if i == 0:
        nl_pos = lines_split[line_number][0].rfind('\n')
        if nl_pos >= 0:
          actual_line_len -= nl_pos + 1
      if longest_word_len[i] < actual_line_len:
        longest_word_len[i] = actual_line_len
  longest_word_len = [ x + 1 for x in longest_word_len]
  line_number = 0
  for line_split in lines_split:
    new_line = ""
    for i in xrange(len(line_split) - 1):
      actual_line_len = len(line_split[i])
      if i == 0:
        if already_indented[line_number]:
          actual_line_len = 0
        else:
          nl_pos = line_split[i].rfind('\n')
          if nl_pos >= 0:
            actual_line_len -= nl_pos + 1
      new_line += line_split[i] + ' ' * (longest_word_len[i] - actual_line_len)
    new_line += line_split[-1]
    nlines.append(new_line)
    line_number += 1
  return nlines

def remove_extra_spaces(line):
  nline = ""
  last_was_space = True
  for i in xrange(len(line)):
    if isspace(line[i]):
      if not last_was_space:
        nline += ' '
      last_was_space = True
    else:
      last_was_space = False
      nline += line[i]
  return nline

def getClassesFromOtherNamespacesDict(lines):
  namespace_to_classes = {}
  for line in lines:
    i = -1
    while i < len(line):
      # look for the scope resolution operator
      i = line.find('::', i + 1)
      if i < 0:
        break
      left_word_end = i - 1
      right_word_start = i + 2
      left_word_start = nextNameStart(line, left_word_end)
      left_word_end += 1
      right_word_end = nextNameEnd(line, right_word_start)
      namespace = line[left_word_start:left_word_end]
      class_name = line[right_word_start:right_word_end]
      i += 2
      if len(class_name) == 0 and len(namespace) == 0:
        # function pointer
        continue
      if namespace not in namespace_to_classes:
        namespace_to_classes[ namespace] = set([class_name])
      else:
        namespace_to_classes[ namespace].add(class_name)
    # check for size_t, which requires the cstddef header
    i = -1
    while i < len(line):
      i = line.find('size_t', i + 1)
      if i < 0:
        break
      if i > 0 and isValidAlpha(line[i - 1]):
        continue
      if i + 6 < len(line) and isValidAlpha(line[i + 6]):
        continue
      if 'std' not in namespace_to_classes:
        namespace_to_classes[ 'std'] = set(['cstddef'])
      else:
        namespace_to_classes[ 'std'].add('cstddef')

  return namespace_to_classes

def dictExtend(x, y):
  for ky, val in y.iteritems():
    if ky in x:
      x[ky].update(val)
    else:
      x[ky] = val

def dictListExtend(x, y):
  for ky, val in y.iteritems():
    if ky in x:
      x[ky].extend(val)
    else:
      x[ky] = val

class BclForwardHeader:

  namespaces = []
  filename = ""
  bcl_includes = set()
  std_includes = set()
  classes = {}
  templates = {}
  single_arg_templates = {}
  typedefs = {}
  external_dependencies = {}
  typedef_to_external_dependencies = {}
  deconvoluted_typedefs = {}
  external_typedef_dependencies = {}
  was_updated = False

  def ExtractFromFile(self, filenm, read_all):
    if filenm.endswith('.h'):
      self.extractFromNormalHeader(filenm)
      return
    self.filename = filenm
    file = open(self.filename, 'r')

    is_primary_fwd_hh = len(os.path.basename(filenm).split('_')) < 3
    in_statement = False
    last_statement = ""
    past_namespace = False
    statement_type_is_typedef = False
    for line in file.xreadlines():
      stripped_line = line.strip()
      if len(stripped_line) == 0:
        continue

      # skip comments
      if stripped_line.startswith('//'):
        continue

      # get includes
      if stripped_line.startswith('#'):
        if stripped_line.startswith('#include '):
          include = stripped_line.split(' ')[ 1]
          inner_include = include[1:-1]
          if include[0] == '"':
            if include.find('/') >= 0:
              self.bcl_includes.add(inner_include)
          else:
            self.std_includes.add(inner_include)
        continue

      if past_namespace:
        if stripped_line[ 0] == '}':
          past_namespace = True
          continue
        print "Illegal statement after close of namespace: " + line

      if not in_statement and stripped_line.startswith('namespace '):
        if len(self.classes) or len(self.typedefs):
          print "Illegal classes or typedefs before innermost namespace in forward header"
        namespace = stripped_line.split(' ')[ 1]
        self.namespaces.append(namespace)
        continue

      if stripped_line[ 0] == '{':
        continue
      elif stripped_line[ 0] == '}':
        past_namespace = True
        continue

      if in_statement:
        last_statement += '\n' + ' ' * get_indent(line) + prune_terminal_comment(stripped_line)
      else:
        statement_type_is_typedef = False
        if stripped_line.startswith('typedef'):
          statement_type_is_typedef = True
        last_statement = prune_terminal_comment(stripped_line)
      if last_statement.endswith(';'):
        in_statement = False
        if len(last_statement) == 1:
          continue
        name = last_statement.split(' ')[-1][:-1]
        full_statement = remove_extra_spaces(last_statement)
        in_templ = False
        if last_statement.startswith('template'):
          in_templ = True
        if not in_templ and len(full_statement) != len(last_statement) and len(full_statement.rsplit('>', 1)[0].rstrip()) <= 120:
          last_statement = full_statement
        if statement_type_is_typedef:
          self.typedefs[name] = last_statement
        elif in_templ and (not is_primary_fwd_hh or last_statement.find('=') >= 0):
          self.templates[name] = last_statement
        elif is_primary_fwd_hh and read_all:
          if in_templ:
            self.templates[name] = last_statement
          elif last_statement.startswith('class') or last_statement.startswith('struct'):
            self.classes[name] = last_statement
      else:
        in_statement = True

  # for speed reasons, this version is not as robust as it could be, ie. code like
  # /*
  # class BCL_API haha
  # */
  # will be registered as a class even though it is not really
  # However, c-style, multiline comments are not supposed to be in bcl code anyway, (the fix obvious bcl guidelines 
  # script should give a suitable warning), and parsing them makes this code many times slower, rendering it less useful
  # for its intended purpose
  def extractFromNormalHeader(self, filenm):
    self.filename = filenm
    file = open(self.filename, 'r')
    # read in the lines
    lines = stripCommentsAndQuotes(file.readlines())

    # remove preprocessor junk
    lines = [ x for x in lines if not x.startswith('#')]

    lines_split = []
    for line in lines:
      if len(line) == 1:
        lines_split.append(line)
        continue

      open_curly_pos = line.find('{')
      close_curly_pos = line.find('}')
      semi_colon_pos = line.find(';')
      if semi_colon_pos == len(line) - 1:
        semi_colon_pos = -1

      # ensure that all {} are alone on their lines
      # ensure that ';' is always followed by a new line
      if open_curly_pos >= 0 or close_curly_pos >= 0 or semi_colon_pos >= 0:
        line_resplit = line.replace('{', '\n{\n').replace('}', '\n}\n').replace(';', ';\n').split()
        line_resplit = [ x.strip() for x in line_resplit]
        line_resplit = [ x for x in line_resplit if len(x) > 0]
        lines_split.extend(line_resplit)
      else:
        lines_split.append(line)

    past_namespace = False
    in_statement = False
    last_statement = ""
    statement_type_is_typedef = False
    statement_type_is_template = False
    statement_type_is_class = False
    curly_depth = 0
    template_scope_depth = 0
    initializer_depth = 0

    for line in lines_split:
      if len(line) == 0:
        continue

      # handle namespaces
      if line.startswith('namespace '):
        if len(self.templates) or len(self.classes):
          past_namespace = True
        elif not past_namespace:
          split_line = line.split(' ')
          if len(split_line) == 2:
            self.namespaces.append(split_line[1])
          elif len(split_line) > 2:
            print "This makes no sense: " + line
            sys.exit(1)
        continue

      if len(line) == 1:
        if line == '{':
          curly_depth += 1
        elif line == '}':
          if curly_depth == len(self.namespaces):
            if len(self.templates) or len(self.classes):
              past_namespace = True
            else:
              del self.namespaces[-1]
          curly_depth -= 1

      if len(self.namespaces) == 0 or self.namespaces[0] != "bcl":
        continue

      in_statement_start = True
      if not in_statement and curly_depth == len(self.namespaces):
        last_statement = ""
        if line.startswith('template'):
          statement_type_is_template = True
          statement_type_is_typedef = False
          statement_type_is_class = False
          template_scope_depth = 0
          in_statement = True
          #print "templ start"
        elif line.startswith('typedef'):
          statement_type_is_typedef = True
          statement_type_is_template = False
          statement_type_is_class = False
          in_statement = True
          #print "typedef start"
        elif line.startswith('class ') or line.startswith('struct '):
          statement_type_is_typedef = False
          statement_type_is_template = False
          statement_type_is_class = True
          in_statement = True
          #print "class start"
        if in_statement and past_namespace:
          print "Illegal statement after close of namespace: " + line + ' in ' + filenm
          sys.exit(1)
      else:
        in_statement_start = False
      if not in_statement:
        continue
      #print "Start? " + str(in_statement_start)
      #print "Statement: " + line

      end_pos_start_looking = 0

      if statement_type_is_template and not statement_type_is_class:
        i = 0
        if in_statement_start:
          i += len('template')
          while i < len(line) and isspace(line[i]):
            i += 1
          end_pos_start_looking = i
        if i < len(line) and (line[i] == '<' or template_scope_depth):
          while i < len(line):
            if line[i] == '<':
              template_scope_depth += 1
            elif line[i] == '>':
              template_scope_depth -= 1
              if template_scope_depth == 0:
                reached_end = True
                i += 1
                break
            i += 1
          end_pos_start_looking = i
        if template_scope_depth == 0 and i < len(line):
          while end_pos_start_looking < len(line) and isspace(line[end_pos_start_looking]):
            end_pos_start_looking += 1
          if end_pos_start_looking < len(line):
            endes = nextNameEnd(line, end_pos_start_looking)
            next_word = line[end_pos_start_looking:endes]
            if next_word == 'class' or next_word == 'struct':
              end_pos_start_looking = endes
              statement_type_is_class = True
            elif len(next_word):
              #print "templated function: " + str(end_pos_start_looking) + ' ' + str(endes) + ' ' + (last_statement + ' ' + line).strip().replace('\n', ' ')
              in_statement = False
              continue
      expected_indent = 2 * (len(self.namespaces) + template_scope_depth)
      if template_scope_depth and (line[0] == '>' or line[0] == '<'):
        expected_indent -= 2

      if in_statement and len(last_statement):
        last_statement += '\n' + ' ' * expected_indent

      if end_pos_start_looking == len(line):
        last_statement += line
        continue

      # look for a semicolon
      semicolon_pos = line.find(';', end_pos_start_looking)

      # check for the end of the typedef
      if statement_type_is_typedef:
        if semicolon_pos < 0:
          last_statement += line
          continue
        else:
          end_pos_start_looking = semi_colon_pos
      elif statement_type_is_class:
        if line.find('<', end_pos_start_looking) >= 0 or line.find('(', end_pos_start_looking) >= 0:
          # partial specialization, no need to insert
          #print "Partial specialization or function: " + last_statement.replace('\n', ' ') + line
          in_statement = False
          continue

        if line == '{':
          end_pos_start_looking = 0
        else:
          # check for a colon, not allowed in class names since iit is used to denote scope
          colon_pos = line.find(':', end_pos_start_looking)

          if colon_pos < 0 and semi_colon_pos < 0:
            # just another continued line
            last_statement += line
            continue
          elif colon_pos >= 0:
            if semi_colon_pos >= 0:
              end_pos_start_looking = min(colon_pos, semi_colon_pos)
            else:
              end_pos_start_looking = colon_pos
          else:
            end_pos_start_looking = semi_colon_pos
      last_statement += line[:end_pos_start_looking]
      #print "Here 3: " + last_statement.replace('\n', ' ') + " " + str(statement_type_is_class)
      in_statement = False
      last_statement = last_statement.rstrip('\n ;:{') + ';'
      if len(last_statement) <= 1:
        continue
      if last_statement[0] == '\n':
        last_statement = last_statement[1:]
      if len(last_statement) <= 1:
        continue
      if statement_type_is_template:
        if last_statement.find('BCL_API') >= 0: 
          print "Unnecessary to have BCL_API on line " + last_statement + " in " + filenm
      elif statement_type_is_class and last_statement.find('BCL_API') < 0:
        print "Missing BCL_API on line " + last_statement + " in " + filenm
      last_statement = last_statement.replace('BCL_API ', '')
      #print "LS: " + last_statement.replace('\n', ' ')
      name = last_statement.split(' ')[-1][:-1]
      full_statement = remove_extra_spaces(last_statement)

      is_unary = full_statement.find(',') <= 0
      if not statement_type_is_template and len(last_statement) != len(full_statement) and len(full_statement.rsplit('>', 1)[0].rstrip()) <= 120:
        last_statement = full_statement
      if statement_type_is_typedef:
        self.typedefs[name] = last_statement
      elif statement_type_is_template :
        self.templates[name] = last_statement
      else:
        self.classes[name] = last_statement

  def combine(self, other):
    if self.namespaces != other.namespaces:
      if len(other.templates) or len(other.classes):
        print "cannot combine forward headers with different namespaces: "
        print "should have: " + str(self.namespaces) + " but received: " + str(other.namespaces) + " from file: " + other.filename
        sys.exit(1)
    else:
      self.bcl_includes = self.bcl_includes.union(other.bcl_includes)
      self.std_includes = self.std_includes.union(other.std_includes)

      for name, val in other.templates.iteritems():
        if name not in self.templates:
          self.templates[ name] = val
        elif val.find('=') >= 0:
          print "WARNING: You will need to manually remove the default template parameter from file " + other.filename
          self.templates[ name] = val
        elif self.templates[ name].find('=') >= 0:
          pass
        else:
          self.templates[ name] = val
      self.classes.update(other.classes)
      old_leng = len(self.typedefs)
      self.typedefs.update(other.typedefs)
      if old_leng < len(self.typedefs):
        print "WARNING: You should remove all namespace-level typedefs from " + other.filename

  def write(self, strm):
    term_filename = self.filename.split('/')[-1]
    header_guard = self.filename.split('/')[-1].upper().replace('.', '_') + '_'
    strm.write(copyright_block)
    strm.write("#ifndef " + header_guard + '\n')
    strm.write("#define " + header_guard + '\n\n')

    # if there are any external dependencies, include that file
    no_depends = False
    if len(self.external_dependencies) or len(self.typedefs):
      strm.write('// include the dependency file for this header\n')
      strm.write('#include "' + self.filename.split('/')[-1].replace('.fwd.hh', '.depends.fwd.hh') + '"\n')
    else:
      no_depends = True
      strm.write("// include bcl_defines.h header\n")
      strm.write('#include "bcl_defines.h"\n')

    if len(self.std_includes) and no_depends:
      strm.write("\n// external includes - sorted alphabetically\n")
      for include in self.std_includes:
        strm.write('#include <' + include + '>\n')

    strm.write('\n// This file contains forward declarations for the ' + self.namespaces[-1] + ' namespace')
    strm.write('\n// This file is mostly automatically generated')
    strm.write('\n// Developers may add typedefs and template default parameters')
    strm.write('\n// all other changes will be removed the next time this file is generated\n')
    indent = 0
    for namespace in self.namespaces:
      strm.write(' ' * indent + 'namespace ' + namespace + '\n' + ' ' * indent + '{\n')
      indent += 2

    strm.write(' ' * (indent - 2) + '/////////////////////\n')
    strm.write(' ' * (indent - 2) + '// regular classes //\n')
    strm.write(' ' * (indent - 2) + '/////////////////////\n')
    strm.write('\n')

    for line in align_trailing_words([ y[1] for y in sorted(self.classes.iteritems())], 1, indent):
      strm.write(' ' * indent + line + '\n')

    if len(self.classes):
      strm.write('\n')

    strm.write(' ' * (indent - 2) + '//////////////////////\n')
    strm.write(' ' * (indent - 2) + '// template classes //\n')
    strm.write(' ' * (indent - 2) + '//////////////////////\n')
    strm.write('\n')

    for line in [ y[1] for y in sorted(self.templates.iteritems())]:
      strm.write(' ' * indent + line + '\n\n')

    strm.write(' ' * (indent - 2) + '//////////////\n')
    strm.write(' ' * (indent - 2) + '// typedefs //\n')
    strm.write(' ' * (indent - 2) + '//////////////\n')
    strm.write('\n')

    for line in align_trailing_words([ y[1] for y in sorted(self.typedefs.iteritems())], 1, indent):
      strm.write(' ' * indent + line + '\n')

    if len(self.typedefs):
      strm.write('\n')

    for x in xrange(len(self.namespaces) - 1, -1, -1):
      indent -= 2
      strm.write(' ' * indent + '} // namespace ' + self.namespaces[x] + '\n')
    strm.write('\n')
    strm.write("#endif // " + header_guard + '\n')

  def updateDependencies(self):
    self.external_dependencies = {}
    if len(self.typedefs) == 0 and len(self.templates) == 0:
      return

    self.external_dependencies = getClassesFromOtherNamespacesDict(self.templates.itervalues())

    for typedef, val in self.typedefs.iteritems():
      new_dependencies = getClassesFromOtherNamespacesDict([val])
      if len(new_dependencies):
        self.typedef_to_external_dependencies[typedef] = new_dependencies
        dictExtend(self.external_dependencies, new_dependencies)

  def deconvolute(self, otherFwds):
    for typedef in self.typedefs.iterkeys():
      self.deconvoluteTypedef(typedef, otherFwds)
    self.std_includes = set()
    if 'std' in self.external_dependencies:
      for x in self.external_dependencies['std']:
        self.std_includes.add(self.getStdHeaders(x))
      self.external_dependencies['std']=[]

    # remove things from the external dependency list that turned out not to be legitimate namespaces
    non_namespaces = []
    for namespace, classes in self.external_dependencies.iteritems():
      if namespace not in otherFwds:
        non_namespaces.append(namespace)
    for namespace in non_namespaces:
      del self.external_dependencies[ namespace]

  def updateDependenciesFinal(self, otherFwds):
    symbols_to_check_dict = self.external_dependencies
    while len(symbols_to_check_dict):
      new_symbols_to_check_dict = {}
      for namespace, classes in symbols_to_check_dict.iteritems():
        if namespace not in otherFwds:
          continue
        fwd_header = otherFwds[namespace]
        for name in classes:
          if name in fwd_header.typedefs:
            dictExtend(new_symbols_to_check_dict, fwd_header.typedef_to_external_dependencies[ name])
            dictExtend(self.external_dependencies, fwd_header.typedef_to_external_dependencies[ name])
            if namespace not in self.external_typedef_dependencies:
              self.external_typedef_dependencies[namespace] = set([name])
            else:
              self.external_typedef_dependencies[namespace].add(name)
      symbols_to_check_dict = new_symbols_to_check_dict
    for namespace, classes in self.external_typedef_dependencies.iteritems():
      self.external_dependencies[namespace] = self.external_dependencies[namespace].difference(classes)

  def get(self, strn):
    if strn in self.classes:
      return self.classes[strn]
    elif strn in self.templates:
      return self.templates[strn]
    return ''

  def getStdHeaders(self, strm):
    if strm == 'string':
      return "string"
    if strm == 'less':
      return "functional"
    if strm.startswith('numeric_limits'):
      return 'limits'
    return strm

  def writeDependencies(self, strm, otherFwds):
    #print "typedefs: " + str(self.typedefs) + " typedefs to extern: " + str(self.typedef_to_external_dependencies)
    #print "external dependencies before: " + str(self.external_dependencies)
    self.updateDependenciesFinal(otherFwds)
    #print "external dependencies after: " + str(self.external_dependencies)
    if len(self.external_dependencies) == 0:
      return False

    term_filename = self.filename.split('/')[-1]
    header_guard = self.filename.split('/')[-1].upper().replace('.', '_').replace('FWD_HH', 'DEPENDS_FWD_HH_')
    strm.write(copyright_block)
    strm.write("#ifndef " + header_guard + '\n')
    strm.write("#define " + header_guard + '\n\n')
    strm.write("// include bcl_defines.h header\n")
    strm.write('#include "bcl_defines.h"\n')
    strm.write("\n// external includes - sorted alphabetically\n")
    for include in self.std_includes:
      strm.write('#include <' + include + '>\n')
    strm.write('\n// AUTOMATICALLY GENERATED FILE -- DO NOT EDIT!')
    strm.write('\n// This file contains forward declarations needed by the forward header for the ' + self.namespaces[-1] + ' namespace\n\n')
    strm.write('namespace bcl\n{')

    strm.write('\n////////////////////////////////')
    strm.write('\n// class forward declarations //')
    strm.write('\n////////////////////////////////\n')
    for namespace, classes in iter(sorted(self.external_dependencies.iteritems())):
      if len(classes) == 0:
        continue
      strm.write('\n  namespace ' + namespace + '\n  {')
      fwd_hdrs = otherFwds[namespace]
      for class_name in classes:
        strm.write('\n    ' + fwd_hdrs.get(class_name) + '\n')
      strm.write('  } // namespace ' + namespace)
    if len(self.external_dependencies):
      strm.write('\n')
    strm.write('\n///////////////////////')
    strm.write('\n// external typedefs //')
    strm.write('\n////////////////////////\n')
    for namespace, classes in iter(sorted(self.external_typedef_dependencies.iteritems())):
      if len(classes) == 0:
        continue
      strm.write('\n  namespace ' + namespace + '\n  {')
      fwd_hdrs = otherFwds[namespace]
      for class_name in classes:
        strm.write('\n    ' + fwd_hdrs.typedefs[class_name] + '\n')
      strm.write('  } // namespace ' + namespace)

    strm.write('\n} // namespace bcl\n\n')
    strm.write("#endif // " + header_guard + '\n')
    return True

  # given a particular typedef to deconvolute, and the set of all forward headers, examine each dependency of the typedef
  # if the dependency is a class, then just add it to the external dependencies and continue.  Otherwise, if it is also a typedef,
  # we need to get the related forward header to deconvolute this typedef first
  def deconvoluteTypedef(self, strn, all_fwd_hdrs):
    if strn not in self.typedef_to_external_dependencies or strn in self.deconvoluted_typedefs:
      return
    symbols_to_check_dict = self.typedef_to_external_dependencies[strn]
    while len(symbols_to_check_dict):
      new_symbols_to_check_dict = {}
      for namespace, class_names in symbols_to_check_dict.iteritems():
        if namespace == 'std':
          continue
        fwd_header = all_fwd_hdrs[namespace]
        for name in class_names:
          if namespace not in self.typedef_to_external_dependencies[strn]:
            self.typedef_to_external_dependencies[strn][namespace] = {}
          if name not in self.typedef_to_external_dependencies[strn][namespace]:
            if namespace not in self.external_dependencies:
              self.external_dependencies[namespace] = set([name])
            else:
              self.external_dependencies[namespace].add(name)
            if name in fwd_header.typedefs:
              if namespace in self.external_typedef_dependencies:
                self.external_typedef_dependencies[namespace].add(name)
              else:
                self.external_typedef_dependencies[namespace] = set([name])
              fwd_header.deconvoluteTypedef(name, all_fwd_hdrs)
              if namespace not in self.typedef_to_external_dependencies[strn]:
                self.typedef_to_external_dependencies[strn][namespace] = set([name])
              else:
                self.typedef_to_external_dependencies[strn][namespace].add(name)
              dictExtend(new_symbols_to_check_dict, fwd_header.typedef_to_external_dependencies[ name])
      symbols_to_check_dict = new_symbols_to_check_dict
    #print strn + " " + str(self.typedef_to_external_dependencies[ strn])
    self.deconvoluted_typedefs.add(strn)


  def __init__(self, fname = "", read_all = False):
    self.namespaces = []
    self.filename = ""
    self.bcl_includes = set()
    self.std_includes = set()
    self.templates = {}
    self.classes = {}
    self.typedefs = {}
    self.external_dependencies = {}
    self.external_typedef_dependencies = {}
    self.was_updated = False
    self.typedef_to_external_dependencies = {}
    self.deconvoluted_typedefs = set()
    if len(fname) and os.path.exists(fname):
      self.ExtractFromFile(fname, read_all)

def usage():
  print "\n" + \
        "usage: CreateNamespaceForwardHeaders.py bcl-path [options]\n" + \
        "options:\n" + \
        "-h/--help print this dialogue\n" + \
        "-o/--output overwrite existing forward headers if necessary (default is to print them to screen)\n" + \
        "-f/--force forces regeneration of all forward headers (e.g. ignore timestamps)\n" + \
        "-e/--existing augments the list of symbols with existing .fwd.hh files\n"
  sys.exit(1)

def main():

  bcl_path = "." + os.sep
  arguments = sys.argv[1::1]
  script_modification_time = os.path.getmtime(sys.argv[0])

  if len(arguments) == 0:
    usage()
  bcl_path = arguments[0]

  output = False
  force = False
  consider_old_fwd_headers = False

  for i in xrange(1, len(arguments)):
    if arguments[i] == '-h' or arguments[i] == '--help':
      usage()
    elif arguments[i] == '-o' or arguments[i] == '--output':
      output = True
    elif arguments[i] == '-f' or arguments[i] == '--force':
      force = True
    elif arguments[i] == '-e' or arguments[i] == '--existing':
      consider_old_fwd_headers = True
    else:
      print arguments[i] + " is not a valid option"
      usage()
  if not bcl_path.endswith('/'):
    bcl_path += '/'
  os.chdir(bcl_path)
  
  # read the copyright block
  global copyright_block
  copyright_block_file = os.path.abspath(bcl_path + "../../documentation/bcl_copyright.txt") if os.path.exists(bcl_path + "../../documentation/bcl_copyright.txt") \
                         else   os.path.abspath(bcl_path + "documentation/bcl_copyright.txt")
  if not os.path.exists(copyright_block_file):
    print "Cannot locate bcl copyright file at " + copyright_block_file + "! Exiting"
    sys.exit(-1)
  else:
    fl = open(copyright_block_file,'r')
    copyright_block = [x.lstrip() for x in fl.readlines() if len(x)]
    if copyright_block[-1][-1] != '\n':
      copyright_block[-1] += '\n'
    if copyright_block[-1] != '\n':
      copyright_block.append('\n') # append an extra new line to separate the copyright block from other code
    copyright_block = ''.join(copyright_block)
    fl.close()
    
  if not os.path.exists('./build'):
    os.mkdir('./build')

  # look for an existing file that denotes the last update times 
  timestamp_file = 'build/.forward_header_update_time.txt'
  last_update_time = 0
  if os.path.exists(timestamp_file):
    last_update_time = os.path.getmtime(timestamp_file)

  generated_fwd_headers = {}
  all_headers_full = getFilesWithSufficesInDirectoryDictionary("include", [".h"])
  if consider_old_fwd_headers:
    dictListExtend(all_headers_full, getFilesWithSufficesInDirectoryDictionary("apps", [".fwd.hh"]))

  # flatten all headers out, so that the forward headers are inserted into the proper place
  all_headers = {}
  for folder, files in all_headers_full.iteritems():
    namespace_folder_depth = 0
    if folder.startswith('include/'):
      namespace_folder_depth += 1
    folder_split = folder.split('/')
    if len(folder_split) > namespace_folder_depth + 1:
      namespace_folder = '/'.join(folder_split[:namespace_folder_depth + 1])
      internal_folder = '/'.join(folder_split[namespace_folder_depth + 1:]) + '/'
      if namespace_folder not in all_headers:
        all_headers[namespace_folder] = [ internal_folder + x for x in files]
      else:
        all_headers[namespace_folder].extend([ internal_folder + x for x in files])
    elif folder not in all_headers:
      all_headers[folder] = files
    else:
      all_headers[folder].extend(files)
  should_write_fwd_hdr = set()
  for folder, files in all_headers.iteritems():
    namespace_folder_depth = 0
    if folder != 'include':
      namespace_folder_depth = 1

    namespace_files = []
    for x in files:
      if x.count('_') <= namespace_folder_depth:
        namespace_files.append(x)
    if len(namespace_files) != 1:
      print "Unclear which file is the namespace header for " + folder + " among: " + str(namespace_files)
      sys.exit(1)
    namespace_file = namespace_files[0][:-2]
    namespace = namespace_file.split('_')[-1]
    namespace_fwd_hdr = folder + '/' + namespace_file
    namespace_depends_fwd_hdr = namespace_fwd_hdr
    namespace_depends_fwd_hdr += '.depends.fwd.hh'
    namespace_fwd_hdr += '.fwd.hh'
    fwd_hdr_mod_time = 0
    must_update = False
    if os.path.exists(namespace_fwd_hdr):
      if os.path.getmtime(namespace_fwd_hdr) > last_update_time:
        must_update = True
      if consider_old_fwd_headers:
        del files[findIndex(files, namespace_fwd_hdr)]
    if os.path.exists(namespace_depends_fwd_hdr):
      if os.path.getmtime(namespace_depends_fwd_hdr) > last_update_time:
        must_update = True
      if consider_old_fwd_headers:
        del files[findIndex(files, namespace_depends_fwd_hdr)]

    # always update if there was no fwd header, or the -f flag was given, or the script was modified more
    # recently than the forward header
    if force or last_update_time < script_modification_time:
      must_update = True
    else:
      if os.path.getmtime(folder) > last_update_time:
        must_update = True
      else:
        for fi in files:
          mod_time = os.path.getmtime(folder + '/' + fi)
          if mod_time > last_update_time:
            must_update = True
            break
    fwd_hdr = BclForwardHeader()
    if not os.path.exists(namespace_fwd_hdr):
      fwd_hdr.namespaces = ['bcl']
      if namespace != 'bcl':
        fwd_hdr.namespaces.append(namespace)
      fwd_hdr.filename = namespace_fwd_hdr
    if must_update:
      should_write_fwd_hdr.add(namespace)
      if os.path.exists(namespace_fwd_hdr):
        fwd_hdr = BclForwardHeader(namespace_fwd_hdr)
      for fi in files:
        fwd_hdr.combine(BclForwardHeader(folder + '/' + fi))
    else:
      if os.path.exists(namespace_fwd_hdr):
        fwd_hdr = BclForwardHeader(namespace_fwd_hdr, True)

    fwd_hdr.updateDependencies()
    generated_fwd_headers[ namespace] = fwd_hdr

  for hdr in generated_fwd_headers.itervalues():
    hdr.deconvolute(generated_fwd_headers)

  for namespace, namespace_fwd in generated_fwd_headers.iteritems():
    if namespace not in should_write_fwd_hdr:
      continue;
    direct = os.path.dirname(namespace_fwd.filename)
    writer = cStringIO.StringIO()
    namespace_fwd.write(writer)
    writerd = cStringIO.StringIO()
    should_writed = namespace_fwd.writeDependencies(writerd, generated_fwd_headers)
    dfilename = namespace_fwd.filename[:-7] + '.depends.fwd.hh'
    if output:
      should_write = True
      st = writer.getvalue()
      if os.path.exists(namespace_fwd.filename):
        ifile = open(namespace_fwd.filename, 'r')
        lines = ifile.read()
        ifile.close()
        if lines == st:
          should_write = False
      if should_write:
        print "Updating " + namespace_fwd.filename
        ofile = open(namespace_fwd.filename, 'w')
        ofile.write(st)
        ofile.close()
      if should_writed:
        should_write = True
        st = writerd.getvalue()
        if os.path.exists(dfilename):
          ifile = open(dfilename, 'r')
          lines = ifile.read()
          ifile.close()
          if st == lines:
            should_write = False
        if should_write:
          print "Updating " + dfilename
          ofile = open(dfilename, 'w')
          ofile.write(writerd.getvalue())
          ofile.close()
      elif os.path.exists(dfilename):
        print "Removing unnecessary dependency file: " + dfilename
        os.remove(dfilename)
    else:
      print namespace_fwd.filename + '\n' + writer.getvalue()
      if should_writed:
        print dfilename + '\n' + writerd.getvalue()

  if output:
    ofile = open(timestamp_file, 'w')
    ofile.write('This file\'s modification time specifies the last time that CreateNamespaceForwardHeaders.py was run\n')
    ofile.write('If files in a given namespace folder (or the folder itself, or the script) have been modified more\n')
    ofile.write('recently, the script is rerun on the given folder (or the whole bcl, if the script itself was modified)\n')
    ofile.close()

if __name__ == '__main__':
  main()
