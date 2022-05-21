import string
import os
import sys
import os.path
import cStringIO

# class represents a single SET(...) command from a CMakeLists file, optionally lower case
class CmakeListSetCommand:
  indent = 0
  variable = ""
  other_variables = []
  sources = set()
  flag_lines = []
  line_number_start = 0
  line_number_end = 0

  # returns line_number_end
  def FromFile(self, lines, line_number_start):
    self.line_number_start = line_number_start
    self.line_number_end = line_number_start
    lstripped_line = lines[ line_number_start].lstrip()
    lstripped_line_upper = lstripped_line.upper()
    self.indent = len(lines[ line_number_start]) - len(lstripped_line)
    if not lstripped_line_upper.startswith("SET("):
      return (-1)
    # get to the 1st line containing a parenthesis
    start_scope_pos = lines[line_number_start].find('(')
    self.line_number_end, end_scope_paren_pos = self.EndOfScope(lines, line_number_start, start_scope_pos)
    if end_scope_paren_pos < 0 or self.line_number_end < 0:
      print "Unclosed parenthesis!"
      return (-1)
    if self.line_number_end >= len(lines):
      self.line_number_end = len(lines) - 1
      end_scope_paren_pos = len(lines[self.line_number_end].rstrip()) - 1

    set_command_lines = lines[self.line_number_start:self.line_number_end + 1:]
    set_command = ''.join(set_command_lines)

    # remove comments from the lines
    set_command_lines_clean = [ line.split('#')[0].split('"')[0].split('\'')[0] for line in set_command_lines]
    set_command_clean = ''.join(set_command_lines_clean)

    if set_command_clean.find('.cpp') < 0 and (set_command_clean.find("SOURCES") < 0 or set_command_clean.find('${') > 0):
      # set command with no sources, not interesting
      return (-1)
    
    if end_scope_paren_pos != len(set_command_lines_clean[-1].rstrip()) - 1:
      print str(line_number_start) + " " + str(start_scope_pos) \
            + " " + str(self.line_number_end) + " " + str(end_scope_paren_pos) \
            + " " + str(len(lines)) + " " + str(len(set_command_lines_clean[-1].rstrip()) - 1)
      print "multiple commands on a single line not allowed, line: " + set_command_clean
      return (-1)

    if set_command_clean != set_command:
      print "could not parse set command: " + ''.join(set_command_lines_clean) + "\nbecause it contained comments or quotes"
      return (-1)

    # now we can finally parse the lines
    set_command_internal = set_command.strip()[start_scope_pos + 1:-1].replace('\n', ' ').replace('\t', ' ')

    #print "set command, internal component: " + set_command_internal

    # make sure there were no comments or quotes in the set command
    if set_command.find('#') > 0 or set_command.find('"') > 0 or set_command.find("'") > 0:
      print "Could not parse cmake set command with internal comments or quotes"
      return (-1)

    set_command_tokens = set_command_internal.split(' ')
    set_command_tokens = [ token for token in set_command_tokens if len(token) > 0]

    self.variable = set_command_tokens[ 0]
    for token_index in xrange(1, len(set_command_tokens)):
      if set_command_tokens[token_index].endswith('.cpp'):
        #print "adding source: " + set_command_tokens[token_index]
        self.sources.add(set_command_tokens[token_index].replace('${CMAKE_CURRENT_SOURCE_DIR}/', ''))
      elif len(self.sources) == 0:
        #print "other variables: " + set_command_tokens[token_index] + " in " + set_command_internal
        self.other_variables.append(set_command_tokens[token_index])
      else:
        #print "flag: " + set_command_tokens[token_index] + " in " + set_command_internal
        self.flag_lines.append(set_command_tokens[token_index])

    if len(self.sources) == 0 and (len(self.other_variables) or len(self.flag_lines)):
      print "no sources in " + set_command_internal
      return (-1)

    #print "sources: " + ','.join(self.sources)
    return self.line_number_end

  def Write(self, x):
    indent_primary_str = " " * self.indent
    indent_internal_str = indent_primary_str + "  "
    indent_source_str = indent_internal_str + '${CMAKE_CURRENT_SOURCE_DIR}/'
    end_internal_line_str = '\n' + indent_internal_str
    end_source_line_str = '\n' + indent_source_str
    if len(self.sources) == 0:
      x.write(indent_primary_str + "SET( " + self.variable + ")\n")
    elif len(self.sources) == 1 and self.line_number_end == self.line_number_start:
      # print on one line
      x.write(indent_primary_str + "SET( " + self.variable + " ${CMAKE_CURRENT_SOURCE_DIR}/" + ''.join([ a for a in self.sources]) + ")\n")
    else:
      x.write(indent_primary_str + "SET(\n")
      x.write(indent_internal_str + self.variable + '\n')
      if len(self.other_variables):
        x.write(indent_internal_str + end_internal_line_str.join(self.other_variables) + '\n')
      sorted_sources = sorted(self.sources)
      x.write(indent_source_str + end_source_line_str.join([ y for y in sorted_sources]) + '\n')
      if len(self.flag_lines):
        x.write(indent_internal_str + end_internal_line_str.join(self.flag_lines) + '\n')
      x.write(indent_primary_str + ')\n')

  # given a set of lines, start position (line_num = line number, pos = position on that line), find the end position of
  # a scope defined by ()
  def EndOfScope(self, lines, line_num, pos):
    if lines[line_num][pos] != '(':
      return line_num, pos

    depth = 1
    new_pos = pos + 1
    new_line = line_num
    while new_pos < len(lines[line_num]):
      if lines[line_num][new_pos] == '(':
        depth += 1
      elif lines[line_num][new_pos] == ')':
        depth -= 1
        if depth == 0:
          break
      new_pos += 1
    if depth > 0:
      new_line += 1
      while new_line < len(lines):
        new_pos = 0
        while new_pos < len(lines[new_line]):
          if lines[new_line][new_pos] == '(':
            depth += 1
          elif lines[new_line][new_pos] == ')':
            depth -= 1
            if depth == 0:
              break
          new_pos += 1
        if depth == 0:
          break
        new_line += 1
    return new_line, new_pos

  def __init__(self):
    self.indent = 0
    self.variable = ""
    self.other_variables = []
    self.sources = set()
    self.flag_lines = []
    self.line_number_start = 0
    self.line_number_end = 0

# check that the cmakelists.txt file in directory contains all the sources in the directory
# sources_dict is a dictionary of directory to source files contained in that directory
# if can_update is true, update the directory/CMakeLists.txt to include all the source
# files in that directory, and no non-existent files
def CheckCmakeLists(directory, sources_dict, can_update):

  # read the cmake lists file
  filename = directory + os.path.sep + "CMakeLists.txt"
  all_set_command_blocks = []
  all_sources = set()
  file = open(filename, 'r')
  lines = file.readlines()
  file.close()

  # get the sources in this directory
  sources_in_dir = sources_dict[ directory]

  # track sources to be added and removed from this file
  removed_sources = set()
  added_sources = set()

  # track the set command with the most source files; new source files will be added to that set command
  largest_set_command_index = -1
  largest_set_command_number_sources = 0

  # initialize line count to zero
  i = 0

  # iterate over lines
  while i < len(lines):

    # skip empty lines
    if len(lines[ i].strip()) == 0:
      i += 1
      continue

    # create a set command
    set_cmd = CmakeListSetCommand()

    # call FromFile to try to parse a set command from line #i in this file
    end_line = set_cmd.FromFile(lines, i)

    if end_line >= 0: # if the line is actually a start of a valid set command

      # copy the original sources for this set command
      original_sources = [ src for src in set_cmd.sources ]

      # iterate over sources, remove non-existent files
      for source in original_sources:
        src_path = directory + os.sep + source
        if source in sources_in_dir or os.path.isfile(src_path) or os.path.isfile(src_path + ".in"):
          all_sources.add(source)
        else:
          set_cmd.sources.remove(source)
          removed_sources.add(source)

      # is the # of sources larger now?
      if len(set_cmd.sources) > largest_set_command_number_sources:
        # yep, so update the largest set command
        largest_set_command_index = len(all_set_command_blocks)
        largest_set_command_number_sources = len(set_cmd.sources)

      # add the set command to the blocks
      all_set_command_blocks.append(set_cmd)

      # continue searching for set commands after the current block
      i = end_line

    # move to the next line
    i += 1

  if largest_set_command_index == -1:
    index = 0
    for set_cmds in all_set_command_blocks:
      if set_cmds.variable.endswith("SOURCES"):
        largest_set_command_index = index
        break
      index += 1

  # check that every source in the sources dictionary was found
  for source in sources_in_dir:
    if source not in all_sources:
      # this source was not found, add it to all_sources, added_sources, and the largest set command block, if applicable
      all_sources.add(source)
      added_sources.add(source)
      if largest_set_command_index != -1:
        all_set_command_blocks[largest_set_command_index].sources.add(source)

  # check if there is anything to remove or add; if not, just return
  if len(removed_sources) + len(added_sources) == 0:
    return

  # warn the user and return if the sources could not be added
  if largest_set_command_index == -1:
    print "Error: " + filename + ": Could not add sources because this file lacks a parseable set command that includes sources"
    return

  # warn the user if there is an ambiguous assignment of sources to a set command
  if len(all_set_command_blocks) > 1 and len(added_sources) > 0:
    print "Warning: ambiguous assignment of new sources in " + filename + " to largest block of sources"

  if can_update:
    # write the modified file to a cstringio object, then write that to a file
    # this is several times faster than writing directly to the file
    out_file = cStringIO.StringIO()
    prev_line_start = 0
    for cmake_set_number in xrange(len(all_set_command_blocks)):
      if(prev_line_start != all_set_command_blocks[cmake_set_number].line_number_start):
        out_file.write(''.join(lines[prev_line_start:all_set_command_blocks[cmake_set_number].line_number_start]))
      all_set_command_blocks[cmake_set_number].Write(out_file)
      prev_line_start = all_set_command_blocks[cmake_set_number].line_number_end + 1
    # write lines following the last set command, if applicable
    if prev_line_start < len(lines):
      out_file.write(''.join(lines[prev_line_start:]))

    # output the actual file
    file = open(filename, 'w')
    file.write(out_file.getvalue())
    file.close()

    # let the user know what was done
    if len(removed_sources):
      print "Removing " + ','.join(removed_sources) + " from " + filename
    if len(added_sources):
      print "Adding " + ','.join(added_sources) + " to " + filename
  else:
    # will not update; just inform the user as to what would be added/removed
    if len(removed_sources):
      print "Running in fix mode would remove " + ','.join(removed_sources) + " from " + filename
    if len(added_sources):
      print "Running in fix mode would add " + ','.join(added_sources) + " to " + filename
