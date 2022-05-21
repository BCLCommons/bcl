#!/usr/bin/python2
'''
Created on Feb 17, 2011

@author: mendenjl
'''

import sys
import os.path
from CodeFileUtils import *
from curses.ascii import isspace, ispunct, isdigit, islower, isalnum, isupper
from copy import deepcopy
from string import lower
import threading, multiprocessing, Queue, time
import cStringIO, commands
from time import gmtime, strftime, sleep

# print the usage for this script
def usage():
  print "BclCleaner.py path-to-bcl [options]"
  print "options: "
  print "-o (--overwrite) automatically overwrite any files with changes"
  print "-s (--stats) write statistics, useful for testing different methods of ambiguous symbol resolution"
  print "   can optionally be followed by a filename to output stats to a file"
  print "-i (--indent) correct indentation (experimental)"
  print "--show-status write out the status bar (primarily for eclipse, whose terminal does not support CR"
  print "--noprogress skip writing even decile percentage completion to the screen"
  print "--nothreads only use one thread, even if many are available"
  print "--directory DIR_NAME only clean files in the given directory"
  sys.exit(1)

# an enum for how spacing should be performed
class Spacing:
  No = -1
  Unk = 0
  Yes = 1

  # these operators usually should have a space before them
  space_lhs_operators = set(['+=', '=', '&=', '^=', '|=', '/=', '<<=', '%=', '*=', '>>=', '-=', '^', '|', '&&', '||', '?', ':', '/', '%', '==', '>=', '<=', '!=', '!', '~', '-' , '+', '{'])

  # these operators usually should not have a space before them
  no_space_lhs_operators = set([')', '}', ';', ']', ',', '.', '->', '.*', '->*', '::*'])

  # these operators usually should have a space after them
  space_rhs_operators = set(['<', '>', ',', ';', '+=', '=', '&=', '^=', '|=', '/=', '<<=', '%=', '*=', '>>=', '-=', '^', '|', '&&', '||', '?', ':', '/', '%', '==', '>=', '<=', '!=', '(', '{', '['])

  # these operators usually should not have a space after them
  no_space_rhs_operators = set(['!', '~', '.', '->', '.*', '->*', '::*'])

# an enum for code token types
class CodeTokenType:
  Space = -1 # a token containing only space
  Comment = 0 # // or /*
  BinaryOperator = 1 # +-*/,etc.  Most operators in C++ fall into this category
  Variable = 2 # a word; could also be a keyword, since we do not handle most keywords currently
  FenceOpen = 3 # ([
  FenceClose = 4 # )]
  BinaryOrUnaryOperator = 5 # &/*, these will be changed, where possible, to BinaryOperator or UnaryLHSOperator
  Function = 7 # Name of a function or constructor call
  FunctionPointerParen = 8 # ( before function pointer, e.g. x ( MyClass::*funcName)(), the ( before MyClass
  SignOrArithmetic = 9 # +/-, these will be changed to either BinaryOperator or UnaryLHSOperator
  UnaryLHSOperator = 10 # !~, also &*+- depending on context
  UnaryRHSOperator = 11 # : inside switches, also ++ and -- when used as postfix
  UnaryUnkOperator = 12 # ++/--; ultimately should be changed to UnaryLHSOperator or UnaryRHSOperator
  Constant = 13 # numbers and strings
  AccessSpec = 14 # public, private, protected
  Return = 15 # currently only return
  ScopeResolution = 16 # :: and ::*
  GTOrTemplate = 17 # >, changes into BinaryOperator or TemplateClose, if it can be determined unambiguously
  LTOrTemplate = 18 # <, changes into BinaryOperator or TemplateOpen, if it can be determined unambiguously
  TemplateOpen = 19 # <
  TemplateClose = 20 # >
  CurlyOpen = 21 # {
  CurlyClose = 22 # }
  SemiColon = 23 # ;
  TypeOrFuncName = 24
  TypeQualifier = 25
  Comma = 26
  Template = 27
  SpecialColon = 28 # colon after switch, case, or access-spec
  CCastFenceClose = 29
  Undefined = 30

  # array of True/False for whether a token, given on the LHS of a possibly unary +-&*, indicates that the token is unary
  is_unary_lhs_type = [ x in set([ SemiColon, Space, Return, BinaryOperator, LTOrTemplate, CCastFenceClose, TemplateOpen, SpecialColon, Comma, FenceOpen, CurlyClose, CurlyOpen, SignOrArithmetic, BinaryOrUnaryOperator, FunctionPointerParen, GTOrTemplate, TemplateClose, UnaryLHSOperator, TypeOrFuncName, TypeQualifier]) for x in xrange(Undefined)]

  # array of True/False for whether a token, given on the RHS of a possibly unary &* indicates that the token is unary
  is_unary_rhs_type = [ x in set([ Comma, FenceClose, TemplateClose, FunctionPointerParen, TypeQualifier]) for x in xrange(Undefined)]

  # array of True/False for whether a token, given on the LHS of a possibly binary +-, indicates that the token is binary
  is_arithmetic_lhs_type = [ x in set([Variable, TypeOrFuncName, FenceClose, UnaryUnkOperator, CCastFenceClose, Constant]) for x in xrange(Undefined)]

  # word types
  word_type = [ x in set([ Variable, Function, TypeOrFuncName, TypeQualifier]) for x in xrange(Undefined)]

  # types which, if they are at the beginning or end of a line, imply a well-defined indent for the next line
  has_well_defined_indent = [ x in set([ AccessSpec, CurlyOpen, CurlyClose, Return, SemiColon]) for x in xrange(Undefined)]

  # namespaces external to the bcl that both 1. have templates and 2. do not conform to our template naming conventions
  # For these scopes, X in X< is always read as a template, unless it is inside a nested scope, e.g.
  # std::vector< will be recognized as a template but std::string::npos< will not
  # is the start of a template type: if X contains only one template decl
  namespaces_non_bcl_template_rules = set(['std', 'cl'])

  # built-in types 
  builtin_types = set(['double', 'long', 'int', 'short', 'unsigned', 'signed', 'char', 'float', 'bool', 'size_t', 'cl_int', 'void', 'uint32_t', 'uint64_t', 'uint16_t', 'int64_t', 'int32_t', 'int16_t'])
  builtin_templated_functions = set([ 'const_cast', 'dynamic_cast', 'reinterpret_cast', 'static_cast', 'hash_map'])
  type_qualifiers = set(['const', 'volatile', 'typedef', 'typename', 'virtual'])
  start_of_line_prevent_continuation_indent = set([ FenceClose, FenceOpen, TemplateOpen, CurlyOpen, TemplateClose, CurlyClose])
  start_of_line_remove_indent = set([ FenceClose, TemplateClose, CurlyClose, AccessSpec])

template_open_types = set([CodeTokenType.LTOrTemplate, CodeTokenType.TemplateOpen])
template_close_types = set([CodeTokenType.GTOrTemplate, CodeTokenType.TemplateClose])

operator_to_type = {}
for op in all_operators:
  if op in set(['+=', '=', '&=', '^=', '|=', '/=', '<<=', '<<', '>>', '%=', '*=', '>>=', '-=', '^', '|', '&&', '||', '?', ':', '/', '%', '.', '->', '==', '>=', '<=', '!=', '->*', '.*']):
    operator_to_type[op] = CodeTokenType.BinaryOperator
  elif op == '!' or op == '~':
    operator_to_type[op] = CodeTokenType.UnaryLHSOperator
  elif op == '++' or op == '--':
    operator_to_type[op] = CodeTokenType.UnaryUnkOperator
  elif op == '*' or op == '&':
    operator_to_type[op] = CodeTokenType.BinaryOrUnaryOperator
  elif op == '+' or op == '-':
    operator_to_type[op] = CodeTokenType.SignOrArithmetic
  elif op == '/*' or op == '//' or op == '*/':
    operator_to_type[op] = CodeTokenType.Comment
  elif op == '(' or op == '[':
    operator_to_type[op] = CodeTokenType.FenceOpen
  elif op == ')' or op == ']':
    operator_to_type[op] = CodeTokenType.FenceClose
  elif op == '>':
    operator_to_type[op] = CodeTokenType.GTOrTemplate
  elif op == '<':
    operator_to_type[op] = CodeTokenType.LTOrTemplate
  elif op == '{':
    operator_to_type[op] = CodeTokenType.CurlyOpen
  elif op == '}':
    operator_to_type[op] = CodeTokenType.CurlyClose
  elif op == ',':
    operator_to_type[op] = CodeTokenType.Comma
  elif op == '::' or op == '::*':
    operator_to_type[op] = CodeTokenType.ScopeResolution
  elif op == ';':
    operator_to_type[op] = CodeTokenType.SemiColon
  else:
    print "unknown operator! " + op
    sys.exit(2)


class CodeStats:
  static_lock = threading.Lock()
  lines_of_code = 0
  gtlt_resolved = 0
  gtlt_unresolved = 0
  binunary_resolved = 0
  binunary_unresolved = 0
  number_errors = 0

  def init(self):
    self.lines_of_code = 0
    self.gtlt_resolved = 0
    self.gtlt_unresolved = 0
    self.binunary_resolved = 0
    self.binunary_unresolved = 0
    self.number_errors = 0

  def add(self, other_stats):
    CodeStats.static_lock.acquire()
    self.lines_of_code += other_stats.lines_of_code
    self.gtlt_resolved += other_stats.gtlt_resolved
    self.gtlt_unresolved += other_stats.gtlt_unresolved
    self.binunary_resolved += other_stats.binunary_resolved
    self.binunary_unresolved += other_stats.binunary_unresolved
    CodeStats.static_lock.release()

CODESTATS = CodeStats()


# classes and functions that are templated in the bcl must conform to the following regex:
# [A-Z]+[0-9]*[a-z][A-Za-z0-9]* 
# This means that they must start with a capital letter, have a lower case letter somewhere in them, and can have 
# any alphanumeric character otherwise
def CouldBeBclTemplatedClassOrFunctionName(s):
  if len(s) == 0 or not s[0] >= 'A' or not s[0] <= 'Z':
    if not s.startswith('t_'):
      return False
    start_pos = 2
  else:
    start_pos = 1

  # check that all characters are alphanumeric
  for i in xrange(start_pos, len(s)):
    if not isalnum(s[i]):
      return False

  # look for a lower-case letter, if the first character was upper case (e.g. everything except t_ types)
  if isupper(s[0]):
    for i in xrange(start_pos, len(s)):
      if islower(s[i]):
        return True
  else:
    for i in xrange(start_pos, len(s)):
      if isupper(s[i]):
        return True
  # no lower-case letter
  return False

class ParenthesisScopeType:
  is_c_cast = False
  type_tokens = [ x in set([ CodeTokenType.Variable, CodeTokenType.ScopeResolution, CodeTokenType.TypeOrFuncName, CodeTokenType.TypeQualifier]) for x in xrange(CodeTokenType.Undefined)]
  type_modifier_tokens = [ x in set([ CodeTokenType.TypeQualifier, CodeTokenType.UnaryLHSOperator, CodeTokenType.BinaryOrUnaryOperator]) for x in xrange(CodeTokenType.Undefined)]

  def __init__(self, tokens, start_token_id):
    template_depth = 0
    is_ambiguous = True
    self.is_c_cast = False
    for i in xrange(start_token_id, len(tokens)):
      typ = tokens[i].typ
      if tokens[i].token == ')':
        if is_ambiguous:
          self.is_c_cast = False
        break
      elif typ == CodeTokenType.TemplateOpen:
        template_depth += 1
      elif typ == CodeTokenType.TemplateClose:
        template_depth -= 1
      elif template_depth == 0:
        if self.is_c_cast:
          if typ == CodeTokenType.ScopeResolution:
            self.is_c_cast = False
          elif not ParenthesisScopeType.type_modifier_tokens[typ]:
            self.is_c_cast = False
            break
          else:
            is_ambiguous = False
        elif not ParenthesisScopeType.type_tokens[typ]:
          self.is_c_cast = False
          break
        elif typ == CodeTokenType.TypeOrFuncName:
          self.is_c_cast = True
          is_ambiguous = False
        elif typ == CodeTokenType.Variable:
          self.is_c_cast = True
        elif typ == CodeTokenType.TypeQualifier:
          is_ambiguous = False

class BclCodeToken:

  token = ""
  lhs_spaces_desired = Spacing.Unk
  rhs_spaces_desired = Spacing.Unk
  typ = -1
  num_spaces = 0

  def fixOperatorSpace(self):
    if self.token.startswith('operator') and ispunct(self.token[-1]):
      space_pos = 8
      space_end_pos = spaceEnd(self.token, space_pos)
      if space_end_pos == space_pos:
        self.token = self.token[:8] + ' ' + self.token[8:]
      elif space_end_pos > space_pos + 1:
        self.token = self.token[:9] + self.token[space_end_pos:]

  def __init__(self, tok, desired_lhs, desired_rhs, typ):
    self.typ = typ
    self.token = tok
    self.lhs_spaces_desired = desired_lhs
    self.rhs_spaces_desired = desired_rhs
    self.num_spaces = 0
    if typ == CodeTokenType.Space:
      self.num_spaces = len(tok)

class CodeState:
  init_list_curly_depth = 0
  func_curly_depth = 0
  paren_depth = 0
  in_c_cast = False
  just_saw_func_decl = False
  in_func_decl = False
  template_depth = 0
  in_typedef = False
  in_func_pointer = False
  in_member_initializer_list = False
  last_was_sequence_point = False
  in_func_call = False
  is_definite_type = False
  last_token = ''
  last_type = CodeTokenType.Space
  saw_type = Spacing.Unk
  namespace_curly_depth = 0

  def registerSequencePoint(self):
    self.just_saw_func_decl = self.in_typedef = self.in_func_pointer = self.in_func_decl = False
    self.saw_type = Spacing.Yes
    self.in_member_initializer_list = False
    self.template_depth = 0
    self.is_definite_type = False
    self.in_c_cast = False
    self.last_was_sequence_point = True

  def GetMinIndent(self, first_token_type):
    indent = self.in_member_initializer_list + self.init_list_curly_depth + self.func_curly_depth + self.paren_depth + self.template_depth + self.namespace_curly_depth

    if indent == 0:
      return 0

    #if self.paren_depth == 0 and self.template_depth == 0 and self.in_member_initializer_list == 0 and not self.last_was_sequence_point:
    #  if first_token_type not in CodeTokenType.start_of_line_prevent_continuation_indent and self.last_type not in CodeTokenType.start_of_line_prevent_continuation_indent:
    #    if self.last_type != CodeTokenType.Comma and self.last_type != CodeTokenType.SpecialColon:
    #      #print "Adding indent because last token was " + self.last_token
    #      indent += 1
    if first_token_type in CodeTokenType.start_of_line_remove_indent:
      indent -= 1
    if first_token_type == CodeTokenType.AccessSpec and ((self.last_token == ':' and self.last_type != CodeTokenType.SpecialColon) or self.last_type == CodeTokenType.Comma) :
      # inheritance requires indent, but previous if statement removes it
      indent += 2
    # special case to handle first { after member initializer list
    if self.in_member_initializer_list and first_token_type == CodeTokenType.CurlyOpen:
      indent -= 1
    return indent * 2

  def IsIndentAbsolute(self, first_token_type):
    indent = self.in_member_initializer_list + self.init_list_curly_depth + self.func_curly_depth + self.paren_depth + self.template_depth + self.namespace_curly_depth

    if indent == 0:
      # probably in preprocessor stuff, which can be indented extra
      return False

    if self.init_list_curly_depth > 0 or (self.last_token == '=' and first_token_type == CodeTokenType.CurlyOpen):
      # initialization lists ambiguate the indentation level
      return False

    if CodeTokenType.has_well_defined_indent[ first_token_type] or CodeTokenType.has_well_defined_indent[ self.last_type]:
      # inheritance and access specifier indentation is well defined
      return True

    # member initializer lists have very well defined indentation rules, because all calls are to constructors, so 
    # there is never an aesthetic resason for indenting the line differents
    if self.in_member_initializer_list:
      return True

    return False

  def ExplainIndent(self, first_token_type):
    print " in_member_initializer_list " + str(self.in_member_initializer_list)
    print " init_list_curly_depth " + str(self.init_list_curly_depth)
    print " func_curly_depth " + str(self.func_curly_depth)
    print " paren_depth " + str(self.paren_depth)
    print " template_depth " + str(self.template_depth)
    print " namespace_curly_depth " + str(self.namespace_curly_depth)
    if first_token_type in CodeTokenType.start_of_line_remove_indent:
      print " first token type required removal"
    if first_token_type == CodeTokenType.AccessSpec and ((self.last_token == ':' and self.last_type != CodeTokenType.SpecialColon) or self.last_type == CodeTokenType.Comma):
      print " inheritance (indent + 2)"
    if self.in_member_initializer_list and first_token_type == CodeTokenType.CurlyOpen:
      print " was in member init list (indent + 1)"


  def update(self, all_tokens, i):
    tok = all_tokens[i]
    typ = tok.typ
    if typ == CodeTokenType.Space or typ == CodeTokenType.Comment:
      return typ
    self.last_was_sequence_point = False
    if typ == CodeTokenType.AccessSpec or typ == CodeTokenType.Template or typ == CodeTokenType.SpecialColon:
      pass
    elif tok.token == ':':
      self.saw_type = Spacing.No
      if self.func_curly_depth == 0 and self.paren_depth == 0 and self.just_saw_func_decl:
        self.in_member_initializer_list = True
    elif typ == CodeTokenType.Return or typ == CodeTokenType.UnaryRHSOperator or typ == CodeTokenType.BinaryOperator or typ == CodeTokenType.Constant:
      self.saw_type = Spacing.No
    elif typ == CodeTokenType.TemplateOpen:
      self.template_depth += 1
      self.saw_type = Spacing.Yes
    elif typ == CodeTokenType.TemplateClose:
      self.template_depth -= 1
      self.saw_type = Spacing.Yes
    elif typ == CodeTokenType.CurlyOpen:
      self.saw_type = Spacing.Unk
      if self.init_list_curly_depth or self.last_token == '=':
        self.init_list_curly_depth += 1
        self.saw_type = Spacing.No
      elif self.func_curly_depth:
        self.func_curly_depth += 1
      elif self.just_saw_func_decl or self.in_member_initializer_list:
        self.just_saw_func_decl = False
        self.in_member_initializer_list = False
        self.in_func_decl = False
        self.func_curly_depth = 1
        self.saw_type = Spacing.Unk
      else:
        self.namespace_curly_depth += 1
    elif typ == CodeTokenType.CurlyClose:
      self.registerSequencePoint()
      if self.init_list_curly_depth:
        self.init_list_curly_depth -= 1
        self.saw_type = Spacing.No
      elif self.func_curly_depth:
        self.func_curly_depth -= 1
      else:
        self.namespace_curly_depth -= 1
    elif typ == CodeTokenType.SemiColon:
      self.registerSequencePoint()
    elif tok.token == '(':
      if typ == CodeTokenType.FunctionPointerParen:
        self.saw_type = Spacing.Yes
        self.in_func_pointer = True
      elif not self.in_member_initializer_list and not self.in_typedef and not self.in_func_decl and self.template_depth == 0 and self.paren_depth == 0 and self.func_curly_depth == 0:
        if self.last_token != ')':
          self.in_func_decl = True
        self.saw_type = Spacing.Yes
        self.in_c_cast = False
      else:
        if self.last_token == 'for':
          self.in_c_cast = False
          self.saw_type = Spacing.Yes
          self.is_definite_type = True
        elif CodeTokenType.word_type[self.last_type] or self.last_type == CodeTokenType.TemplateClose or self.last_token == 'sizeof':
          self.in_c_cast = False
          if self.saw_type == Spacing.Unk:
            self.saw_type = Spacing.No
        else:
          if self.last_type == CodeTokenType.FenceClose and self.last_token == ')':
            self.in_c_cast = False
          else:
            self.in_c_cast = ParenthesisScopeType(all_tokens, i + 1).is_c_cast
          if self.in_c_cast:
            self.saw_type = Spacing.Yes
          else:
            self.saw_type = Spacing.Unk
      self.paren_depth += 1
    elif tok.token == ')':
      if self.in_func_pointer:
        self.in_func_pointer = False
        if i >= len(all_tokens) - 1 or all_tokens[i + 1].token != '(':
          self.in_func_pointer = False
        else:
          typ = CodeTokenType.CCastFenceClose
        self.saw_type = Spacing.Yes
      else:
        if self.in_c_cast:
          typ = CodeTokenType.CCastFenceClose
          self.in_c_cast = False
        elif self.in_func_decl and self.paren_depth == 1:
          self.just_saw_func_decl = True
          self.in_func_decl = False
        else:
          self.saw_type = Spacing.No
      self.paren_depth -= 1
    elif typ == CodeTokenType.UnaryLHSOperator:
      if tok.token != '*' and tok.token != '&':
        self.saw_type = Spacing.No
    elif typ == CodeTokenType.TypeOrFuncName or typ == CodeTokenType.Variable:
      if self.is_definite_type:
        typ = CodeTokenType.TypeOrFuncName
        self.is_definite_type = False
      if typ == CodeTokenType.TypeOrFuncName and self.saw_type == Spacing.Unk:
        self.saw_type = Spacing.Yes
      elif self.func_curly_depth == 0 and tok.token.startswith('s_'):
        self.saw_type = Spacing.No
    elif typ == CodeTokenType.TypeQualifier:
      if tok.token == 'typedef':
        self.in_typedef = True
      self.saw_type = Spacing.Yes
    elif typ == CodeTokenType.Comma:
      if not self.init_list_curly_depth:
        if self.in_func_decl and self.paren_depth == 1:
          self.saw_type = Spacing.Yes
    elif tok.token == '[' or tok.token == ']':
      self.saw_type = Spacing.No
    elif typ == CodeTokenType.BinaryOrUnaryOperator:
      if CodeTokenType.is_unary_lhs_type[ self.last_type] or self.template_depth or self.in_typedef or self.in_c_cast or self.in_func_pointer:
        typ = CodeTokenType.UnaryLHSOperator
      elif self.saw_type == Spacing.No:
        typ = CodeTokenType.BinaryOperator
      elif self.paren_depth == 1 and self.in_func_decl:
        typ = CodeTokenType.UnaryLHSOperator
      elif self.paren_depth == 0:
        typ = CodeTokenType.UnaryLHSOperator
      else:
        typ = CodeTokenType.BinaryOperator
    self.last_token = tok.token
    self.last_type = typ
    return typ

  def __init__(self, x = None):
    if x is None:
      self.paren_depth = 0
      self.init_list_curly_depth = 0
      self.func_curly_depth = 0
      self.just_saw_func_decl = False
      self.in_func_decl = False
      self.template_depth = 0
      self.in_typedef = False
      self.in_func_pointer = False
      self.in_member_initializer_list = False
      self.last_token = ''
      self.last_type = CodeTokenType.Space
      self.saw_type = Spacing.Unk
      self.is_definite_type = False
      self.in_c_cast = False
      self.last_was_sequence_point = True
      self.namespace_curly_depth = 0
    else:
      self.paren_depth = x.paren_depth
      self.init_list_curly_depth = x.init_list_curly_depth
      self.func_curly_depth = x.func_curly_depth
      self.just_saw_func_decl = x.just_saw_func_decl
      self.in_func_decl = x.in_func_decl
      self.template_depth = x.template_depth
      self.in_typedef = x.in_typedef
      self.in_func_pointer = x.in_func_pointer
      self.in_member_initializer_list = x.in_member_initializer_list
      self.last_token = x.last_token
      self.last_type = x.last_type
      self.saw_type = x.saw_type
      self.is_definite_type = x.is_definite_type
      self.in_c_cast = x.in_c_cast
      self.last_was_sequence_point = x.last_was_sequence_point
      self.namespace_curly_depth = x.namespace_curly_depth

op_lhs_space_pref = {}
op_rhs_space_pref = {}
for op in all_operators:
  if op in Spacing.space_lhs_operators:
    op_lhs_space_pref[ op] = Spacing.Yes
  elif op in Spacing.no_space_lhs_operators:
    op_lhs_space_pref[ op] = Spacing.No
  else:
    op_lhs_space_pref[ op] = Spacing.Unk
  if op in Spacing.space_rhs_operators:
    op_rhs_space_pref[ op] = Spacing.Yes
  elif op in Spacing.no_space_rhs_operators:
    op_rhs_space_pref[ op] = Spacing.No
  else:
    op_rhs_space_pref[ op] = Spacing.Unk

def IsNotFunctionPointerStart(line, i):
  if line[i] != '(':
    return True
  en = line.find(')', i)
  if en < 0:
    return True
  middle = line[i + 1:en:1].replace(' ', '')
  if en + 1 >= len(line):
    return True
  star_pos = middle.find('*')
  if star_pos < 0 or star_pos == len(middle) - 1:
    return True;

  if middle.find('(') >= 0:
    return True

  next_word_end = nextNameEnd(middle, star_pos + 1)
  if next_word_end != len(middle):
    return True
  if star_pos != 0 and (star_pos == 1 or middle[star_pos - 2:star_pos] != '::'):
    return True

  en += 1
  while en < len(line) and isspace(line[en]):
    en += 1
  if en < len(line):
    if line[en] != '(':
      return True
  else:
    return True
  if line.find('.') >= 0 or line.find('+') >= 0 or line.find('-') >= 0 or line.find('/') >= 0:
    return True
  if star_pos == 0 and lower(middle) == middle:
    return True
  return False

class BclCodeLine:

  code_tokens = []
  code_spaces = []
  ending_template_paren_depths = []
  state = CodeState()
  min_indent = 0
  is_indent_definite = False
  next_paren_is_func_pointer = False

  def GetNextToken(self, lin, st):
    pos = st
    ch = lin[pos]
    if isspace(ch):
      while pos < len(lin) and isspace(lin[pos]):
        pos += 1
      return BclCodeToken(lin[st:pos], Spacing.Unk, Spacing.Unk, CodeTokenType.Space)
    elif isValidAlphaStart(ch):
      pos = nextNameEnd(lin, st)
      name = lin[st:pos]
      if name == "operator":
        spos = spaceEnd(lin, pos)
        npos = walkToEndOfOperator(lin, spos)
        if npos == spos and spos > pos:
          npos = lin.find('(', spos)
        if npos > spos:
          op = lin[spos:npos]
          if op == '(' or op == '[':
            npos += 1
          name = lin[st:npos]
        return BclCodeToken(name, Spacing.Unk, Spacing.No, CodeTokenType.Function)
      elif name == "template":
        return BclCodeToken(name, Spacing.Unk, Spacing.Unk, CodeTokenType.Template)
      elif name == 'return':
        return BclCodeToken(name, Spacing.Unk, Spacing.Yes, CodeTokenType.Return)
      elif name == 'public' or name == 'private' or name == 'protected':
        return BclCodeToken(name, Spacing.Unk, Spacing.Unk, CodeTokenType.AccessSpec)
      elif name == 'delete' or name == 'sizeof':
        return BclCodeToken(name, Spacing.Yes, Spacing.Unk, CodeTokenType.UnaryLHSOperator)
      spos = spaceEnd(lin, pos)
      if spos == pos and pos < len(lin) - 1 and lin[pos] == ':' and lin[pos + 1] == ':':
        spos = walkToEndOfOperator(lin, pos)
        #print "St: " + str(st) + " end:  " + str(pos)
        op = lin[pos:spos]
        if op == '::' and spos < len(lin) and isValidAlphaStart(lin[spos]):
          pos = nextNameEnd(lin, spos)
          class_or_func_name = lin[spos:pos]
          if name in CodeTokenType.namespaces_non_bcl_template_rules or CouldBeBclTemplatedClassOrFunctionName(class_or_func_name):
            return BclCodeToken(lin[st:pos], Spacing.Unk, Spacing.Unk, CodeTokenType.TypeOrFuncName)
      could_be_name = CouldBeBclTemplatedClassOrFunctionName(name) or name in CodeTokenType.builtin_types or name in CodeTokenType.builtin_templated_functions
      if spos < len(lin) and could_be_name and (lin[spos] == '&' or lin[spos] == '*'):
        old_spos = spos
        spos = spaceEnd(lin, spos + 1)
      if spos < len(lin) and lin[spos] == '(':
        if IsNotFunctionPointerStart(lin, spos):
          return BclCodeToken(name, Spacing.Unk, Spacing.No, CodeTokenType.Function)
        else:
          self.next_paren_is_func_pointer = True
          return BclCodeToken(name, Spacing.Unk, Spacing.Yes, CodeTokenType.TypeOrFuncName)
      elif could_be_name:
        return BclCodeToken(name, Spacing.Unk, Spacing.Unk, CodeTokenType.TypeOrFuncName)
      elif name in CodeTokenType.type_qualifiers:
        return BclCodeToken(name, Spacing.Unk, Spacing.Unk, CodeTokenType.TypeQualifier)
      return BclCodeToken(name, Spacing.Unk, Spacing.Unk, CodeTokenType.Variable)
    elif isdigit(ch) or ch == '.' or ch == '-' or ch == '+':
      pos = walkToEndOfCPPNumber(lin, st)
      if pos != st:
        if not isdigit(ch):
          return BclCodeToken(lin[st:pos], Spacing.Yes, Spacing.Yes, CodeTokenType.Constant)
        return BclCodeToken(lin[st:pos], Spacing.Yes, Spacing.Yes, CodeTokenType.Constant)
    elif ch == "'" or ch == '"':
      pos += 1
      while pos < len(lin) and lin[pos] != ch:
        if lin[pos] == '\\':
          pos += 1
        pos += 1
      if pos < len(lin) and lin[pos] == ch:
        pos += 1
      return BclCodeToken(lin[st:pos], Spacing.Yes, Spacing.Yes, CodeTokenType.Constant)
    elif ch == '#':
      return BclCodeToken(lin[st:], Spacing.Unk, Spacing.Unk, CodeTokenType.Comment)
    elif ch == '\\':
      return BclCodeToken(lin[st:], Spacing.Unk, Spacing.No, CodeTokenType.Comment)
    pos = walkToEndOfOperator(lin, st)
    #print "St: " + str(st) + " end:  " + str(pos)
    op = lin[st:pos]
    if op == '//':
      return BclCodeToken(lin[st:], Spacing.Unk, Spacing.Unk, CodeTokenType.Comment)
    elif op == '/*':
      pos = lin.find('*/', st + 1)
      if pos == -1:
        pos = len(lin)
      else:
        pos += 2
      return BclCodeToken(lin[st:pos], Spacing.Unk, Spacing.Unk, CodeTokenType.Comment)
    elif op == '(' and self.next_paren_is_func_pointer:
      self.next_paren_is_func_pointer = False
      return BclCodeToken(op, Spacing.Unk, Spacing.Yes, CodeTokenType.FunctionPointerParen)
    if (op == '&' or op == '*') and self.next_paren_is_func_pointer:
      return BclCodeToken(op, Spacing.Yes, Spacing.No, CodeTokenType.UnaryLHSOperator)
    return BclCodeToken(op, op_lhs_space_pref[op], op_rhs_space_pref[op], operator_to_type[op])

  def OptimizeSpacing(self):
    if len(self.code_spaces) == 0:
      return

    for i in xrange(len(self.code_tokens)):
      if self.code_tokens[i].rhs_spaces_desired == Spacing.No:
        self.code_spaces[i + 1] = 0
      elif self.code_tokens[i].rhs_spaces_desired == Spacing.Yes:
        self.code_spaces[i + 1] = 1
      else:
        self.code_spaces[i + 1] = min(self.code_spaces[i + 1], 1)
    self.code_spaces[-1] = 0

  def CountBinaryOps(self):
    ops = 0
    for tok in self.code_tokens:
      if tok.typ == CodeTokenType.BinaryOperator or tok.typ == CodeTokenType.BinaryOrUnaryOperator:
        ops += 1
    return ops

  def AddSpaces(self):
    for i in xrange(len(self.code_tokens)):
      if self.code_tokens[i].rhs_spaces_desired == Spacing.Yes:
        self.code_spaces[i + 1] = max(self.code_spaces[i + 1], 1)

  def GetString(self):
    if len(self.code_spaces) == 0:
      return '\n'
    if len(self.code_tokens) == 0:
      return ' ' * self.code_spaces[0] + '\n'
    st = ""
    for i in xrange(len(self.code_tokens)):
      st += self.code_spaces[i] * ' ' + self.code_tokens[i].token
    return st + '\n'

  def getLastCodeTokenType(self):
    if len(self.code_tokens) == 0:
      return CodeTokenType.Space
    if self.code_tokens[-1].typ == CodeTokenType.Comment:
      if len(self.code_tokens) == 1:
        return CodeTokenType.Space
      else:
        return self.code_tokens[-2].typ
    return self.code_tokens[-1].typ

  def LabelObviousTemplates(self, prev_token_typ, previous_depths):
    next_to_last_token = len(self.code_tokens) - 1
    template_paren_depths = []
    paren_depth = 0
    prev_type = prev_token_typ
    for i in xrange(len(self.code_tokens)):
      # template< 
      if prev_type == CodeTokenType.Template and self.code_tokens[i].typ == CodeTokenType.LTOrTemplate:
        self.code_tokens[i].typ = CodeTokenType.TemplateOpen
        self.code_tokens[i].lhs_spaces_desired = Spacing.No
        template_paren_depths.append(paren_depth)
      elif self.code_tokens[i].typ == CodeTokenType.FenceOpen or self.code_tokens[i].typ == CodeTokenType.FunctionPointerParen:
        paren_depth += 1
      elif self.code_tokens[i].typ == CodeTokenType.FenceClose:
        paren_depth -= 1
      elif self.code_tokens[i].typ == CodeTokenType.GTOrTemplate and len(template_paren_depths):
        if template_paren_depths[-1] == paren_depth:
          self.code_tokens[i].typ = CodeTokenType.TemplateClose
          if i and prev_type != CodeTokenType.TemplateClose:
            self.code_tokens[i].lhs_spaces_desired = Spacing.No
          del template_paren_depths[-1]
      prev_type = self.code_tokens[i].typ

    template_terminal_chars = set(['>', ',', '::', '::*', ';'])
    if prev_token_typ in template_close_types and self.code_tokens[0].token in template_terminal_chars:
      prev_token_typ = CodeTokenType.TemplateClose
      if self.code_tokens[0].token == '>':
        self.code_tokens[0].typ = CodeTokenType.TemplateClose
    elif prev_token_typ in template_open_types:
      if self.code_tokens[0].token == '>':
        prev_token_typ = CodeTokenType.TemplateOpen
        self.code_tokens[0].typ = CodeTokenType.TemplateClose
      elif len(self.code_tokens) > 1 and self.code_tokens[1].token == '>':
        if self.code_tokens[0].typ == CodeTokenType.Variable or self.code_tokens[0].typ == CodeTokenType.TypeOrFuncName or self.code_tokens[0].typ == CodeTokenType.Constant:
          prev_token_typ = CodeTokenType.TemplateOpen
          self.code_tokens[1].typ = CodeTokenType.TemplateClose

    # change <X> into templates. strictly speaking, x < y > z is legitimate C++, but in practice is never used because 
    # the x < y evaluates to true or false, so there is no need to compare >
    prev_type = self.code_tokens[0].typ
    for i in xrange(1, len(self.code_tokens) - 1):
      # template< 
      t = self.code_tokens[i].typ
      if prev_type in template_open_types and self.code_tokens[i + 1].typ in template_close_types:
        if t == CodeTokenType.Variable or t == CodeTokenType.TypeOrFuncName or t == CodeTokenType.Constant:
          if i > 1:
            if self.code_tokens[i - 2].typ == CodeTokenType.Variable:
              self.code_tokens[i - 2].typ = CodeTokenType.TypeOrFuncName
          elif prev_token_typ == CodeTokenType.Variable:
            prev_token_typ = CodeTokenType.TypeOrFuncName
          self.code_tokens[i - 1].typ = CodeTokenType.TemplateOpen
          self.code_tokens[i - 1].lhs_spaces_desired = Spacing.No
          self.code_tokens[i + 1].typ = CodeTokenType.TemplateClose
          self.code_tokens[i + 1].lhs_spaces_desired = Spacing.No
      prev_type = self.code_tokens[i].typ

    prev_type = prev_token_typ
    for i in xrange(len(self.code_tokens)):
      if self.code_tokens[i].typ == CodeTokenType.LTOrTemplate:
        if prev_type == CodeTokenType.UnaryUnkOperator or prev_type == CodeTokenType.Constant:
          self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
          self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
          self.code_tokens[i].typ = CodeTokenType.BinaryOperator
        elif i < next_to_last_token:
          if self.code_tokens[i + 1].typ == CodeTokenType.UnaryUnkOperator:
            self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
            self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
            self.code_tokens[i].typ = CodeTokenType.BinaryOperator
          elif self.code_tokens[i + 1].typ in template_close_types:
            self.code_tokens[i].lhs_spaces_desired = self.code_tokens[i + 1].lhs_spaces_desired = Spacing.No
            self.code_tokens[i].typ = CodeTokenType.TemplateOpen
            self.code_tokens[i + 1].typ = CodeTokenType.TemplateClose
      elif self.code_tokens[i].typ == CodeTokenType.TemplateOpen:
        if i < next_to_last_token and self.code_tokens[i + 1].typ in template_close_types:
            self.code_tokens[i + 1].lhs_spaces_desired = Spacing.No
            self.code_tokens[i + 1].typ = CodeTokenType.TemplateClose
      elif self.code_tokens[i].typ == CodeTokenType.GTOrTemplate:
        if prev_type == CodeTokenType.UnaryUnkOperator:
          self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
          self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
          self.code_tokens[i].typ = CodeTokenType.BinaryOperator
        elif prev_type == CodeTokenType.TemplateClose:
          self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
          self.code_tokens[i].typ = CodeTokenType.TemplateClose
        elif prev_type == CodeTokenType.TypeOrFuncName or prev_type == CodeTokenType.TypeQualifier:
          self.code_tokens[i].lhs_spaces_desired = Spacing.No
          self.code_tokens[i].typ = CodeTokenType.TemplateClose
        elif i < next_to_last_token:
          if self.code_tokens[i + 1].typ in set([CodeTokenType.UnaryUnkOperator, CodeTokenType.Constant, CodeTokenType.UnaryLHSOperator]):
            self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
            self.code_tokens[i].typ = CodeTokenType.BinaryOperator
          elif self.code_tokens[i + 1].typ in template_close_types:
            self.code_tokens[i].typ = CodeTokenType.TemplateClose
            self.code_tokens[i + 1].typ = CodeTokenType.TemplateClose
            self.code_tokens[i + 1].lhs_spaces_desired = Spacing.Yes
          elif self.code_tokens[i + 1].typ == CodeTokenType.SemiColon:
            self.code_tokens[i].typ = CodeTokenType.TemplateClose
            self.code_tokens[i].lhs_spaces_desired = Spacing.No
      prev_type = self.code_tokens[i].typ

    prev_type = self.code_tokens[0].typ
    for i in xrange(1, len(self.code_tokens)):
      if prev_type in template_close_types and self.code_tokens[i].token in template_terminal_chars:
        #if self.code_tokens[i].token == '(' or self.code_tokens[i].typ == CodeTokenType.ScopeResolution:

        if prev_type == CodeTokenType.GTOrTemplate:
          self.code_tokens[i - 1].typ = CodeTokenType.TemplateClose
          self.code_tokens[i - 1].lhs_spaces_desired = Spacing.No
        if self.code_tokens[i].typ == CodeTokenType.ScopeResolution:
          self.code_tokens[i].lhs_spaces_desired = Spacing.No
      prev_type = self.code_tokens[i].typ

    prev_type = prev_token_typ
    # try using the names to identify template start characters
    for i in xrange(len(self.code_tokens)):
      if self.code_tokens[i].typ == CodeTokenType.LTOrTemplate:
        if prev_type == CodeTokenType.TypeOrFuncName:
          self.code_tokens[i].typ = CodeTokenType.TemplateOpen
          self.code_tokens[i].lhs_spaces_desired = Spacing.No
          self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
        else:
          self.code_tokens[i].typ = CodeTokenType.BinaryOperator
          self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
          self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
      prev_type = self.code_tokens[i].typ

    prev_type = prev_token_typ
    template_paren_depths = previous_depths
    paren_depth = 0
    hit_sequence_point = False
    if prev_token_typ == CodeTokenType.SemiColon or prev_token_typ == CodeTokenType.CurlyClose or prev_token_typ == CodeTokenType.CurlyOpen:
      hit_sequence_point = True
    for i in xrange(len(self.code_tokens)):
      if self.code_tokens[i].typ == CodeTokenType.TemplateOpen:
        template_paren_depths.append(paren_depth)
      elif self.code_tokens[i].typ == CodeTokenType.FenceOpen or self.code_tokens[i].typ == CodeTokenType.FunctionPointerParen:
        paren_depth += 1
      elif self.code_tokens[i].typ == CodeTokenType.FenceClose:
        paren_depth -= 1
      elif self.code_tokens[i].typ == CodeTokenType.CurlyClose or self.code_tokens[i].typ == CodeTokenType.CurlyOpen or self.code_tokens[i].typ == CodeTokenType.SemiColon:
        hit_sequence_point = True
        template_paren_depths = []
      elif self.code_tokens[i].typ == CodeTokenType.GTOrTemplate or self.code_tokens[i].typ == CodeTokenType.TemplateClose:
        if len(template_paren_depths):
          if template_paren_depths[-1] == paren_depth:
            self.code_tokens[i].typ = CodeTokenType.TemplateClose
            if i and prev_type != CodeTokenType.TemplateClose:
              self.code_tokens[i].lhs_spaces_desired = Spacing.No
            del template_paren_depths[-1]
          elif paren_depth > template_paren_depths[-1] or hit_sequence_point:
            self.code_tokens[i].typ = CodeTokenType.BinaryOperator
            self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
          else:
            print "Bad parens for template? In " + self.GetString() + " Template depths;  " + str(template_paren_depths)
        elif self.code_tokens[i].typ == CodeTokenType.GTOrTemplate:
          self.code_tokens[i].typ = CodeTokenType.BinaryOperator
          self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
        else:
          print "Overclosed template? In " + self.GetString() + " Template depths;  " + str(template_paren_depths)
      prev_type = self.code_tokens[i].typ
    for i in xrange(len(self.code_tokens)):
      if self.code_tokens[i].typ in template_close_types and i < next_to_last_token and self.code_tokens[i + 1].token == '(':
        self.code_tokens[i].rhs_spaces_desired = Spacing.No
    self.ending_template_paren_depths = [ x - paren_depth for x in template_paren_depths]
    return prev_token_typ

  def combineBuiltinMultiWordTypes(self):
    i = 1
    while i < len(self.code_tokens):
      if self.code_tokens[i].typ == CodeTokenType.TypeOrFuncName and self.code_tokens[i].token in CodeTokenType.builtin_types and self.code_tokens[i - 1].typ == CodeTokenType.TypeOrFuncName and self.code_tokens[i - 1].token in CodeTokenType.builtin_types:
        self.code_tokens[i - 1].token += ' ' + self.code_tokens[i].token
        del self.code_tokens[i]
        del self.code_spaces[i]
        if i < len(self.code_tokens) and self.code_tokens[i].typ == CodeTokenType.TypeOrFuncName and self.code_tokens[i].token in CodeTokenType.builtin_types:
          self.code_tokens[i - 1].token += ' ' + self.code_tokens[i].token
          del self.code_tokens[i]
          del self.code_spaces[i]
      else:
        i += 1

  def __init__(self, line, filename, prev_token_typ = CodeTokenType.Space, prev_template_depths = [], state = CodeState()):
    self.code_tokens = []
    self.code_spaces = []
    self.ending_template_paren_depths = []
    self.state = state
    self.min_indent = 0
    self.is_indent_definite = False
    self.next_paren_is_func_pointer = False

    #print "Tokenizing " + line
    just_added_token = True
    i = 0
    while i < len(line):
      #print "Next: sz: " + str(len(line) - i) + " ln: " + line[i:]
      next_tok = self.GetNextToken(line, i)
      if not next_tok.num_spaces:
        if just_added_token:
          self.code_spaces.append(0)
        self.code_tokens.append(next_tok)
        just_added_token = True
      else:
        self.code_spaces.append(next_tok.num_spaces)
        just_added_token = False
      i += len(next_tok.token)
    if just_added_token:
      self.code_spaces.append(0)
    else:
      self.code_spaces[-1] = 0

    if len(self.code_tokens) == 0:
      return
    self.code_tokens[0].lhs_spaces_desired = Spacing.Unk
    self.code_tokens[-1].rhs_spaces_desired = Spacing.Unk
    self.code_tokens[-1].token = self.code_tokens[-1].token.rstrip('\n\r')
    self.combineBuiltinMultiWordTypes()

    prev_type = prev_token_typ
    for i in xrange(len(self.code_tokens)):
      if prev_type == CodeTokenType.TypeQualifier:
        if self.code_tokens[i].typ == CodeTokenType.Variable:
          self.code_tokens[i].typ = CodeTokenType.TypeOrFuncName
          if i + 2 < len(self.code_tokens) and self.code_tokens[i + 1].typ == CodeTokenType.ScopeResolution:
            self.code_tokens[i + 2].typ = CodeTokenType.TypeOrFuncName
        elif self.code_tokens[i].typ == CodeTokenType.BinaryOrUnaryOperator:
          i_back = i - 2
          if i_back >= 0:
            if self.code_tokens[i_back].typ == CodeTokenType.Variable or self.code_tokens[i_back].typ == CodeTokenType.TypeOrFuncName:
              self.code_tokens[i_back].typ = CodeTokenType.TypeOrFuncName
      prev_type = self.code_tokens[i].typ

    prev_type = prev_token_typ

    # switch / case / default modifcations
    if len(self.code_tokens) > 1:
      if self.code_tokens[1].token == ':' and (self.code_tokens[0].typ == CodeTokenType.AccessSpec or self.code_tokens[0].token == "default"):
        self.code_tokens[1].typ = CodeTokenType.SpecialColon
        self.code_tokens[1].lhs_spaces_desired = Spacing.No
      elif self.code_tokens[0].token == "case":
        for i in xrange(2, len(self.code_tokens)):
          if self.code_tokens[i].token == ':':
            self.code_tokens[i].typ = CodeTokenType.SpecialColon
            self.code_tokens[i].lhs_spaces_desired = Spacing.No
            break

    for i in xrange(len(self.code_tokens)):
      # Destructors, function pointers
      if prev_type == CodeTokenType.ScopeResolution and (self.code_tokens[i].token == '~' or self.code_tokens[i].token == '*'):
        self.code_tokens[i].lhs_spaces_desired = self.code_tokens[i].rhs_spaces_desired = Spacing.No
        self.code_tokens[i].tok = CodeTokenType.UnaryLHSOperator
      prev_type = self.code_tokens[i].typ

    # constants where the + or - is a binary operator instead of a unary operator
    i = 0
    prev_type = prev_token_typ
    while i < len(self.code_tokens):
      if self.code_tokens[i].typ == CodeTokenType.Constant and (self.code_tokens[i].token[0] == '+' or self.code_tokens[i].token[0] == '-'):
        if not CodeTokenType.is_unary_lhs_type[prev_type]:
          sig = self.code_tokens[i].token[0]
          self.code_tokens[i].token = self.code_tokens[i].token[1:]
          self.code_tokens.insert(i, BclCodeToken(sig, Spacing.Yes, Spacing.Yes, CodeTokenType.BinaryOperator))
          self.code_spaces.insert(i + 1, 0)
          self.code_tokens[i].lhs_spaces_desired = self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
          i += 1
      prev_type = self.code_tokens[i].typ
      i += 1
    if len(self.code_tokens) > 1:
      # skip before any trailing comment
      i = len(self.code_tokens) - 1
      if self.code_tokens[i].typ == CodeTokenType.Comment:
        i -= 1
      if self.code_tokens[i].token == ';':
        i -= 1
        while i >= 0 and self.code_tokens[i].token == ';':
          del self.code_tokens[i]
          del self.code_spaces[i]
          i -= 1

    next_to_last_token = len(self.code_tokens) - 1

    prev_type = prev_token_typ
    for i in xrange(len(self.code_tokens)):
      if self.code_tokens[i].typ == CodeTokenType.SignOrArithmetic:
        if CodeTokenType.is_arithmetic_lhs_type[ prev_type]:
          self.code_tokens[i].typ = CodeTokenType.BinaryOperator
          self.code_tokens[i].lhs_spaces_desired = self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
        else:
          self.code_tokens[i].typ = CodeTokenType.UnaryLHSOperator
          self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
          self.code_tokens[i].rhs_spaces_desired = Spacing.No
      prev_type = self.code_tokens[i].typ

    if self.code_tokens[-1].typ == CodeTokenType.Comment:
      if self.code_tokens[-1].token.startswith('/*') and self.code_tokens[-1].token.endswith('*/'):
        self.code_tokens[-1].token = '//' + self.code_tokens[-1].token[2:-2]

    number_template_or_gtlt = 0
    number_unary_or_binary = 0
    for tok in self.code_tokens:
      if tok.typ == CodeTokenType.GTOrTemplate or tok.typ == CodeTokenType.LTOrTemplate:
        number_template_or_gtlt += 1
      elif tok.typ == CodeTokenType.BinaryOrUnaryOperator:
        number_unary_or_binary += 1

    prev_token_typ = self.LabelObviousTemplates(prev_token_typ, prev_template_depths)
    number_template_or_gtlt_b = 0
    for tok in self.code_tokens:
      if tok.typ == CodeTokenType.GTOrTemplate or tok.typ == CodeTokenType.LTOrTemplate:
        number_template_or_gtlt_b += 1
        print "Unresolved token: " + tok.token + " in " + line
    CODESTATS.gtlt_resolved += number_template_or_gtlt - number_template_or_gtlt_b
    CODESTATS.gtlt_unresolved += number_template_or_gtlt_b

    prev_type = prev_token_typ
    for i in xrange(len(self.code_tokens) - 1):
      if self.code_tokens[i].typ == CodeTokenType.UnaryUnkOperator:
        if prev_type == CodeTokenType.UnaryRHSOperator or prev_type == CodeTokenType.Variable or self.code_tokens[i + 1].typ == CodeTokenType.SemiColon:
          self.code_tokens[i].typ = CodeTokenType.UnaryRHSOperator
          self.code_tokens[i].lhs_spaces_desired = Spacing.No
          self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
        elif CodeTokenType.is_unary_lhs_type[prev_type] or self.code_tokens[i + 1].typ == CodeTokenType.Variable or self.code_tokens[i + 1].typ == CodeTokenType.Function or self.code_tokens[i + 1].typ == CodeTokenType.UnaryUnkOperator:
          self.code_tokens[i].typ = CodeTokenType.UnaryLHSOperator
          self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
          self.code_tokens[i].rhs_spaces_desired = Spacing.No
        else:
          print "Ambiguous LHS/RHS line: " + line
      prev_type = self.code_tokens[i].typ

    if self.code_tokens[-1].typ == CodeTokenType.UnaryUnkOperator:
      if prev_type == CodeTokenType.UnaryRHSOperator or prev_type == CodeTokenType.Variable:
        self.code_tokens[-1].typ = CodeTokenType.UnaryRHSOperator
        self.code_tokens[-1].lhs_spaces_desired = Spacing.No
        self.code_tokens[-1].rhs_spaces_desired = Spacing.Yes
      elif CodeTokenType.is_unary_lhs_type[prev_type]:
        self.code_tokens[-1].typ = CodeTokenType.UnaryLHSOperator
        self.code_tokens[-1].lhs_spaces_desired = Spacing.Yes
        self.code_tokens[-1].rhs_spaces_desired = Spacing.No
      else:
        print "Ambiguous LHS/RHS line: " + line

    self.min_indent = self.state.GetMinIndent(self.code_tokens[ 0].typ)
    self.is_indent_definite = self.state.IsIndentAbsolute(self.code_tokens[ 0].typ)
    for i in xrange(len(self.code_tokens)):
       old_token_type = self.code_tokens[i].typ
       new_token_type = self.state.update(self.code_tokens, i)
       if old_token_type != new_token_type:
         self.code_tokens[i].typ = new_token_type
         if new_token_type == CodeTokenType.BinaryOperator:
           self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
           self.code_tokens[i].rhs_spaces_desired = Spacing.Yes
           #print "Token " + str(i) + " was a binary operator"
         elif new_token_type == CodeTokenType.UnaryLHSOperator:
           self.code_tokens[i].lhs_spaces_desired = Spacing.Yes
           self.code_tokens[i].rhs_spaces_desired = Spacing.No
           #print "Token " + str(i) + " was a unary operator"
         elif new_token_type == CodeTokenType.CCastFenceClose:
           #print "Token " + str(i) + " was a c-cast close parens"
           if i + 1 < len(self.code_tokens) and (self.code_tokens[i + 1].token == '&' or self.code_tokens[i + 1] == '*'):
             self.code_tokens[i].rhs_spaces_desired = Spacing.No

    for i in xrange(len(self.code_tokens)):
      if self.code_tokens[i].typ == CodeTokenType.BinaryOrUnaryOperator:
        print "Unresolved token in " + line
        break

    # recount the # of unresolved symbols

    number_unary_or_binary_b = 0
    for tok in self.code_tokens:
      if tok.typ == CodeTokenType.BinaryOrUnaryOperator:
        number_unary_or_binary_b += 1
    CODESTATS.binunary_resolved += number_unary_or_binary - number_unary_or_binary_b
    CODESTATS.binunary_unresolved += number_unary_or_binary_b

    for i in xrange(1, len(self.code_tokens)):
      lhs_desired_desired = self.code_tokens[i - 1].rhs_spaces_desired
      if lhs_desired_desired == Spacing.Unk:
        self.code_tokens[i - 1].rhs_spaces_desired = self.code_tokens[i].lhs_spaces_desired
      elif self.code_tokens[i].lhs_spaces_desired == Spacing.Unk:
        self.code_tokens[i].lhs_spaces_desired = lhs_desired_desired
      elif self.code_tokens[i].lhs_spaces_desired != lhs_desired_desired:
        if (self.code_tokens[ i - 1].typ != CodeTokenType.FenceOpen and self.code_tokens[i - 1].typ != CodeTokenType.FunctionPointerParen) or self.code_tokens[ i].token != ';':
          self.code_tokens[i - 1].rhs_spaces_desired = self.code_tokens[i].lhs_spaces_desired = Spacing.No
        else:
          self.code_tokens[i - 1].rhs_spaces_desired = self.code_tokens[i].lhs_spaces_desired = Spacing.Yes

# given a line containing #include, get the file that follows, e.g.
# #include "vector.h" -> vector.h
# #include <vector> -> vector
# there may be space before the #include and after the " or >, so don't rely on absolute positions
def get_header_part(x):
  include_pos = x.find("#include ")
  if include_pos < 0:
    return ""

  first_char = x[ include_pos + len("#include ")]
  last_char = ""
  if first_char == '"':
    last_char = '"'
  elif first_char == '<':
    last_char = '>'
  else:
    return ""
  end_pos = x.find(last_char, include_pos + len("#include ") + 1)
  if end_pos < 0:
    end_pos = len(x)

  return x[ include_pos + len("#include ") + 1: end_pos:1]

# split a file name into a file and a suffix part
# e.g. vector.h -> vector , h
# e.g. vector.fwd.hh -> vector , fwd.hh
# e.g. vector -> vector, (empty string)
# e.g. (empty string) -> (empty string), (empty string)
def split_file_suffix(x):
  if len(x) == 0:
    return "", ""
  suffix_pos = x.find('.')
  if x < 0:
    return x, ""
  if suffix_pos == len(x) - 1:
    return x[ 0:suffix_pos:1], ""
  return x[ 0:suffix_pos:1], x[ suffix_pos + 1::1]

# given two filenames, return -1,0,1 depending on whether the filenames are <,=,> than each other, respectively
# comparison is lexicographically, first by file suffix, then by file name
def compare_header_strings(str_a, str_b):
  header_part_a = get_header_part(str_a)
  header_part_b = get_header_part(str_b)
  #print "a: " + header_part_a + " b: " + header_part_b
  if header_part_a == header_part_b:
    return 0
  a_suffix_pos = header_part_a.find('.')
  b_suffix_pos = header_part_b.find('.')
  suffix_part_a = ""
  suffix_part_b = ""
  if a_suffix_pos >= 0 and b_suffix_pos >= 0:
    suffix_part_a = header_part_a[a_suffix_pos + 1::]
    suffix_part_b = header_part_b[b_suffix_pos + 1::]
    header_part_a = header_part_a[0:a_suffix_pos:1]
    header_part_b = header_part_b[0:b_suffix_pos:1]
  elif a_suffix_pos >= 0:
    return 1
  elif b_suffix_pos >= 0:
    return (-1)
  if header_part_a > header_part_b:
    return 1
  elif header_part_a < header_part_b:
    return (-1)
  elif suffix_part_a > suffix_part_b:
    return 1
  elif suffix_part_a < suffix_part_b:
    return (-1)
  return 0

# given a list of lines (e.g. a header block), check whether they are already sorted
def HeadersAreSorted(lines):
  if len(lines) < 2:
    return True
  for i in xrange(1, len(lines)):
    if compare_header_strings(lines[i - 1], lines[i]) > 0:
      return False
  return True

# given a line of code and a set of namespaces that we are currently in, remove
# the unnecessary namespace qualifier from all functions, e.g. 
# x::DoOtherStuff(); 
# in namespace x, this function will remove the redundant namespace qualifier:
# DoStuff();
def removeUnneededNamespaces(line, namespaces):
  if len(namespaces) == 0 or line.count(namespaces[-1] + "::") == 0:
    return line

  end_line_position = line.find('//')
  quote_position = line.find('"')
  if quote_position >= 0:
    if end_line_position < 0:
      end_line_position = quote_position
    elif end_line_position > quote_position:
      end_line_position = quote_position
  if end_line_position < 0:
    end_line_position = len(line)
  #print "namespaces: " + '::'.join(namespaces) + "::"
  new_line = ""
  last_scope_end_pos = 0
  next_scope_end_pos = line.find('::')
  while next_scope_end_pos > 0 and next_scope_end_pos < end_line_position:
    next_scope_pos = next_scope_end_pos - 1
    while next_scope_pos > 0 and isValidAlpha(line[next_scope_pos - 1]):
      next_scope_pos -= 1
    new_line += line[last_scope_end_pos:next_scope_pos:1]
    last_scope_end_pos = next_scope_pos
    namespace_name = line[ next_scope_pos:next_scope_end_pos:1]
    namespace_name_pos = findIndex(namespaces, namespace_name)
    name_was_overqualified = 0
    if namespace_name_pos >= 0:
      overqualified_name = '::'.join(namespaces[namespace_name_pos:]) + "::"
      #print "okay"
      if len(overqualified_name) + next_scope_pos <= len(line):
        #print line[next_scope_pos:len(overqualified_name) + next_scope_pos:1]
        if line.find(overqualified_name, next_scope_pos, len(overqualified_name) + next_scope_pos) == next_scope_pos:
          next_scope_end_pos = next_scope_pos + len(overqualified_name)
          last_scope_end_pos = next_scope_end_pos
          name_was_overqualified = 1
          #print "overqualified name: "
    if name_was_overqualified == 0:
      next_scope_end_pos += 2
      new_line += line[last_scope_end_pos:next_scope_end_pos:1]
      last_scope_end_pos = next_scope_end_pos
    if next_scope_end_pos >= end_line_position:
      next_scope_end_pos = -1
    if next_scope_end_pos > 0:
      next_scope_end_pos = line.find('::', last_scope_end_pos)
  new_line += line[last_scope_end_pos::1]
  #if new_line != line:
    #print "Line was: " + line + " now:     " + new_line
  return new_line

def GetCppStringRange(line):
  first_non_space = 0
  while first_non_space < len(line) and isspace(line[first_non_space]):
    first_non_space += 1
  if first_non_space == len(line) or line[first_non_space] == '#':
    return -1, -1
  comment_pos = line.find('//', first_non_space)
  if comment_pos < 0:
    comment_pos = len(line)
  #if line[first_non_space:].startswith('typedef') or line[first_non_space:].startswith('e_'):
  #  return -1, -1
  return first_non_space, comment_pos

def findAll(line, x):
  p = []
  i = line.find(x)
  while i >= 0:
    p.append(i)
    if i + 1 < len(line):
      i = line.find(x, i + 1)
    else:
      break
  return p

def largeSpacesAreAligned(line_a, line_b):

  line_a = line_a.strip()
  line_b = line_b.strip()
  ln = min(len(line_a), len(line_b))

  if ln <= 5 or (line_a.find('  ') < 0 or line_b.find('  ') < 0):
    return False

  # find all spaces in a and b
  spaces_a = findAll(line_a, ' ')
  spaces_b = findAll(line_b, ' ')

  space_ranges_a = [ [ spaces_a[ 0], spaces_a[ 0]]]
  i = 1
  for i in xrange(1, len(spaces_a)):
    if spaces_a[i] == space_ranges_a[-1][1] + 1:
      space_ranges_a[-1][1] += 1
    else:
      space_ranges_a.append([ spaces_a[ i], spaces_a[ i]])

  # remove trivial space ranges (size == 1)
  space_ranges_a = [ x for x in space_ranges_a if x[0] != x[1] and x[1] < ln]

  space_ranges_b = [ [ spaces_b[ 0], spaces_b[ 0]]]
  i = 1
  for i in xrange(1, len(spaces_b)):
    if spaces_b[i] == space_ranges_b[-1][1] + 1:
      space_ranges_b[-1][1] += 1
    else:
      space_ranges_b.append([ spaces_b[ i], spaces_b[ i]])
  # remove trivial space ranges (size == 1)
  space_ranges_b = [ x for x in space_ranges_b if x[0] != x[1] and x[1] < ln]

  # At least one non-trivial space in A must start at a non-trivial space in B
  non_trivial_spaces_starts = set([x[0] for x in space_ranges_a])
  non_trivial_spaces_starts.intersection_update(set([x[0] for x in space_ranges_b]))

  if len(non_trivial_spaces_starts):
    #print "Line " + line_a + " and " + line_b + " are aligned by space starts"
    return True
  non_trivial_spaces_ends = set([x[1] for x in space_ranges_a])
  non_trivial_spaces_ends.intersection_update(set([x[1] for x in space_ranges_b]))
  if len(non_trivial_spaces_ends):
    #print "Line " + line_a + " and " + line_b + " are aligned by space ends"
    return True
  #print "Line " + line_a + " and " + line_b + " are not aligned by long spaces"
  return False


def linesAreAligned(line_a, line_b):
  # if both have any of the following, and they are all aligned, but not aligned in the properly spaced string, then the alignment was deliberate:
  # // comments.  Aligned comments, by themselves, are enough to consider lines aligned
  # ,'s Spaces can also match commas, and commas preceded by )}], need not match
  # ([  
  # +-*/&|%=!^ For this match, all characters in the set (including space) are considered equivalent
  # Nonessential, but sufficient is matching ' m_', ' s_', ' g_', all of which are considered equivalent, provided there
  # is exactly one on each line

  if len(line_a) == 0 or len(line_b) == 0:
    return False

  basic_punct = '+-*/&|%=!^'

  # first, X out everything in quotes
  line_a_copy = deepcopy(line_a.rstrip())
  line_b_copy = deepcopy(line_b.rstrip())
  if len(line_a_copy) > len(line_b_copy):
    line_a_copy = line_a[:len(line_b_copy)]
  elif len(line_a_copy) < len(line_b_copy):
    line_b_copy = line_b[:len(line_a_copy)]
  ln = len(line_a_copy)
  if HasQuotes(line_a_copy):
    line_a_copy = OverWriteQuotes(line_a_copy, 'X')
  if HasQuotes(line_b_copy):
    line_b_copy = OverWriteQuotes(line_b_copy, 'X')

  # check for comments
  comment_a_pos = line_a_copy.find('//')
  if comment_a_pos < 0:
    comment_a_pos = line_a_copy.find('/*')
  comment_b_pos = line_b_copy.find('//')
  if comment_b_pos < 0:
    comment_b_pos = line_b_copy.find('/*')
  if comment_a_pos >= 0 and comment_b_pos >= 0:
    #print "Line " + line_a_copy + " is comment-aligned with " + line_b_copy
    return comment_a_pos == comment_b_pos
  if comment_a_pos >= 0:
    line_a_copy = line_a_copy[:comment_a_pos].rstrip()
    ln = len(line_a_copy)
    line_b_copy = line_b_copy[:ln]
  elif comment_b_pos >= 0:
    line_b_copy = line_b_copy[:comment_b_pos].rstrip()
    ln = len(line_b_copy)
    line_a_copy = line_a_copy[:ln]

  if ln < 5:
    return False

  # check whether they match simply
  if largeSpacesAreAligned(line_a_copy, line_b_copy):
    return True

  shift_pos_a = findAll(line_a_copy, '{')
  shift_pos_a.extend(findAll(line_a_copy, '}'))
  shift_pos_b = findAll(line_b_copy, '{')
  shift_pos_b.extend(findAll(line_b_copy, '}'))
  if len(shift_pos_a) or len(shift_pos_b):
    matches = shift_pos_a == shift_pos_b
    #if matches:
    #  print "Line " + line_a_copy + " matches { in " + line_b_copy
    #else:
    #  print "Line " + line_a_copy + " failed match on { in " + line_b_copy
    return matches

  trimmed_char = "X"
  if line_a_copy[-1] != line_b_copy[-1]:
    if (line_a_copy[-1] == ',' or line_a_copy[-1] == ';'):
      if not (line_b_copy[-1] == ',' or line_b_copy[-1] == ';'):
        ln -= 1
        trimmed_char = line_a_copy[-1]
        line_a_copy = line_a_copy[:-1]
        line_b_copy = line_b_copy[:-1]
    elif (line_b_copy[-1] == ',' or line_b_copy[-1] == ';'):
      ln -= 1
      trimmed_char = line_b_copy[-1]
      line_a_copy = line_a_copy[:-1]
      line_b_copy = line_b_copy[:-1]

  # comments were of no help.  Try matching commas
  found_matching = False
  paren_depth_a = 0
  paren_depth_b = 0
  for i in xrange(ln):
    delta_paren_a = 0
    if line_a_copy[i] == '(':
      delta_paren_a = 1
    elif line_a_copy[i] == ')':
      delta_paren_a = -1
    delta_paren_b = 0
    if line_b_copy[i] == '(':
      delta_paren_b = 1
    elif line_b_copy[ i] == ')':
      delta_paren_b = -1
    if delta_paren_a * delta_paren_b == -1:
      #print "Line " + line_a_copy + " is parenthesis-misaligned with " + line_b_copy
      return False
    paren_depth_a += delta_paren_a
    paren_depth_b += delta_paren_b
    if paren_depth_a != paren_depth_b:
      continue
    if line_a_copy[i] == ',':
      if line_b_copy[i] == ',':
        found_matching = True
      elif line_b_copy[i] == ' ':
        continue
      elif i + 1 < ln and line_a_copy[i + 1] == ' ':
        space_end = spaceEnd(line_a_copy, i + 1)
        comma_pos = line_b_copy.find(',', i)
        if comma_pos < 0:
          return False
        elif comma_pos < space_end and line_b_copy[comma_pos - 1] != ' ':
          continue
        else:
          #print "Line " + line_a_copy + " is comma-misaligned with " + line_b_copy
          return False
    elif line_b_copy[i] == ',':
      if line_a_copy[i] == ' ':
        continue
      elif i + 1 < ln and line_b_copy[i + 1] == ' ':
        space_end = spaceEnd(line_b_copy, i + 1)
        comma_pos = line_a_copy.find(',', i)
        if comma_pos < 0:
          return False
        elif comma_pos < space_end and line_a_copy[comma_pos - 1] != ' ':
          continue
        else:
          #print "Line " + line_a_copy + " is comma-misaligned with " + line_b_copy
          return False
  if paren_depth_a != paren_depth_b and min(paren_depth_a, paren_depth_b) != 0:
    #print "Line " + line_a_copy + " is parenthesis-different from " + line_b_copy
    return False
  if found_matching:
    #print "Line " + line_a_copy + " matches (, in " + line_b_copy
    return True

  for i in xrange(ln):
    if basic_punct.find(line_a_copy[i]) >= 0:
      if basic_punct.find(line_b_copy[i]) >= 0:
        found_matching = True
        continue
      elif line_b_copy[i] == ' ':
        continue
      #print "Line " + line_a_copy + " is different in basic punct from " + line_b_copy
      return False
    elif basic_punct.find(line_b_copy[i]) >= 0:
      if line_a_copy[i] == ' ':
        continue
      #print "Line " + line_a_copy + " is different in basic punct from " + line_b_copy
      return False
  if found_matching:
    #print "Line " + line_a_copy + " matches basic punct from " + line_b_copy
    return True

  shift_pos_a = findAll(line_a_copy, '<<')
  shift_pos_b = findAll(line_b_copy, '<<')
  if len(shift_pos_a) or len(shift_pos_b):
    matches = shift_pos_a == shift_pos_b
    #if matches:
    #  print "Line " + line_a_copy + " matches << in " + line_b_copy
    #else:
    #  print "Line " + line_a_copy + " failed match on << in " + line_b_copy
    return matches
  shift_pos_a = findAll(line_a_copy, '>>')
  shift_pos_b = findAll(line_b_copy, '>>')
  if len(shift_pos_a) or len(shift_pos_b):
    matches = shift_pos_a == shift_pos_b
    #if matches:
    #  print "Line " + line_a_copy + " matches >> in " + line_b_copy
    #else:
    #  print "Line " + line_a_copy + " failed match on >> in " + line_b_copy
    return matches

  shift_pos_a = findAll(line_a_copy, '(')
  shift_pos_b = findAll(line_b_copy, '(')
  if len(shift_pos_a) or len(shift_pos_b):
    matches = shift_pos_a == shift_pos_b
    #if matches:
    #  print "Line " + line_a_copy + " matches ( in " + line_b_copy
    #else:
    #  print "Line " + line_a_copy + " failed match on ( in " + line_b_copy
    return matches

  shift_pos_a = findAll(line_a_copy, '[')
  shift_pos_b = findAll(line_b_copy, '[')
  if len(shift_pos_a) or len(shift_pos_b):
    matches = shift_pos_a == shift_pos_b
    #if matches:
    #  print "Line " + line_a_copy + " matches [ in " + line_b_copy
    #else:
    #  print "Line " + line_a_copy + " failed match on [ in " + line_b_copy
    return matches

  # first, X out everything in quotes
  line_a_copy = deepcopy(line_a.rstrip())
  line_b_copy = deepcopy(line_b.rstrip())
  if HasQuotes(line_a_copy):
    line_a_copy = OverWriteQuotes(line_a_copy, 'X')
  if HasQuotes(line_b_copy):
    line_b_copy = OverWriteQuotes(line_b_copy, 'X')
  if comment_a_pos >= 0:
    line_a_copy = line_a_copy[:comment_a_pos].rstrip()
  if comment_b_pos >= 0:
    line_b_copy = line_b_copy[:comment_b_pos].rstrip()

  if line_a_copy[-1] == ';' or line_a_copy[-1] == ',':
    line_a_copy = line_a_copy[:-1].rstrip()
  if line_b_copy[-1] == ';' or line_b_copy[-1] == ',':
    line_b_copy = line_b_copy[:-1].rstrip()

  xa = line_a_copy.rfind(' ')
  xb = line_b_copy.rfind(' ')
  if xa != xb:
    #print "Line " + line_a_copy + " no match in " + line_b_copy
    return False
  var_a_end = nextNameEnd(line_a_copy, xa + 1)
  var_b_end = nextNameEnd(line_b_copy, xb + 1)
  matches = var_a_end == var_b_end and min(var_a_end, var_b_end) == min(len(line_a_copy), len(line_b_copy))
  #if matches:
  #  print "Line " + line_a_copy + " matches variable in " + line_b_copy
  #else:
  #  print "Line " + line_a_copy + " no match in " + line_b_copy
  return matches

copyright_block = []
max_window_sz = 2

class BclCleanerStatus:
  def __init__(self, total_files, show_status, show_progress):
    self.total_files = total_files
    self.file_count = 0
    self.show_status = show_status
    self.last_percent = 0
    self.last_line_length = 0
    self.update_lock = threading.Lock()
    self.show_progress = show_progress

  def update(self, number_errors, err_string, filename, nfiles = 1):
    self.update_lock.acquire()
    if len(err_string):
      print err_string
    if number_errors:
      self.last_line_length = 0
      global CODESTATS
      CODESTATS.number_errors += number_errors
    self.file_count += nfiles
    if self.show_status:
      percent = 100.0 * self.file_count / float(self.total_files)
      stars = int(percent / 5)
      star_str = '*' * stars
      space_str = ' ' * (20 - stars)
      next_line = "\rBcl Cleaner progress [" + star_str + space_str + '] ' + str(int(percent)) + '% '
      next_line += filename
      additional_spaces = ""
      if len(next_line) < self.last_line_length:
        additional_spaces = ' ' * (self.last_line_length - len(next_line))
      sys.stdout.write(next_line + additional_spaces)
      self.last_line_length = len(next_line)
      sys.stdout.flush()
    elif self.show_progress:
      percent = int(100.0 * self.file_count / float(self.total_files))
      if percent >= self.last_percent + 10:
        sys.stdout.write(str(percent) + "% ")
        sys.stdout.flush()
        self.last_percent = percent
    self.update_lock.release()

status_bar = BclCleanerStatus(0, False, False)

# clean a file with a given namespace and filename, and whether or not to overwrite
def CleanFile(namespace, filename, should_overwrite):
  spacing_str = "Should replace "
  if should_overwrite:
    spacing_str = "Replaced "
  real_namespace = namespace[2::1]
  file_namespace = real_namespace
  output_stream = cStringIO.StringIO()
  if real_namespace.startswith('include/'):
    real_namespace = real_namespace[len('include/'):]

  # open the file, examine the includes
  ifile = open(namespace + os.sep + filename, 'r')
  morelines = ifile.readlines()
  ifile.close()
  
  is_example = real_namespace.startswith('example/') and len(real_namespace.split('/')) > 2
  must_use_example_macro = is_example and filename != 'example_io_file.cpp'
  log_prefix = "File: " + filename + ' '
  
  can_clean_namespaces = not filename.endswith("bcl_linal_operations_cpu.cpp")
  
  if morelines == None or len(morelines) == 0:
    output_stream.write(log_prefix + " is empty and should be removed!")

  last_line_was_empty = False
  namespace_includes = []
  bcl_includes = []
  std_includes = []
  header_to_suffices = {}

  included_headers = set()
  have_includes = 0
  n_empty_lines_removed = 0
  n_sorted_header_blocks = 0
  n_stripped_header_names = 0
  n_duplicate_headers = 0
  n_lines_with_redundant_namespaces = 0
  n_lines_wrong_macro = 0
  n_bad_spacing = 0
  n_long_comments = 0
  n_bad_header_guards = 0
  n_tabs = 0
  n_fixed_copyright = 0
  for line in morelines:
    n_tabs += line.count('\t')
  if n_tabs > 0:
    morelines = [ line.replace('\t', '  ') for line in morelines]
  current_namespaces = []
  current_namespace_line = []
  bracket_depth = 0
  need_namespace_bracket = 0
  forward_header_count = 0
  preprocessor_if_depth = 0
  header_guard_lines = [ -1, -1, -1 ]
  bcl_copyright_line_start = -1
  bcl_copyright_line_end = -1
  bcl_copyright_prefix = "// (c)"
  original_lines = []
  prev_include_block_size = 0
  for line in morelines:
    lstripped_line = line.lstrip()

    n_tabs += line.count('\t')
    # test for bcl copyright
    if header_guard_lines[0] < 0:
      if ((len(lstripped_line) == 0 and bcl_copyright_line_start >= 0) or lstripped_line.startswith(bcl_copyright_prefix)) and ( bcl_copyright_line_start < 0 or bcl_copyright_line_end == len(original_lines) - 1):
        if bcl_copyright_line_start < 0:
          bcl_copyright_line_start = len(original_lines)
        bcl_copyright_line_end = len(original_lines)
        bcl_copyright_line_pos = bcl_copyright_line_end - bcl_copyright_line_start
        if bcl_copyright_line_pos >= len(copyright_block):
          n_fixed_copyright += 1
        else:
          if copyright_block[bcl_copyright_line_pos].rstrip() != lstripped_line.rstrip() and (len(copyright_block[bcl_copyright_line_pos]) > 1 or len(lstripped_line) > 0):
            n_fixed_copyright += 1
          original_lines.append(copyright_block[bcl_copyright_line_pos])
        continue
      bcl_copyright_line_pos = bcl_copyright_line_end - bcl_copyright_line_start
      if bcl_copyright_line_pos+1 < len(copyright_block):
        output_stream.write("Copyright too short; was " + str(bcl_copyright_line_pos) + " in length but should have had " + str(len(copyright_block)) + " lines")
        if bcl_copyright_line_end < 0:
          original_lines.extend(copyright_block)
        else:
          original_lines.extend(copyright_block[bcl_copyright_line_pos+1:])
        n_fixed_copyright += 1
        bcl_copyright_line_end = len(copyright_block) + bcl_copyright_line_start
      
    # check for header guard lines
    if lstripped_line.startswith('#'):
      # pre-processor commands
      comment_pos = line.find('/*')
      if comment_pos >= 0 and line.rstrip().endswith('*/'):
        old_line = line
        comment = line.rstrip()[ comment_pos + 2:-2].rstrip()
        line = line[:comment_pos] + '//' + comment + '\n'
        lstripped_line = line.lstrip()
        n_bad_spacing += 1
        output_stream.write("Replaced " + old_line + "with " + line)

      if lstripped_line.startswith('#if'):
        preprocessor_if_depth += 1
        included_headers = set()
        if header_guard_lines[ 0] < 0 and preprocessor_if_depth == 1 and lstripped_line.startswith('#ifndef'):
          header_guard_lines[ 0] = len(original_lines)
      elif lstripped_line.startswith('#define'):
        if preprocessor_if_depth == 1 and header_guard_lines[ 1] < 0:
          header_guard_lines[ 1] = len(original_lines)
      elif lstripped_line.startswith('#endif'):
        preprocessor_if_depth -= 1
        included_headers = set()
        if preprocessor_if_depth == 0 and header_guard_lines[ 2] < 0:
          header_guard_lines[ 2] = len(original_lines)
      elif lstripped_line.startswith("#else") or lstripped_line.startswith("#elif"):
        included_headers = set()
      elif lstripped_line.startswith("#include "):
        if last_line_was_empty:
          if len(original_lines) > 1 and ( original_lines[-2].startswith('#include') or original_lines[-2].startswith('//')):
            verb = " Should delete" if not should_overwrite else " Deleted"
            output_stream.write( filename + verb + " blank line between includes\n")
            if should_overwrite:
              n_empty_lines_removed += 1
              original_lines.pop()
          else:
            output_stream.write( filename + " Should add a comment before include block starting with " + lstripped_line)
        last_line_was_empty = False
        num_spaces_start = len(line) - len(lstripped_line)
        stripped_line = lstripped_line.rstrip()[ len("#include")::1].strip()

        if len(stripped_line) == 0:
          original_lines.append(line)
          continue

        first_char = stripped_line[ 0]

        if first_char != '"' and first_char != '<':
          original_lines.append(line)
          continue

        last_char = ''
        if first_char == '"':
          last_char = '"'
        else:
          last_char = '>'

        post_header_component = ""
        end_pos = stripped_line.find(last_char, 1)
        if end_pos < 0:
          end_pos = len(stripped_line)
        else:
          post_header_component = stripped_line[ end_pos + 1::1]
        header_component = stripped_line[ 1 : end_pos : 1]

        if len(real_namespace) and header_component.startswith(real_namespace):
          header_component = header_component[ len(real_namespace)::1]
          n_stripped_header_names += 1
        have_includes = 1
        if header_component not in included_headers:
          revised_header_line = line[ 0:num_spaces_start:1] + "#include " + first_char + header_component + last_char + post_header_component + '\n'
          if first_char == '<':
            std_includes.append(revised_header_line)
          elif header_component.find('/') < 0:
            # if there are already any non-namespace includes, then the headers will end up sorted
            if len(std_includes) + len(bcl_includes) > 0:
              n_sorted_header_blocks += 1
            namespace_includes.append(revised_header_line)
          else:
            if len(std_includes) > 0:
              n_sorted_header_blocks += 1
            bcl_includes.append(revised_header_line)
          included_headers.add(header_component)
          file_name, file_suffix = split_file_suffix(header_component)
          #print "filename: " + file_name + " suffix: " + file_suffix
          if file_name not in header_to_suffices:
            header_to_suffices[ file_name] = set()
            header_to_suffices[ file_name].add(file_suffix)
          elif file_suffix not in header_to_suffices[ file_name] and header_component != filename + "pp":
            header_to_suffices[ file_name].add(file_suffix)
            output_stream.write("File: " + filename \
                  + " Warning: (likely) redundant include: it should only be necessary to include one of " \
                  + " and ".join([ file_name + '.' + suffix for suffix in header_to_suffices[ file_name]]))
          if file_suffix == 'fwd.hh':
            forward_header_count += 1
        else:
          if should_overwrite:
            output_stream.write("File: " + filename + " Log: Removed duplicate header " + header_component)
          else:
            output_stream.write("File: " + filename + " Log: Should remove duplicate header " + header_component)
          n_duplicate_headers += 1
        continue

    if have_includes == 1:
      have_includes = 0
      # write out the header lines, sorted
      #print ''.join(namespace_includes)
      if len(namespace_includes) > 0:
        if not HeadersAreSorted(namespace_includes):
          namespace_includes = sorted(namespace_includes, compare_header_strings)
          n_sorted_header_blocks += 1
        original_lines.extend(namespace_includes)

      #print ''.join(bcl_includes)
      if len(bcl_includes) > 0:
        if not HeadersAreSorted(bcl_includes):
          bcl_includes = sorted(bcl_includes, compare_header_strings)
          n_sorted_header_blocks += 1
        original_lines.extend(bcl_includes)

      #print ''.join(std_includes)
      if len(std_includes) > 0:
        if not HeadersAreSorted(std_includes):
          std_includes = sorted(std_includes, compare_header_strings)
          n_sorted_header_blocks += 1
        original_lines.extend(std_includes)

      if len(std_includes) > 0 and len(bcl_includes) > 0:
          output_stream.write("File: " + filename \
                + " Warning: external library includes mixed with bcl includes; please separate them with a comment")

      prev_include_block_size = len(namespace_includes) + len(bcl_includes)

      file_name, file_suffix = split_file_suffix(filename)
      if forward_header_count != 0:
        if file_suffix.endswith('pp'):
          output_stream.write("File: " + filename + " Warning: implementation files should not need forward headers")
        elif len(bcl_includes) + len(namespace_includes) != forward_header_count:
          output_stream.write(\
          "File: " + filename + " Warning: forward headers mixed with bcl headers; they should be separated by a comment")

      allowed_app_includes = set(['bcl_app_interface', 'bcl_app_apps', 'bcl_app_groups', 'bcl_app_group_handler', 'bcl_app_interface_release'])
      # check for bcl_app_ inclusions outside the cpp and the example
      for header in included_headers:
        if header.startswith('bcl_app_') or header.find('/bcl_app_') >= 0:
          header_name, header_suffix = split_file_suffix(header)
          header_name = header_name.split('/')[-1]
          if file_name != header_name and not file_name.startswith('example') and header_name not in allowed_app_includes:
            output_stream.write(\
              "ERROR: File: " + filename + " includes " + header \
              + " but only examples and the associated source file are allowed to include application headers");

      namespace_includes = []
      bcl_includes = []
      std_includes = []
      forward_header_count = 0

    if lstripped_line.rstrip().startswith("namespace "):
      prev_include_block_size = 0
      need_namespace_bracket = 1
      end_name_pos = nextVarNameEnd(lstripped_line, len("namespace "))
      namespace_name = lstripped_line[len("namespace "):end_name_pos:1]
      current_namespaces.append(namespace_name)
      current_namespace_line.append(len(original_lines))

    line_sans_quotes = stripComments([ OverWriteQuotes(line, 'X') ])
    if len(line_sans_quotes) == 0:
      line_sans_quotes = ""
    else:
      line_sans_quotes = line_sans_quotes[0]
    open_bracket_count = line_sans_quotes.count('{')
    close_bracket_count = line_sans_quotes.count('}')
    if need_namespace_bracket == 1 and open_bracket_count > 0:
      need_namespace_bracket = 0
    bracket_depth += open_bracket_count - close_bracket_count
    if bracket_depth < 0:
      print '\n' + filename + " has a problem with brackets closure! On line: " + line
      bracket_depth = 0

    if bracket_depth + need_namespace_bracket < len(current_namespaces):
      leaving_namespace = current_namespaces[bracket_depth + need_namespace_bracket]
      leaving_namespace_line = current_namespace_line[bracket_depth + need_namespace_bracket]
      current_namespaces = current_namespaces[:bracket_depth + need_namespace_bracket]

      old_line = line
      line = ' ' * (2 * len(current_namespaces)) + '} // namespace ' + leaving_namespace + '\n'
      if old_line.rstrip() != line.rstrip():
        if len(original_lines) - leaving_namespace_line >= 10 or old_line.find('//') > 0:
          lstripped_line = line.lstrip()
          output_stream.write("\nBad namespace comment, was " + old_line + " now: " + line)
          n_bad_spacing += 1
    if last_line_was_empty and len(lstripped_line) == 0:
      n_empty_lines_removed += 1
      continue
    if len(lstripped_line) == 0:
      last_line_was_empty = True
    else:
      last_line_was_empty = False

    namespace_cleaned_line = removeUnneededNamespaces(line, current_namespaces) if can_clean_namespaces else line

    if namespace_cleaned_line != line:
      n_lines_with_redundant_namespaces += 1

    if must_use_example_macro and \
      ( \
        namespace_cleaned_line.find('io::File::MustOpenIFStream') >= 0 \
        or namespace_cleaned_line.find('io::File::MustOpenOFStream') >= 0 \
      ):
      n_lines_wrong_macro += 1
      original_lines.append \
      ( \
        namespace_cleaned_line.replace('io::File::MustOpenIFStream','BCL_ExampleMustOpenInputFile') \
        .replace('io::File::MustOpenOFStream','BCL_ExampleMustOpenOutputFile') \
        .replace('\r', '')
      )
    else:
      original_lines.append(namespace_cleaned_line.replace('\r', ''))

  StatusSkipped = 0
  StatusReplaced = max_window_sz + 1

  # line statuses:
  #   0 = skipped, ignored, or line was same w/ vs w/o spacing
  #   -max_window_sz <= x <= max_window_sz, x != 0, aligned with line x away
  #   max_window_sz+1 = was replaced
  line_status = [ StatusSkipped ] * len(original_lines)

  prev_token = CodeTokenType.Space
  prev_template_depths = []
  prev_stack = CodeState()
  freely_remove_spaces = [] # number of spaces that can be freely removed
  freely_remove_space_indent_depth = [] # indent depth for freely removable spaces
  prev_indentation_level = 0
  in_continued_line = False
  comment_block_line = -1
  for i in xrange(len(original_lines)):
    # line devoid of spaces
    stripped_line = original_lines[ i].strip()

    current_line_continues = False
    if len(stripped_line) and stripped_line[-1] == '\\':
      current_line_continues = True

    # skip empty lines, also preprocessor lines (begin with '#') and commented lines
    if len(stripped_line) == 0 or stripped_line[0] == '#':
      last_was_comment_block_start = False
      last_was_comment_block_mid_line = False
      in_comment_block = False
      in_continued_line = current_line_continues
      continue

    indent_change = 0
    actual_indent = original_lines[ i].find(stripped_line[0])

    if stripped_line.startswith('//'):
      prev_indentation_level = prev_stack.GetMinIndent(CodeTokenType.Comment)
      number_slashes = len(stripped_line) - len(stripped_line.strip('/'))
      if comment_block_line < 0 and number_slashes >= 4 and len(stripped_line) >= 7 and len(original_lines[i]) < 100 and i < len(original_lines) - 2:
        if number_slashes == len(stripped_line):
          next_line = original_lines[i + 2].strip()
          if len(next_line) >= 7 and next_line.count('/') == len(next_line):
            next_line = original_lines[i + 1].strip()
            if len(next_line) >= 7 and next_line.startswith('//') and next_line.endswith('//') and next_line.count('/') < len(next_line) and not next_line.find('//!') >= 0:
              number_slashes_last_line = len(next_line) - len(next_line.strip('/'))
              if number_slashes_last_line == 4 or number_slashes_last_line == 6:
                comment_block_line = 0
      if actual_indent or comment_block_line >= 0:
        while len(freely_remove_space_indent_depth) and freely_remove_space_indent_depth[-1] >= actual_indent:
          del freely_remove_space_indent_depth[-1]
          del freely_remove_spaces[-1]
        # comment with preceeding spaces, indent may need fixing
        expected_indent = prev_stack.GetMinIndent(CodeTokenType.Comment)
        prev_indentation_level = expected_indent
        if comment_block_line >= 0:
          comment_block_line += 1
          if comment_block_line == 3:
            comment_block_line = -1
          # comment block line, actual indent is 2 less
          expected_indent -= 2
          if expected_indent < 0:
            expected_indent = 0
          if actual_indent != expected_indent:
            if not len(freely_remove_space_indent_depth) or freely_remove_space_indent_depth[ -1] < actual_indent:
              freely_remove_spaces.append(actual_indent - expected_indent)
              freely_remove_space_indent_depth.append(actual_indent)
            space_added_line = ' ' * expected_indent + stripped_line
            print spacing_str + " line:\n" + original_lines[i] + "with:\n" + space_added_line
            original_lines[ i] = space_added_line + '\n'
            n_bad_spacing += 1
          in_continued_line = current_line_continues
          continue
        comment_block_line = -1
        if actual_indent < expected_indent:
          space_added_line = ' ' * expected_indent + stripped_line
          print spacing_str + " line:\n" + original_lines[i] + "with:\n" + space_added_line
          original_lines[ i] = space_added_line + '\n'
          n_bad_spacing += 1
        elif actual_indent > expected_indent and len(freely_remove_space_indent_depth):
          indent_after_removal = actual_indent - freely_remove_spaces[-1]
          if indent_after_removal >= expected_indent:
            space_added_line = ' ' * indent_after_removal + stripped_line
            print spacing_str + " line:\n" + original_lines[i] + "with:\n" + space_added_line
            original_lines[ i] = space_added_line + '\n'
            n_bad_spacing += 1
      in_continued_line = current_line_continues
      continue;

    while len(freely_remove_space_indent_depth) and freely_remove_space_indent_depth[-1] >= actual_indent:
      expected_indent = prev_stack.GetMinIndent(CodeTokenType.CurlyClose)
      if freely_remove_space_indent_depth[-1] == actual_indent and stripped_line[0] == '}':
        indent_change = actual_indent - expected_indent
      del freely_remove_space_indent_depth[-1]
      del freely_remove_spaces[-1]
    comment_block_line = -1
    # generate the optimally-spaced string
    bcl_line = BclCodeLine(original_lines[i], filename, prev_token, prev_template_depths, prev_stack)
    # update the previous token
    last_token_type = bcl_line.getLastCodeTokenType()
    if last_token_type != CodeTokenType.Space:
      prev_token = last_token_type

    if in_continued_line and actual_indent < 2:
      bcl_line.min_indent = max(bcl_line.min_indent, 2)

    if len(bcl_line.code_tokens) > 1 and bcl_line.code_tokens[0].typ == CodeTokenType.AccessSpec and bcl_line.code_tokens[1].typ == CodeTokenType.SpecialColon:
      expected_indent = bcl_line.min_indent
      if actual_indent != expected_indent:
        if not len(freely_remove_space_indent_depth) or freely_remove_space_indent_depth[ -1] < actual_indent:
          freely_remove_spaces.append(actual_indent - expected_indent)
          freely_remove_space_indent_depth.append(actual_indent)

    if bcl_line.code_spaces[ 0] < bcl_line.min_indent:
      indent_change = bcl_line.code_spaces[ 0] - bcl_line.min_indent
      bcl_line.code_spaces[ 0] = bcl_line.min_indent
    elif bcl_line.code_spaces[ 0] == bcl_line.min_indent and len(freely_remove_space_indent_depth):
      del freely_remove_space_indent_depth[-1]
      del freely_remove_spaces[-1]
    elif bcl_line.code_spaces[ 0] > bcl_line.min_indent and len(freely_remove_space_indent_depth):
      indent_after_removal = bcl_line.code_spaces[ 0] - freely_remove_spaces[-1]
      if indent_after_removal >= bcl_line.min_indent:
        bcl_line.code_spaces[ 0] = indent_after_removal
        indent_change = freely_remove_spaces[-1]
      else:
        del freely_remove_space_indent_depth[-1]
        del freely_remove_spaces[-1]
    if abs(bcl_line.min_indent - prev_indentation_level) > 2 and len(original_lines[i]) < 180:
      output_stream.write(log_prefix + " statement at line #" + str(i + 1) + " needs to be split up: " + original_lines[i])
    prev_indentation_level = bcl_line.min_indent
    prev_template_depths = bcl_line.ending_template_paren_depths
    prev_stack = bcl_line.state
    in_continued_line = current_line_continues

    global CODESTATS
    CODESTATS.lines_of_code += 1

    # check whether the line is already aligned or is aligned by the >= 3 spaces heuristic, e.g.
    # >= 3 spaces anywhere after the start of the line suggests it was deliberately aligned, so do not change anything
    optimized_spacing = False
    if line_status[i] == 0 and stripped_line.find('   ') < 0:
      optimized_spacing = True
      bcl_line.OptimizeSpacing()

    space_added_line = bcl_line.GetString()

    # check whether the strings differ, if not, no point in checking alignment
    # also skip lines that increased in size more than 7, because most likely these lines had omitted all the spaces for a reason
    if space_added_line.rstrip() == original_lines[i].rstrip() or (len(space_added_line) > len(original_lines[ i]) + 7 and space_added_line.count(';') == 0):
      continue

    if optimized_spacing:
      # check for alignment with other strings
      window_sz = 1
      while window_sz <= max_window_sz:
        if i >= window_sz and line_status[ i - window_sz] != StatusReplaced and linesAreAligned(original_lines[i], original_lines[i - window_sz]):
          line_status[ i] = -window_sz
          line_status[ i - window_sz] = window_sz
          break
        if i + window_sz < len(original_lines) and linesAreAligned(original_lines[i], original_lines[i + window_sz]):
          line_status[ i] = window_sz
          line_status[ i + window_sz] = -window_sz
          break
        window_sz += 1
      if line_status[ i] == StatusSkipped:
        output_stream.write(spacing_str + " line:\n" + original_lines[i] + "with:\n" + space_added_line)
        n_bad_spacing += 1
        line_status[ i] = StatusReplaced
        original_lines[ i] = space_added_line
      elif indent_change:
        space_added_line = ""
        if indent_change > 0:
          space_added_line = original_lines[i][ indent_change:]
        else:
          space_added_line = ' ' * -indent_change + original_lines[i]
        output_stream.write(spacing_str + " line:\n" + original_lines[i] + "with:\n" + space_added_line)
        n_bad_spacing += 1
        line_status[ i] = StatusReplaced
        original_lines[ i] = space_added_line
    else:
      output_stream.write(spacing_str + " line:\n" + original_lines[i] + "with:\n" + space_added_line)
      n_bad_spacing += 1
      line_status[ i] = StatusReplaced
      original_lines[ i] = space_added_line

  if len(original_lines) > 1 and filename.endswith('h'):
    log_str_add = log_prefix
    log_str_fix = log_prefix
    if should_overwrite:
      log_str_add += "Added "
      log_str_fix += "Fix "
    else:
      log_str_add += "Should add "
      log_str_fix += "Should fix "
    # determine what the header guard should be
    desired_guard_no_end_underscore = filename.split('/')[-1].upper().replace('.', '_')
    desired_header_guard = desired_guard_no_end_underscore + '_'

    if header_guard_lines[2] < 0:
      n_bad_header_guards += 1
      print log_str_add + "header guard"
      if should_overwrite:
        original_lines[ 0] = "#ifndef " + desired_header_guard + "\n#define " + desired_header_guard + '\n\n' + original_lines[0].lstrip()
        original_lines[-1] += '#endif // ' + desired_header_guard + '\n'
    else:

      # get the actual header guard variable
      actual_guard = original_lines[header_guard_lines[0]].split(' ')[-1].strip()
      actual_define = original_lines[header_guard_lines[1]].split(' ')[-1].strip()

      bad_ifndef = desired_guard_no_end_underscore != actual_guard.strip('_')
      check_comment = True
      if actual_guard == actual_define:
        if bad_ifndef:
          check_comment = False
          n_bad_header_guards += 1
          output_stream.write(log_str_fix + "header guard, was " + actual_guard + ', correct: ' + desired_header_guard)
          actual_guard = actual_define = desired_header_guard
          original_lines[ header_guard_lines[0]] = "#ifndef " + desired_header_guard + "\n"
          original_lines[ header_guard_lines[1]] = "#define " + desired_header_guard + '\n'
          original_lines[ header_guard_lines[2]] = '#endif // ' + desired_header_guard + '\n'
      elif not bad_ifndef and header_guard_lines[0] == header_guard_lines[1] - 1:
        original_lines[ header_guard_lines[1]] = "#define " + actual_guard + '\n'
        actual_define = actual_guard
        n_bad_header_guards += 1
        output_stream.write(log_str_fix + "header guard define, was " + actual_define + ', correct: ' + desired_header_guard)
      else:
        # the ifndef is wrong and the define did not immediately follow, this probably means
        n_bad_header_guards += 1
        output_stream.write(log_prefix + "Warning: Mismatched ifndef / define in header guard; cannot be fixed automatically")
        check_comment = False
      if check_comment:
        comment_pos = original_lines[ header_guard_lines[2]].find('//')
        if comment_pos >= 0:
          comment = original_lines[ header_guard_lines[2]][ comment_pos + 2:].rstrip()
          if actual_guard != comment.strip():
            original_lines[ header_guard_lines[2]] = '#endif // ' + actual_guard + '\n'
            n_bad_header_guards += 1
            output_stream.write(log_str_fix + "incorrect header guard endif comment")
        else:
          n_bad_header_guards += 1
          output_stream.write(log_str_fix + "missing header guard endif comment")
          original_lines[ header_guard_lines[2]] = '#endif // ' + actual_guard + '\n'
  if n_tabs > 0:
    log_str = log_prefix
    if should_overwrite:
      log_str += "Replaced "
    else:
      log_str += "Should replace "
    output_stream.write(log_str + str(n_tabs) + " tabs with spaces")
  if n_empty_lines_removed > 0:
    log_str = log_prefix
    if should_overwrite:
      log_str += "Removed "
    else:
      log_str += "Should remove "
    output_stream.write(log_str + str(n_empty_lines_removed) + " consecutive empty lines")
  if n_sorted_header_blocks > 0:
    log_str = log_prefix
    if should_overwrite:
      log_str += "Sorted"
    else:
      log_str += "Should sort"
    output_stream.write(log_str + " header blocks")
  if n_stripped_header_names > 0:
    log_str = log_prefix
    if should_overwrite:
      log_str += "Removed"
    else:
      log_str += "Should remove"
    output_stream.write(log_str + " unnecessary namespace-folder component from " + str(n_stripped_header_names) + " includes")
  if n_fixed_copyright:
    log_str = log_prefix
    if should_overwrite:
      log_str += "Fixed"
    else:
      log_str += "Should fixe"
    output_stream.write(log_str + " bcl copyright block. ")
  if n_lines_with_redundant_namespaces > 0:
    log_str = log_prefix
    if should_overwrite:
      log_str += "Removed"
    else:
      log_str += "Should remove"
    output_stream.write(log_str + " unnecessary namespace scope resolution component from " + str(n_lines_with_redundant_namespaces) + " lines")
  
  if n_lines_wrong_macro > 0:
    log_str = log_prefix
    if should_overwrite:
      log_str += "Switch"
    else:
      log_str += "Should switch"
    output_stream.write(log_str + " io::File::MustOpen[IO]Fstream commands in example to BCL_ExampleMustOpen{Input,Output}File on " + str(n_lines_wrong_macro) + " lines")

  if n_bad_spacing > 0:
    log_str = log_prefix
    if should_overwrite:
      log_str += "Fixed "
    else:
      log_str += "Should fix "
    output_stream.write(log_str + str(n_bad_spacing) + " lines with bad spacing")
#      if n_long_comments:
#        log_str = log_prefix
#        if should_overwrite:
#          log_str += "Changed "
#        else:
#          log_str += "Should change "
#        file_errs += n_long_comments
#        print log_str + str(n_long_comments) + " C-style comments (/**/) into C++ style comments (//)"
#      if len(original_lines[-1].strip()):
#        if should_overwrite:
#          log_str += "Added "
#        else:
#          log_str += "Should add "
#        file_errs += 1
#        print log_str + " an empty line at the end of the file"
  file_errs = n_tabs + n_lines_with_redundant_namespaces + n_empty_lines_removed + n_sorted_header_blocks + n_stripped_header_names + n_duplicate_headers + n_bad_spacing + n_bad_header_guards + n_fixed_copyright + n_lines_wrong_macro
  if should_overwrite and file_errs > 0:
    file_out = open(namespace + os.sep + filename, 'w')
    file_out.write(''.join(original_lines))
    file_out.close()
    output_stream.write("overwriting " + namespace + os.sep + filename)
  output_str = output_stream.getvalue()
  global status_bar
  status_bar.update(file_errs, output_str, filename, 1)

# actual path to this script, relative to the bcl path
script_name = "scripts/code/BclCleaner.py"

class ThreadCleanFile(threading.Thread):
    def __init__(self, queue, id_num):
        threading.Thread.__init__(self)
        self.queue = queue
        self.id_num = id_num

    def run(self):
      while True:
        #grabs info from the queue
        info = self.queue.get()
        directory = info[0]
        should_output = info[1]
        number_files = info[2]
        pickle_file = "BclCleaner." + str(self.id_num) + ".txt"
        cmd = script_name + " " + os.getcwd() + " --directory " + directory + " --nothreads --noprogress --stats " + pickle_file
        if should_output:
          cmd += ' -o '
        (status, output) = commands.getstatusoutput(cmd)

        if status == 0:
          # get the statistics too
          new_stats = CodeStats()
          while not os.path.exists(pickle_file):
            sleep(0.5)
          infile = open(pickle_file, 'r')
          line = infile.readline().strip().split(' ')
          try:
            new_stats.number_errors = int(line[0])
            new_stats.lines_of_code = int(line[1])
            new_stats.binunary_resolved = int(line[2])
            new_stats.binunary_unresolved = int(line[3])
            new_stats.gtlt_resolved = int(line[4])
            new_stats.gtlt_unresolved = int(line[5])
            os.remove(pickle_file)
          except:
            print "Bad output file at " + pickle_file
          global CODESTATS
          global status_bar
          CODESTATS.add(new_stats)
          status_bar.update(new_stats.number_errors, output, directory, number_files)
        else:
          status_bar.update(1, output, directory, number_files)

        #signals to queue job is done
        self.queue.task_done()

def CleanFiles(namespace_file_dict, should_overwrite, show_status, show_progress, threaded):
  global status_bar

  total_files = 0
  for dir, namespace_files in namespace_file_dict.items():
    for namespace, files in namespace_files.items():
      total_files += len(files)

  status_bar = BclCleanerStatus(total_files, show_status, show_progress)

  if threaded:
    #spawn a pool of threads, and pass them queue instance
    ncpus = multiprocessing.cpu_count()
    queue = Queue.Queue()
    jobs = [ThreadCleanFile(queue, i) for i in range(ncpus)]
    for job in jobs:
        job.setDaemon(True)
        job.start()

    for dir, namespace_files in namespace_file_dict.items():
      for namespace, files in namespace_files.items():
        if len(files):
          queue.put([namespace, should_overwrite, len(files)])
    queue.join()
  else:
    for dir, namespace_files in namespace_file_dict.items():
      for namespace, files in namespace_files.items():
        for filename in files:
          CleanFile(namespace, filename, should_overwrite)
  print ""

def is_source_file(x):
  return x.endswith('.cpp') or x.endswith('.fwd.hh') or x.endswith('.h') or x.endswith('.hpp')

def main():

  global script_name, copyright_block
  script_name = os.path.abspath(sys.argv[0])
  copyright_block_file = os.path.abspath(os.path.dirname(script_name) + "/../../documentation/bcl_copyright.txt")
  if not os.path.exists(copyright_block_file):
    print "Cannot locate bcl copyright file! Exiting"
    sys.exit(-1)
  else:
    fl = open(copyright_block_file,'r')
    copyright_block = [x.strip() + '\n' for x in fl.readlines() if len(x)]
    if copyright_block[-1] != '\n':
      copyright_block.append('\n') # append an extra new line to separate the copyright block from other code
    fl.close()
  should_overwrite = False
  show_stats = False
  show_status = False
  show_progress = True
  args = sys.argv[1::1]
  subdir = ""
  next_is_subdir = False
  threaded = True
  stats_file = ""
  last_was_stats = False
  if len(args) == 0 or lower(args[0].lstrip('-')) in set(['h', 'help']):
    usage()
  for arg in args[1:]:
    if arg == "h" or arg == "help" or arg == '-h' or arg == '--help':
      usage()
    elif arg == "-o" or arg == "--overwrite":
      should_overwrite = True
    elif arg == "-s" or arg == "--stats":
      show_stats = True
      last_was_stats = True
    elif arg == "--show-status":
      show_status = True
    elif arg == "--directory":
      next_is_subdir = True
    elif arg == "--nothreads":
      threaded = False
    elif arg == "--noprogress":
      show_progress = False
    elif next_is_subdir:
      next_is_subdir = False
      subdir = arg
      if not subdir.endswith(os.sep):
        subdir += os.sep
    elif last_was_stats:
      stats_file = arg
      last_was_stats = False
    else:
      print "Unrecognized flag: " + arg
      usage()

  # save the original working directory
  original_directory = os.getcwd()

  os.chdir(args[0])


  global status_bar
  directories_to_files = {}
  number_errors = 0
  if subdir == "":
    directories_to_files["./source"] = getFilesWithSufficesInDirectoryDictionary("./source", [ ".cpp" ])
    directories_to_files["./apps"] = getFilesWithSufficesInDirectoryDictionary("./apps", [ ".h", ".hpp", ".fwd.hh", ".cpp" ])
    directories_to_files["./example"] = getFilesWithSufficesInDirectoryDictionary("./example", [ ".h", ".hpp", ".fwd.hh", ".cpp" ])
    directories_to_files["./include"] = getFilesWithSufficesInDirectoryDictionary("./include", [ ".h", ".hpp", ".fwd.hh"])
    CleanFiles(directories_to_files, should_overwrite, show_status, show_progress, threaded)
  else:
    files = [ str(file) for file in os.listdir(subdir) if os.path.isfile(subdir + file) and is_source_file(file)]
    status_bar = BclCleanerStatus(len(files), show_status, show_progress)
    for filename in files:
      try:
        CleanFile(subdir, filename, should_overwrite)
      except:
        print "Could not clean file: ",filename," check for bad c++ syntax or >> used to close templates, or long c-style comments /* spread over multiple lines */"

  # return to the original directory
  os.chdir(original_directory)
  if stats_file == "":
    print "BclCleaner finished, " + str(CODESTATS.number_errors) + " errors found"
    print "Lines of code: " + str(CODESTATS.lines_of_code)

    if show_stats:
      print "Statistics:"
      total_binunary = max(CODESTATS.binunary_resolved + CODESTATS.binunary_unresolved, 1)
      total_gtlt = max(CODESTATS.gtlt_resolved + CODESTATS.gtlt_unresolved, 1)
      print "*/& resolved directly: " + str(CODESTATS.binunary_resolved) + " = " + strWithPrecision(100.0 * CODESTATS.binunary_resolved / float(total_binunary), 2) + "%"
      print "*/& unresolved: " + str(CODESTATS.binunary_unresolved) + " = " + strWithPrecision(100.0 * CODESTATS.binunary_unresolved / float(total_binunary), 2) + "%"
      print "<> resolved: " + str(CODESTATS.gtlt_resolved) + " = " + strWithPrecision(100.0 * CODESTATS.gtlt_resolved / float(total_gtlt), 2) + "%"
      print "<> unresolved: " + str(CODESTATS.gtlt_unresolved) + " = " + strWithPrecision(100.0 * CODESTATS.gtlt_unresolved / float(total_gtlt), 2) + "%"
  else:
    fil = open(stats_file, 'w')
    fil.write(\
      str(CODESTATS.number_errors) + " " + str(CODESTATS.lines_of_code) + " "\
      + str(CODESTATS.binunary_resolved) + " " + str(CODESTATS.binunary_unresolved) + " "\
      + str(CODESTATS.gtlt_resolved) + " " + str(CODESTATS.gtlt_unresolved))
    fil.close()

if __name__ == '__main__':
  main()
