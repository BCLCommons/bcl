'''
Created on Apr 14, 2010
@brief Functions used in parsing code-files, some of which are bcl-specific
@author: mendenjl
'''

import sys
import os
import cStringIO
from curses.ascii import isspace, isalpha, isdigit
from string import lower

# write a string with a certain precision
def strWithPrecision(number, precision):
  number_str = str(number)
  decimal_point = number_str.find('.')
  if decimal_point < 0 or decimal_point + precision + 1 >= len(number_str):
    return number_str
  return number_str[:decimal_point + precision + 1]

def Contains(listy, stringy):
  i = 0
  for strins in listy:
    if strins == stringy and stringy == strins:
      i = 1
  return i

def findIndexStringStartingWith(listy, stringy, startPos = 0):
  for i in xrange(max(startPos, 0), len(listy) , 1):
    if listy[i].startswith(stringy):
      return i
  return (-1)

def rfindIndexStringStartingWith(listy, stringy, startPos = sys.maxint):
  for i in xrange(min(len(listy), startPos), 0, -1):
    if listy[i].startswith(stringy):
      return i
  return (-1)

def maxLeng(listy):
  i = 0
  for a in listy:
    if len(str(a)) > i:
      i = len(str(a))
  return i

def findIndex(listy, stringy):
  i = -1
  j = 0
  for strins in listy:
    if strins == stringy and stringy == strins:
      i = j
      break
    j += 1
  return i

def count(listy, stringy):
  i = 0
  for strins in listy:
    if strins == stringy and stringy == strins:
      i += 1
  return i

def StartOfFenceScope(lines, line_num, pos, startScope, endScope):
  if lines[line_num][pos] != endScope:
    return line_num, pos
  else:
    depth = 1
    new_pos = pos - 1
    new_line = line_num
    while new_pos >= 0:
      if lines[line_num][new_pos] == startScope:
        depth += 1
      elif lines[line_num][new_pos] == endScope:
        depth -= 1
        if depth == 0:
          break
      new_pos -= 1
    if depth > 0:
      while new_line >= 0:
        new_pos = 0
        while new_pos >= 0:
          if lines[new_line][new_pos] == startScope:
            depth += 1
          elif lines[new_line][new_pos] == endScope:
            depth -= 1
            if depth == 0:
              break
          new_pos -= 1
        if depth == 0:
          break
        new_line -= 1
    return new_line, new_pos

def EndOfFenceScope(lines, line_num, startScope, endScope, pos):
  if lines[line_num][pos] != startScope:
    return line_num, pos
  else:
    depth = 1
    new_pos = pos + 1
    new_line = line_num
    while new_pos < len(lines[line_num]):
      if lines[line_num][new_pos] == startScope:
        depth += 1
      elif lines[line_num][new_pos] == endScope:
        depth -= 1
        if depth == 0:
          break
      new_pos += 1
    if depth > 0:
      new_line += 1
      while new_line < len(lines):
        new_pos = 0
        while new_pos < len(lines[new_line]):
          if lines[new_line][new_pos] == startScope:
            depth += 1
          elif lines[new_line][new_pos] == endScope:
            depth -= 1
            if depth == 0:
              break
          new_pos += 1
        if depth == 0:
          break
        new_line += 1
    return new_line, new_pos

# put start iter at the first char, endIter at last char, or endIter at -1 if no variable name
def nextVarNameStart(originalString, startIter):
  #print "var: " + str(startIter) + " " + str(len(originalString))
  while startIter < len(originalString) and (isValidAlpha(originalString[startIter]) == 0 or originalString[startIter].isdigit()):
    startIter += 1
  if startIter >= len(originalString):
    startIter = -1
  elif nextVarNameEnd(originalString, startIter) == -1:
    while startIter < len(originalString) and isValidAlpha(originalString[startIter]):
      startIter += 1
    startIter = nextVarNameStart(originalString, startIter)
  return startIter

  # put start iter at the first char, endIter at last char, or endIter at -1 if no variable name
def nextNonArrayVarNameStart(originalString, startIter):
  #print "var: " + str(startIter) + " " + str(len(originalString))
  while startIter < len(originalString) and (isValidAlpha(originalString[startIter]) == 0 or originalString[startIter].isdigit()):
    startIter += 1
  if startIter >= len(originalString):
    startIter = -1
  elif nextNonArrayVarNameEnd(originalString, startIter) == -1:
    while startIter < len(originalString) and isValidAlpha(originalString[startIter]):
      startIter += 1
    startIter = nextNonArrayVarNameStart(originalString, startIter)
  return startIter

def nextFuncNameStart(originalString, startIter):
  while startIter < len(originalString) and (isValidAlpha(originalString[startIter]) == 0 or originalString[startIter].isdigit()):
    startIter += 1
  if startIter >= len(originalString):
    startIter = -1
  elif nextFuncNameEnd(originalString, startIter) == -1:
    while startIter < len(originalString) and isValidAlpha(originalString[startIter]):
      startIter += 1
    return nextFuncNameStart(originalString, startIter)
  return startIter

def nextNonArrayVarNameEnd(originalString, startItera):
  startIter = startItera
  endIter = -1
  if startIter < len(originalString):
    endIter = startIter
    while endIter < len(originalString) and isValidAlpha(originalString[endIter]):
      endIter += 1
    basicEndIter = endIter
    while endIter < len(originalString) and originalString[endIter] == ' ':
      endIter += 1
    if endIter < len(originalString) and (originalString[endIter] == '(' or originalString[endIter] == '['):
      endIter = -1
    else:
      endIter = basicEndIter
  else:
    endIter = -1
  return endIter

def nextNameEnd(originalString, startItera):
  startIter = startItera
  endIter = -1
  #print "func: " + str(startIter)
  if startIter < len(originalString):
    endIter = startIter
    while endIter < len(originalString) and isValidAlpha(originalString[endIter]):
      endIter += 1
  return endIter

def nextNameStart(originalString, startItera):
  startIter = startItera
  endIter = -1
  #print "func: " + str(startIter)
  if startIter >= 0 and isValidAlpha(originalString[startIter]):
    endIter = startIter
    while endIter >= 0 and isValidAlpha(originalString[endIter]):
      endIter -= 1
    endIter += 1
  return endIter

def spaceEnd(originalString, i):
  j = i
  while j < len(originalString) and isspace(originalString[j]):
    j += 1
  return j

def nextVarNameEnd(originalString, startItera):
  basicEndIter = nextNameEnd(originalString, startItera)
  endIter = -1
  if basicEndIter != -1:
    nextChPos = spaceEnd(originalString, basicEndIter)
    if nextChPos == len(originalString) or originalString[nextChPos] != '(':
      endIter = basicEndIter
  return endIter

def nextFuncNameEnd(originalString, startItera):
  basicEndIter = nextNameEnd(originalString, startItera)
  endIter = -1
  if basicEndIter != -1:
    nextChPos = spaceEnd(originalString, basicEndIter)
    if nextChPos < len(originalString) and originalString[nextChPos] == '(':
      endIter = basicEndIter
  return endIter

def findReplaceAllVariableNames(originalString, listOfVarsToFind, listOfVarsToReplace):
  startIter = 0
  endIter = 0
  startIter = nextVarNameStart(originalString, startIter)

  while startIter > -1:
    endIter = nextVarNameEnd(originalString, startIter)

    localPos = findIndex(listOfVarsToFind, originalString[startIter:endIter])
    if localPos > -1:
      originalString = originalString[0:startIter:1] + '(' + str(listOfVarsToReplace[localPos]) + ')' + originalString[endIter::1]
    else:
      startIter = endIter + 1
    startIter = nextVarNameStart(originalString, startIter)
  return originalString

def HasQuotes(line):
  return line.find("'") >= 0 or line.find('"') >= 0

def OverWriteQuotes(originalString, fill_ch):
  line = ""
  i = 0
  ln = len(originalString)
  while i < ln:
    ch = originalString[i]
    line += ch
    if ch == "'" or ch == '"':
      i += 1
      if i == ln:
        return line
      while i < ln - 1 and originalString[i] != ch:
        line += fill_ch
        if originalString[i] == '\\':
          line += fill_ch
          i += 1
        i += 1
      if i < ln:
        if originalString[i] == ch:
          line += ch
        else:
          line += fill_ch
    i += 1
  return line

def getAllArguments(originalString, startPosa):
  startPos = startPosa
  argList = []
  while originalString[startPos] != '(':
    startPos += 1
  parenDepth = 1
  #print "Original: " + originalString
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
  #print "now arglist: " + str(argList)
  return argList

def functionArgumentsLastParen(originalString, startPosa):
  startPos = startPosa
  while originalString[startPos] != '(':
    startPos += 1
  parenDepth = 1

  startPos += 1
  while parenDepth > 0:
    if originalString[startPos] == '(':
      parenDepth += 1
    elif originalString[startPos] == ')':
      parenDepth -= 1
    startPos += 1
  return startPos - 1

def isMathOp(aCh):
  i = (aCh == '+' or aCh == '-' or aCh == '/' or aCh == '%' or aCh == '*' or aCh == '^' or aCh == '&' or aCh == '|' or aCh == '>' or aCh == '<' or aCh == '!' or aCh == '=')
  return i

def isExtendedMathOp(aCh):
  i = (aCh == '+' or aCh == '-' or aCh == '/' or aCh == '%' or aCh == '*' or aCh == '^' or aCh == '&' or aCh == '|' or aCh == '>' or aCh == '<' or aCh == '!' or aCh == '=' or aCh == '(' or aCh == ')')
  return i

def endOfThisInteger(strn, startPos):
  i = startPos
  if strn[i] == '-' or strn[i] == '+':
    i += 1
  if i < len(strn) and strn[i].isdigit():
    i += 1
    while i < len(strn) and strn[i].isdigit():
      i += 1
    return i
  elif i < len(strn):
    return -1
  else:
    return i

def startOfNextInteger(strn, startPos):
  i = startPos
  last = '+'
  if startPos > 0 and startPos < len(strn):
    last = strn[startPos - 1]
  goodHere = 0
  while goodHere != 1 and i < len(strn):
    # if the last thing we saw was anything other than a math operator
    # and the current thing is anything other than a sign or digit, then we aren't at the right point
    while i < len(strn):
      if isValidAlpha(last) == 0 and (strn[i].isdigit() == 1 or ((strn[i] == '-' or strn[i] == '+') and i + 1 < len(strn) and strn[i + 1].isdigit())):
        break
      if strn[i] != ' ':
        last = strn[i]
      i += 1
    if i < len(strn):
      endPos = endOfThisInteger(strn, i)
      if endPos == len(strn) or (strn[endPos] != '.' and strn[endPos] != 'E' and strn[endPos] != 'e'):
        goodHere = 1
      else:
        i = walkToEndOfNumber(strn, i) + 1
        if i < len(strn):
          last = strn[i - 1]
  if goodHere == 1:
    return i
  else:
    return -1

# find all things that look like variables, also return a list of indices of the starting / ending characters in this version
def findAllNonArrayVariablesAndIndices(line):
  tokens = []
  listOIndices = []
  inToken = 0

  # eliminate all array variables, function calls, math operations, etc.
  ite = len(line) - 1
  while ite > -1:
    if line[ite] == '(' or line[ite] == '[':
      inToken = 0
      ite -= 1
      while ite > -1 and isValidAlpha(line[ite]):
        ite -= 1
      ite += 1
    elif isExtendedMathOp(line[ite]) or line[ite] == ']' or line[ite] == ',' or line[ite] == '?' or line[ite] == ':' or line[ite] == ' ':
      inToken = 0
    else:
      if inToken == 0:
        tokens.append("")
        listOIndices.append(ite)
      inToken = 1
      tokens[len(tokens) - 1] += line[ite]
      listOIndices[len(tokens) - 1] = ite
    ite -= 1

  i = 0
  while i < len(tokens):
    # Get rid of any tokens that were blank, a single space, or that started off with a digit or anything else that isn't a valid starting character
    tokens[i] = tokens[i].strip()
    if len(tokens[i]) == 0 or (len(tokens[i]) == 1 and tokens[i][0] == ' ') or not isValidAlphaStart(tokens[i][len(tokens[i]) - 1]):
      del tokens[i]
      del listOIndices[i]
    else:
      # reverse the string
      tokens[i] = tokens[i][::-1]
      i += 1
  return tokens, listOIndices

# find all things that look like variables  
def findAllNonArrayVariables(line):
  return findAllNonArrayVariablesAndIndices(line)

def deltaBracketDepth(aLine):
  diffInDepth = 0
  pos = 0
  while pos < len(aLine):
    if aLine[pos] == '{':
      diffInDepth += 1
      pos += 1
    elif aLine[pos] == '}':
      diffInDepth -= 1
      pos += 1
    else:
      pos += 1
  return diffInDepth

def findEndOfUncommentedCPPFunction(ourLines, startLine):
  bracketDepth = 0
  while startLine < len(ourLines):
    bracketDepth += ourLines[startLine].count('{') - ourLines[startLine].count('}')
    if bracketDepth <= 0:
      break
    startLine += 1
  if startLine == len(ourLines):
    print "Unbalanced brackets in " + '\n'.join(ourLines) + "terminating "
    sys.exit(1)
  return startLine

def isInteger(strn): # assumes strn is a spaceless number
  retVal = len(strn) > 0 and ((strn[0] == '-' and len(strn) > 1 and strn[1].isdigit()) or strn[0].isdigit())
  i = 1
  while i < len(strn) and retVal == 1:
    retVal = strn[i].isdigit()
    i += 1
  return retVal

def walkToEndOfNumber(strn, startPos):
  i = startPos
  if len(strn) > i:
    gotP = 0
    if strn[i].isdigit():
      i += 1
    elif len(strn) > i + 1:
      if (strn[i] == '-' or strn[i] == '+') and strn[i + 1] == '.' and len(strn) > i + 2 and strn[i + 2].isdigit():
        i += 3
        gotP = 1
      elif strn[i] == '.' and strn[i + 1].isdigit():
        i += 2
        gotP = 1
      elif (strn[i] == '-' or strn[i] == '+') and strn[i + 1].isdigit():
        i += 2
    if i > startPos:
      if gotP == 0:
        while i < len(strn) and strn[i].isdigit():
          i += 1
        if i < len(strn) and strn[i] == '.':
          i += 1
      while i < len(strn) and strn[i].isdigit():
        i += 1
      if i < len(strn) and (strn[i] == 'e' or strn[i] == 'E'):
        last = i
        i += 1
        if i < len(strn) and (strn[i] == '+' or strn[i] == '-'):
          i += 1
        if i < len(strn) and strn[i].isdigit():
          i += 1
          while i < len(strn) and strn[i].isdigit():
            i += 1
        else:
          i = last
  return i

ishex = [ isdigit(chr(x)) or (chr(x) >= 'a' and chr(x) <= 'f') or (chr(x) >= 'A' and chr(x) <= 'F') for x in xrange(256)]

def walkToEndOfCPPNumber(strn, startPos):
  en = walkToEndOfNumber(strn, startPos)
  if en == startPos or en == len(strn) or not isalpha(strn[en]):
    return en
  low = lower(strn[en])

  # handle floating point suffices
  if low == 'f' or low == 'd':
    return en + 1

  # handle hexidecimal constants
  if low == 'x':
    en += 1
    while en < len(strn) and ishex[ord(strn[en])]:
      en += 1
    if en == len(strn) or not isalpha(strn[en]):
      return en
    low = lower(strn[en])

  # handle integral suffices
  if low == 'u':
    en += 1
    if en == len(strn):
      return en
    low = lower(strn[en])
  if low == 'l':
    en += 1
    if en < len(strn) and lower(strn[en]) == 'l':
      return en + 1
  return en

def walkToStartOfNumber(strn, startPos):
  i = startPos
  gotSign = 0
  gotP = 0
  if i >= 0 and (strn[startPos].isdigit() or (strn[startPos] == '.' and startPos != 0 and strn[startPos - 1].isdigit())):
    while i >= 0 and strn[i].isdigit():
      i -= 1
    if i >= 0 and strn[i] == '-' or strn[i] == '+':
      i -= 1
      gotSign = 1
    last = i
    if i >= 0 and (strn[i] == 'e' or strn[i] == 'E'):
      i -= 1
      if i >= 0 and strn[i] == '.':
        i -= 1
        gotP = 1
      if i >= 0 and strn[i].isdigit():
        i -= 1
        while i >= 0 and strn[i].isdigit():
          i -= 1
        if gotP == 0 and i >= 0 and strn[i] == '.':
          i -= 1
          while i >= 0 and strn[i].isdigit():
            i -= 1
        if i >= 0 and strn[i] == '-' or strn[i] == '+':
          i -= 1
      else:
        i = last
        gotP = 0
    elif i >= 0 and strn[i] == '.' and gotSign == 0:
      i -= 1
      while i >= 0 and strn[i].isdigit():
        i -= 1
      if i >= 0 and strn[i] == '-' or strn[i] == '+':
        i -= 1
    if i >= 0 and (strn[i + 1] == '-' or strn[i + 1] == '+') and not (isMathOp(strn[i]) or strn[i] == '('):
      i += 1

  return i

all_operators = set(['+', '-', '+=', '=', '&=', '^=', '|=', '/=', '<<=', '<<', '>>', '%=', '*=', '>>=', '-=', '&', '^', '|', '&&', '||', ',', '?', ':', '::', '/', '%', '*', '(', ')', '[', ']', '{', '}', '.', '->', '--', '++', '==', '>=', '<=', '<', '>', '!=', '!', '~', ';', '//', '/*', '*/', '->*', '.*', '::*'])
def walkToEndOfOperator(strn, start_pos):
  if start_pos >= len(strn):
    return start_pos
  i = start_pos + 1
  ch = strn[start_pos]
  while i < len(strn) and ch in all_operators:
    ch += strn[i]
    i += 1
  if ch not in all_operators:
    i -= 1
  return i

def isAFloat(strn):
  return walkToEndOfNumber(strn, 0) == len(strn)

def walkToEndOfNextScope(strn, startPos):# assumes no comments or spaces
  #print "walk to end of next scope: "+strn+"\n"+str(startPos)+"\n"
  endPos = startPos
  if strn[startPos] == '-':
    endPos += 1
  if strn[endPos].isdigit() == 1:
    endPos = walkToEndOfNumber(strn, endPos)
  elif(strn[endPos] == '('):
    # print "Started off at paren"
    parenDepth = 1
    endPos += 1
    while parenDepth > 0:
      if strn[endPos] == '(':
        parenDepth += 1
      elif strn[endPos] == ')':
        parenDepth -= 1
      endPos += 1

  #  print "finished at " + str(endPos)
  elif isValidAlphaStart(strn[endPos]):
    while len(strn) > endPos and isValidAlpha(strn[endPos]):
      endPos += 1
    if endPos < len(strn) and strn[endPos] == '(':
      parenDepth = 1
      endPos += 1
      while parenDepth > 0:
        if strn[endPos] == '(':
          parenDepth += 1
        elif strn[endPos] == ')':
          parenDepth -= 1
        endPos += 1
  else:
    print "no valid scope!: " + strn
    sys.exit(1)
  return endPos

def walkToStartOfScope(strn, endPos):
  startPos = endPos
  if strn[startPos] == ')': #this'll be easy
    parenDepth = -1
    startPos -= 1
    while parenDepth < 0:
      if strn[startPos] == '(':
        parenDepth += 1
      elif strn[startPos] == ')':
        parenDepth -= 1
      startPos -= 1
    while startPos >= 0 and isValidAlpha(strn[startPos]):
      startPos -= 1
    if startPos >= 0 and strn[startPos] == '-':
      if startPos >= 1:
        if isMathOp(strn[startPos - 1]) or strn[startPos - 1] == '(':
          startPos -= 1
      else:
        startPos -= 1
  else: # the more difficult case; first we'll try assuming it's a variable and see how far that gets us back.  Then try assuming it's a number.  Whichever gives us the lowest result is correct.  Note that if it is a number and the first letter is '-', we should not look for another negative sign
    while startPos >= 0 and isValidAlpha(strn[startPos]):
      startPos -= 1
    variablePos = startPos
    numberPos = endPos
    if strn[endPos].isdigit() and startPos > -1:
      numberPos = walkToStartOfNumber(strn, endPos)
    if variablePos < numberPos:
      startPos = variablePos
      if startPos >= 0 and strn[startPos] == '-':
        if startPos >= 1:
          if isMathOp(strn[startPos - 1]) or strn[startPos - 1] == '(':
            startPos -= 1
        else:
          startPos -= 1
    else:
      startPos = numberPos

  return startPos + 1

def rpartition(anstr, sep):
  pos = anstr.rfind(sep)
  if pos >= 0:
    return [anstr[0:pos:1], sep, anstr[(pos + len(sep))::1]]
  else:
    return ['', '', anstr]

def partition(anstr, sep):
  pos = anstr.find(sep)
  if pos >= 0:
    return [anstr[0:pos:1], sep, anstr[(pos + len(sep))::1]]
  else:
    return [anstr, '', '']

def isValidAlpha(character):
  val = 0
  if character.isalpha() and character != ' ':
    val = 1
  elif character.isdigit():
    val = 1
  elif character == '_':
    val = 1
  return val

def isValidAlphaStart(character):
  val = 0
  if character.isalpha() and character.isdigit() == 0 and character != ' ':
    val = 1
  elif character == '_':
    val = 1
  return val

def stripComments(ourLines):
  buffer = ""
  for strings in ourLines:
    if strings != None:
      strings = strings.strip()
      if(len(strings)):
        buffer += strings + '\n'
  walker = 0
  placer = cStringIO.StringIO()
  inEscape = 0
  if len(buffer):
    while walker + 1 < len(buffer):
      if buffer[walker] == '/' and buffer[walker + 1] == '/':
        #print "short comment"
        walker += 2
        while walker < len(buffer) and buffer[walker] != '\n':
          walker += 1
      elif buffer[walker] == '/' and buffer[walker + 1] == '*':
        #print "long comment"
        walker += 2
        while (walker + 1 < len(buffer) and (buffer[walker] != '*' or buffer[walker + 1] != '/')):
          walker += 1
        walker += 1
        if walker < len(buffer):
          walker += 1
      else:
        if inEscape == 0:
          if(buffer[walker] == '\\'):
            #   print "escape"
            if buffer[walker + 1] != '\n':
              inEscape = 1
              placer.write(buffer[walker])
              walker += 1
            else:
              walker += 2
          elif(buffer[walker] == '\''):
            #print "short quote"
            placer.write(buffer[walker]) # this will be the '
            walker += 1
            if walker < len(buffer):
              placer.write(buffer[walker]) # this will be the character
              if buffer[walker] == '\\':
                walker += 1
                placer.write(buffer[walker]) # this will be the optional post-escape escape character
                walker += 1
              elif buffer[walker] != '\'':
                walker += 1
              if walker < len(buffer) and buffer[walker] == '\'': # if there's not a \' then it must be a derivative
                placer.write(buffer[walker])
                walker += 1
          elif(buffer[walker] == '"'):
            #print "long quote"
            placer.write(buffer[walker])
            walker += 1
            while(walker < len(buffer) and (buffer[walker] != '"' or inEscape == 1)):
              inEscape = (inEscape == 0 and buffer[walker] == '\\')
              placer.write(buffer[walker])
              walker += 1
            if walker < len(buffer):
              placer.write(buffer[walker])
              walker += 1
          else:
            #print "standard"
            if buffer[walker] == '\t':
              placer.write(' ')
            else:
              placer.write(buffer[walker])
            walker += 1
        else:
          inEscape = 0
          if buffer[walker] == '\t':
            placer.write(' ')
          else:
            placer.write(buffer[walker])
          walker += 1

  if walker < len(buffer):
    placer.write(buffer[walker])
  buffer = placer.getvalue()
  ourLines = buffer.split('\n')
  ourLines = [ourLines[i].strip() for i in xrange(len(ourLines)) if ourLines[i] != None]
  ourLines = [ourLines[i] for i in xrange(len(ourLines)) if len(ourLines[i]) > 0]
  return(ourLines)

def stripCommentsAndQuotes(ourLines):
  buffer = ""
  for strings in ourLines:
    if strings != None:
      strings = strings.strip()
      if(len(strings)):
        buffer += strings + '\n'
  walker = 0
  placer = cStringIO.StringIO()
  inEscape = 0
  if len(buffer):
    while walker + 1 < len(buffer):
      if buffer[walker] == '/' and buffer[walker + 1] == '/':
        #print "short comment"
        walker += 2
        while walker < len(buffer) and buffer[walker] != '\n':
          walker += 1
      elif buffer[walker] == '/' and buffer[walker + 1] == '*':
        #print "long comment"
        walker += 2
        while (walker + 1 < len(buffer) and (buffer[walker] != '*' or buffer[walker + 1] != '/')):
          walker += 1
        walker += 1
        if walker < len(buffer):
          walker += 1
      else:
        if inEscape == 0:
          if(buffer[walker] == '\\'):
            #   print "escape"
            if buffer[walker + 1] != '\n':
              inEscape = 1
              placer.write(buffer[walker])
              walker += 1
            else:
              walker += 2
          elif(buffer[walker] == '\''):
            #print "short quote"
            placer.write(buffer[walker]) # this will be the '
            walker += 1
            if walker < len(buffer):
              if buffer[walker] == '\\':
                walker += 2
              elif buffer[walker] != '\'':
                walker += 1
              if walker < len(buffer) and buffer[walker] == '\'': # if there's not a \' then it must be a derivative
                placer.write(buffer[walker])
                walker += 1
          elif(buffer[walker] == '"'):
            #print "long quote"
            placer.write(buffer[walker])
            walker += 1
            while(walker < len(buffer) and (buffer[walker] != '"' or inEscape == 1)):
              inEscape = (inEscape == 0 and buffer[walker] == '\\')
              walker += 1
            if walker < len(buffer):
              placer.write(buffer[walker])
              walker += 1
          else:
            #print "standard"
            if buffer[walker] == '\t':
              placer.write(' ')
            else:
              placer.write(buffer[walker])
            walker += 1
        else:
          inEscape = 0
          if buffer[walker] == '\t':
            placer.write(' ')
          else:
            placer.write(buffer[walker])
          walker += 1

  if walker < len(buffer):
    placer.write(buffer[walker])
  buffer = placer.getvalue()
  ourLines = buffer.split('\n')
  ourLines = [ourLines[i].strip() for i in xrange(len(ourLines)) if ourLines[i] != None]
  ourLines = [ourLines[i] for i in xrange(len(ourLines)) if len(ourLines[i]) > 0]
  return(ourLines)


def stripMultilineComments(ourLines):
  buffer = ""
  for strings in ourLines:
    if strings != None:
      strings = strings.strip()
      if(len(strings)):
        buffer += strings + '\n'
  walker = 0
  placer = cStringIO.StringIO()
  inEscape = 0
  if len(buffer):
    while walker + 1 < len(buffer):
      if buffer[walker] == '/' and buffer[walker + 1] == '/':
        #print "short comment"
        placer.write("//")
        walker += 2
        while walker < len(buffer) and buffer[walker] != '\n':
          placer.write(buffer[walker])
          walker += 1
      elif buffer[walker] == '/' and buffer[walker + 1] == '*':
        #print "long comment"
        walker += 2
        while (walker + 1 < len(buffer) and (buffer[walker] != '*' or buffer[walker + 1] != '/')):
          walker += 1
        walker += 1
        if walker < len(buffer):
          walker += 1
      else:
        if inEscape == 0:
          if(buffer[walker] == '\\'):
            #   print "escape"
            if buffer[walker + 1] != '\n':
              inEscape = 1
              placer.write(buffer[walker])
              walker += 1
            else:
              walker += 2
          elif(buffer[walker] == '\''):
            #print "short quote"
            placer.write(buffer[walker]) # this will be the '
            walker += 1
            if walker < len(buffer):
              placer.write(buffer[walker]) # this will be the character
              if buffer[walker] == '\\':
                walker += 1
                placer.write(buffer[walker]) # this will be the optional post-escape escape character
                walker += 1
              elif buffer[walker] != '\'':
                walker += 1
              if walker < len(buffer) and buffer[walker] == '\'': # if there's not a \' then it must be a derivative
                placer.write(buffer[walker])
                walker += 1
          elif(buffer[walker] == '"'):
            #print "long quote"
            placer.write(buffer[walker])
            walker += 1
            while(walker < len(buffer) and (buffer[walker] != '"' or inEscape == 1)):
              inEscape = (inEscape == 0 and buffer[walker] == '\\')
              placer.write(buffer[walker])
              walker += 1
            if walker < len(buffer):
              placer.write(buffer[walker])
              walker += 1
          else:
            #print "standard"
            if buffer[walker] == '\t':
              placer.write(' ')
            else:
              placer.write(buffer[walker])
            walker += 1
        else:
          inEscape = 0
          if buffer[walker] == '\t':
            placer.write(' ')
          else:
            placer.write(buffer[walker])
          walker += 1

  if walker < len(buffer):
    placer.write(buffer[walker])
  buffer = placer.getvalue()
  ourLines = buffer.split('\n')
  ourLines = [ourLines[i].strip() for i in xrange(len(ourLines)) if ourLines[i] != None]
  ourLines = [ourLines[i] for i in xrange(len(ourLines)) if len(ourLines[i]) > 0]
  return(ourLines)

def isCommentBlockDelimitingLine(line):
  for character in line:
    if character != '/' and character != ' ':
      return False
  return True

def getDoxyCommentBlocks(ourLines):
  buffer = ""
  for strings in ourLines:
    if strings != None:
      strings = strings.strip()
      if len(strings) and '/' in strings:
        buffer += strings + '\n'
      else:
        strings += '\n'
  if len(buffer) == 0:
    return []
  if buffer[-1] != '\n':
    buffer += '\n'
  walker = 0
  ourBlocks = []
  while walker + 2 < len(buffer):
    if buffer[walker] == '/' and buffer[walker + 1] == '/' and buffer[walker + 2] == '!':
      blocks_with_tags = {}
      while walker + 2 < len(buffer) and buffer[walker] == '/' and buffer[walker + 1] == '/' and buffer[walker + 2] == '!':
        #print "short comment"
        walker += 3
        next_new_line = buffer.find('\n', walker)
        next_line = buffer[ walker:next_new_line].lstrip()
        walker = next_new_line + 1
        if len(next_line) > 2 and next_line.startswith('@'):
          next_line = next_line[1:].strip()
          tag_and_rest_of_line = partition(next_line, ' ')
          if len(tag_and_rest_of_line) == 3:
            if tag_and_rest_of_line[0] in blocks_with_tags:
              blocks_with_tags[ tag_and_rest_of_line[0]] += '\n' + tag_and_rest_of_line[2]
            else:
              blocks_with_tags[ tag_and_rest_of_line[0]] = tag_and_rest_of_line[2]
      if len(blocks_with_tags) > 0:
        ourBlocks.append(blocks_with_tags)
    elif buffer[walker] == '/' and buffer[walker + 1] == '/':
      walker = buffer.find('\n', walker + 1) + 1
    elif buffer[walker] == '/' and buffer[walker + 1] == '*':
      #print "long comment"
      walker += 2
      while (walker + 1 < len(buffer) and (buffer[walker] != '*' or buffer[walker + 1] != '/')):
        walker += 1
      walker += 1
      if walker < len(buffer):
        walker += 1
    else:
      old_char = buffer[walker]
      walker += 1
      if(old_char == '\\'):
        #   print "escape"
        walker += 1
      elif(old_char == '\''):
        #print "short quote"
        if buffer[walker] == '\\':
          walker += 2
        elif buffer[walker] != '\'':
          walker += 1
        if walker < len(buffer) and buffer[walker] == '\'': # if there's not a \' then it must be a derivative
          walker += 1
      elif(old_char == '"'):
        #print "long quote"
        while(walker < len(buffer) and buffer[walker] != '"'):
          if buffer[walker] == '\\':
            walker += 1
          walker += 1
        if walker < len(buffer):
          walker += 1

  return ourBlocks

def getDoxyBlocksContainingTag(ourLines, tag):
  new_blocks = []
  for block in getDoxyCommentBlocks(ourLines):
    if tag in block:
      new_blocks.append(block)

  return new_blocks

def getClassesAndStructs(ourLines):
  lines = stripComments(ourLines)
  lines = [ line.strip() for line in lines if (line.strip().startswith('struct ') or line.strip().startswith('class ')) and not line.strip().endswith(';')]
  for i in xrange(len(lines)):
    if lines[i].startswith('class '):
      lines[i] = lines[i][6:]
    else:
      lines[i] = lines[i][7:]
    if lines[i].startswith('BCL_API '):
      lines[i] = lines[i][8:]
    if lines[i].find(' ') > 0:
      lines[i] = lines[i][:lines[i].find(' ')]
    if lines[i].find('<') > 0:
      lines[i] = lines[i][:lines[i].find('<')]
    if lines[i].find(':') > 0:
      lines[i] = lines[i][:lines[i].find(':')]
  return lines

def extractStrings(ourLines, prep):
  ourStrings = []
  for i in xrange(len(ourLines)):
    if '\'' in ourLines[i] or '"' in ourLines[i]:
      walker = 0
      lastStrEnd = 0
      strParts = []
      while walker < len(ourLines[i]):
        if ourLines[i][walker] == '\'' or ourLines[i][walker] == '"':
          term = ourLines[i][walker]
          strParts.append(ourLines[i][lastStrEnd:walker:1] + prep + str(len(ourStrings) + 1))
          strStart = walker
          walker += 1
          while walker < len(ourLines[i]) and ourLines[i][walker] != term:
            if ourLines[i][walker] == '\\':
              walker += 1
            walker += 1
          walker = min(walker + 1, len(ourLines[i]))
          lastStrEnd = walker
          ourStrings.append(ourLines[i][strStart:lastStrEnd:1])
        else:
          walker += 1
      ourLines[i] = ''.join(strParts) + ourLines[i][lastStrEnd::1]
  return ourLines, ourStrings

def getAdditionalHeadersFromLines(lines, moreHeaders):
  for line in lines:
    if line.find("#include") >= 0:
      line = line.replace("#include", "").strip()
      line = line.replace('"', "").strip()
      if (line.endswith(".h") or line.endswith(".fwd.h")) and moreHeaders.find("include" + os.sep + line + " ") < 0 and line.find('/') < 0 and os.path.exists("include" + os.sep + line):
        moreHeaders += "include" + os.sep + line + " "
        ifile = open("include" + os.sep + line, 'r')
        morelines = ifile.readlines()
        ifile.close()
        moreHeaders = getAdditionalHeadersFromLines(morelines, moreHeaders)
  return moreHeaders

def getFilesInDirectoryDictionary(directory, suffix):
  all_files = {}
  for root, dirs, files in os.walk(directory):
    all_files[ root] = [ str(file) for file in files if os.path.isfile(root + os.sep + file) and file.endswith(suffix)]
    if '.svn' in dirs:
        dirs.remove('.svn')  # don't visit svn directories
  return all_files

def getFileInDirectoryDictionary(directory, filename):
  all_files = {}
  for root, dirs, files in os.walk(directory):
    all_files[ root] = [ str(file) for file in files if os.path.isfile(root + os.sep + file) and file == filename]
    if '.svn' in dirs:
        dirs.remove('.svn')  # don't visit svn directories
  return all_files

def getFilesWithSufficesInDirectoryDictionary(directory, suffices):
  all_files = {}
  for root, dirs, files in os.walk(directory):
    all_files[ root] = []
    for suffix in suffices:
      new_files = [ str(file) for file in files if os.path.isfile(root + os.sep + file) and file.endswith(suffix)]
      if len(new_files) > 0:
        all_files[ root].extend(new_files)
    if '.svn' in dirs:
        dirs.remove('.svn')  # don't visit svn directories
  return all_files

def getFilesFromDirectory(directory, suffix):
  all_files = []
  for root, dirs, files in os.walk(directory):
    all_files += [ str(root + os.sep + file) for file in files if (os.path.isfile(root + os.sep + file) and file.endswith(suffix))]
    if '.svn' in dirs:
        dirs.remove('.svn')  # don't visit svn directories
  return all_files

def writeTwoLists(afile, prefix, vars, joiner, RHSs, suffix):
  if len(vars):
    endStr = suffix + prefix
    afile.write(prefix + endStr.join([(str(vars[i]) + joiner + str(RHSs[i])) for i in xrange(len(vars))]) + suffix)

def writeOneList(afile, prefix, vars, suffix):
  if len(vars):
    endStr = suffix + prefix
    afile.write(prefix + endStr.join([str(vars[i]) for i in xrange(len(vars))]) + suffix)

def writeOneSet(afile, prefix, vars, suffix):
  if len(vars):
    endStr = suffix + prefix
    afile.write(prefix + endStr.join([str(var) for var in vars]) + suffix)
