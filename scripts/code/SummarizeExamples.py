#!/usr/bin/python2
'''
Created on Nov 5, 2010

Functions:
Analyzes examples, classes, and their connections
Takes statistics on # of example checks, # of publicly accessible functions, by namespace, developer, and by file

@author: mendenjl
'''

import os
import sys
import os.path
import datetime
from operator import *
import subprocess
from symbol import exec_stmt
from CodeFileUtils import *
from BclExampleFile import *
from BclClassFile import BCLClassFile
from Ctag import *
from HtmlTableHelperFuncs import *

# directories in which app files live

def DeveloperStatsSqlStatements(lines):
  format_str = "ds_date,ds_name,ds_n_public_functions,ds_n_example_checks,ds_max_percent_checked,ds_n_classes,ds_n_examples_complete,ds_n_examples_incomplete"
  rows = [ ','.join([ item.strip(' ') for item in lines[i].split('|') if len(item) > 0]) for i in xrange(len(lines)) if len(lines[i]) > 0]
  now = datetime.datetime.now()
  date = str(now.year) + '-' + str(now.month) + '-' + str(now.day);
  sql_liner = cStringIO.StringIO()
  for row in rows:
    sql_liner.write('INSERT into developer_statistics(' + format_str + ") values(" + date + "," + row + ");\n")
  return sql_liner.getvalue()

def usage():
  print "\n" + \
        "usage: SummarizeExamples.py bcl-path [options]\n" + \
        "options:\n" + \
        "-h/--help print this dialogue\n" + \
        "-o/--output-directory print all results to html files in this directory\n" + \
        "-s/--screen print all results to screen (default is to print them to html files)\n" + \
        "-d/--developer-stats compute developer stats\n" + \
        "-e/--errors find errors in example/class comment linkages\n" + \
        "-f/--functions print information about each function\n" + \
        "-n/--namespaces print information about each namespace\n" + \
        "-x/--examples print information about each example\n" + \
        "--sql print sql statements to insert information into the database\n" + \
        "If none of -d, -e, -x, -f , -n are specified, all information will be printed"
  sys.exit(1)

def main():

  bcl_path = "." + os.sep
  usr_dir = os.getenv("HOME")
  html_path = usr_dir + os.sep + "tmp" + os.sep
  arguments = sys.argv[1::1]

  if len(arguments) == 0:
    usage()
  bcl_path = arguments[0]

  output_to_screen = False
  output_developer_stats = False
  output_errors = False
  output_functions = False
  output_namespace_info = False
  output_example_info = False
  output_sql_files = False
  last_was_output_dir = False

  for i in xrange(1, len(arguments)):
    if last_was_output_dir:
      pass
    elif arguments[i] == '-h' or arguments[i] == '--help':
      usage()
    elif arguments[i] == '-o' or arguments[i] == '--output-directory':
      if i < len(arguments) - 1:
        last_was_output_dir = True
        html_path = arguments[i + 1]
        continue
      else:
        print "-o/--output-directory requires an output directory!"
        sys.exit(1)
    elif arguments[i] == '-s' or arguments[i] == '--screen':
      output_to_screen = True
    elif arguments[i] == '-d' or arguments[i] == '--developer-stats':
      output_developer_stats = True
    elif arguments[i] == '-e' or arguments[i] == '--errors':
      output_errors = True
    elif arguments[i] == '-f' or arguments[i] == '--functions':
      output_functions = True
    elif arguments[i] == '-n' or arguments[i] == '--namespaces':
      output_namespace_info = True
    elif arguments[i] == '-x' or arguments[i] == '--examples':
      output_example_info = True
    elif arguments[i] == '--sql':
      output_sql_files = True;
    else:
      print arguments[i] + " is not a valid option"
      usage()
    last_was_output_dir = False

  exit_after_errors = False
  if not (output_developer_stats or output_errors or output_functions or output_namespace_info or output_example_info):
    output_developer_stats = output_errors = output_functions = output_namespace_info = output_example_info = True
  elif not (output_developer_stats or output_functions or output_namespace_info or output_example_info):
    exit_after_errors = True

  if not bcl_path.endswith('/'):
    bcl_path += '/'
  if not html_path.endswith('/'):
    html_path += '/'

  # get a list of all the example files
  examples_with_paths = getFilesFromDirectory(bcl_path + "example", "cpp")
  #print "reading examples"
  example_dictionary = {}
  for ex in examples_with_paths:
    exampl = BCLExampleFile()
    exampl.ExtractFromFile(ex)
    if len(exampl.short_filename.strip()):
      example_dictionary[ exampl.short_filename] = exampl


  # get a list of the headers and sources
  sources = getFilesFromDirectory(bcl_path + "source", ".cpp")
  sources.extend(getFilesFromDirectory(bcl_path + "apps", ".cpp"))

  headers = getFilesFromDirectory(bcl_path + "include", ".h")
  headers.extend(getFilesFromDirectory(bcl_path + "apps", ".h"))

  unique_class_id = {}
  for source in sources:
    unique_class_id[ source[source.rfind('/') + 1:-4]] = source
  for header in headers:
    unique_class_id[ header[header.rfind('/') + 1:-2]] = header

  #print "reading source files"
  class_dict = {}
  for class_name in unique_class_id.values():
    new_class = BCLClassFile()
    new_class.ExtractFromFile(class_name, example_dictionary)
    if len(new_class.short_filename):
      if new_class.namespace != 'bcl':
        class_dict[ new_class.namespace + os.sep + new_class.short_filename] = new_class
      else:
        class_dict[ new_class.short_filename] = new_class

  #print "finding errors"
  errors = []
  for class_type in class_dict.values():
    error_strs = class_type.GetErrorString()
    if len(error_strs):
      errors.extend([ "| " + class_type.short_filename + " | " + error + " |" for error in error_strs])
  for example in example_dictionary.values():
    error_strs = example.GetErrorString(class_dict)
    if len(error_strs):
      errors.extend([ "| " + example.short_filename + " | " + error + " |" for error in error_strs])

  if not os.path.exists(html_path) and not output_to_screen:
    os.mkdir(html_path)

  if len(errors) and output_errors:
    heading = str(len(errors)) + " errors found in examples and class comments, links, and statuses"
    header_row = "| Filename | Error | Detail |"
    if not output_to_screen:
      error_file = open(html_path + "ExampleAndCommentErrors.html", "w")
      error_file.write(GetHTMLHeader() + '\n')
      error_file.write(HtmlTableFromTwikiTable(heading, [], header_row, sorted(errors)))
      error_file.write(GetHTMLFooter() + '\n')
      error_file.close()
    else:
      print heading + '\n' + header_row + '\n' + '\n'.join(sorted(errors))
  if exit_after_errors:
    sys.exit(0)

  # cmd1 contains the command to get all the ctags from the include and source directories and print them to stdout
  cmd1 = "ctags -R --c-kinds=f+p --fields=+f-i+k-K-l+m-n+s+S-z+t+a  -f '-' " + bcl_path + "include " + bcl_path + "source " + bcl_path + "apps "
  for line in os.popen4(cmd1)[1].readlines():
    ctag = Ctag(line)

    # ignore header lines  and other lines that do not seem to declare a valid function
    if len(ctag.function_name) == 0:
      #print "non-function " + line.strip()
      continue

    # ignore pure virtual functions
    if ctag.implementation_type == 2:
      #print "pure virtual function " + line.strip()
      continue

    # ignore destructors
    if ctag.function_name[0] == '~':
      #print "destructor" + line.strip()
      continue

    # ignore protected and private functions
    if ctag.access_type > 1:
      #print "protected or private " + line.strip()
      continue

    # ignore class function that are not accompanied by an access scope (e.g. implementations outside the class defs)
    if ctag.holder_type == 0 and ctag.access_type == 0:
      #print "class function implementation outside class def " + line
      continue;

    # ignore implementations of non-class functions in cpp's
    if ctag.function_type == 'f' and ctag.holder_type == 1 and ctag.file_name.endswith(".cpp"):
      #print "function implementation outside header file " + line
      continue

    #ignore clone and get class identifier
    if ctag.holder_type == 0 and ctag.signature == "() const\n" and (ctag.function_name == "Clone" or ctag.function_name == "GetClassIdentifier" or ctag.function_name == "GetAlias"):
      continue

    #ignore sizeof operator as well as any function which, by name, are clearly variables or template parameters
    if ctag.function_name == "sizeof" or ctag.function_name.startswith("m_") \
       or ctag.function_name.startswith("g_") or ctag.function_name.startswith("s_") \
       or ctag.function_name.startswith("t_"):
      continue

    # get the file name, but consider .fwd.hh's to be the same as .h
    file = ctag.file_name
    if file.endswith('.fwd.hh'):
      file = file[:-len('.fwd.hh')] + ".h"

    # actual applications are in the apps directory rather than app. 
    # The intermediate path is not needed, so discard it.  From here on, all file paths follow the convention of
    # namespace/file
    if file.find('/bcl_app_') > 0:
      file = 'app/' + file.split('/')[-1]

    if file in class_dict.keys():
      #print "found " + file + " for line " + line.strip()
      class_dict[ file].ConsiderTag(ctag)
    #else:
      #print "couldn't find " + file + " when considering " + line.strip()

  for class_type in class_dict.values():
    class_type.DistributeFunctions()

  if output_functions:
    function_list = []
    num_functions = 0
    for class_type in class_dict.values():
      funcs = class_type.GetAllFunctions()
      if len(funcs):
        num_functions += len(funcs)
        function_list.extend([ class_type.short_filename + " | " + func for func in funcs ])

    if len(function_list):
      heading = str(num_functions) + " externally accessible functions founds in include and source directories"
      subheadings = []
      subheadings.append("Clone and GetClassIdentifier are ignored")
      header_row = "| Filename | Class | Function | Author(s) |"
      if not output_to_screen:
        func_file = open(html_path + "Functions.html", "w")
        func_file.write(GetHTMLHeader() + '\n')
        subheadings.append("Sorting may take >10 seconds due to the size of the table")
        func_file.write(HtmlTableFromTwikiTable(heading, subheadings, header_row, sorted(function_list)))
        func_file.write(GetHTMLFooter() + '\n')
        func_file.close()
      else:
        print heading + '\n' + '\n'.join(subheadings) + '\n' + header_row + '\n' + '\n'.join(sorted(function_list))

  developer_names = set([])
  namespace_names = set([])
  for class_type in class_dict.values():
    namespace_names.add(class_type.namespace)
    for class_comment in class_type.classes.values():
      if class_comment.example_unnecessary == 0:
        for author in class_comment.authors:
          developer_names.add(author)
  for example in example_dictionary.values():
    for author in example.authors:
      developer_names.add(author)

  developer_stats = {}
  for developer in developer_names:
    new_developer = BCLDeveloperStatistic()
    new_developer.name = developer
    for class_type in class_dict.values():
      new_developer.ConsiderFile(class_type)
    developer_stats[ developer] = new_developer

  namespace_stats = {}
  for namespace in namespace_names:
    new_namespace = BCLNamespaceStatistic()
    new_namespace.name = namespace
    for class_type in class_dict.values():
      new_namespace.ConsiderFile(class_type)
    namespace_stats[ namespace] = new_namespace

  for example in example_dictionary.values():
    for author in example.authors:
      developer_stats[ author].ConsiderExample(example)
    namespace_stats[ example.namespace].ConsiderExample(example)

  total_example_checks = 0
  total_functions = 0
  total_classes = 0
  total_examples_complete = 0
  total_examples_incomplete = 0

  if len(developer_stats):
    developer_stats_list = []
    for developer, stats in sorted(developer_stats.items()):
      #max_percent = 100.0 * float(stats.number_of_example_checks) / float(max(stats.number_of_example_checks, stats.number_of_functions))
      max_percent_distributed = 0.0
      total_example_checks += stats.number_of_example_checks_distributed
      total_functions += stats.number_of_functions_distributed
      total_classes += stats.number_of_classes_distributed
      total_examples_complete += stats.number_of_examples_complete_distributed
      total_examples_incomplete += stats.number_of_examples_incomplete_distributed
      if stats.number_of_example_checks_distributed > 0 or stats.number_of_functions_distributed > 0:
        max_percent_distributed = 100.0 * float(stats.number_of_example_checks_distributed) / float(max(stats.number_of_example_checks_distributed, stats.number_of_functions_distributed))
      if output_developer_stats:
        developer_stats_list.append('| ' + developer + ' | ' \
            + strWithPrecision(stats.number_of_functions_distributed, 1) + ' | ' \
            + strWithPrecision(stats.number_of_example_checks_distributed, 1) + ' | '\
            + strWithPrecision(max_percent_distributed, 1) + ' | ' \
            + strWithPrecision(stats.number_of_classes_distributed, 1) + ' | ' \
            + strWithPrecision(stats.number_of_examples_complete_distributed, 1) + ' | ' \
            + strWithPrecision(stats.number_of_examples_incomplete_distributed, 1))
    max_percent = 100.0 * float(total_example_checks) / float(max(1, max(total_example_checks, total_functions)))
    if output_developer_stats:
      developer_stats_list.append('| bcl-all | ' + str(int(total_functions)) + ' | ' \
            + str(int(total_example_checks)) + ' | ' + strWithPrecision(max_percent, 1) + ' | '\
            + str(int(total_classes)) + ' | ' + str(int(total_examples_complete)) + ' | '\
            + str(int(total_examples_incomplete)))
    if output_developer_stats:
      heading = "Developer Statistics"
      header_row = "| developer | ~# public functions written | # example checks | max possible % of public functions checked | # classes | # examples complete | # examples incomplete |"
      developer_stats_list = sorted(developer_stats_list)
      if not output_to_screen:
        developer_stats_file = open(html_path + "DeveloperStats.html", "w")
        developer_stats_file.write(GetHTMLHeader() + '\n')
        developer_stats_file.write(HtmlTableFromTwikiTable(heading, [], header_row, developer_stats_list))
        developer_stats_file.write(GetHTMLFooter() + '\n')
        developer_stats_file.close()
      else:
        print heading + '\n' + header_row + '\n' + '\n'.join(developer_stats_list)
      if output_sql_files:
        sql_file = open(html_path + "DeveloperStats.sql", "w")
        sql_file.write(DeveloperStatsSqlStatements(developer_stats_list))
        sql_file.close()

  if output_example_info:
    header_row = "| example | *status* | # example checks | # functions that should be tested | max possible % of public functions checked |"
    heading = "Example information"
    example_table_rows = []
    for example in example_dictionary.values():
      if example.number_of_checks != 0:
        max_percent = 100.0 * float(example.number_of_checks) / float(max(example.number_of_checks, example.number_of_functions))
        example_table_rows.append("| " + example.short_filename[8:-4] + " | " + example.status + " | " + str(example.number_of_checks) + " | " \
              + str(example.number_of_functions) + " | " + strWithPrecision(max_percent, 1) + " |")
    for class_type in class_dict.values():
      if class_type.NeedsAnExample() == 1:
        hypothetical_example_name = class_type.short_filename[4:].split('.')[0]

        status = "non-existent"
        if len(class_type.examples):
          status = "empty"
        example_table_rows.append("| " + hypothetical_example_name + " | " + status + " | 0 | " + str(class_type.TotalNumberFunctions()) + " | 0.0 |")

    if not output_to_screen:
      example_table_file = open(html_path + "ExamplesData.html", "w")
      example_table_file.write(GetHTMLHeader() + '\n')
      example_table_file.write(HtmlTableFromTwikiTable(heading, [], header_row, sorted(example_table_rows)))
      example_table_file.write(GetHTMLFooter() + '\n')
      example_table_file.close()
    else:
      print heading + '\n' + header_row + '\n' + '\n'.join(sorted(example_table_rows))

  if output_namespace_info:
    example_table_rows = []
    total_examples_complete = 0
    total_examples_incomplete = 0
    namespace_stats_list = []
    for namespace, stats in sorted(namespace_stats.items()):
      #max_percent = 100.0 * float(stats.number_of_example_checks) / float(max(stats.number_of_example_checks, stats.number_of_functions))
      max_percent = 0.0
      total_examples_complete += stats.number_of_examples_complete
      total_examples_incomplete += stats.number_of_examples_incomplete
      if stats.number_of_example_checks > 0 or stats.number_of_functions > 0:
        max_percent = 100.0 * float(stats.number_of_example_checks) / float(max(stats.number_of_example_checks, stats.number_of_functions, 1))
        # add offset if there are no classes to prevent division-by-zero errors
        if stats.number_of_classes == 0:
          stats.number_of_classes = 1
        namespace_stats_list.append('| ' + namespace + ' | ' \
            + str(stats.number_of_classes) + ' | '
            + str(stats.number_of_examples_complete) + ' | '
            + str(stats.number_of_examples_incomplete) + ' | '
            + strWithPrecision(100.0 * float(stats.number_of_examples_complete) / float(stats.number_of_classes), 1) + ' | '
            + strWithPrecision(100.0 * float(stats.number_of_examples_incomplete) / float(stats.number_of_classes), 1) + ' | '
            + strWithPrecision(100.0 - 100.0 * float(stats.number_of_examples_incomplete + stats.number_of_examples_complete) / float(stats.number_of_classes), 1) + ' | '
            + strWithPrecision(stats.number_of_example_checks, 1) + ' | '\
            + strWithPrecision(stats.number_of_functions, 1) + ' | ' \
            + strWithPrecision(max_percent, 1))
    max_percent = 100.0
    if total_functions > 0:
      max_percent = 100.0 * float(total_example_checks) / float(max(total_example_checks, total_functions))
    # add offset if there are no classes to prevent division-by-zero errors
    if total_classes == 0:
      total_classes = 1
    namespace_stats_list.append('| total | ' + str(int(total_classes)) + ' | ' \
            + str(int(total_examples_complete)) + ' | ' + str(int(total_examples_incomplete)) + ' | ' \
            + strWithPrecision(100.0 * float(total_examples_complete) / float(total_classes), 1) + ' | '
            + strWithPrecision(100.0 * float(total_examples_incomplete) / float(total_classes), 1) + ' | '
            + strWithPrecision(100.0 - 100 * float(total_examples_incomplete + total_examples_complete) / float(total_classes), 1) + ' | '
            + strWithPrecision(total_example_checks, 1) + ' | '\
            + strWithPrecision(total_functions, 1) + ' | ' \
            + strWithPrecision(max_percent, 1))
    header_row = "| namespace | #classes | #examples complete | #examples incomplete | %examples complete | %incomplete | %unwritten | # example checks | # functions that should be tested | max possible % of public functions checked |"
    heading = "Example data summarized by namespace"
    if not output_to_screen:
      example_namespace_file = open(html_path + "NamespaceExampleData.html", "w")
      example_namespace_file.write(GetHTMLHeader() + '\n')
      example_namespace_file.write(HtmlTableFromTwikiTable(heading, [], header_row, sorted(namespace_stats_list)))
      example_namespace_file.write(GetHTMLFooter() + '\n')
      example_namespace_file.close()
    else:
      print heading + '\n' + header_row + '\n' + '\n'.join(sorted(namespace_stats_list))

if __name__ == '__main__':
  main()
