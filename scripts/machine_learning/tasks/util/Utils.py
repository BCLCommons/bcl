'''
File system utilities
Created April 2013
@author: mendenjl
'''

import os, sys, time, commands, shutil, os.path
from time import gmtime, strftime, sleep

# python equivalent of mkdir -p $dr
def mkdirmp( dr):
  if not os.path.exists( dr):
    current_path_component = ''
    for p in dr.split(os.sep):
      current_path_component += p + os.sep
      if p != '.' and not os.path.exists(current_path_component):
        os.mkdir(current_path_component)
        
# python equivalent of rm -rf
def rm_rf(path):
  # have to handle files that might spring up while deleting the contents of the folder, so try to do it twice; if
  # it still fails, use the system-level fallback
  try:
    if os.path.isdir(path):
      shutil.rmtree(path)
    elif os.path.exists(path):
      os.remove(path)
  except:
    try:
      if os.path.isdir(path):
        shutil.rmtree(path)
      elif os.path.exists(path):
        os.remove(path)
    except:
      os.system('rm -rf ' + path)

# try to execute a command up to N times
def tryExecute(com, max_retries, msg, seconds_between_tries):
  retries = 0
  (status, output) = commands.getstatusoutput(com)
  if status == 0:
    return (True, output)
  elif status != -1 and status != 255:
    return (False, output)
  err_out = output
  while retries < max_retries:
    print msg + ', retrying up to ' + str(max_retries - retries) + ', last error was: ' + err_out
    retries += 1
    sleep(seconds_between_tries)
    (status, output) = commands.getstatusoutput(com)
    if status != 0:
      err_out = output
    else:
      return (True, output)
  return (False, err_out)
