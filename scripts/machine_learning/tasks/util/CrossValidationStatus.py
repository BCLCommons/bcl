'''
Status bar and cross validation information printer
Created April 2013
@author: mendenjl
'''
import os, sys, os.path
from time import gmtime, strftime, sleep
from Utils import *

# class that handles updating the status bar
class BclCVStatus:

  # enum for overall status
  ERRORS = -1
  FINISHED = 1
  RUNNING = 2
  NEVER_RAN = 3

  # enum for improvement direction
  SMALLER_IS_BETTER = -1
  LARGER_IS_BETTER = 1

  def __init__(self, show_status, show_progress, status_filename, files_remove_on_success, result_file, abort_command, expect_blind_dataset):
    self.show_status = show_status
    self.last_percent = 0
    self.last_line_length = 0
    self.show_progress = show_progress
    self.last_line = ''
    self.status_filename = status_filename
    self.jobs_in_queue = 0
    self.jobs_running = 0
    self.jobs_finished = 0
    self.jobs_error = 0
    self.jobs_pbs_error = 0
    self.jobs_merge_status = 0
    self.overall_status = BclCVStatus.RUNNING
    self.files_remove_on_success = files_remove_on_success
    self.result_file = result_file
    self.wasFinalized = False
    self.abort_command = abort_command
    self.final_result = 0
    self.improvement_type = BclCVStatus.SMALLER_IS_BETTER
    self.independent_ave_file = status_filename + '.raw_ind_ave.txt'
    self.independent_average = 0
    self.final_blind_result = 0
    self.final_consensus_result = 0
    self.expect_blind_dataset = expect_blind_dataset

  # just update the status variables without writing anything
  def updateStatus(self):
    if not os.path.exists(os.path.dirname(self.status_filename)):
      self.jobs_in_queue = 0
      self.jobs_running = 0
      self.jobs_finished = 1
      self.jobs_error = 0
      self.jobs_pbs_error = 0
      self.jobs_merge_status = 1
      self.overall_status = BclCVStatus.FINISHED
      if len(self.result_file) and not os.path.exists(self.result_file):
        self.overall_status = BclCVStatus.NEVER_RAN
      return self.overall_status
    # try opening the status file.  if this fails, it probably just means that another job is writing to it
    maxsleep=30
    while not os.path.exists(self.status_filename) and maxsleep > 0:
      maxsleep-=2
      sleep(2)
    statuses = []
    try:
      fil = open(self.status_filename, 'r')
      try:
        statuses = [ int(x) for x in fil.readline().strip().split(' ')]
      except:
        pass
      fil.close()
    except:
      pass
    if len(statuses) == 6:
      self.jobs_in_queue = statuses[0]
      self.jobs_running = statuses[1]
      self.jobs_finished = statuses[2] - statuses[3] - statuses[4]
      self.jobs_error = statuses[3]
      self.jobs_pbs_error = statuses[4]
      self.jobs_merge_status = statuses[5]
      self.overall_status = BclCVStatus.RUNNING
      if self.jobs_error or self.jobs_pbs_error:
        self.overall_status = BclCVStatus.ERRORS
      elif not self.jobs_in_queue and not self.jobs_running and self.jobs_merge_status == BclCVStatus.FINISHED:
        self.overall_status = BclCVStatus.FINISHED

    if self.jobs_finished and os.path.exists(self.independent_ave_file):
      maxsleep=30
      while not os.path.exists(self.independent_ave_file) and maxsleep > 0:
        sleep(2)
        maxsleep -= 2

      try:
        fil = open(self.independent_ave_file, 'r')
        try:
          self.independent_average = float(fil.readline().strip())
        except:
          pass
        fil.close()
      except:
        pass

    return self.overall_status

  # write out any error files
  def writeErrors(self):
    if self.jobs_error:
      print "\nErrors occurred while running the bcl: tails of log files with errors:"
      got_errs = False
      while not got_errs:
        try:
          fil = open(self.status_filename + '.bclerr', 'r')
          try:
            print ''.join(fil.readlines())
            got_errs = True
          except:
            pass
          fil.close()
        except:
          pass
    if self.jobs_pbs_error:
      print "\nErrors occurred due to job requirements, details:\n"
      got_errs = False
      err_str = ''
      while not got_errs:
        if os.path.exists(self.status_filename + '.pbserr'):
          try:
            fil = open(self.status_filename + '.pbserr', 'r')
            try:
              err_str = ''.join(fil.readlines())
              got_errs = len(err_str.strip()) > 0
            except:
              pass
            fil.close()
          except:
            pass
        else:
          sleep(1)
      if got_errs:
        print err_str

  # write the current status
  def writeStatus(self):
    tail = ''
    if self.jobs_error:
      tail += 'Jobs failed due to errors while running bcl model:Train!\n'
    elif self.jobs_pbs_error:
      tail += 'Jobs failed due to errors in job requirements (memory, walltime, gpu, or ppn)\n'
    elif self.jobs_merge_status == BclCVStatus.ERRORS:
      tail += 'Jobs failed due to errors while merging predictions with bcl model:MergePredictions and model:ComputeJuryStatistics!\n'
    else:
      tail = ''
      if self.jobs_finished and self.independent_average:
        tail += ' Raw independent average: ' + str(self.independent_average)[:5:1] + ' '
      tail += ' Jobs status: '
      if self.jobs_merge_status > 0:
        if self.jobs_merge_status == BclCVStatus.RUNNING:
          tail += 'All CV jobs finished; merging results and calculating final objective function on consensus predictor'
        elif self.jobs_merge_status == BclCVStatus.FINISHED:
          tail += 'Reading final consensus predictor result'
      else:
        if self.jobs_finished:
          tail += str(self.jobs_finished) + ' finished '
        if self.jobs_running:
          tail += str(self.jobs_running) + ' running '
        if self.jobs_in_queue:
          tail += str(self.jobs_in_queue) + ' enqueued'
    total_jobs = self.jobs_in_queue + self.jobs_running + self.jobs_finished + self.jobs_error
    percent = int(100.0 * self.jobs_finished / float(max(total_jobs, 1)))
    if self.show_status:
      stars = int(percent / 5)
      star_str = '*' * stars
      space_str = ' ' * (20 - stars)
      next_line = "\rBCL model:Train Cross Validation [" + star_str + space_str + '] ' + str(int(percent)) + '% ' + tail
      additional_spaces = ""
      if len(next_line) < self.last_line_length:
        additional_spaces = ' ' * (self.last_line_length - len(next_line))
      if self.last_line.strip() != next_line.strip():
        sys.stdout.write(next_line + additional_spaces)
        self.last_line_length = len(next_line)
        sys.stdout.flush()
        self.last_line = next_line
    elif self.show_progress:
      if percent >= self.last_percent + 10:
        sys.stdout.write(str(percent) + "% ")
        sys.stdout.flush()
        self.last_percent = percent

  # update and write the current status
  def updateWriteStatus(self):
    self.updateStatus()
    self.writeStatus()

  # try to finalize this run; e.g. removing unnecessary files, displaying the final result, and writing any errors
  def tryFinalize(self, final_result_message = 'Consensus Result:'):
    if self.wasFinalized:
      return True
    self.updateStatus()
    if self.overall_status != BclCVStatus.RUNNING:
      self.wasFinalized = True
      if self.overall_status == BclCVStatus.ERRORS:
        self.writeErrors()
        if len(self.abort_command):
          (status, output) = tryExecute(self.abort_command, 2, 'aborting jobs', 1)
          if not status:
            print "Failed to abort jobs, message: " + output
      elif self.overall_status == BclCVStatus.FINISHED:
        results = self.getDatasetResults(self.result_file)
        split_result = results.split('\t')
        if split_result[0] != 'nan':
          self.final_result = float(split_result[0])
        else:
          self.final_result = None
        if len(split_result) > 3:
          direction = split_result[2].strip()
          if direction.startswith('Larger'):
            self.improvement_type = BclCVStatus.LARGER_IS_BETTER
          else:
            self.improvement_type = BclCVStatus.SMALLER_IS_BETTER
        print '\n' + final_result_message + ' ' + ''.join(results).strip()
        results = []
        if self.expect_blind_dataset or os.path.exists(self.result_file.replace('ind_merged', 'blind')):
          results = self.getDatasetResults(self.result_file.replace('ind_merged', 'blind'))
          split_result = results.split('\t')
          if split_result[0] != 'nan':
            self.final_blind_result = float(split_result[0])
            print 'Blind dataset result: ' + str(self.final_blind_result)
          else:
            self.final_blind_result = None
          results = self.getDatasetResults(self.result_file.replace('ind_merged', 'consensus'))
          split_result = results.split('\t')
          if split_result[0] != 'nan':
            self.final_consensus_result = float(split_result[0])
            print 'Visible dataset result: ' + str(self.final_consensus_result)
          else:
            self.final_consensus_result = None
        for fil in self.files_remove_on_success:
          rm_rf(fil)
      return True
    return False

  def getDatasetResults(self, filename):
    results = []
    for i in xrange(60):
      if os.path.exists(filename):
        res_file = open(filename, 'r')
        results = res_file.readlines()
        res_file.close()
        break
      else:
        sleep(3)
    if len(results) == 0:
      print "Could not read results file at " + str(self.result_file)
      sys.exit(1)
    return results[0]
