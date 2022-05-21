'''
Run a command given in a queue in a separate thread represented by this class.
Class also supports a status bar
Created April 2013
@author: mendenjl
'''
import commands, threading, Queue

# runs a job locally; requires a status bar variable with an updateWriteStatus method
class ThreadLocalRun(threading.Thread):

  def __init__(self, queue, status_bar):
    threading.Thread.__init__(self)
    self.queue = queue
    self.status_bar = status_bar

  def run(self):
    while True:
      #grabs info from the queue
      try:
        command_to_run = self.queue.get(True, 5)
      except Queue.Empty:
        # exit; this prevents threads from living on forever when nothing remains in the queue
        exit()
      if self.status_bar != None:
        self.status_bar.updateWriteStatus()

      # insert into the out queue after calling deheader
      (status, output) = commands.getstatusoutput(command_to_run)

      if status != 0:
        print "error during job submission: " + output + ' of ' + command_to_run
      if self.status_bar != None:
        self.status_bar.updateWriteStatus()

      #signals to queue job is done
      self.queue.task_done()
      
      
