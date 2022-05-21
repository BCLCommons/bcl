#!/usr/bin/env python2.7
'''
Script to perform descriptor selection based on dataset scoring with methods like InformationGain, FScore, Input sensitivity
based on cross-validation

@author Mariusz Butkiewicz
@date 04/10/2013
'''
import BclModelCrossValidation as bcl_cv
import argparse
import ConfigParser
import os, sys, time, operator

class BclModelDescriptorSelection:

  def __init__(self, selection_type = None):
    self.selection_type = selection_type

  def addSelectionMethodOptions(self, parser):
    '''
    add descriptor selection flags to parser
    '''
    selection_group = parser.add_argument_group('selection options', 'Descriptor selection options')
    selection_group.add_argument('-r', '--range', nargs = 3, metavar = ('<min>', '<max>', '<step-size>'), default = [10, 20, 10], type = int, help = 'sampling range', dest = 'sample_range')

    return parser

  def run(self, given_flags = []):
    '''
    '''
    parser = bcl_cv.BclCommandCreator.getParser()
    parser = self.addSelectionMethodOptions(parser)

    if len(given_flags) == 0:
      (option_args, remaining_args) = parser.parse_known_args()
    else:
      (option_args, remaining_args) = parser.parse_known_args(given_flags)

    name = option_args.name
    if self.selection_type != None:
      name += '_' + self.selection_type

    gather_status = []
    rounds = []
    results = {}

    for round in xrange(option_args.sample_range[0], option_args.sample_range[1] + 1, option_args.sample_range[2]):

      print "\nTraining with top ", round, " feature columns: "
      option_args.name = name + '_top_' + str(round) + '_features'
      option_args.top_n_features = round

      if not option_args.local_jobs > 0:
        option_args.just_submit = True
        option_args.show_status = False

      commander = bcl_cv.BclCommandCreator(option_args, remaining_args)
      outputf = open(commander.log_files_path + os.sep + 'command.txt', 'w')
      outputf.write(' '.join(sys.argv))
      outputf.close()
      commander.run()

      gather_status.append(commander.status_bar)
      rounds.append(round)

    had_errors = False
    all_finished = False
    length_last_line = 0
    while not had_errors and not all_finished:
      jobs_in_queue = jobs_running = jobs_finished = jobs_error = jobs_pbs_error = jobs_merge_status = 0

      had_errors = False
      all_finished = True
      for i in xrange(len(rounds)):
        cmd = gather_status[i]
        if not cmd.tryFinalize('Result with ' + str(rounds[i]) + ' descriptors: '):
          all_finished = False
        elif cmd.overall_status == cmd.ERRORS:
          had_errors = True
        else:
          results[rounds[i]] = cmd.final_result

        jobs_in_queue += cmd.jobs_in_queue
        jobs_running += cmd.jobs_running
        jobs_finished += cmd.jobs_finished
        jobs_error += cmd.jobs_error
        jobs_pbs_error += cmd.jobs_pbs_error
        if cmd.jobs_merge_status == bcl_cv.BclCVStatus.RUNNING:
          jobs_merge_status += 1

      if not option_args.local_jobs > 0:
        line = '\rPBS Status:  q:' + str(jobs_in_queue)
        line += ' r:' + str(jobs_running)
        line += ' f:' + str(jobs_finished)
        line += ' e:' + str(jobs_error)
        line += ' pe:' + str(jobs_pbs_error)
        line += ' m:' + str(jobs_merge_status)

        if length_last_line > len(line):
          line += ' ' * (length_last_line - len(line))
        sys.stdout.write(line)
        sys.stdout.flush()
        length_last_line = len(line)

        time.sleep(3)

    if had_errors:
      print "Had errors, exiting"
      sys.exit(1)

    print "\nCombine final results .."

    final_result_filename = 'results/final_cv_' + name + '_result.txt'

    smaller_is_better = gather_status[0].improvement_type == bcl_cv.BclCVStatus.SMALLER_IS_BETTER
    sorted_results = sorted(results.iteritems(), key = operator.itemgetter(1), reverse = smaller_is_better)

    print 'results', results
    print 'sorted res:', sorted_results

    final_resultf = open(final_result_filename, 'w')
    final_resultf.write('# Descriptors\tResult\n')
    final_resultf.write('\n'.join([str(x[0]) + '\t' + str(x[1]) for x in sorted_results]))
    final_resultf.close()

    print 'Best result: Round: ' + str(sorted_results[-1][0]) + " Value: " + str(sorted_results[-1][1])
    print 'Combine final results .. done.'
    option_args.name = name

def main():
 '''
 '''
 BclModelDescriptorSelection().run()

if __name__ == '__main__':
  main()
