#!/usr/bin/env python2.7
'''
Script to score a dataset using scoring methods like InformationGain, FScore, Input sensitivity

@author Mariusz Butkiewicz
@date 04/10/2013
'''

import argparse
import ConfigParser
import os, sys
import types

class BclModelDatasetScore:

  def __init__(self):
    pass

  @staticmethod
  def getParser():
    '''
    get an instance of ArgumentParser for dataset score
    '''
    parser = argparse.ArgumentParser()
    return parser


  def getSelectionMethodOptions(self, parser):
    '''
    add descriptor selection flags to parser 
    '''
    selection_group = parser.add_argument_group('selection options', 'These options include descriptor selection options to chose')
    choices = ['InformationGain', 'FScore', 'Combine']
    help = 'InformationGain\n'
    help += 'FScore\n'
    help += 'InformationGain * FScore\n'
    selection_group.add_argument('-s', '--scoring-type', choices = choices, help = help, dest = 'scoring_type', required = True)

    selection_group = parser.add_argument_group('mandatory flags', '')
    selection_group.add_argument('-d', '--datasets', nargs = '*', action = 'append', help = 'datasets the descriptor scores are based on', dest = 'datasets', required = True)
    selection_group.add_argument('-b', '--bcl', help = 'bcl executable, default: bcl.exe', dest = 'bcl', default = 'bcl.exe')
    selection_group.add_argument('-o', '--output_score_file', help = 'score output file containing the scores for every feature column, default: score.out', dest = 'score_file', required = True)
    selection_group.add_argument('-c', '--cutoff', help = 'cutoff between categories of data', dest = 'cutoff', required = True)
    selection_group.add_argument('-f', '--features', help = 'features file name', dest = 'features', required = True)
    parser.add_argument \
    (\
      "--blind", \
      help = "Optional partition/chunk id of the blind chunk, which is never evaluated by the model. " \
             + " The first parameter is the zero-indexed blind id, the second is the number of partitions for the blind. " \
             + " E.g. to ignore the second third of the dataset completely during this cross-validation, --blind 1 3", \
      nargs = 2, \
      dest = 'blind_id_nchunks', \
      type = int \
    )

    return parser

  def getGnuplotPng(self, score_file, scoring_type):
    '''
    parse dataset score file and plot results with gnuplot
    '''

    if os.path.exists(score_file + '.png'):
      return

    main_cmd = 'grep -A2 \'bcl::linal::Vector<float>\' ' + score_file + ' | tail -n1 | tr \'\t\' \'\n\' | grep -v \'e\' | sort -nr '
    main_cmd += '> ' + score_file + '.rawplot'

    if os.system(main_cmd) != 0:
      print "Could not execute BCL model:Score with: " + main_cmd
      sys.exit(1)

    gnuplot_file = open(score_file + '.gnuplot', 'w')
    gnuplot_file.write('set term png size 600,400\n')
    gnuplot_file.write('set output "' + score_file + '.png"\n')
    gnuplot_file.write('set grid\n')
    gnuplot_file.write('set title "' + scoring_type + '"\n')
    gnuplot_file.write('plot "' + score_file + '.rawplot"\n')
    gnuplot_file.close()

    main_cmd = 'gnuplot ' + score_file + '.gnuplot >& /dev/null'

    if os.system(main_cmd) != 0:
      print "Could not execute gnuplot with: " + main_cmd
      sys.exit(1)

    print '\nGenerate GnuPlot png file with sorted dataset scores. You can look at the png file with:'
    print 'gthumb ' + score_file + '.png'

    # clean up intermediate files
    try:
      os.remove(score_file + '.rawplot')
      os.remove(score_file + '.gnuplot')
    except OSError:
      print 'Could not clean up intermediate gnuplot files!'


  def getBclScoreString(self, selection, cutoff):
    '''
    '''
    sdict = {}
    sdict['InformationGain'] = str('Partition(partitioner=InformationGain,cutoff=' + str(cutoff) + ')')
    sdict['FScore'] = 'FScore(cutoff=' + str(cutoff) + ')'
    sdict['Combine'] = 'Multiply(FScore(cutoff=' + str(cutoff) + '),Partition(partitioner=InformationGain,cutoff=' + str(cutoff) + '))'

    return sdict[selection]


  def run(self, given_flags = []):
    '''
    '''
    parser = self.getParser()
    parser = self.getSelectionMethodOptions(parser)

    if len(given_flags) == 0:
      (opt, remaining_args) = parser.parse_known_args()
    else:
      (opt, remaining_args) = parser.parse_known_args(given_flags)

    bcl = str(opt.bcl)
    cutoff = str(opt.cutoff)
    self.has_blind_set = False
    if opt.blind_id_nchunks and len(opt.blind_id_nchunks) == 1:
      opt.blind_id_nchunks.append(opt.chunks)
    elif opt.blind_id_nchunks and len(opt.blind_id_nchunks) == 0:
      opt.blind_id_nchunks = None

    if opt.blind_id_nchunks:
      opt.name += '_Blind' + str(opt.blind_id_nchunks[0])
      self.has_blind_set = True
    blinder = ''
    if self.has_blind_set:
      blinder = ',number chunks=' + str(opt.blind_id_nchunks[1]) + ',chunks="[0,' + str(opt.blind_id_nchunks[1]) + ') - [' + str(opt.blind_id_nchunks[0]) + ']"'

    datasets = "Combined( " + ",".join([ str("Subset(filename=" + blinder + "".join(ds) + ")") for ds in opt.datasets[0]]) + ")"

    score_file = str(opt.score_file)
    scoring_type = self.getBclScoreString(str(opt.scoring_type), cutoff)
    descriptor_labels_file = str(opt.features)

    if os.path.exists(score_file):
      print 'Score file ', score_file, ' exists already! Skipping dataset scoring!'
      self.getGnuplotPng(score_file, str(opt.scoring_type))
      return

    print 'Chosen score type: ', scoring_type

    main_cmd = bcl + " descriptor:ScoreDataset -source '" + datasets + "' -output " + score_file + " -score '" + scoring_type + "' -feature_labels " + descriptor_labels_file + " -opencl Disable"

    if os.system(main_cmd) != 0:
      print "Could not execute BCL model:Score with: " + main_cmd
      sys.exit(1)

    self.getGnuplotPng(score_file, str(opt.scoring_type))


def main():
 '''
 '''
 BclModelDatasetScore().run()

if __name__ == '__main__':
  main()
