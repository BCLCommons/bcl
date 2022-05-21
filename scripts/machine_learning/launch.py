#!/usr/bin/env python2.7
'''
Script to select one of the following applications
  - dataset scoring using methods like InformationGain, FScore, InputSensitivity
  - descriptor selection based on a scored dataset
  - n x m - fold cross-validation

@author Mariusz Butkiewicz
@date 04/10/2013
'''

from tasks import BclModelCrossValidation as bcl_cv
from tasks import BclModelDatasetScore as bcl_score
from tasks import BclModelDescriptorSelection as bcl_desc_sel
from tasks import BclModelDependentDescriptorSelection as bcl_model_desc_sel
import argparse
import ConfigParser
import os, sys
import types

def get_parser():
  parser = argparse.ArgumentParser\
  (\
    epilog = 'Additional flags with arguments will be passed on to the bcl unchanged. '\
             'Flags shown in the above help (except --config-file and --type) can also be stored in a configuration file, '
             '--datasets x.bin y.bin --max-minutes 60 in the config file would be written as\n' \
             'datasets: [x.bin,y.bin]\n' \
             'max-minutes: 60', \
    add_help = False \
  )
  return parser


def get_flags(parser):
  '''
  add descriptor selection flags to parser
  '''
  selection_group = parser.add_argument_group('selection options', 'select whether to perform a simple cross validation or descriptor selection method')
  choices = ['dataset_scoring', 'feature_selection', 'cross_validation', 'feature_selection_by_model']
  help = 'dataset_scoring -- Score a dataset \n'
  help += 'feature_selection -- Perform feature optimization given a score file \n'
  help += 'cross_validation -- Cross validate a model \n'
  help += 'feature_selection_by_model -- Uses a model-dependent iterative scoring function '
  selection_group.add_argument('-t', '--type', choices = choices, help = help, dest = 'selection_type', required = True)
  parser.add_argument(\
    '--config-file', \
    dest = 'config_file', \
    metavar = '<config_f>', \
    default = 'config.ini', \
    help = 'Configuration file to use, if not config.ini' \
  )
  parser.add_argument(\
    '-h', '--help', \
    help = 'show this help message and then exit. If --type is set, also show the help for the selection method before exiting', \
    dest = 'help', \
    action = 'store_true', \
    default = False \
  )
  return parser


def parse_config_file(opts):
  ''' 
  parse the config file for training / crossvalidating models
  '''
  config = ConfigParser.SafeConfigParser()
  if os.path.exists(opts.config_file):
    config.read(opts.config_file)
  return config


def get_default_flags(opts, remaining_args):
  '''
  helper function to reconstruct flags as a dictionary from the 'variables' section in the config file.
  '''
  flags = config_section_map('variables', parse_config_file(opts), {}, False)
  flags = dict(zip(flags[::2], flags[1::2]))
  ra_dict = {}

  ra_index = 0
  while ra_index < len(remaining_args):
    if remaining_args[ra_index].startswith('--'):
      tag = remaining_args[ra_index].lstrip('-')
      ra_index += 1
      values = []
      while ra_index < len(remaining_args) and not remaining_args[ra_index].startswith('-'):
        values.append(remaining_args[ra_index])
        ra_index += 1
      if len(values) == 1:
        ra_dict[tag] = values[0]
      elif len(values) == 0:
        ra_dict[tag] = ''
      else:
        ra_dict[tag] = values
    else:
      ra_index += 1
  used_tags = set()
  for tag, values in flags.iteritems():
    if tag in ra_dict:
      flags[tag] = ra_dict[tag]
      used_tags.add(tag)

  remaining_args_reconstructed = []
  ra_index = 0
  while ra_index < len(remaining_args):
    if remaining_args[ra_index].startswith('--'):
      tag = remaining_args[ra_index].lstrip('-')
      ra_index += 1
      if tag not in used_tags:
        remaining_args_reconstructed.append(remaining_args[ra_index - 1])
        continue
      while ra_index < len(remaining_args) and not remaining_args[ra_index].startswith('-'):
        ra_index += 1
    else:
      remaining_args_reconstructed.append(remaining_args[ra_index])
      ra_index += 1
  remaining_args = remaining_args_reconstructed
  return flags, remaining_args_reconstructed

def config_section_map(section, config, default_flags, key_as_flag = True):
  '''
  helper function to exctract the flag information of a specific section in the config file
  '''
  sdict = {}
  options = []
  if config.has_section(section):
    options = config.options(section)

  for option in options:
    try:
      sdict[option] = config.get(section, option, 0, default_flags)
      if sdict[option] == -1:
        DebugPrint("skip: %s" % option)
    except:
      print("exception on %s!" % option)
      sdict[option] = None

  flags = []
  for key, value in sdict.iteritems():
    if key_as_flag == True:
      flags.append("--" + key)
    else:
      flags.append(key)

    if value == None or len(value) == 0:
      continue

    try:
      if isinstance(eval(value), (list, tuple)):
        for item in eval(value):
          flags.append(str(item).strip('"\''))
      else:
        flags.append(str(value).strip('"\''))
    except:
      flags.append(str(value).strip('"\''))

  return flags

def run_cross_validation(opts, remaining_args):
  '''
  run cross-validation with parameters read from the config file.
  The final model will be evaluated by a ROC curve type plot or a correlation plot.
  '''
  parser = bcl_cv.BclCommandCreator.getParser()

  config = parse_config_file(opts)
  default_flags, remaining_args = get_default_flags(opts, remaining_args)
  flags = []
  flags.extend(config_section_map('bcl', config, default_flags))
  flags.extend(config_section_map('main', config, default_flags))
  flags.extend(config_section_map('learning', config, default_flags))
  flags.extend(config_section_map('cv', config, default_flags))
  flags.extend(config_section_map('dataset', config, default_flags))
  flags.extend(remaining_args)

  (option_args, remaining_args) = parser.parse_known_args(flags)

  commander = bcl_cv.BclCommandCreator(option_args, remaining_args)

  outputf = open(commander.log_files_path + os.sep + 'command.txt', 'w')
  outputf.write(' '.join(sys.argv))
  reader = open(opts.config_file, 'r')
  outputf.write(''.join(reader.readlines()))
  reader.close()
  outputf.close()

  commander.run()

def run_ann_connection_weights(opts, remaining_args):
  '''
  run ann-based connection weights
  '''
  parser = bcl_ann_cw.BclANNDescriptorSelectionByConnectionWeights.getParser()

  config = parse_config_file(opts)
  default_flags, remaining_args = get_default_flags(opts, remaining_args)
  flags = []
  flags.extend(config_section_map('bcl', config, default_flags))
  flags.extend(config_section_map('main', config, default_flags))
  flags.extend(config_section_map('learning', config, default_flags))
  flags.extend(config_section_map('cv', config, default_flags))
  flags.extend(config_section_map('dataset', config, default_flags))
  flags.extend(config_section_map('ann_cw', config, default_flags))
  flags.extend(remaining_args)

  (option_args, remaining_args) = parser.parse_known_args(flags)

  print "option_args.ds_rounds: ", option_args.ds_rounds

  commander = bcl_ann_cw.BclANNDescriptorSelectionByConnectionWeights(option_args, remaining_args)
  outputf = open(commander.cross_validation_runner.log_files_path + os.sep + 'command.txt', 'w')
  outputf.write(' '.join(sys.argv))
  reader = open(opts.config_file, 'r')
  outputf.write(''.join(reader.readlines()))
  reader.close()
  outputf.close()

  commander.run()

def run_md_descriptor_selection(opts, remaining_args):
  '''
  run model-dependent feature selection
  '''
  parser = bcl_model_desc_sel.BclModelDependentDescriptorSelection.getParser()

  config = parse_config_file(opts)
  default_flags, remaining_args = get_default_flags(opts, remaining_args)
  flags = []
  flags.extend(config_section_map('bcl', config, default_flags))
  flags.extend(config_section_map('main', config, default_flags))
  flags.extend(config_section_map('learning', config, default_flags))
  flags.extend(config_section_map('cv', config, default_flags))
  flags.extend(config_section_map('dataset', config, default_flags))
  flags.extend(config_section_map('descriptor-selection-model-dependent', config, default_flags))
  flags.extend(remaining_args)

  (option_args, remaining_args) = parser.parse_known_args(flags)

  print "option_args.ds_rounds: ", option_args.ds_rounds

  commander = bcl_model_desc_sel.BclModelDependentDescriptorSelection(option_args, remaining_args)
  outputf = open(commander.cross_validation_runner.log_files_path + os.sep + 'command.txt', 'w')
  outputf.write(' '.join(sys.argv))
  reader = open(opts.config_file, 'r')
  outputf.write(''.join(reader.readlines()))
  reader.close()
  outputf.close()

  commander.run()


def run_dataset_scoring(opts, remaining_args):
  '''
  run score evaluation of one or many datasets given in the config file
  '''
  config = parse_config_file(opts)
  default_flags, remaining_args = get_default_flags(opts, remaining_args)
  flags = []
  flags.extend(config_section_map('bcl', config, default_flags))
  flags.extend(config_section_map('variables', config, default_flags))
  flags.extend(config_section_map('score', config, default_flags))
  flags.extend(config_section_map('dataset', config, default_flags))
  flags.extend(remaining_args)

  bcl_score.BclModelDatasetScore().run(flags)


def run_descriptor_selection(opts, remaining_args):
  '''
  run descriptor selection based on a given dataset score file and cross-validation information from the config file 
  '''
  config = parse_config_file(opts)
  default_flags, remaining_args = get_default_flags(opts, remaining_args)
  flags = []
  flags.extend(config_section_map('bcl', config, default_flags))
  flags.extend(config_section_map('main', config, default_flags))
  flags.extend(config_section_map('learning', config, default_flags))
  flags.extend(config_section_map('cv', config, default_flags))
  flags.extend(config_section_map('dataset', config, default_flags))
  flags.extend(config_section_map('descriptor-selection', config, default_flags))
  flags.extend(['--score-file', str(config.get('score', 'output_score_file'))])
  flags.extend(remaining_args)

  bcl_desc_sel.BclModelDescriptorSelection(config.get('score', 'scoring-type')).run(flags)

def main():
  '''
  main method
  '''
  parser = get_parser()
  parser = get_flags(parser)

  try:
    (opts, remaining_args) = parser.parse_known_args()
  except SystemExit:
    # There may have been a validation error, but the user may have also specified help, so be sure to honor that request
    # for example submit.py -t not_an_option --help
    if '-h' in os.sys.argv or '--help' in os.sys.argv or len(os.sys.argv) == 1:
      parser.print_help()
    raise
  
  # test whether help was requested
  if opts.help:
    parser.print_help()
    remaining_args.append('-h')
    if not opts.selection_type:
      # print help and exit; case where no selection method was specified, 
      sys.exit(0)
    sys.stderr.write('---------------------------------------------------------------------------\n')
    sys.stderr.write('Flags and options available for selection mode: ' + opts.selection_type + '\n')
    sys.stderr.write('---------------------------------------------------------------------------\n')
    sys.stderr.write('\n')
  
  if len(opts.config_file) and not os.path.exists(opts.config_file) and not opts.config_file == 'config.ini':
    print "No configuration file found at " + opts.config_file
    sys.exit(-1)

  # if help was requested, it will be handled below
  if opts.selection_type == 'cross_validation':
    run_cross_validation(opts, remaining_args)
  elif opts.selection_type == 'dataset_scoring':
    run_dataset_scoring(opts, remaining_args)
  elif opts.selection_type == 'feature_selection':
    run_dataset_scoring(opts, remaining_args)
    run_descriptor_selection(opts, remaining_args)
  elif opts.selection_type == 'feature_selection_by_model':
    run_md_descriptor_selection(opts, remaining_args)
  else:
    print "No implementation found for ", opts.selection_type

if __name__ == '__main__':
    main()
