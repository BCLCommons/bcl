#!/usr/bin/env python2.7
'''
Script to perform descriptor selection based on a model-dependent scoring functions

@author Jeffrey Mendenhall
@date 05/07/2013
'''
import BclModelCrossValidation as bcl_cv
import argparse
import ConfigParser
import os, sys, time, operator, multiprocessing
from curses.ascii import isalpha

class BclModelDependentDescriptorSelection:

  @staticmethod
  def getParser():
    parser = bcl_cv.BclCommandCreator.getParser()
    ds_weights = parser.add_argument_group('descriptor selection weights', 'Weights for different components of the descriptor selection score')
    ds_weights.add_argument \
    (\
      '--weight-consistency', \
      type = float, \
      default = 0.0, \
      dest = 'consistency', \
      help = 'the influence of sign consistency on the final result.'\
            ' Higher values place more emphasis on columns that have the same effect on the output across the networks'\
    )
    ds_weights.add_argument \
    (\
      '--weight-consistency-best', \
      default = 0.0, \
      type = float, \
      dest = 'consistency_best', \
      help = 'the influence of sign consistency of the best models on the final result. See note on --weight-cohen about model performance assessment' \
    )
    ds_weights.add_argument \
    (\
      '--weight-sqr', \
      default = 0.0, \
      type = float, \
      dest = 'sqrweight', \
      help = 'influence of the average (weight influence ^ 2) (and rescaled) on the score output.'\
            ' Higher values emphasize columns with the largest weights, independently of sign'\
    )
    ds_weights.add_argument \
    (\
      '--weight-abs', \
      default = 0.0, \
      type = float, \
      dest = 'absweight', \
      help = 'influence of the average (weight) (and rescaled) on the score output.'\
            ' Higher values emphasize columns with the largest weights in the same direction'\
    )
    ds_weights.add_argument \
    (\
      '--weight-utility', \
      default = 0.0, \
      type = float, \
      dest = 'utilityweight', \
      help = "Change the weight for the utility function of input sensitivity.  The utility function considers whether models " \
             "that gave incorrect output for a particular feature consistently used the feature more than the models that performed " \
              "well on that particular feature"\
    )
    selection_group = parser.add_argument_group('selection options', 'Descriptor selection options')
    selection_group.add_argument \
    (\
      '--min-features', \
      default = 1, \
      type = int, \
      dest = 'min_features', \
      help = 'minimum number of features to allow'\
    )
    selection_group.add_argument \
    (\
      '--ds-rounds', \
      default = 10, \
      type = int, \
      dest = 'ds_rounds', \
      help = 'number of descriptor selection rounds to perform'\
    )
    selection_group.add_argument\
    (\
       "--ds-best-models", \
       help = "Set this flag to use only the best model (based on final-objective-function from each independent chunk " \
              "for computing the model-dependent descriptor selection", \
       dest = 'ds_only_best_models', \
       action = 'store_true', \
       default = False \
    )
    selection_group.add_argument \
    (\
      '--min-features-removed-per-round', \
      default = 1, \
      type = int, \
      dest = 'min_features_removed', \
      help = 'minimum number of features to remove per round'\
    )
    selection_group.add_argument \
    (\
      '--max-features-removed-per-round', \
      default = 100, \
      type = int, \
      dest = 'max_features_removed', \
      help = 'maximum number of features to remove per round'\
    )
    #method = selection_group.add_mutually_exclusive_group()
    selection_group.add_argument \
    (\
      '--attrition-rate', \
      default = None, \
      type = float, \
      dest = 'attrition', \
      help = 'desired fraction of features removed per round; bounded by --min-features-removed-per-round ' \
             'and --max-features-removed-per-round'\
    )
    sensitivity = parser.add_argument_group\
    (\
      'Sensitivity-related options', \
      'Flags related to the use of input sensitivity (IS). IS is generally more accurate than looking at neural network '\
      'weights, but is much slower to compute.  IS is more general, however, and will work with any machine learning method'\
    )
    sensitivity.add_argument \
    (\
      '--sensitivity', \
      dest = 'sensitivity', \
      default = False, \
      action = 'store_true', \
      help = 'Give this flag to use input sensitivity rather than neural network weights (the default) to comptue pseudo-derivative. '\
             'This enables the use of all --is-... flags (input-sensitivity related flags)'
    )
    sensitivity.add_argument \
    (\
      '--scoring-host', \
      dest = 'score_host', \
      default = os.getenv('HOSTNAME'), \
      help = "The host to login to for scoring models" \
    )
    sensitivity.add_argument \
    (\
      '--is-threads', \
      dest = 'score_host_threads', \
      default = multiprocessing.cpu_count(), \
      type = int, \
      help = "The # of threads to use for computing input sensitivity" \
    )
    sensitivity.add_argument \
    (\
      '--is-feature-count', \
      dest = 'is_feature_count', \
      default = 200, \
      type = int, \
      help = "The # of features from the datasets to use for input sensitivity. " \
             "Unless --is-use-static-features, selects the best features for IS, which are usually features that are "\
             "predicted rather differently across the models (e.g. high standard deviation), or in the case of classification " \
             " targets, as close as possible to the same # of models with correct predictions as incorrect" \
    )
    sensitivity.add_argument \
    (\
      '--is-use-static-features', \
      dest = 'is_static_features', \
      default = False, \
      action = 'store_true', \
      help = "Causes a new dataset to be created (if it doesn't already exist) that will be used with input sensitivity.  "\
             "rather than selecting features that should be ideal for computing IS, uses a static set, which may be "\
             "beneficial in some cases, particularly if prediction is slow" \
    )
#    method.add_argument \
#    (\
#      '--score-cutoff', \
#      default = None, \
#      type = float, \
#      dest = 'score_cutoff', \
#      help = 'scores below this value will be removed each round. ' + \
#      ' NOTE: This feature does not currently respect --min-features-removed-per-round and --max-features-removed-per-round'\
#    )
    selection_group.add_argument \
    (\
      '--no-random-seed', \
      default = False, \
      dest = 'no_random_seed', \
      action = 'store_true', \
      help = 'set this flag if you do not want a random seed provided to each network'
    )
    selection_group.add_argument \
    (\
      '--continue', \
      dest = 'continued', \
      help = "Set to continue the descriptor selection run with the same id/name. If the round number is omitted, " \
      + "continues exactly where the last run left off. If the round number is given, proceeds using the models " \
      + " generated in the given round, rescores them, and redetermines the features to use", \
      metavar = ('<optional round number>'), \
      nargs = '?', \
      # the spaces is important for const, otherwise checking options.continued 
      # returns false because python evaluates empty strings as false! 
      const = -1, \
      default = -2 \
    )
    selection_group.add_argument \
    (\
      '--continue-after-scoring', \
      dest = 'continue_post_scoring', \
      help = "Set to continue the descriptor selection run with the same id/name with the round number after the scoring " \
             + "phase for that round.  This is appropriate if the score weights were left the same", \
      action = 'store_true', \
      default = False
    )
    return parser

  RND_STR = 'rnd$$$ROUND_N$$$'
  LASTRND_STR = 'rnd$$$LASTROUND_N$$$'
  N_FEATURES_STR = '$$$N_FEATURES$$$'
  INPUT_SENS_FEATURES = '$$$INPUT_SENS_FEATURES$$$'

  def __init__(self, option_args, remaining_args):
    # ensure reasonable feature values
    opt_scores = [option_args.consistency, option_args.consistency_best, option_args.sqrweight, option_args.absweight, option_args.utilityweight]

    if min(opt_scores) < 0.0:
      print "All weights must be >= 0 (minimum weight was: " + str(min(opt_scores))
      sys.exit(-1)
    if max(opt_scores) <= 0.0:
      print "At least one weight must be non-zero!"
      sys.exit(-1)
    if option_args.min_features < 1:
      print "--min_features must be given a number > 0, automatically setting to 1"
      option_args.min_features = 1
    if option_args.ds_rounds < 1:
      print "--ds-rounds < 1, no descriptor selection will be performed"
    if option_args.min_features_removed < 1:
      print "--min-features-removed changed to minimum value of 1"
      option_args.min_features_removed = 1
    if option_args.max_features_removed < option_args.min_features_removed:
      print "--max-features-removed must be given a value > --min-features-removed, autosetting max = min"
      option_args.max_features_removed = option_args.min_features_removed
    if option_args.iterate.find('NeuralNetwork') < 0 and not option_args.sensitivity:
      print "This descriptor selection method only works with neural networks, use --sensitivty to use input sensitivity instead"
      sys.exit(1)
    if option_args.attrition == None:
      print "No attrition rate given assuming it is equal to max-features-removed / number-features"

    self.options = option_args
    self.remaining_args = remaining_args
    if not self.options.no_random_seed and '-random_seed' not in self.remaining_args and '--random_seed' not in self.remaining_args:
      self.remaining_args.append('-random_seed')

    if self.options.blind_id_nchunks and len(self.options.blind_id_nchunks) == 1:
      self.options.blind_id_nchunks.append(option_args.chunks)
    elif self.options.blind_id_nchunks and len(self.options.blind_id_nchunks) == 0:
      self.options.blind_id_nchunks = None
    self.has_blind_set = False
    if self.options.blind_id_nchunks:
      self.options.name += '_Blind' + str(self.options.blind_id_nchunks[0])
      self.has_blind_set = True

    self.base_name = self.options.name
    rnd_str = BclModelDependentDescriptorSelection.RND_STR
    self.name_tmpl = self.base_name + "_" + rnd_str
    self.score_files = []
    self.options.model_storage = 'File'
    self.options.just_submit = False
    self.options.show_status = True
    self.options.print_ind = True
    self.current_round = 0
    self.round_results = {}
    self.round_results_raw = {}
    self.blind_results = {}
    self.visible_results = {}
    self.improvement_type = -1
    self.method = 'NeuralNetworkWeights'
    self.method_alias = 'ANN-Weight'
    self.uses_number_features = False
    if self.options.sensitivity:
      self.method = 'InputSensitivity'
      if self.options.iterate.strip().startswith('NeuralNetwork'):
        self.method = 'InputSensitivityNeuralNetwork'
      else:
        self.uses_number_features = True
      self.method_alias = 'Input-Sensitivity'

    original_training_size_limit = self.options.train_size_limit
    original_no_write_scripts = self.options.no_write_scripts
    self.options.train_size_limit = 2
    self.options.no_write_scripts = True
    self.cross_validation_runner = bcl_cv.BclCommandCreator(self.options, self.remaining_args)
    self.options.train_size_limit = original_training_size_limit
    self.options.no_write_scripts = original_no_write_scripts
    self.training_dataset_small = self.cross_validation_runner.replaceDynamicVariables(self.cross_validation_runner.training_dataset_tmpl, 0, 1, 0, 0)

    continuation_round = self.options.continued
    if continuation_round == -2:
      self.options.continued = False
    else:
      self.options.continued = True

    self.feature_counts = [self.cross_validation_runner.dataset_feature_cols]
    self.score_file_path = './feature-scores/' + self.base_name
    if not os.path.exists(self.score_file_path) or self.options.continued == False:
      bcl_cv.BclCommandCreator.mkdirmp(self.score_file_path)
    self.result_labels_command = ''
    if self.options.results:
      self.result_labels_command = ' --result_labels ' + self.options.results + ' '
    self.score_file_name_tmpl = self.score_file_path + os.sep + rnd_str + '.score'
    self.descriptor_file_name_tmpl = self.score_file_path + os.sep + rnd_str + '_cols' + BclModelDependentDescriptorSelection.N_FEATURES_STR + '.obj'
    self.select_log = self.cross_validation_runner.log_files_path + os.sep + "select_" + rnd_str + '.log'
    self.score_log = self.cross_validation_runner.log_files_path + os.sep + "score_" + rnd_str + '.log'
    self.select_command = self.options.bcl + " descriptor:RefineByScore -select '"

    self.last_round_features = self.options.features
    self.last_round_labels_command = ''
    if self.options.features:
      self.last_round_labels_command = ' -feature_labels ' + self.options.features + ' '

    # create the randomized dataset for input sensitivity targets
    self.score_host_access_cmd = ''
    if self.options.score_host != os.getenv('HOSTNAME'):
      self.score_host_access_cmd = 'ssh ' + self.options.score_host + ' cd ' + os.getcwd() + ' '
      full_cmd = self.score_host_access_cmd
      (status, output) = bcl_cv.tryExecute(full_cmd, 2, 'sshing into ' + self.options.score_host, 3)
      if not status:
        print "Could not ssh onto " + self.options.score_host + " with " + full_cmd + ' aborting'
        sys.exit(-1)
      self.score_host_access_cmd += ' \; '

    if self.options.sensitivity:
      n_features = self.options.is_feature_count
      n_datasets = len(self.options.datasets)
      features_per_dataset = n_features / n_datasets
      if self.options.is_static_features:
        rand_dataset_name = os.path.splitext(self.options.datasets[0])[0] + '.rand' + str(n_features) + '.bin'
        if not os.path.exists(rand_dataset_name):
          print "Creating dataset at " + rand_dataset_name
          complete_dataset = ''
          if n_datasets > 1:
            complete_dataset = 'Balanced('
          for data in self.options.datasets:
            complete_dataset += 'Chunk(chunks="[0, ' + str(features_per_dataset) + ')",dataset=Randomize(Subset(filename=' + data + '))),'
          complete_dataset = complete_dataset[:-1]
          if n_datasets > 1:
            complete_dataset += ')'
          complete_dataset = "'" + complete_dataset + "'"
          generate_script_name = self.cross_validation_runner.log_files_path + os.sep + "randomize_dataset.sh"
          script = open(generate_script_name, 'w')
          script.write("#!/bin/sh\n")
          full_cmd = self.options.bcl + ' descriptor:GenerateDataset -source ' + complete_dataset + ' -opencl Disable '
          full_cmd += self.last_round_labels_command + self.result_labels_command + ' -output ' + rand_dataset_name
          script.write(full_cmd + '\n')
          script.close()
          os.chmod(generate_script_name, 0755)
          (status, output) = bcl_cv.tryExecute(self.score_host_access_cmd + ' ' + generate_script_name, 2, 'Creating the randomized dataset', 3)
          if not status:
            print "Could not create the randomized dataset; error: " + output
            print "Command was: " + full_cmd
            sys.exit(-1)
          os.remove(generate_script_name)
        self.training_dataset_small = 'Subset(filename=' + rand_dataset_name + ')'
      else:
        complete_dataset = ''
        if n_datasets > 1:
          complete_dataset = 'Combined('
        blinder = ''
        if self.has_blind_set:
          blinder = ',number chunks=' + str(self.options.blind_id_nchunks[1]) + ',chunks="[0,' + str(self.options.blind_id_nchunks[1]) + ') - [' + str(self.options.blind_id_nchunks[0]) + ']"'
        for data in self.options.datasets:
          complete_dataset += 'Subset(filename=' + data + blinder + '),'
        complete_dataset = complete_dataset[:-1]
        if n_datasets > 1:
          complete_dataset += ')'
        self.training_dataset_small = complete_dataset
    #blind_id_nchunks
    if self.options.attrition == None:
      self.options.attrition = float(self.options.max_features_removed) / float(self.feature_counts[0])
    self.select_command += "Top(" + BclModelDependentDescriptorSelection.N_FEATURES_STR + ")' "
    self.select_command += " -score_file " + self.score_file_name_tmpl + " -output " + self.descriptor_file_name_tmpl + \
                           ' -logger File ' + self.select_log
    self.score_command = self.options.bcl + " descriptor:ScoreDataset -source '" + self.training_dataset_small + "' " \
      + "-output " + self.score_file_name_tmpl + " -opencl Disable -score '" \
      + self.method + "("
    if self.options.sensitivity and self.uses_number_features:
      self.score_command += 'delta=1.0,'
      if not self.options.is_static_features:
         self.score_command += 'feature limit=' + BclModelDependentDescriptorSelection.INPUT_SENS_FEATURES + ','
    balancing_str = 'balance=True,categorical=False'
    if self.options.final_obj:
      if self.options.final_obj.find('CategoricalMax') >= 0:
        balancing_str = 'balance=False,categorical=True'
      elif self.options.final_obj.find('Accuracy') >= 0:
        balancing_str = 'balance=False,categorical=False'
    self.score_command += \
      "storage=File(directory=./models/" + self.name_tmpl.replace(self.RND_STR, self.LASTRND_STR) + ",prefix=model" \
      + ("" if not self.options.ds_only_best_models else ",pick best=True") + "),"\
      + "weights=(consistency=" + str(self.options.consistency) + ",square=" + str(self.options.sqrweight) \
      + ",absolute=" + str(self.options.absweight) + ",utility=" + str(self.options.utilityweight) \
       + ",consistency best=" + str(self.options.consistency_best) + ',' + balancing_str \
      + "))' " + self.result_labels_command + ' -logger File ' + self.score_log
    if self.options.score_host_threads:
      self.score_command += ' -scheduler PThread ' + str(self.options.score_host_threads) + ' '

    self.labels_command = ' -feature_labels ' + self.descriptor_file_name_tmpl + ' '
    self.final_output_file = self.cross_validation_runner.independent_files_path + os.sep + 'final_descriptor_selection.txt'
    self.continuing = self.options.continued
    self.old_feature_counts = []
    if self.continuing:
      existing_obj_rnd_cols = [ str(x) for x in os.listdir(self.score_file_path) if str(x).endswith('.obj')]
      alpha = ''.join([str(chr(x)) for x in xrange(256) if isalpha(chr(x))]) + '.'
      feature_cols_dict = {}
      continuation_round = int(continuation_round)
      self.max_round = continuation_round
      if continuation_round == -1:
        self.max_round = None
      for filename in existing_obj_rnd_cols:
        rnd_cols = filename.translate(None, alpha).split('_')
        round = int(rnd_cols[0])
        cols = int(rnd_cols[1])
        if self.max_round == None or round < self.max_round:
          feature_cols_dict[round] = cols
        else:
          os.remove(self.score_file_path + os.sep + filename)
      last_rnd = 0
      for rnd, cols in sorted(feature_cols_dict.iteritems(), key = operator.itemgetter(0)):
        if last_rnd + 1 == rnd:
          self.old_feature_counts.append(cols)
          last_rnd += 1
        else:
          break

  def performIteration(self):
    # update name & features flag
    self.options.name = self.replaceDynamicVariables(self.name_tmpl)

    if self.current_round:
      self.options.features = self.replaceDynamicVariables(self.descriptor_file_name_tmpl)

    bcl_status = None
    if self.continuing:
      bcl_status = bcl_cv.BclCommandCreator.CreateStatus(self.options)
      if bcl_status.updateStatus() != bcl_status.FINISHED:
        self.continuing = False
      else:
        bcl_status.tryFinalize()
    if not self.continuing:
      # perform cross-validation
      print "\nTraining models for round " + str(self.current_round) + " of " + self.method_alias + \
            " based feature minimization with " + str(self.feature_counts[-1]) + " features "
      commander = bcl_cv.BclCommandCreator(self.options, self.remaining_args)
      commander.run()
      bcl_status = commander.status_bar
      if bcl_status.overall_status == bcl_status.ERRORS:
        print "stopping due to errors"
        sys.exit(-1)
    self.round_results[self.current_round] = bcl_status.final_result
    self.round_results_raw[self.current_round] = bcl_status.independent_average
    if self.has_blind_set:
      self.blind_results[self.current_round] = bcl_status.final_blind_result
      self.visible_results[self.current_round] = bcl_status.final_consensus_result
    self.improvement_type = bcl_status.improvement_type

    # update round number
    self.current_round += 1

    # update the # of features
    current_feature_count = self.feature_counts[-1]
    if current_feature_count == self.options.min_features:
      return
    # calculate the # of features desired using the attrition ratio
    features_using_attrition = int(current_feature_count * (1.0 - self.options.attrition))
    target_features_this_rnd = min(current_feature_count - self.options.min_features_removed, features_using_attrition)
    min_features_bnd_this_rnd = max(current_feature_count - self.options.max_features_removed, self.options.min_features)
    new_feature_count = max(min_features_bnd_this_rnd, target_features_this_rnd)

    # append # features to list
    self.feature_counts.append(new_feature_count)

    if not self.options.continue_post_scoring:
      if self.continuing and self.max_round != None and self.current_round == self.max_round:
        self.continuing = False
      if self.continuing:
        if len(self.old_feature_counts) < self.current_round:
          self.continuing = False
        else:
          self.feature_counts[-1] = self.old_feature_counts[self.current_round - 1]

    # perform scoring

    # add --feature_labels if not on round 1 or if feature labels was given
    full_command = self.replaceDynamicVariables(self.score_command + self.last_round_labels_command)
    if self.options.dry_run:
      print "Would execute: " + str(full_command)
    else:
      if self.continuing and not os.path.exists(self.replaceDynamicVariables(self.score_file_name_tmpl)):
        self.continuing = False
      if not self.continuing:
        print "Scoring models from round " + str(self.current_round - 1) + " to yield " + self.replaceDynamicVariables(self.score_file_name_tmpl)
        if len(self.score_host_access_cmd) == 0:
          (status, output) = bcl_cv.tryExecute(full_command, 2, "Scoring models", 10)
        else:
          score_script_name = self.cross_validation_runner.log_files_path + os.sep + "score.sh"
          script = open(score_script_name, 'w')
          script.write("#!/bin/sh\n")
          script.write(full_command)
          script.close()
          os.chmod(score_script_name, 0755)
          (status, output) = bcl_cv.tryExecute(self.score_host_access_cmd + score_script_name, 2, 'Scoring the models', 3)
        if not status:
          print "Error while scoring models using command: "
          print full_command

          sys.exit(-1)
    self.last_round_labels_command = self.replaceDynamicVariables(self.labels_command)

    if self.options.continue_post_scoring:
      if self.continuing and self.max_round != None and self.current_round == self.max_round:
        self.continuing = False
      if self.continuing:
        if len(self.old_feature_counts) < self.current_round:
          self.continuing = False
        else:
          self.feature_counts[-1] = self.old_feature_counts[self.current_round - 1]

    # perform selection for next round
    full_command = self.replaceDynamicVariables(self.select_command)
    if self.options.dry_run:
      print "Would execute: " + str(full_command)
    else:
      if not self.continuing:
        print "Selecting Features for round " + str(self.current_round)
        (status, output) = bcl_cv.tryExecute(full_command, 2, "Selecting features", 10)
        if not status:
          print "Error while selecting features"
          sys.exit(-1)
    return

  def replaceDynamicVariables(self, strn):
    strx = strn.replace(self.RND_STR, 'rnd' + str(self.current_round)) \
               .replace(self.N_FEATURES_STR, str(self.feature_counts[-1])) \
               .replace(self.LASTRND_STR, 'rnd' + str(self.current_round - 1))
    if strx.find(self.INPUT_SENS_FEATURES) >= 0:
      strx = strx.replace(self.INPUT_SENS_FEATURES, str(int(float(self.feature_counts[ 0]) / float(self.feature_counts[-2]) * self.options.is_feature_count)))
    return strx

  def writeResults(self):
    smaller_is_better = self.improvement_type == bcl_cv.BclCVStatus.SMALLER_IS_BETTER
    sorted_results = sorted(self.round_results.iteritems(), key = operator.itemgetter(1), reverse = smaller_is_better)
    final_resultf = open(self.final_output_file, 'w')
    final_resultf.write('Round\t# Descriptors\tConsensus Result\tRaw Result')
    base_results = [str(x[0]) + '\t' + str(self.feature_counts[x[0]]) + '\t' + str(x[1]) + '\t' + str(self.round_results_raw[x[0]]) for x in sorted_results]
    if self.has_blind_set:
      final_resultf.write('\tBlind Result\tVisible Result\n')
      base_results = [ base_results[x] + '\t' + str(self.blind_results[sorted_results[x][0]]) + '\t' + str(self.visible_results[sorted_results[x][0]]) for x in range(len(base_results)) ]
    else:
      final_resultf.write('\n')
    final_resultf.write('\n'.join(base_results))
    final_resultf.close()

  def run(self):
    while self.feature_counts[-1] >= self.options.min_features and self.current_round < self.options.ds_rounds:
      print "Iteration ", self.current_round
      self.performIteration()
      self.writeResults()
    print "Stopping " + ' '.join([str(x) for x in [self.feature_counts[-1], self.options.min_features, self.current_round, self.options.ds_rounds]])

    smaller_is_better = self.improvement_type == bcl_cv.BclCVStatus.SMALLER_IS_BETTER
    sorted_results = sorted(self.round_results.iteritems(), key = operator.itemgetter(1), reverse = smaller_is_better)

    print 'results', self.round_results
    print 'sorted res:', sorted_results

    self.writeResults()

    print 'Best result: Round: ' + str(sorted_results[-1][0]) + " Value: " + str(sorted_results[-1][1])
    print 'Combine final results .. done.'

def main():
  '''
  '''
  parser = BclModelDependentDescriptorSelection.getParser()
  (options, remaining_args) = parser.parse_known_args()
  commander = BclModelDependentDescriptorSelection(options, remaining_args)
  outputf = open(commander.cross_validation_runner.log_files_path + os.sep + 'command.txt', 'w')
  outputf.write(' '.join(sys.argv))
  outputf.close()
  commander.run()

if __name__ == '__main__':
  main()
