#!/usr/bin/env python2.7
'''
Script to perform a single cross-validation run and collect results using bcl model:Train, model:PredictionMerge, and model:ComputeJuryStatistics
Created April 2013
@author: mendenjl

'''

import os
import sys
import os.path
import multiprocessing, commands
import argparse, datetime
import uuid, hashlib
from curses.ascii import isdigit
from time import gmtime, strftime, sleep
import shutil, threading, Queue, random
from util.CrossValidationStatus import *
from util.Utils import rm_rf, tryExecute
import util.ThreadLocalRun 

# class that handles running a single cross validation; everything from memory estimation, PBS submission, merging results
# and evaluating the objective function on the final merged independent set goes here
class BclCommandCreator:

  #statics
  ind_chunk_tmpl = '$$$INDEPENDENT$$$'
  mon_chunk_tmpl = '$$$MONITORING$$$'
  key_chunk_tmpl = '$$$KEY$$$'
  repeat_chunk_tmpl = '$$$REPEAT$$$'
  job_file_tmpl = 'independent' + ind_chunk_tmpl + '_monitoring' + mon_chunk_tmpl + '_number' + repeat_chunk_tmpl
  final_indmerged_filename = 'final_objective.ind_merged.txt'

  # map from iterate primary type to memory multiplier type
  iterate_to_memory_multiplier = {}
  iterate_to_memory_multiplier['DecisionTree'.lower()] = 1.1
  iterate_to_memory_multiplier['KappaNearestNeighbor'.lower()] = 2.2
  iterate_to_memory_multiplier['Kohonen'.lower()] = 1.1
  iterate_to_memory_multiplier['LinearRegression'.lower()] = 2.2
  iterate_to_memory_multiplier['NeuralNetwork'.lower()] = 1.1
  iterate_to_memory_multiplier['OpenCLResilientPropagation'.lower()] = 2.2
  iterate_to_memory_multiplier['OpenCLSimplePropagation'.lower()] = 2.2
  iterate_to_memory_multiplier['OpenclIterateSequentialMinimalOptimization'.lower()] = 3.2
  iterate_to_memory_multiplier['SupportVectorMachine'.lower()] = 2.2

  # map from gpu selection choice to actual bcl flag
  gpu_map = {}
  gpu_map['Intel'.lower()] = '-opencl "Intel(R)_OpenCL" TYPE_CPU '
  gpu_map['ATI'.lower()] = '-opencl ATI_Stream TYPE_GPU '
  gpu_map['NVIDIA'.lower()] = '-opencl NVIDIA_CUDA TYPE_GPU '
  gpu_map['AMD'.lower()] = '-opencl AMD_Accelerated_Parallel_Processing TYPE_GPU '
  gpu_map['Disable'.lower()] = '-opencl Disable'
  gpu_map['None'.lower()] = ''

  @staticmethod
  def getParser():
    parser = argparse.ArgumentParser\
    (\
      epilog = "Additional flags with arguments will be passed on to the bcl unchanged\n" \
    )

    parser.add_argument('--bcl', help = "BCL executable to use", default = '/blue/meilerlab/apps/Linux2/x86_64/bin/bcl.exe', dest = 'bcl')

    training_group = parser.add_argument_group('training', 'These options directly affect training that is performed')
    dataset_group = parser.add_argument_group('dataset', "These options affect which dataset file(s) and which feature/result columns from those datasets are used for training")
    output_group = parser.add_argument_group('output', "These options affect what data is output and where")
    resources_group = parser.add_argument_group('resources', "These options determine which computational resources are used")
    cv_group = parser.add_argument_group('cross-validation', "These options determine the splitting used for cross-validation")

    training_group.add_argument\
    (\
       "-l" , "--learning-method", \
       help = "BCL machine learning method", \
       dest = 'iterate', \
       required = True \
    )
    dataset_group.add_argument\
    (\
       "-d", "--datasets", \
       help = "dataset files (must be .bin files)", \
       dest = 'datasets', \
       nargs = '+', \
       metavar = '<dataset>', \
       required = True \
    )
    dataset_group.add_argument\
    (\
       "--training-size-limit", \
       help = "limit the training dataset size to this many rows", \
       dest = 'train_size_limit', \
       type = int
    )
    dataset_group.add_argument\
    (\
       "--monitoring-size-limit", \
       help = "limit the monitoring dataset size to this many rows", \
       dest = 'mon_size_limit', \
       type = int
    )
    dataset_group.add_argument\
    (\
       "--include-monitoring-in-training", \
       help = "set to true to include the monitoring dataset in the training", \
       dest = 'include_m_in_t', \
       action = 'store_true', \
       default = False \
    )
    dataset_group.add_argument\
    (\
       "--swap-training-and-monitoring", \
       help = "set to true to include one chunk of the training dataset in the monitoring", \
       dest = 'swap_t_and_m', \
       action = 'store_true', \
       default = False \
    )
    dataset_group.add_argument\
    (\
       "--no-cross-validation", \
       help = "Set to true to merge the training, monitoring, and independent set. " \
              + " This can be useful to produce a final set of models on all the data", \
       dest = 'no_cv', \
       action = 'store_true', \
       default = False \
    )
    dataset_group.add_argument\
    (\
       "--watch-blind", \
       help = "Monitoring and independent evaluation will take place on the blind set. Implies -result_averaging_window 0.", \
       dest = 'watch_blind', \
       action = 'store_true', \
       default = False \
    )
    # general cross-validation options
    cv_group.add_argument\
    (\
       "-m", "--monitoring-id-range", \
       help = "Min Max value for monitoring id range", \
       dest = 'mon_id_range', \
       nargs = 2, \
       metavar = ('<min>', '<max>'), \
       default = [0, 4], \
       type = int \
    )
    cv_group.add_argument\
    (\
       "-i", "--independent-id-range", \
       help = "Min Max value for independent id range", \
       dest = 'ind_id_range', \
       nargs = 2, \
       metavar = ('<min>', '<max>'), \
       default = [0, 4], \
       type = int \
    )
    cv_group.add_argument\
    (\
      "--monitor-independent-set", \
      help = "If set, monitoring and independent datasets will be one in the same",
      dest = 'mon_ind_same',
      action = 'store_true' \
    )
    cv_group.add_argument\
    (\
       "-n", "--cross-validations", \
       help = "Number of cross-validation chunks", \
       dest = 'chunks', \
       metavar = '<#CV>', \
       default = 5, \
       type = int \
    )
    cv_group.add_argument \
    (\
      "--blind", \
      help = "Optional partition/chunk index of the blind chunk, which is never evaluated by the model. " \
             + " The first parameter is the zero-indexed blind index, the second is the number of partitions for the blind. " \
             + " E.g. to ignore the second third of the dataset completely during this cross-validation, --blind 1 3", \
      nargs = 2, \
      dest = 'blind_id_nchunks', \
      type = int \
    )
    cv_group.add_argument \
    (\
      "--blind-descriptor", \
      help = "ID descriptor used to select the blind set by dataset IDs, e.g. in a protein dataset, PdbID can be used" \
             + "to select individual proteins for the blind set", \
      dest = 'blind_descriptor', \
    )
    cv_group.add_argument \
    (\
      "--blind-ids", \
      help = "Optional ID(s) to exclude , which is never evaluated by the model. " \
             + " The first parameter is the zero-indexed blind id, the second is the number of partitions for the blind. " \
             + " E.g. to ignore the second third of the dataset completely during this cross-validation, --blind 1 3", \
      dest = 'blind_descriptor_ids', \
      nargs = '+'
    )
    cv_group.add_argument \
    (\
      "--blind-features", \
      help = "Feature descriptor (s) used to select the blind set by feature range", \
      dest = 'blind_features', \
      nargs = '+' 
    )
    cv_group.add_argument \
    (\
      "--blind-feature-ranges", \
      help = "Optional ID(s) to exclude , which is never evaluated by the model. " \
             + " The first parameter is the zero-indexed blind id, the second is the number of partitions for the blind. " \
             + " E.g. to ignore the second third of the dataset completely during this cross-validation, --blind 1 3", \
      dest = 'blind_feature_ranges', \
      nargs = '+'
    )
    cv_group.add_argument \
    (\
      "--cv-repeats", \
      help = "Number of models to train on each cross-validation set", \
      dest = 'cv_repeats', \
      type = int, \
      default = 1 \
    )
    cv_group.add_argument \
    (\
     "--complete", \
     help = "Run any failed jobs. Use this if the cross-validation is stopped by jobs that failed to avoid redoing the " \
            + "whole cross-validation. Using this option when the cross-validation is still running is not recommended " \
            + "and will result in termination of the cross-validation.", \
     dest = 'complete', \
     action = 'store_true' \
    )
    cv_group.add_argument \
    (\
     "--no-flock-submit", \
     help = "Set this flag to not use flock to prevent multiple instances of this script from launching jobs simultaneously. " \
            " flock is problematic when the file system is behaving poorly on a cluster head node. Note that if the " \
            " cluster is generally misbehaving, it will be impossible to detect when all jobs are done, preventing " \
            "results from being calculated. Likewise, only use this option if only the head node is experiencing issues with file locking", \
     dest = 'no_flock', \
     action = 'store_true' \
    )
    cv_group.add_argument\
    (\
      "--flock-max-wait", \
      help = "sets a timeout for waiting for flock for updating the status file. Updating the status file " \
      "takes a few microseconds, so ordinarily these locks are instantaneously acquired. On NFS shares, " \
      " file locking doesn't always work reliably, so the flock call sometimes doesn't terminate without a timeout, but " \
      " beware that this may lead to status file corruption in exceptionally rare circumstances where two processes then " \
      " edit the status file simultaneously.",
      dest = "flock_max_wait", \
      type = int, \
      default = 30 \
    )
    balanced_combined = dataset_group.add_mutually_exclusive_group()
    balanced_combined.add_argument\
    (\
       "--balanced", \
       help = "Set to balance input from all datasets", \
       dest = 'balanced', \
       action = 'store_true' \
    )
    balanced_combined.add_argument\
    (\
       "--combined", \
       help = "Set to combine (not balance) input from all datasets.  This is the default behavior", \
       dest = 'balanced', \
       action = 'store_false', \
       default = False \
    )
    dataset_group.add_argument\
    (\
      "--features", \
      help = "File (or string) containing feature descriptors for use in training", \
      dest = "features"\
    )
    dataset_group.add_argument\
    (\
      "--score-file", \
      help = "(strictly optional) features score file, must be used with --select-top-features flag", \
      dest = "score_file"\
    )
    dataset_group.add_argument\
    (\
      "--top-features", \
      help = "chooses the highest scoring <N> features from the --score-file to use in place of --features", \
      dest = "top_n_features", \
      type = int \
    )
    dataset_group.add_argument\
    (\
      "--results", \
      help = "File (or string) containing result descriptors for use in training", \
      dest = "results"\
    )
    dataset_group.add_argument\
    (\
      "--filter-by-descriptor", \
      help = "Descriptor (which must be in the features file) to filter the training data set. " \
             "Specify desired range with --filter-by-descriptor-range", \
      dest = "filter_descriptor", \
      nargs = '+', \
      metavar = '<filters>' \
    )
    dataset_group.add_argument\
    (\
      "--filter-by-descriptor-range", \
      help = "Numeric range for the descriptor specified by --filter-by-descriptor, syntax is like [0,5)+[7,9]", \
      dest = "filter_descriptor_range", \
      nargs = '+', \
      metavar = '<filters>' \
    )
    dataset_group.add_argument\
    (\
      "--select-id-descriptors", \
      help = "ID descriptor (which must be in the features file) to filter the training data set. " \
             "Specify desired values as comma separated lists in --select-ids", \
      dest = "filter_id_descriptor", \
      nargs = '+', \
      metavar = '<filters>' \
    )
    dataset_group.add_argument\
    (\
      "--select-ids", \
      help = "Comma-separated strings; to select for training range for the descriptor specified by --filter-by-ids, " \
      "syntax is like \"App,Ban,Chr\" \"0,1\"", \
      dest = "filter_ids", \
      nargs = '+', \
      metavar = '<filters>' \
    )
    dataset_group.add_argument\
    (\
      "--remove-id-descriptors", \
      help = "ID descriptor (which must be in the features file) to exclude from the training data set. " \
             "Specify desired values as comma separated lists in --remove-ids", \
      dest = "remove_id_descriptor", \
      nargs = '+', \
      metavar = '<filters>' \
    )
    dataset_group.add_argument\
    (\
      "--remove-ids", \
      help = "Comma-separated strings; to select for training range for the descriptor specified by --filter-by-ids, " \
      "syntax is like \"App,Ban,Chr\" \"0,1\"", \
      dest = "remove_ids", \
      nargs = '+', \
      metavar = '<filters>' \
    )
    dataset_group.add_argument\
    (\
      "--filter-by-result", \
      help = "Descriptor (which must be in the results file) to filter the training data set. " \
             "Specify desired range with --filter-by-result-range", \
      dest = "filter_result", \
      nargs = '+', \
      metavar = '<filters>' \
    )
    dataset_group.add_argument\
    (\
      "--filter-by-result-range", \
      help = "Numeric range for the descriptor specified by --filter-by-result, syntax is like [0,5)+[7,9]", \
      dest = "filter_result_range", \
      nargs = '+', \
      metavar = '<filters>' \
    )
    dataset_group.add_argument\
    (\
      "--filter-all", \
      help = "If any filters are specified, they are applied to training, monitoring, and independent sets (by default, they are applied only to training)", \
      dest = "filter_all", \
      default = False, \
      action = 'store_true' \
    )
    dataset_group.add_argument\
    (\
      "--filter-evaluation", \
      help = "If any filters are specified, they are applied to monitoring, and independent sets only (by default, they are applied only to training)", \
      dest = "filter_eval", \
      default = False, \
      action = 'store_true' \
    )
    output_group.add_argument\
    (\
       "--final-objective-function", \
       help = "objective function to call and merge on independent datasets", \
       dest = 'final_obj', \
       metavar = '<objective-function>' \
    )
    output_group.add_argument\
    (\
       "--id", \
       help = "identifier; models will be placed in models/{id}/{id}_i{independent}_m{monitoring}.model," + \
              " log files in log_files/{id}/{id}_i{independent}_m{monitoring}.log.txt", \
       dest = 'name', \
       default = strftime("%Y_%m_%d_%H_%M_%S", gmtime()), \
       metavar = '<name>' \
    )
    output_group.add_argument\
    (\
       "--print-independent-predictions", \
       help = "If set, merges all independent predictions into results/{id}/{id}.txt.gz " + \
       "and writes independent predictions from each training run are saved in results/{id}/{id}_i{independent}_m{monitoring}.txt.gz", \
       dest = 'print_ind', \
       default = False, \
       action = 'store_true' \
    )
    output_group.add_argument\
    (\
       "--print-training-predictions", \
       help = "If set, merges all training predictions into results/{id}/{id}.txt.gz " + \
       "and writes training predictions from each training run are saved in results/{id}/{id}_i{independent}_m{monitoring}.txt.gz", \
       dest = 'print_train', \
       default = False, \
       action = 'store_true' \
    )
    output_group.add_argument\
    (\
       "--use-median", \
       help = "Used only if --print-independent-predictions was given.  Causes median prediction to be computed for each CV, rather than the average", \
       dest = 'use_median', \
       default = False, \
       action = 'store_true' \
    )
    output_group.add_argument\
    (\
       "--use-jury", \
       help = "Used only if --print-independent-predictions was given.  Causes jury prediction to be computed for each CV, rather than the average", \
       dest = 'use_jury', \
       default = False, \
       action = 'store_true' \
    )
    output_group.add_argument\
    (\
       "--use-min", \
       help = "Used only if --print-independent-predictions was given.  Causes min prediction to be computed for each CV, rather than the average", \
       dest = 'use_min', \
       default = False, \
       action = 'store_true' \
    )
    output_group.add_argument\
    (\
       "--use-max", \
       help = "Used only if --print-independent-predictions was given.  Causes max prediction to be computed for each CV, rather than the average", \
       dest = 'use_max', \
       default = False, \
       action = 'store_true' \
    )
    training_group.add_argument\
    (\
       "--max-iterations", \
       help = "# of steps to train the learning method ", \
       dest = 'iterations', \
       default = 100, \
       type = int \
    )
    training_group.add_argument\
    (\
       "--max-minutes", \
       help = "Maximum # of minutes to allow for training." + \
       "Walltime will be this time + time to allow for dataset loading + some allowance", \
       dest = 'minutes', \
       default = 120, \
       type = int\
    )
    output_group.add_argument\
    (\
       "--store-model", \
       help = "Where to store model (by default, no model is stored)", \
       dest = 'model_storage', \
       choices = ['Db', 'File'] \
    )
    output_group.add_argument\
    (\
       "--store-metadata", \
       help = "Where to store metadata (objective function, descriptors, etc.) (by default, same location as --store-model, if given)", \
       dest = 'meta_storage', \
       choices = ['Db', 'File'] \
    )

    pbs_local = resources_group.add_mutually_exclusive_group()
    ncpus = multiprocessing.cpu_count()
    pbs_local.add_argument\
    (\
       "--local", \
       help = "Set to run jobs locally", \
       dest = 'local_jobs', \
       type = int, \
       metavar = '<# of jobs to run simultaneously, default= ' + str(ncpus) + '>', \
       nargs = '?', \
       default = None, \
       const = ncpus \
    )
    pbs_local.add_argument\
    (\
       "--pbs", \
       help = "Set to run jobs over pbs with qsub", \
       dest = 'pbs_jobs', \
       metavar = ('<extra resource args to qsub>'), \
       nargs = '?', \
       # the spaces is important for const, otherwise checking options.pbs_jobs 
       # returns false because python evaluates empty strings as false! 
       const = ' ', \
       default = False \
    )
    pbs_local.add_argument\
    (\
       "--slurm", \
       help = "Set to run jobs over slurm with sbatch", \
       dest = 'slurm_jobs', \
       metavar = ('<extra resource args to qsub>'), \
       nargs = '?', \
       # the spaces is important for const, otherwise checking options.slurm_jobs 
       # returns false because python evaluates empty strings as false! 
       const = ' ', \
       default = False \
    )
    pbs_local.add_argument\
    (\
       "--gparallel", \
       help = "Set to run jobs over gnu parallels.  Note that host is ignored when using this option", \
       dest = 'gparallel', \
       metavar = ('<arguments to GNU parallels>'), \
       nargs = '?', \
       # the spaces is important for const, otherwise checking options.pbs_jobs 
       # returns false because python evaluates empty strings as false! 
       const = ' ', \
       default = False \
    )
    resources_group.add_argument\
    (\
       "--pbs-resources", \
       help = "additional resources to request explicitly when requesting ppn and nodes with -l flag", \
       dest = 'pbs_resources', \
       metavar = ('<extra resource args to qsub>'), \
       nargs = '?', \
       # the spaces is important for const, otherwise checking options.pbs_jobs 
       # returns false because python evaluates empty strings as false! 
       const = ' ', \
       default = False \
    )
    resources_group.add_argument\
    (\
       "--slurm-resources", \
       help = "additional resources to request explicitly when from slurm using the gres flag", \
       dest = 'slurm_resources', \
       metavar = ('<extra resource args to slurm>'), \
       nargs = '?', \
       # the spaces is important for const, otherwise checking options.pbs_jobs 
       # returns false because python evaluates empty strings as false! 
       const = ' ', \
       default = False \
    )
    resources_group.add_argument\
    (\
       "--bundle-jobs", \
       help = "(PBS only) If specified, combine multiple jobs to improve job priority on pbs clusters that are not fond of serial jobs. " \
              + "Currently this is implemented by distributing the jobs to different nodes, so do not set this number higher than there are "
              + "typically nodes on the cluster", \
       dest = 'job_bundle_size', \
       default = 1, \
       type = int \
    )
    resources_group.add_argument\
    (\
       "--host", \
       help = "The host to login to for job submission", \
       dest = 'host', \
       metavar = ('<hostname>'), \
       default = os.getenv('HOSTNAME') \
    )
    resources_group.add_argument\
    (\
       "--just-submit", \
       help = "Just submit jobs and return", \
       dest = 'just_submit', \
       action = 'store_true', \
       default = False \
    )
    resources_group.add_argument\
    (\
       "--max-reruns", \
       help = "Maximum # of times to attempt re-running a job if it fails", \
       dest = 'max_reruns', \
       type = int, \
       metavar = '<RERUNS>', \
       default = 1 \
    )
    resources_group.add_argument\
    (\
       "--opencl", \
       help = "Opencl driver to use; set to None if you receive errors like opencl no such flag available", \
       dest = 'opencl_driver', \
       default = 'Disable', \
       choices = [ 'Intel', 'ATI', 'NVIDIA', 'AMD', 'Disable', 'None'] \
    )
    resources_group.add_argument\
    (\
       "--opencl-gpu-id", \
       help = "Opencl gpu id; implies TYPE_GPU; but opencl flag still must be set", \
       dest = 'opencl_gpu_id', \
       type = int, \
       default = -1, \
       metavar = '<gpu_id>' \
    )
    resources_group.add_argument\
    (\
       "--dry-run", \
       help = "Set to not actually run any training commands, just output the files", \
       dest = 'dry_run', \
       default = False, \
       action = 'store_true' \
    )
    resources_group.add_argument\
    (\
       "--no-write-scripts", \
       help = "Implies --dry-run, nothing is run and no scripts are written", \
       dest = 'no_write_scripts', \
       default = False, \
       action = 'store_true' \
    )
    resources_group.add_argument\
    (\
       "--max-pbs-jobs", \
       help = "maximum number of jobs to allow for running on pbs at a time", \
       dest = 'max_pbs_jobs', \
       default = -1, \
       type = int \
    )
    resources_group.add_argument\
    (\
       "--max-slurm-jobs", \
       help = "maximum number of jobs to allow for running on SLURM at a time", \
       dest = 'max_slurm_jobs', \
       default = -1, \
       type = int \
    )
    resources_group.add_argument\
    (\
       "--slurm-partition", \
       help = "If given, will be passed to the --partition flag for slurm", \
       dest = 'slurm_partition', \
       default = ''
    )
    resources_group.add_argument\
    (\
       "--memory-offset", \
       help = "Always add this much memory (in units of MB) to job requirements; this is necessary when using a relatively large blind dataset", \
       dest = 'memory_offset', \
       default = 0, \
       type = int \
    )
    resources_group.add_argument\
    (\
       "--no-failure-propagation", \
       help = "Normally if any job fails during cross validation, all remaining jobs are stopped. Setting this flag prevents that from happening", \
       dest = 'no_fail_prop', \
       action = 'store_true', \
       default = False \
    )

    output_group.add_argument\
    (\
       "--round", \
       help = "Used to control stored round #; only used if --store-model was given", \
       dest = 'round', \
       type = int, \
       default = 0 \
    )
    resources_group.add_argument\
    (\
       "--pthreads", \
       help = "# of pthreads to run the bcl with", \
       dest = 'pthreads', \
       type = int \
    )
    output_group.add_argument\
    (\
       "--override-memory-multiplier", \
       help = "If you're getting out-of-memory errors on the cluster, set this to something larger than the default for the current iterate", \
       dest = 'memory_multiplier', \
       type = float \
    )
    output_group.add_argument\
    (\
       "--remove-files-on-success", \
       help = "Setting this flag causes removal of all log, job, and non-merged independent predictions upon successful CV completion", \
       dest = 'remove_files', \
       action = 'store_true', \
       default = False \
    )
    output_group.add_argument\
    (\
       "--show-status", \
       help = "Set this flag to display status lines whenever new jobs are started running or complete", \
       dest = 'show_status', \
       action = 'store_true', \
       default = False \
    )
    output_group.add_argument\
    (\
       "--no-plot", \
       help = "If set, no plots will be generated", \
       dest = 'no_plot', \
       action = 'store_true', \
       default = False \
    )
    output_group.add_argument\
    (\
       "--predict-with-best-models", \
       help = "Set this flag to use only the best model (based on final-objective-function from each independent chunk " \
              "for computing results on blind and visible datasets.  All models will still be used when computing the " \
              "consensus result to avoid biasing.  This option is currently only valid if --store-model is set to File ", \
       dest = 'only_best_models', \
       action = 'store_true', \
       default = False \
    )

    return parser

  # method to retrieve what the final result output would be given the name of the cv
  @staticmethod
  def CreateStatus(options):
    log_files_path = './log_files' + os.sep + options.name
    scripts_path = log_files_path + os.sep + 'scripts'
    status_filename = log_files_path + os.sep + 'status.txt'
    sha1obj=hashlib.sha1()
    sha1obj.update(status_filename)
    status_lockfilename = ( "/tmp/bcl_cv_" + sha1obj.hexdigest() if options.local_jobs else status_filename) + '.lock'
    
    files_to_remove = [ scripts_path + os.sep + 'wakeup', status_lockfilename]
    if options.remove_files:
      files_to_remove = [ log_files_path ]
    results_file = ''
    if options.print_ind:
      results_file = 'results' + os.sep + options.name + os.sep + BclCommandCreator.final_indmerged_filename
    expect_blind_dataset = False
    if options.blind_id_nchunks and len(options.blind_id_nchunks) == 2:
      expect_blind_dataset = True
    return BclCVStatus(options.show_status, False, status_filename, files_to_remove, results_file, '', expect_blind_dataset)

  # python equivalent to removing a directory if it exists, then calling mkdir -p on it, e.g. rm -rf x/y; mkdir -p x/y
  @staticmethod
  def mkdirmp(directory_name):
    # remove files from any existing directory
    if os.path.exists(directory_name):
      rm_rf(directory_name)
    paths = directory_name.split(os.sep)
    paths = [x for x in paths if len(x) > 0 and x != '.']
    current_path_component = ''
    for p in paths:
      current_path_component += p + os.sep
      if not os.path.exists(current_path_component):
        os.mkdir(current_path_component)
    return directory_name
  
  def mkdirmpWithSymlinks(self,directory_name):
    if not self.options.complete or not os.path.exists(directory_name):
      self.mkdirmp(directory_name)
    nm = directory_name
    paths = [x for x in nm.split(os.sep) if len(x) > 0 and x != '.']
    if len(paths) == 2:
      primary_path = paths[0]
      abs_path = os.path.abspath(directory_name)
      for p in self.output_directories:
        dest = directory_name + os.sep + p.rstrip('/').split(os.sep)[-2]
        try:
          os.symlink(p, dest)
        except:
          pass
        try:
          os.symlink(abs_path, p + os.sep + primary_path)
        except:
          pass
      self.output_directories.append(os.path.abspath(directory_name))
    return directory_name

  # replaces all occurrences of job-dependent variables like ind_chunk_templ with the given integer
  def replaceDynamicVariables(self, strn, ind, mon, key, cv_repeat_id):
    return strn.replace(self.ind_chunk_tmpl, str(ind))\
               .replace(self.mon_chunk_tmpl, str(mon))\
               .replace(self.key_chunk_tmpl, str(key))\
               .replace(self.repeat_chunk_tmpl, str(cv_repeat_id))

  # write out shell commands used during job running and interaction with queue and epilogue scripts
  def writeCommonCommands(self):
    outputf = open(self.proepi_filename, 'w')

    outputf.write( """#!/bin/bash
export BCL_MAX_RERUNS="""+str(self.options.max_reruns) + """
export BCL_CV_STATUS_FILENAME="""+self.status_filename + """
export BCL_CV_STATUS_LOCK_FILENAME="""+self.status_lockfilename + """
export BCL_CV_WAKEUP_FILENAME="""+self.wakeup_file + """
export BCL_CV_LOGFILES_PATH="""+self.log_files_path+ """
export PATH=""" + os.environ['PATH'] + """
export LD_LIBRARY_PATH=""" + os.environ['LD_LIBRARY_PATH'] + """
export PBS_O_WORKDIR=""" + os.getcwd() + """
export BCL_CV_SCRIPTS_PATH="""+self.scripts_path + """
export FLOCK_MAX_WAIT=""" + str(self.options.flock_max_wait) + "\n" + 
( 'export BCL_MERGE_SCRIPT=' + self.merge_script + '\n' if self.merge_script else ""))
    if self.options.local_jobs and len(self.job_keys_to_run) == 1:
      outputf.write('alias flock=echo\n')
    script_dirname = os.path.dirname(os.path.abspath(__file__))
    if script_dirname.startswith('/ssd') or script_dirname.startswith('/hd'):
      hostname=os.getenv('HOST', "")
      if len(hostname) == 0:
        hostname=os.getenv('HOSTNAME', "")
        if hostname.find('.') > 0:
          hostname = hostname[:hostname.find('.')]
      script_dirname = '/net/' + hostname + script_dirname
    outputf.write(". " + script_dirname + "/util/functions.sh\n")
    outputf.close()
    os.chmod(self.proepi_filename, 0755)

  def checkLogFileFinished(self, log_file):
    if not os.path.exists(log_file):
      return False
    ifile = open(log_file,'r')
    lines = [x.strip().lower() for x in ifile.readlines()]
    ifile.close()
    if len(lines) == 0 or len([x for x in lines if x.find('bcl has run for') >= 0]) == 0:
      return False
    test1 = len([x for x in lines if x.find('independent predictions written') >= 0]) > 0 if self.options.print_ind else True
    test2 = len([x for x in lines if x.find('final model written') >= 0]) > 0 if self.options.model_storage else True
    test3 = len([x for x in lines if x.find('independent data: final score') >= 0]) > 0 if self.options.final_obj else True
    return test1 and test2 and test3
  
  def writePbsScripts(self):
    local_cwd = os.getcwd()
    effective_host = self.options.host
    short_host_name = effective_host
    if self.options.host and len(self.options.host):
      effective_host = self.options.host

    ppos = effective_host.find('.')
    if effective_host.find('.') >= 0:
      short_host_name = effective_host[:ppos]

    primary_lock_file = '~/.BclModelCrossValidation.' + effective_host + '.lock'

    outputf = open(self.abort_script, 'w')
    outputf.write('#!/bin/bash\n')
    outputf.write('if [ "$HOSTNAME" != "' + effective_host + '" ] && [ `echo $HOSTNAME | sed \'s/\..*//g\'` != "' + effective_host + '" ]; then\n')
    outputf.write('ssh ' + short_host_name + ' ' + os.path.abspath(self.abort_script) + '\n')
    outputf.write('exit_status=$?\n')
    if effective_host != short_host_name:
      outputf.write('if [[ $exit_status -ne 0 ]]; then\n')
      outputf.write('ssh ' + effective_host + ' ' + os.path.abspath(self.abort_script) + '\nfi\n')
    outputf.write('exit\nfi\n\n')
    outputf.close()
    os.chmod(self.abort_script, 0755)
    self.all_submission_commands.append("echo 'kill -9 '$$' 2>&1 >& /dev/null' >> " + self.abort_script + '\n')
    self.n_pbs_jobs = len(self.job_keys_to_run)
    if len(self.job_keys_to_run):
      self.all_submission_commands.append("touch " + primary_lock_file + '\n')
      self.all_submission_commands.append("sleep 1\n")
      self.all_submission_commands.append('(\n')
      if not self.options.no_flock:
        self.all_submission_commands.append('  flock -x -w 60 300 \n')
        self.all_submission_commands.append('  [[ $? -eq 0 ]] || echo "Failed to lock ' + primary_lock_file + \
                                            ' after 60 seconds. Proceeding anyway, assuming file locking is broken. ' \
                                            'If multiple submission scripts are running, you will need to verify that ' \
                                            'all jobs have been accepted by the scheduler. Set --no-flock-submit if you ' \
                                            'are sure this is the only process submitting jobs for you on this cluster gateway" \n')
      bundle_size = self.options.job_bundle_size
      bundles = []
      bundles_filename=""
      outputb = None
      if bundle_size > 1:
        bundles_filename = self.scripts_path + os.sep + "bundles.txt"
        outputb = open(bundles_filename,'w')
      for (ind,mon,key,repeat,realid) in self.job_keys_to_run:
        script_name = self.scripts_path + os.sep + self.options.name + '.' + str(realid) + '.sh'
        script_name = script_name.replace('/./', '/')
        ljob_file_path = self.replaceDynamicVariables(self.job_file_path, ind, mon, realid, repeat) if bundle_size == 1 else self.job_file_path.replace('$$$BUNDLE$$$',str(len(bundles)-1))
        outputf = open(script_name, 'w')
        outputf.write('#!/bin/bash\n')
        if bundle_size == 1:
          outputf.write('echo "Starting job ' + str(realid) + ' on "`cat $PBS_NODEFILE`\n')
        else:
          outputf.write('export PATH=' + os.environ['PATH'] + '\n')
          outputf.write('export LD_LIBRARY_PATH=''' + os.environ['LD_LIBRARY_PATH'] + '\n')
          outputf.write('cd ' + local_cwd + '\n')

        outputf.write(". " + self.proepi_filename + '\n')
        outputf.write('PrologueStatus ' + self.status_filename + '\n')
        outputf.write('run_command_check_status ' + self.const_args + ' ')
        outputf.write(self.replaceDynamicVariables(self.variable_args, ind, mon, realid, repeat) + '\n')
        if bundle_size > 1:
          if len(bundles) == 0 or len(bundles[-1]) == bundle_size:
            bundles.append([])
          bundles[-1].append(realid)
          outputb.write(script_name + '\n')
          outputf.write('DecrementBundleCount ''' + str(len(bundles)-1) + '\n')
        outputf.write('exit 0')
        outputf.close()
        os.chmod(script_name, 0755)
        if bundle_size == 1:
          var_sub = self.replaceDynamicVariables(self.variable_submission_cmd, ind, mon, realid, repeat)
          epilogue_script_name = script_name + '.epilogue.sh'
          base_sub = self.base_submission_cmd + ',epilogue=' + epilogue_script_name + ' '
          outputf = open(epilogue_script_name, 'w')
          outputf.write(\
           "#!/bin/bash\ncd " + local_cwd + "\n. \"" + self.proepi_filename + '\"\nPbsEpilogue ' \
           + ljob_file_path + ' $@\nexit 0')
          outputf.close()
          os.chmod(epilogue_script_name, 0755)
          var_sub += " " + script_name + " "
    
          # touch the output file; necessary due to NFS syncing
          self.all_submission_commands.append("touch  " + self.replaceDynamicVariables(self.job_file_path, ind, mon, realid, repeat))
          self.all_submission_commands.append(base_sub + var_sub + " | sed 's/^\([0-9]*\).*/qdel \\1\\n/' >> " + self.abort_script)
          self.all_submission_commands.append('sleep 1')
  
      if bundle_size > 1:
        outputb.close()
        self.n_pbs_jobs = len(bundles)
        bundles_script=self.scripts_path + os.sep + "run_bundle.sh"
        outputb = open( bundles_script,'w')
        outputb.write('#!/bin/bash\n. "' + self.proepi_filename + '"\nRunPBSBundle ' + str(bundle_size) + "\nexit 0")
        outputb.close()
        os.chmod(bundles_script, 0755)
        for bundleid in range(len(bundles)):
          job_path = self.job_file_path.replace('$$$BUNDLE$$$',str(bundleid))
          var_sub = self.variable_submission_cmd.replace('$$$BUNDLE$$$',str(bundleid))
          epilogue_script_name = self.scripts_path + os.sep + os.path.basename(job_path.replace('.job.txt','.epilogue.sh'))
          base_sub = self.base_submission_cmd + ',epilogue=' + epilogue_script_name + ' '
          outputf = open(epilogue_script_name, 'w')
          outputf.write(\
           "#!/bin/bash\ncd " + local_cwd + "\n. \"" + self.proepi_filename + '\"\nPbsEpilogue ' \
           + job_path + ' $@\nexit 0')
          outputf.close()
          os.chmod(epilogue_script_name, 0755)
          var_sub += " -vBUNDLE=" + str(bundleid) + " " + bundles_script + " "
          self.all_submission_commands.append("touch  " + job_path)
          self.all_submission_commands.append(base_sub + var_sub + " | sed 's/^\([0-9]*\).*/qdel \\1\\n/' >> " + self.abort_script)
          self.all_submission_commands.append('sleep 1')
  
      self.all_submission_commands.append(') 300> ' + primary_lock_file)
    output_main = open(self.main_submission_script, 'w')
    output_main.write('#!/bin/bash\n')
    output_main.write('set -e\n')
    output_main.write('echo ' + str(self.number_jobs) + ' 0 ' + str(self.original_number_jobs-self.number_jobs) + ' 0 0 0 > ' + self.status_filename + "\n")
    output_main.write('touch ' + self.status_lockfilename + ".; sleep 1\n")
    output_main.write('\n'.join(self.all_submission_commands))
    output_main.write('\nchmod a+x ' + self.abort_script)
    output_main.close()
    os.chmod(self.main_submission_script, 0755)
    if not self.options.dry_run:
      outputf = open(self.status_filename, 'w')
      outputf.write(str(self.number_jobs) + ' 0 ' + str(self.original_number_jobs-self.number_jobs) + ' 0 0 0')
      outputf.close()

  def writeSlurmScripts(self):
    local_cwd = os.getcwd()
    effective_host = self.options.host
    short_host_name = effective_host
    if self.options.host and len(self.options.host):
      effective_host = self.options.host

    ppos = effective_host.find('.')
    if effective_host.find('.') >= 0:
      short_host_name = effective_host[:ppos]

    # the only way to force slurm to indicate the job id at submission time (as qsub always does) is with this non-sense 
    # statement, which ordinarily only works when srun is run with -vv. Mysteriously it works sometimes under -v as well
    # but that's because the job id is reported differently in the verbose output depending on whether the job could be
    # distributed at submission time or not. 
    # A better way would be to setup a one-line prologue script to print out the job id, however, the only prologue script under
    # the user's control isn't run until the beginning of job execution, which is long after we need it. The SLURM
    # slurmctld prolog is run at submission time and could output this information, however it would require a 
    # cooperative scheduler admin to add this simple fix, which we sadly do not have. 
    # TL;DR : this is fragile and will probably break whenever they mess with slurm's verbose output
    cmd_jobid_from_slurm_output_to_abort_script = " 2>&1 | egrep -m1 '^srun: (job|jobid) [0-9]' | cut -f3 -d\  | sed 's/^\([0-9]*\).*/scancel \\1\\n/' >> " + self.abort_script + " &"

    primary_lock_file = '~/.BclModelCrossValidation.' + effective_host + '.lock'

    outputf = open(self.abort_script, 'w')
    outputf.write('#!/bin/bash\n')
    outputf.write('if [ "$HOSTNAME" != "' + effective_host + '" ] && [ `echo $HOSTNAME | sed \'s/\..*//g\'` != "' + effective_host + '" ]; then\n')
    outputf.write('ssh ' + short_host_name + ' ' + os.path.abspath(self.abort_script) + '\n')
    outputf.write('exit_status=$?\n')
    if effective_host != short_host_name:
      outputf.write('if [[ $exit_status -ne 0 ]]; then\n')
      outputf.write('ssh ' + effective_host + ' ' + os.path.abspath(self.abort_script) + '\nfi\n')
    outputf.write('exit\nfi\n\n')
    outputf.close()
    os.chmod(self.abort_script, 0755)
    self.all_submission_commands.append("echo 'kill -9 '$$' 2>&1 >& /dev/null' >> " + self.abort_script + '\n')
    self.n_pbs_jobs = len(self.job_keys_to_run)
    if len(self.job_keys_to_run):
      self.all_submission_commands.append("touch " + primary_lock_file + '\n')
      self.all_submission_commands.append("sleep 1\n")
      self.all_submission_commands.append('(\n')
      if not self.options.no_flock:
        self.all_submission_commands.append('  flock -x -w 60 300 \n')
        self.all_submission_commands.append('  [[ $? -eq 0 ]] || echo "Failed to lock ' + primary_lock_file + \
                                            ' after 60 seconds. Proceeding anyway, assuming file locking is broken. ' \
                                            'If multiple submission scripts are running, you will need to verify that ' \
                                            'all jobs have been accepted by the scheduler. Set --no-flock-submit if you ' \
                                            'are sure this is the only process submitting jobs for you on this cluster gateway" \n')
      bundle_size = self.options.job_bundle_size
      bundles = []
      bundles_filename=""
      outputb = None
      if bundle_size > 1:
        bundles_filename = self.scripts_path + os.sep + "bundles.txt"
        outputb = open(bundles_filename,'w')
      for (ind,mon,key,repeat,realid) in self.job_keys_to_run:
        script_name = self.scripts_path + os.sep + self.options.name + '.' + str(realid) + '.sh'
        script_name = script_name.replace('/./', '/')
        ljob_file_path = self.replaceDynamicVariables(self.job_file_path, ind, mon, realid, repeat) if bundle_size == 1 else self.job_file_path.replace('$$$BUNDLE$$$',str(len(bundles)-1))
        outputf = open(script_name, 'w')
        outputf.write('#!/bin/bash\n')
        if bundle_size == 1:
          outputf.write('echo "Starting job ' + str(realid) + ' on $SLURM_NODELIST"\n')
        else:
          outputf.write('export PATH=' + os.environ['PATH'] + '\n')
          outputf.write('export LD_LIBRARY_PATH=''' + os.environ['LD_LIBRARY_PATH'] + '\n')
          outputf.write('cd ' + local_cwd + '\n')

        outputf.write(". " + self.proepi_filename + '\n')
        outputf.write('PrologueStatus ' + self.status_filename + '\n')
        outputf.write('run_command_check_status_slurm ' + self.const_args + ' ')
        outputf.write(self.replaceDynamicVariables(self.variable_args, ind, mon, realid, repeat) + '\n')
        if bundle_size > 1:
          if len(bundles) == 0 or len(bundles[-1]) == bundle_size:
            bundles.append([])
          bundles[-1].append(realid)
          outputb.write(script_name + '\n')
          outputf.write('DecrementBundleCount ''' + str(len(bundles)-1) + '\n')
        outputf.write('exit 0')
        outputf.close()
        os.chmod(script_name, 0755)
        if bundle_size == 1:
          var_sub = self.replaceDynamicVariables(self.variable_submission_cmd, ind, mon, realid, repeat)
          epilogue_script_name = script_name + '.epilogue.sh'
          base_sub = self.base_submission_cmd + ' --epilog=' + epilogue_script_name + ' '
          outputf = open(epilogue_script_name, 'w')
          outputf.write(\
           "#!/bin/bash\ncd " + local_cwd + "\n. \"" + self.proepi_filename + '\"\nSLURMEpilogue ' \
           + ljob_file_path + ' $@\nexit 0')
          outputf.close()
          os.chmod(epilogue_script_name, 0755)
          var_sub += " " + script_name + " "
    
          # touch the output file; necessary due to NFS syncing
          self.all_submission_commands.append("touch  " + self.replaceDynamicVariables(self.job_file_path, ind, mon, realid, repeat))
          self.all_submission_commands.append(base_sub + var_sub + cmd_jobid_from_slurm_output_to_abort_script)
          self.all_submission_commands.append('sleep 1')
  
      if bundle_size > 1:
        outputb.close()
        self.n_pbs_jobs = len(bundles)
        bundles_script=self.scripts_path + os.sep + "run_bundle.sh"
        outputb = open( bundles_script,'w')
        outputb.write('#!/bin/bash\n. "' + self.proepi_filename + '"\nRunSlurmBundle ' + str(bundle_size) + "\nexit 0")
        outputb.close()
        os.chmod(bundles_script, 0755)
        for bundleid in range(len(bundles)):
          job_path = self.job_file_path.replace('$$$BUNDLE$$$',str(bundleid))
          var_sub = self.variable_submission_cmd.replace('$$$BUNDLE$$$',str(bundleid))
          epilogue_script_name = self.scripts_path + os.sep + os.path.basename(job_path.replace('.job.txt','.epilogue.sh'))
          base_sub = self.base_submission_cmd + ' --epilog=' + epilogue_script_name + ' '
          outputf = open(epilogue_script_name, 'w')
          outputf.write(\
           "#!/bin/bash\ncd " + local_cwd + "\n. \"" + self.proepi_filename + '\"\nSLURMEpilogue ' \
           + job_path + ' $@\nexit 0')
          outputf.close()
          os.chmod(epilogue_script_name, 0755)
          var_sub += " -vBUNDLE=" + str(bundleid) + " " + bundles_script + " "
          self.all_submission_commands.append("touch  " + job_path)
          self.all_submission_commands.append(base_sub + var_sub + cmd_jobid_from_slurm_output_to_abort_script)
          self.all_submission_commands.append('sleep 1')
  
      self.all_submission_commands.append(') 300> ' + primary_lock_file)
    output_main = open(self.main_submission_script, 'w')
    output_main.write('#!/bin/bash\n')
    output_main.write('set -e\n')
    output_main.write('echo ' + str(self.number_jobs) + ' 0 ' + str(self.original_number_jobs-self.number_jobs) + ' 0 0 0 > ' + self.status_filename + "\n")
    output_main.write('touch ' + self.status_lockfilename + "; sleep 1\n")
    output_main.write('\n'.join(self.all_submission_commands))
    output_main.write('\nchmod a+x ' + self.abort_script)
    output_main.close()
    os.chmod(self.main_submission_script, 0755)
    if not self.options.dry_run:
      outputf = open(self.status_filename, 'w')
      outputf.write(str(self.number_jobs) + ' 0 ' + str(self.original_number_jobs-self.number_jobs) + ' 0 0 0')
      outputf.close()

  def run(self):
    if not self.options.dry_run and not self.options.just_submit:
      self.status_bar.updateWriteStatus()
    queue = None
    local_job_list = []
    if self.options.pbs_jobs:
      if self.options.dry_run:
        print "Would execute: " + self.host_access_cmd + self.main_submission_script
      else:
        if self.options.max_pbs_jobs and self.options.max_pbs_jobs > 0:
          number_jobs_currently_running = self.options.max_pbs_jobs + 1
          max_jobs_that_can_still_be_running = max(self.options.max_pbs_jobs, self.n_pbs_jobs) - self.n_pbs_jobs
          while number_jobs_currently_running > max_jobs_that_can_still_be_running:
            username = os.getenv('USER', 'user')
            (success, output) = commands.getstatusoutput(self.host_access_cmd + 'qselect -u ' + username + ' -s QR | wc -l')
            new_number_jobs_currently_running = int(output.strip(' \n'))
            if new_number_jobs_currently_running <= max_jobs_that_can_still_be_running:
              break
            if new_number_jobs_currently_running < number_jobs_currently_running:
              print str(new_number_jobs_currently_running) + " jobs still running; waiting until this reaches " + str(max_jobs_that_can_still_be_running)
            sleep(30)
            number_jobs_currently_running = new_number_jobs_currently_running
        nsecs = 10
        success = False
        output = 'No locks available'
        (success, output) = tryExecute(self.host_access_cmd + self.main_submission_script, 0, 'pbs host connection for job submission', 10)
        while not success and output.find('No locks available') >= 0 and nsecs < 71:
          print "NFS file system locking issue. Waiting an additional ",nsecs," seconds and will try again"
          sleep(nsecs)
          nsecs += 60
          (success, output) = tryExecute(self.host_access_cmd + self.main_submission_script, 0, 'pbs host connection for job submission', nsecs)
        if not success and not self.options.no_flock:
          if output.find('No locks available') >= 0:
            tryExecute('sed -i \'s/ flock /# flock /\' ' + self.main_submission_script, 0, 'Turning off file locking', nsecs)
            print "NFS file system locking issue. Attempting to run without locks. This may result in jobs not being " \
                  "scheduled by the scheduler. If this is happening often, write trouble ticket to restart lockd, and statd on this machine"
            nsecs = 10
            (success, output) = tryExecute(self.host_access_cmd + self.main_submission_script, 0, 'pbs host connection for job submission', 10)
          if not success:
            print "\nFailed to submit commands over pbs, error was: " + output + "\nLikely you cannot qsub from " + self.options.host + "; use --host with a machine name that you can qsub from"
            sys.exit(1)
        if self.options.just_submit:
          print "\nJobs submitted"
          self.success = True
          return
    elif self.options.slurm_jobs:
      if self.options.dry_run:
        print "Would execute: " + self.host_access_cmd + self.main_submission_script
      else:
        if self.options.max_slurm_jobs and self.options.max_slurm_jobs > 0:
          number_jobs_currently_running = self.options.max_slurm_jobs + 1
          max_jobs_that_can_still_be_running = max(self.options.max_slurm_jobs, self.n_pbs_jobs) - self.n_pbs_jobs
          while number_jobs_currently_running > max_jobs_that_can_still_be_running:
            username = os.getenv('USER', 'user')
            (success, output) = commands.getstatusoutput(self.host_access_cmd + 'squeue -u ' + username + ' -t PD,R | wc -l')
            new_number_jobs_currently_running = int(output.strip(' \n'))
            if new_number_jobs_currently_running <= max_jobs_that_can_still_be_running:
              break
            if new_number_jobs_currently_running < number_jobs_currently_running:
              print str(new_number_jobs_currently_running) + " jobs still running; waiting until this reaches " + str(max_jobs_that_can_still_be_running)
            sleep(30)
            number_jobs_currently_running = new_number_jobs_currently_running
        nsecs = 10
        success = False
        output = 'No locks available'
        (success, output) = tryExecute(self.host_access_cmd + self.main_submission_script, 0, 'SLURM host connection for job submission', 10)
        while not success and output.find('No locks available') >= 0 and nsecs < 71:
          print "NFS file system locking issue. Waiting an additional ",nsecs," seconds and will try again"
          sleep(nsecs)
          nsecs += 60
          (success, output) = tryExecute(self.host_access_cmd + self.main_submission_script, 0, 'SLURM host connection for job submission', nsecs)
        if not success and not self.options.no_flock:
          if output.find('No locks available') >= 0:
            tryExecute('sed -i \'s/ flock /# flock /\' ' + self.main_submission_script, 0, 'Turning off file locking', nsecs)
            print "NFS file system locking issue. Attempting to run without locks. This may result in jobs not being " \
                  "scheduled by the scheduler. If this is happening often, write trouble ticket to restart lockd, and statd on this machine"
            nsecs = 10
            (success, output) = tryExecute(self.host_access_cmd + self.main_submission_script, 0, 'SLURM host connection for job submission', 10)
          if not success:
            print "\nFailed to submit commands over SLURM, error was: " + output + "\nLikely you cannot srun from " + self.options.host + "; use --host with a machine name that you can srun from"
            sys.exit(1)
        if self.options.just_submit:
          print "\nJobs submitted"
          self.success = True
          return
    elif self.options.local_jobs:
      if self.options.dry_run:
        print "Would execute:\n" + '\n'.join([ self.host_access_cmd + ' ' + cm for cm in self.all_submission_commands])
      else:
        if not self.options.just_submit:
          #spawn a pool of threads, and pass them queue instance
          ncpus = self.options.local_jobs
          queue = Queue.Queue()
          local_job_list = [util.ThreadLocalRun.ThreadLocalRun(queue, None if not self.options.show_status else self.status_bar) for i in range(ncpus)]
          for job in local_job_list:
            job.daemon = True
          for job in local_job_list:
            job.start()
          for cm in self.all_submission_commands:
            queue.put(cm)
        else:
          if self.options.local_jobs < len(self.all_submission_commands):
            sys.exit( """Error: cannot use just-submit to run jobs locally unless all jobs can be run simultaneously. 
Either remove --just-submit, use a different scheduler, or increase the number of local cores to match the jobs, e.g. 
--local """ + str( len(self.all_submission_commands)))
          for cm in self.all_submission_commands:
            print cm
            (status, output) = commands.getstatusoutput(cm)
            if status != 0:
              print "error during job submission: " + output + ' of ' + cm
          return
        
    elif self.options.gparallel:
      # GNU parallels command
      main_cmd = 'echo -e "' + '\\n'.join([str(x[4]) for x in self.job_keys_to_run]) + "\" | parallel -u " + self.options.gparallel + ' ' + str(os.getcwd()) + self.scripts_path[1:] + os.sep + self.options.name + '.{.}.sh &'
      if self.options.dry_run:
        print "Would execute:\n" + main_cmd
      else:
        if os.system(main_cmd) != 0:
          print "Could not execute gnu-parallels with " + main_cmd
          sys.exit(1)

    if self.options.dry_run:
      return

    while not os.path.exists(self.wakeup_file):
      sleep(5)
      self.status_bar.updateWriteStatus()

    self.status_bar.updateWriteStatus()
    if self.options.local_jobs > 0:
      queue.join()
    if not self.status_bar.tryFinalize():
      sleep(1)
      self.status_bar.tryFinalize()
    self.success = self.status_bar.overall_status != BclCVStatus.ERRORS

  def writeWalltime(self):
    localstr = str(int(self.walltime / 60)) + ":"
    mint = int(self.walltime % 60)
    if mint < 10:
      localstr += '0'
    return localstr + str(mint) + ':00'

  def writeLocalScripts(self):
    local_cwd = str(os.getcwd())
    optional_nohup = ' screen -d -m ' if self.options.just_submit else ' '
    for (ind,mon,key,repeat,realid) in self.job_keys_to_run:
      script_name = self.scripts_path + os.sep + self.options.name + '.' + str(realid) + '.sh'
      self.all_submission_commands.append(self.host_access_cmd + self.base_submission_cmd + optional_nohup + script_name)
      job_file = self.replaceDynamicVariables(self.job_file_path, ind, mon, realid, repeat)
      outputf = open(script_name, 'w')
      outputf.write('#!/bin/bash\ncd ' + local_cwd + '\n')
      outputf.write('. ' + self.proepi_filename + '\n')
      outputf.write('PrologueStatus ' + self.status_filename + '\n')
      outputf.write('run_command_check_status ' + self.const_args + ' ' + self.replaceDynamicVariables(self.variable_args, ind, mon, realid, repeat) + ' &> ' + job_file)
      outputf.close()
      os.chmod(script_name, 0755)
    if self.options.dry_run:
      main_submission_script = self.scripts_path + os.sep + self.options.name + "_main.sh"
      output_main = open(main_submission_script, 'w')
      output_main.write('#!/bin/bash\n')
      output_main.write('echo ' + str(self.number_jobs) + ' 0 ' + str(self.original_number_jobs-self.number_jobs) + ' 0 0 0 > ' + self.status_filename + "\n")
      output_main.write('\n'.join(self.all_submission_commands))
      output_main.close()
      os.chmod(main_submission_script, 0755)
    else:
      outputf = open(self.status_filename, 'w')
      outputf.write(str(self.number_jobs) + ' 0 ' + str(self.original_number_jobs-self.number_jobs) + ' 0 0 0')
      outputf.close()
      
  def getDatasetInfo(self, dataset):
    dataset_nm = self.replaceDynamicVariables(dataset, 0, 1, 0, 0).strip('"\'')
    unique_filename = str(uuid.uuid4())
    full_command = self.options.bcl + ' descriptor:GenerateDataset -info -opencl Disable -source \'' + dataset_nm + '\' -output ' + unique_filename + ' ' + self.labels_command
    if commands.getstatusoutput(full_command)[0] != 0:
      if commands.getstatusoutput(full_command.replace('-opencl Disable', ''))[0] != 0:
        print "could not access a dataset to retrieve its size!"
        print full_command.replace('-opencl Disable', '')
        sys.exit(-1)
    file_to_read = file(unique_filename, 'r')
    file_lines = file_to_read.readlines()
    file_to_read.close()
    os.remove(unique_filename)
    file_lines = [ x.strip() for x in file_lines if len(x) > 0 and isdigit(x[0])]
    if len(file_lines) < 5:
      print "GenerateDataset -info returned " + str(len(file_lines)) + " lines, but expected at least 5!"
      sys.exit(1)
    if not file_lines[0].lower().endswith('mb'):
      print "Expected 1st line from GenerateDataset -info -output ... to end with MB"
      sys.exit(1)

    file_tokens = [x.split(' ')[0] for x in file_lines]
    # MB, feature cols, result cols, id cols
    return [ int(file_tokens[0][:-2]), int(file_tokens[1]), int(file_tokens[2]), int(file_tokens[3]) ]

  def writeMergeScript(self):
    if not self.options.print_ind:
      return
    script = "#!/bin/bash\n"
    script += ". " + self.proepi_filename + "\n"
    if not self.options.local_jobs:
      script += "sleep 5 \n"
    # if several descriptor files exist, and this script is being used, then its impossible that they are different
    # descriptor files (since we rm -rf the models directory before starting this), so condense them down into one 
    # descriptor file
    script += "remove_duplicate_descriptor_files " + self.model_files_path + "\n"
    self.merge_log_file = self.log_files_path + os.sep + 'log_merge.txt'
    base_flags = ' -logger File ' + self.merge_log_file + ' '
    if self.options.opencl_driver != 'None':
      base_flags += ' -opencl Disable '
    key = 0
    ind_path_templ = self.independent_files_path + os.sep + self.job_file_tmpl + '.gz'
    ind_path_norepeat_templ = ind_path_templ.replace(' _number' + self.repeat_chunk_tmpl, '')
    merged_independent_files = []
    monitoring_merged_str = str(self.options.mon_id_range[0]) + '-' + str(self.options.mon_id_range[1])
    ind_merged_str = str(self.options.ind_id_range[0]) + '-' + str(self.options.ind_id_range[1])
    all_ind_files = ''
    final_merged_file = self.replaceDynamicVariables(ind_path_templ + '.txt', ind_merged_str, monitoring_merged_str, -1, 0)
    self.primary_result_file = self.independent_files_path + os.sep + BclCommandCreator.final_indmerged_filename
    self.blind_result_file = self.independent_files_path + os.sep + 'final_objective.blind.txt'
    self.consensus_result_file = self.independent_files_path + os.sep + 'final_objective.consensus.txt'
    merge_suffix = ''
    if self.options.use_median:
      merge_suffix = ' -median '
    elif self.options.use_jury:
      merge_suffix = " -jury '" + self.options.final_obj + "' "
    elif self.options.use_max:
      merge_suffix = ' -max '
    elif self.options.use_min:
      merge_suffix = ' -min '
    statistics_suffix =''
    if self.options.no_plot:
      statistics_suffix = ' -no_plot ' 
    if len(self.model_retrieval_string):
      script += 'run_merge_command ' + self.options.bcl + ' model:PredictionMerge ' + \
                self.model_retrieval_string.replace(' -storage_model ',' -input_model_storage ') \
                + ' -output ' + final_merged_file + ' ' + merge_suffix + base_flags + '\n'
    else:
      merge_command_prefix = 'run_merge_command ' + self.options.bcl + ' model:PredictionMerge -input '
      for ind in xrange(self.options.ind_id_range[0], self.options.ind_id_range[1] + 1):
        script += merge_command_prefix
        files_this_monitoring = ''
        for mon in xrange(self.options.mon_id_range[0], self.options.mon_id_range[1] + 1):
          if (mon == ind) != self.options.mon_ind_same:
            continue
          for repeat in xrange(self.options.cv_repeats):
            files_this_monitoring += self.replaceDynamicVariables(ind_path_templ, ind, mon, key, repeat) + ' '
            key += 1
        script += files_this_monitoring
        all_ind_files += files_this_monitoring
        merged_independent_files.append(self.replaceDynamicVariables(ind_path_norepeat_templ, ind, monitoring_merged_str, key, 0))
        script += ' -output ' + merged_independent_files[-1] + ' ' + merge_suffix + ' ' + base_flags + '\n'
        all_ind_files += merged_independent_files[-1] + ' '
    if not self.options.local_jobs:
      script += "sleep 2 \n"
    script += 'run_merge_command ' + self.options.bcl + ' model:ComputeStatistics ' + ' -message_level Verbose -input '\
              + final_merged_file + base_flags + statistics_suffix
    if self.options.final_obj:
      script += " -filename_obj_function " + self.primary_result_file + " -obj_function '" + self.options.final_obj + "' "
    else:
      script += " -table_name " + self.primary_result_file
    script += '\n'
    if self.have_blind and self.options.model_storage:
      self.blind_output_file = self.independent_files_path + os.sep + 'blind_predictions.txt'
      script += 'run_merge_command ' + self.options.bcl + ' model:Test -average -output ' + self.blind_output_file
      script += ' -retrieve_dataset \'' + self.blind_dataset_label + '\' ' + self.model_retrieval_string
      script += base_flags + '\n'
      script += 'run_merge_command ' + self.options.bcl + ' model:ComputeStatistics ' + base_flags + ' -message_level Verbose -input '\
                + self.blind_output_file + statistics_suffix
      if self.options.final_obj:
        script += " -filename_obj_function " + self.blind_result_file + " -obj_function '" + self.options.final_obj + "' "
      else:
        script += " -table_name " + self.blind_result_file
      script += '\n'
      self.consensus_output_file = self.independent_files_path + os.sep + 'consensus_predictions.txt'
      script += 'run_merge_command ' + self.options.bcl + ' model:Test -average -output ' + self.consensus_output_file
      script += ' -retrieve_dataset \'' + self.consensus_dataset_label + '\' ' + self.model_retrieval_string
      script += base_flags + '\n'
      script += 'run_merge_command ' + self.options.bcl + ' model:ComputeStatistics ' + base_flags + ' -message_level Verbose -input '\
                + self.consensus_output_file + statistics_suffix
      if self.options.final_obj:
        script += " -filename_obj_function " + self.consensus_result_file + " -obj_function '" + self.options.final_obj + "' "
      else:
        script += " -table_name " + self.consensus_result_file
    script += '\n'
    if self.options.remove_files:
      script += 'rm ' + all_ind_files + '\n'
    self.merge_script = self.scripts_path + os.sep + 'merge_ind_and_compute_results.sh'
    outputf = open(self.merge_script, 'w')
    outputf.write(script)
    outputf.close()
    os.chmod(self.merge_script, 0755)
    if not self.number_jobs:
      just_run_merge = not os.path.exists(self.primary_result_file)
      print just_run_merge," ",os.path.exists(self.primary_result_file)," ",self.primary_result_file
      if not just_run_merge and self.have_blind:
        print just_run_merge," ",os.path.exists(self.blind_result_file)," ",os.path.exists(self.consensus_result_file)
        if not os.path.exists(self.blind_result_file) or not os.path.exists(self.consensus_result_file):
          just_run_merge = True
      if just_run_merge:
        if self.options.pbs_jobs or self.options.slurm_jobs:
          output_main = open(self.main_submission_script, 'a')
          output_main.write('\n' + self.merge_script + '\n')
        else:
          self.all_submission_commands.append(self.host_access_cmd + self.merge_script)
      else:
        if self.options.pbs_jobs or self.options.slurm_jobs:
          output_main = open(self.main_submission_script, 'a')
          output_main.write('\necho 0 0 ' + str(self.original_number_jobs) + ' 0 0 1 > ' + self.status_filename + '\n')
          output_main.write('\ntouch ' + self.wakeup_file + '\n')
          output_main.close()
        else:
          self.all_submission_commands.append('echo 0 0 ' + str(self.original_number_jobs) + ' 0 0 1 > ' + self.status_filename)
          self.all_submission_commands.append('touch ' + self.wakeup_file)

  def __init__(self, option_args, remaining_args):
    # short circuit; if the --complete flag was given, and results are already available, then skip all the rest of the steps
    self.have_blind = ( option_args.blind_descriptor_ids and len(option_args.blind_descriptor_ids)) \
                      or ( option_args.blind_id_nchunks and len(option_args.blind_id_nchunks)) \
                      or ( option_args.blind_features and len(option_args.blind_feature_ranges))
    if option_args.complete and option_args.print_ind:
      resfile='./results' + os.sep + option_args.name + os.sep + BclCommandCreator.final_indmerged_filename
      if os.path.exists( resfile):
        status_filename = 'log_files' + os.sep + option_args.name + os.sep + 'status.txt'
        stat = BclCVStatus(False, False, status_filename, [], resfile, '', self.have_blind)
        stat.updateWriteStatus()
        sys.exit(0)
    
    # validation
    if option_args.chunks < 3:
      print "# chunks must be at least 3 to allow for unique training, monitoring, and independent datasets"
      sys.exit(-1)
    if option_args.mon_size_limit != None and option_args.mon_size_limit < 1:
      print "At least monitoring data point is required"
      sys.exit(-1)
    if option_args.train_size_limit != None and option_args.train_size_limit < 1:
      print "At least training data point is required"
      sys.exit(-1)
    if option_args.ind_id_range[0] < 0 or option_args.mon_id_range[0] < 0:
      print "Independent and monitoring id range mins must be at least 0"
      sys.exit(-1)
    if option_args.ind_id_range[0] > option_args.ind_id_range[1] or option_args.mon_id_range[0] > option_args.mon_id_range[1]:
      print "Independent and monitoring id range maxs must be >= range mins!"
      sys.exit(-1)
    if option_args.ind_id_range[1] >= option_args.chunks or option_args.mon_id_range[1] >= option_args.chunks:
      print "Independent and monitoring id range maxs must be < # chunks!"
      sys.exit(-1)
    if sum([int(x.endswith('.bin')) for x in option_args.datasets]) != len(option_args.datasets):
      print "All datasets must be bin files!"
      sys.exit(-1)
    if option_args.host != os.getenv("HOSTNAME") and option_args.gparallel:
      print "Host flag cannot be given with parallels flag"
      sys.exit(-1)
    if (option_args.score_file != None) != (option_args.top_n_features != None):
      print "--score-file must be given with --top-features"
      sys.exit(-1)
    elif option_args.score_file != None and not os.path.exists(option_args.score_file):
      print "Given score file: " + option_args.score_file + " does not exist!"
      sys.exit(-1)
    if option_args.use_jury and not option_args.final_obj:
      print "Jury can only be given if a final objective function is given!"
      sys.exit(-1)
    if sum([option_args.use_jury, option_args.use_min, option_args.use_max, option_args.use_median]) > 1:
      print "Can only combine independent predictions via one of use-jury, use-max, use-min, use-median"
      sys.exit(-1)
    if option_args.model_storage != 'File' and option_args.only_best_models:
      print '--best-model can only be used when storing models in files currently (e.g. --model_storage File must be given)'
      sys.exit(-1)
    if option_args.filter_descriptor and not option_args.filter_descriptor_range:
      print "--filter-descriptor requires that --filter-descriptor-range also be given"
      sys.exit(-1)
    if option_args.filter_result and not option_args.filter_result_range:
      print "--filter-result-range must be provided whenever --filter-result is provided"
      sys.exit(-1)
    if option_args.filter_result and len(option_args.filter_result) != len(option_args.filter_result_range):
      print "Must have the same number of result filter ranges as result filters"
      sys.exit(-1)
    if option_args.filter_descriptor and len(option_args.filter_descriptor) != len(option_args.filter_descriptor_range):
      print "Must have the same number of descriptor filter ranges as descriptor filters"
      sys.exit(-1)
    if option_args.filter_id_descriptor and len(option_args.filter_id_descriptor) != len(option_args.filter_ids):
      print "Must have the same number of id select ranges as id select descriptors"
      sys.exit(-1)
    if option_args.remove_ids and len(option_args.remove_id_descriptor) != len(option_args.remove_ids):
      print "Must have the same number of id exclude ranges as id filters"
      sys.exit(-1)
    if option_args.blind_descriptor_ids is None:
      option_args.blind_descriptor_ids = []
    if option_args.blind_feature_ranges is None:
      option_args.blind_features = []
      option_args.blind_feature_ranges = []
    if len(option_args.blind_descriptor_ids) and option_args.blind_id_nchunks and len(option_args.blind_id_nchunks):
      print "Blind descriptor ids cannot be used with blind chunks at this time!"
      sys.exit(-1)
    if len(option_args.blind_feature_ranges) and option_args.blind_id_nchunks and len(option_args.blind_id_nchunks):
      print "Blind descriptor ranges cannot be used with blind chunks at this time!"
      sys.exit(-1)
    if len(option_args.blind_feature_ranges) and len(option_args.blind_descriptor_ids):
      print "Blind descriptor ranges cannot be used with blind desriptor ids at this time!"
      sys.exit(-1)
    if len(option_args.blind_feature_ranges) and len(option_args.blind_feature_ranges) != len(option_args.blind_features):
      print "Must have the same # of blind feature ranges as blind features!"
      sys.exit(-1)
    self.have_blind = len(option_args.blind_descriptor_ids) or len(option_args.blind_feature_ranges) \
                      or ( option_args.blind_id_nchunks and len(option_args.blind_id_nchunks))
    if not '-random_seed' not in remaining_args and '--random_seed' not in remaining_args and options.cv_repeats > 1:
      print "Repeated cross validation runs would be the same since --random_seed was not given; auto-setting random-seed"
      remaining_args.append('-random_seed')
    # if jobs are started locally, it is likely that all the jobs will start at exactly the same time; likewise, all of
    # them will receive the same random seed, which is not good.  Here we instead use /dev/urandom to get random numbers
    if '-random_seed' in remaining_args or '--random_seed' in remaining_args:
      remaining_args = [ x for x in remaining_args if x.lstrip('-') != 'random_seed']
      remaining_args.append('-random_seed')
      remaining_args.append('`od -A n -t dI -N 4 /dev/urandom | tr -d -`')

    self.output_directories = []
    processor_count = 0
    if option_args.pbs_jobs:
      processor_count += 1
    if option_args.local_jobs:
      processor_count += 1
    if option_args.gparallel:
      processor_count += 1
    if option_args.slurm_jobs:
      processor_count += 1
    if processor_count != 1:
      print "Job run location must be given (i.e. one of --local|--gparallel|--pbs|--slurm) must be given"
      sys.exit(-1)
    # ensure that the bcl is accessible
    if commands.getstatusoutput(option_args.bcl + " -help")[0] != 0:
      print option_args.bcl + " was not accessible!"
      sys.exit(-1)
    if option_args.host.split('.')[0] != os.getenv('HOSTNAME',option_args.host).split('.')[0]:
      host_access = "ssh " + option_args.host + ' cd ' + os.getcwd()
      host_bcl = host_access + ' \; ' + option_args.bcl + " -help"
      (was_okay, msg) = tryExecute(host_bcl, 10, option_args.bcl + " was not accessible from " + option_args.host, 10)
      if not was_okay:
        if commands.getstatusoutput(host_access)[0] != 0:
          print "Could not access host: " + option_args.host + " output: " + msg
          sys.exit(-1)
        else:
          print option_args.bcl + " was not accessible from " + option_args.host + " output: " + msg
    if option_args.name.find(' ') >= 0:
      print "Cannot have spaces in id!"
      sys.exit(-1)
    elif str(os.getcwd()).find(' ') >= 0:
      print "cannot have spaces in current working directory"
      sys.exit(-1)

    # handle non-local host
    self.host_access_cmd = ""
    hostname=os.getenv('HOSTNAME',option_args.host)
    if option_args.host.split('.')[0] != hostname.split('.')[0]:
      self.host_access_cmd = "ssh " + option_args.host + ' cd ' + os.getcwd() + ' \; '

    self.options = option_args
    self.abort_script = ''
    this_id = option_args.name
    self.log_files_path = './log_files/' + this_id
    self.scripts_path = self.log_files_path + os.sep + 'scripts'
    if self.options.complete and not os.path.exists(self.scripts_path):
      print "Job was not started previously, ignoring complete flag"
      self.options.complete = option_args.complete = False 
    self.independent_files_path = ''
    self.model_files_path = ''
    self.proepi_filename = self.scripts_path + os.sep + 'functions.sh'
    self.status_filename = self.log_files_path + os.sep + "status.txt"
    
    sha1obj=hashlib.sha1()
    sha1obj.update(self.status_filename)
    self.status_lockfilename = ( "/tmp/bcl_cv_" + sha1obj.hexdigest() if self.options.local_jobs else self.status_filename) + '.lock'
    
    self.main_submission_script = self.scripts_path + os.sep + "main.sh"
    self.wakeup_file = self.scripts_path + os.sep + 'wakeup'

    putative_abort_script = self.scripts_path + os.sep + 'abort.sh'
    if os.path.exists(putative_abort_script) and os.path.exists(self.status_filename) and os.path.exists(self.status_lockfilename):
      # check for the lock file; if it is absent, then the job is finished
      try:
        temp_status = BclCVStatus(False, False, self.status_filename, [], '', '', False)
        temp_status.updateStatus()
        if temp_status.jobs_in_queue or temp_status.jobs_running:
          print "Warning: cross validation with ID " + self.options.name + " exists and is still running; stopping"
          (status, output) = commands.getstatusoutput(putative_abort_script)
          if status != 0:
            print "Failed to stop previous cross validation run. Resubmitting jobs anyway"
          else:
            # sleep for a few seconds to allow pbs time to really kill all the jobs
            sleep(10)
          os.unlink(self.status_lockfilename)
      except:
        pass
    if self.options.pbs_jobs or self.options.slurm_jobs:
      self.abort_script = putative_abort_script

    self.mkdirmpWithSymlinks(self.log_files_path)
    if self.options.complete:
      # manually remove any problematic files
      if os.path.exists(self.status_lockfilename):
        os.unlink( self.status_lockfilename)
      if os.path.exists(self.status_filename + '.bclerr'):
        os.unlink( self.status_filename + '.bclerr')
      if os.path.exists(self.status_filename + '.pbserr'):
        os.unlink( self.status_filename + '.pbserr')
      if os.path.exists(self.status_filename + '.slurmerr'):
        os.unlink( self.status_filename + '.slurmerr')
      if os.path.exists(self.abort_script):
        os.unlink( self.abort_script)
      if os.path.exists(self.wakeup_file):
        os.unlink( self.wakeup_file)

    if not self.options.no_write_scripts:
      self.mkdirmpWithSymlinks(self.scripts_path)

    if self.options.score_file != None:
      base_features_file = self.options.score_file[self.options.score_file.rfind(os.sep) + 1:]
      base_features_file = str(os.path.splitext(base_features_file)[0])

      new_features_file = self.log_files_path + os.sep + base_features_file + '.top' + str(self.options.top_n_features) + '.txt'
      create_features_command = self.options.bcl + ' descriptor:RefineByScore -select "Top(' \
                                + str(self.options.top_n_features) + ')" -output ' + new_features_file + ' -score_file ' \
                                + self.options.score_file
      (status, output) = commands.getstatusoutput(create_features_command)
      if status != 0:
        print "Failed to create new feature set using command:\m" + create_features_command + "\nproblem: " + output
        sys.exit(-1)
      self.options.features = new_features_file

    self.labels_command = ""
    if option_args.results != None and option_args.results != '':
      self.labels_command += '-result_labels "' + option_args.results + '" '
    if option_args.features != None and option_args.features != '':
      self.labels_command += '-feature_labels "' + option_args.features + '" '

    self.training_dataset_tmpl = ""
    self.monitoring_dataset_tmpl = ""
    self.independent_dataset_tmpl = ""
    self.variable_args = ""
    self.merge_script = ""
    self.using_opencl = option_args.iterate.lower().startswith('opencl')
    if self.using_opencl and self.options.opencl_driver == 'Disable':
      print "Using an opencl iterate; automatically switiching to NVIDIA "
      self.options.opencl_driver = 'NVIDIA'

    self.const_args = self.options.bcl + " model:Train '" + self.options.iterate + "' "

    train_chunks_prefix = ''
    train_chunks_postfix = ''
    mon_chunks_prefix = ''
    mon_chunks_postfix = ''
    if option_args.train_size_limit != None:
      train_chunks_prefix = 'Rows(rows="[0, ' + str(option_args.train_size_limit) + ')",dataset='
      train_chunks_postfix = ')'
    if option_args.mon_size_limit != None:
      mon_chunks_prefix = 'Rows(rows="[0, ' + str(option_args.mon_size_limit) + ')",dataset='
      mon_chunks_postfix = ')'

    train_chunk = mon_chunk = '0, ' + str(option_args.chunks) + ') '
    if not self.options.no_cv:
      train_chunk += '- [' + self.ind_chunk_tmpl + ']'
      if not self.options.include_m_in_t:
        train_chunk += ' - [' + self.mon_chunk_tmpl + ']'
      mon_chunk = self.mon_chunk_tmpl + ']'
      if self.options.swap_t_and_m:
        x = mon_chunk
        mon_chunk = train_chunk
        train_chunk = x

    if self.options.blind_id_nchunks and len(self.options.blind_id_nchunks) == 1:
      self.options.blind_id_nchunks.append(option_args.chunks)

    dataset_var = "$$$DATASET$$$"
    training_dataset_filter = dataset_var
    if self.options.filter_descriptor:
      for i in range(len(self.options.filter_descriptor)):
        training_dataset_filter = 'FeatureRange(dataset=' + training_dataset_filter + ',range="' + self.options.filter_descriptor_range[i] + '",feature="' + self.options.filter_descriptor[i] + '")'
    if self.options.filter_result:
      for i in range(len(self.options.filter_result_range)):
        training_dataset_filter = 'ResultRange(dataset=' + training_dataset_filter + ',range="' + self.options.filter_result_range[i] + '",result="' + self.options.filter_result[i] + '")'
    if self.options.filter_ids:
      for i in range(len(self.options.filter_ids)):
        ids_quoted = '"' + '","'.join(self.options.filter_ids[i].split(',')) + '"'
        training_dataset_filter = 'SelectIDs(invert=False,dataset=' + training_dataset_filter + ',id type="' + self.options.filter_id_descriptor[i] + '",ids(' + ids_quoted + '))'
    if self.options.remove_ids:
      for i in range(len(self.options.remove_ids)):
        ids_quoted = '"' + '","'.join(self.options.remove_ids[i].split(',')) + '"'
        training_dataset_filter = 'SelectIDs(invert=True,dataset=' + training_dataset_filter + ',id type="' + self.options.remove_id_descriptor[i] + '",ids(' + ids_quoted + '))'
    
    mon_ind_dataset_filter = training_dataset_filter if self.options.filter_all or self.options.filter_eval else dataset_var
    if self.options.filter_eval:
      training_dataset_filter = dataset_var
    subset_prefix = ""
    if not self.have_blind:
      subset_prefix = 'Subset(number chunks=' + str(option_args.chunks) + ',chunks="[$$$CHUNKS$$$",filename="$$$FILENAME$$$")'
    elif self.options.blind_id_nchunks and len(self.options.blind_id_nchunks):
      blind_chunks = self.options.blind_id_nchunks[1]
      blind_id = self.options.blind_id_nchunks[0]
      if blind_id >= blind_chunks:
        print "Bad blind id #: " + str(blind_id) + " should be less than # chunks: " + str(blind_chunks)
        sys.exit(-1)
      full_blind_subset = 'Subset(number chunks=' + str(blind_chunks) + ',chunks="[0,' + str(blind_chunks) + ')-[' + str(blind_id) + ']",filename="$$$FILENAME$$$")'
      subset_prefix = 'Chunks(number chunks=' + str(option_args.chunks) + ',chunks="[$$$CHUNKS$$$",dataset=' + full_blind_subset + ')'
    elif self.options.blind_descriptor:
      blind_descriptor = self.options.blind_descriptor
      blind_descriptor_strings = self.options.blind_descriptor_ids
      ids_quoted = '"' + '","'.join(blind_descriptor_strings) + '"'
      print ids_quoted
      print blind_descriptor
      print str(option_args.chunks)
      subset_prefix = 'SelectIDs(invert=True,id type="' + blind_descriptor + '",ids('+ids_quoted+'),dataset=Subset(number chunks=' + str(option_args.chunks) + ',chunks="[$$$CHUNKS$$$",filename="$$$FILENAME$$$"))'
    elif self.options.blind_features:
      subset_prefix = 'Subset(number chunks=' + str(option_args.chunks) + ',chunks="[$$$CHUNKS$$$",filename="$$$FILENAME$$$")'
      for x in range(len(self.options.blind_features)):
        blind_feature = self.options.blind_features[x]
        blind_feature_range = self.options.blind_feature_ranges[x]
        subset_prefix = 'FeatureRange(invert=True,dataset=' + subset_prefix + ',range="' + blind_feature_range + '",feature="' + blind_feature + '")'
    independent_prefix = mon_ind_dataset_filter.replace(dataset_var,subset_prefix.replace('$$$CHUNKS$$$', self.ind_chunk_tmpl + ']'))
    monitoring_prefix = mon_chunks_prefix + mon_ind_dataset_filter.replace(dataset_var,subset_prefix.replace('$$$CHUNKS$$$', mon_chunk))
    monitoring_postfix = mon_chunks_postfix
    training_prefix = train_chunks_prefix + training_dataset_filter.replace(dataset_var,subset_prefix.replace('$$$CHUNKS$$$', train_chunk))
    training_postfix = train_chunks_postfix

    independent_outer_prefix = train_outer_prefix = mon_outer_prefix = ""
    independent_outer_postfix = train_outer_postfix = mon_outer_postfix = " "
    if len(option_args.datasets) > 1:
      chosen_type = 'Balanced(' if option_args.balanced else 'Combined('
      independent_outer_prefix = independent_outer_prefix + 'Combined('
      independent_outer_postfix = ')' + independent_outer_postfix
      train_outer_prefix = train_outer_prefix + train_chunks_prefix + chosen_type
      train_outer_postfix = ')' + train_chunks_postfix + train_outer_postfix
      mon_outer_prefix = mon_outer_prefix + mon_chunks_prefix + 'Combined('
      mon_outer_postfix = ')' + mon_chunks_postfix + mon_outer_postfix

    independent_dataset_str = independent_outer_prefix
    monitoring_dataset_str = mon_outer_prefix
    training_dataset_str = train_outer_prefix

    total_dataset_size = 0
    for dataset in option_args.datasets:
      independent_dataset_str += independent_prefix.replace('$$$FILENAME$$$', dataset) + ','
      monitoring_dataset_str += monitoring_prefix.replace('$$$FILENAME$$$', dataset) + monitoring_postfix + ','
      training_dataset_str += training_prefix.replace('$$$FILENAME$$$', dataset) + training_postfix + ','
      total_dataset_size += os.path.getsize(dataset)

    # remove trailing comma, add suffices 
    self.independent_dataset_tmpl = independent_dataset_str[:-1] + independent_outer_postfix
    self.monitoring_dataset_tmpl = monitoring_dataset_str[:-1] + mon_outer_postfix
    self.training_dataset_tmpl = training_dataset_str[:-1] + train_outer_postfix
    if self.options.swap_t_and_m and self.options.mon_ind_same:
      self.independent_dataset_tmpl = self.monitoring_dataset_tmpl

    if self.have_blind:
      blind_str = ''
      consensus_str = ''
      if self.options.blind_id_nchunks and len(self.options.blind_id_nchunks) == 2:
        blind_chunks = self.options.blind_id_nchunks[1]
        blind_id = self.options.blind_id_nchunks[0]
        blind_str = 'Subset(number chunks=' + str(blind_chunks) + ',chunks=[' + str(blind_id) + '],filename="$$$FILENAME$$$")'
        consensus_str = 'Subset(number chunks=' + str(blind_chunks) + ',chunks="[0,' + str(blind_chunks) + ') - [' + str(blind_id) + ']",filename="$$$FILENAME$$$")'
      elif len(self.options.blind_descriptor_ids):
        blind_descriptor = self.options.blind_descriptor
        blind_descriptor_strings = self.options.blind_descriptor_ids
        ids_quoted = '"' + '","'.join(blind_descriptor_strings) + '"'
        blind_str = 'SelectIDs(invert=False,id type="' + blind_descriptor + '",ids('+ids_quoted+'),dataset=Subset(filename="$$$FILENAME$$$"))'
        consensus_str = blind_str.replace('invert=False','invert=True',1)
      else:
        blind_str = 'Subset(filename="$$$FILENAME$$$")'
        for x in range(len(self.options.blind_features)):
          blind_feature = self.options.blind_features[x]
          blind_feature_range = self.options.blind_feature_ranges[x]
          blind_str = 'FeatureRange(dataset=' + blind_str + ',invert=False,range="' + blind_feature_range + '",feature="' + blind_feature + '")'
        consensus_str = blind_str.replace('invert=False','invert=True')
         
      if len(option_args.datasets) > 1:
        self.blind_dataset_label = 'Combine('
        self.consensus_dataset_label = 'Combine('
        for dataset in option_args.datasets:
          self.blind_dataset_label += mon_ind_dataset_filter.replace(dataset_var,blind_str.replace('$$$FILENAME$$$', dataset)) + ','
          self.consensus_dataset_label += mon_ind_dataset_filter.replace(dataset_var,consensus_str.replace('$$$FILENAME$$$', dataset)) + ','
        self.blind_dataset_label = self.blind_dataset_label[:-1] + ')'
        self.consensus_dataset_label = self.consensus_dataset_label[:-1] + ')'
      else:
        self.blind_dataset_label = mon_ind_dataset_filter.replace(dataset_var,blind_str.replace('$$$FILENAME$$$', option_args.datasets[0]))
        self.consensus_dataset_label = mon_ind_dataset_filter.replace(dataset_var,consensus_str.replace('$$$FILENAME$$$', option_args.datasets[0]))
    
    if self.options.watch_blind:
      self.monitoring_dataset_tmpl = self.independent_dataset_tmpl = self.blind_dataset_label
    
    total_dataset_size = total_dataset_size >> 20
    print datetime.datetime.today()
    print "Total bin file sizes: " + str(total_dataset_size) + "MB"
    mb_per_sec = 4
    estimated_min_to_load = total_dataset_size / mb_per_sec / 60
    print "Max minutes to load dataset (@" + str(mb_per_sec) + "MB/s): " + str(estimated_min_to_load)

    # add an extra 10 minutes just to be safe
    self.walltime = option_args.minutes + estimated_min_to_load + min(option_args.minutes + 2, 20);
    print "Estimated max walltime per job: " + self.writeWalltime()
    self.const_args += ' -max_minutes ' + str(option_args.minutes)
    self.const_args += ' -max_iterations ' + str(option_args.iterations)
    self.const_args += ' ' + self.labels_command + ' '

    # disable opencl unless it is really needed
    self.const_args += ' ' + self.gpu_map[self.options.opencl_driver.lower()] + ' '
    if self.options.opencl_driver != 'Disable' and self.options.opencl_driver != 'None':
      # set the gpu id, if it was given; or use the one given in the pbs gpu file
      if self.options.opencl_gpu_id >= 0:
        self.const_args += ' ' + str(self.options.opencl_gpu_id) + ' '
      elif self.options.pbs_jobs:
        self.const_args += ' `awk \'{ if( FNR != 1) printf( " % s",", ") ; printf( " % i", substr( $1, length( $1), 1))}\' $PBS_GPUFILE` '
      elif self.options.slurm_jobs:
        self.const_args += ' $SLURM_JOB_GPUS ' 
        
    # threads
    if option_args.pthreads != None:
      self.const_args += ' -scheduler PThread ' + str(option_args.pthreads)

    if remaining_args != None and len(remaining_args):
      self.const_args += ' ' + ' '.join(remaining_args) + ' '

    if option_args.final_obj != None:
      self.const_args += " -final_objective_function '" + option_args.final_obj + "' "

    # determine dataset sizes
    self.dataset_mb = 0
    self.independent_mb = 0
    self.dataset_feature_cols = 0
    self.dataset_result_cols = 0
    self.dataset_rows = 0
    self.training_rows = 0

    # get monitoring dataset info
    monitoring_dataset_info = self.getDatasetInfo(self.monitoring_dataset_tmpl)
    self.dataset_mb = monitoring_dataset_info[0]
    self.dataset_feature_cols = monitoring_dataset_info[2]
    self.dataset_result_cols = monitoring_dataset_info[3]

    # compute the ratio between result columns and total columns; result columns require more memory due to prediction
    # and objective functions
    result_col_multiplier = (self.dataset_result_cols - 1) * 4 + 18
    self.result_column_ratio = float(result_col_multiplier) / float(self.dataset_feature_cols + result_col_multiplier)

    self.independent_mb = self.getDatasetInfo(self.independent_dataset_tmpl)[0]
    self.dataset_mb += self.independent_mb

    #training dataset info
    training_dataset_info = self.getDatasetInfo(self.training_dataset_tmpl)
    training_dataset_mb = training_dataset_info[0]

    # expected memory for the bcl executable and reading monitoring and independent datasets
    base_mem = 150 + (1.01 + self.result_column_ratio) * self.dataset_mb + self.options.memory_offset
    multiplier = 0
    if option_args.memory_multiplier != None:
      multiplier = option_args.memory_multiplier + self.result_column_ratio
      self.training_mb = multiplier * training_dataset_mb + base_mem
    else:
      iterate_decomposed = option_args.iterate.replace('(', ' ').replace(')', ' ').replace('=', ' ').lower().split(' ')
      multiplier = 0
      for iterate, it_multiplier in self.iterate_to_memory_multiplier.iteritems():
        if iterate in iterate_decomposed:
          multiplier = it_multiplier
          break
      if multiplier == 0:
        print "Iterate: " + option_args.iterate + " requires an unknown multiple of dataset memory for training."
        print "Assuming 3x for now; if this is wrong, update BclModelCrossValidation.py script or set flag "
        print " -override - memory - multiplier with the proper value"
        multiplier = 3.03 + self.result_column_ratio

      # balanced datasets require up to twice the memory because the core datasets are loaded first, then balanced and discared
      if multiplier < 1.5 and self.options.balanced:
        multiplier = 2.02
      multiplier += self.result_column_ratio
      self.training_mb = multiplier * training_dataset_mb + base_mem
    self.training_mb = int(self.training_mb + 1) * self.options.job_bundle_size
    if self.options.show_status:
      print "Given iterate/dataset combination has memory multiplier of " + str(int(multiplier * 100) / 100.0)
      print "Training is estimated to consume up to " + str(self.training_mb) + " MB " + ( " per CV run " if self.options.job_bundle_size == 1 else " per bundle of " + str(self.options.job_bundle_size) + " jobs")

    # create variable parts of the string
    self.variable_args += " -training '" + self.training_dataset_tmpl + "' " \
                          + " -monitoring '" + self.monitoring_dataset_tmpl + "' " \
                          + " -independent '" + self.independent_dataset_tmpl + "' "

    if option_args.print_ind:
      self.independent_files_path = "./results/" + this_id
      self.mkdirmpWithSymlinks(self.independent_files_path)
      self.variable_args += " --print_independent_predictions " + self.independent_files_path + os.sep + self.job_file_tmpl + '.gz '
      if option_args.print_train:
        self.variable_args += " --print_training_predictions " + self.independent_files_path + os.sep + 'train_' + self.job_file_tmpl + '.gz '

    if option_args.model_storage:
      self.model_storage_string = ''
      if option_args.model_storage == 'Db':
        self.model_retrieval_string = " -storage_model 'Db(session=" + this_id + ")' "
        self.model_storage_string = self.model_retrieval_string[:-3] + ",cv_ind_id=" + self.ind_chunk_tmpl + ",cv_mon_id=" + self.mon_chunk_tmpl + ",cv_total=" + str(option_args.chunks) + ")' "
      else:
        self.model_files_path = './models/' + this_id
        if not self.options.no_write_scripts:
          self.mkdirmpWithSymlinks(self.model_files_path)
        self.model_retrieval_string = " -storage_model 'File(directory=" + self.model_files_path + ",prefix=\"model\")' "
        self.model_storage_string = self.model_retrieval_string[:-3] + ",write_descriptors=1,key=" + self.key_chunk_tmpl + ")' "
        if self.options.only_best_models:
          self.model_retrieval_string = self.model_retrieval_string[:-3] + ",pick best=True)' "
      self.variable_args += self.model_storage_string

    if option_args.meta_storage:
      rond = str(option_args.round)
      if option_args.meta_storage == 'Db':
        self.variable_args += " -storage_descriptor_selection 'Db(session=" + this_id + ",round=" + rond + ",cv_ind_id=" + self.ind_chunk_tmpl + ",cv_mon_id=" + self.mon_chunk_tmpl + ",cv_total=" + str(option_args.chunks) + ")' "
      else:
        if not self.options.no_write_scripts:
          self.model_files_path = self.mkdirmpWithSymlinks('./models/' + this_id)
        self.variable_args += " -storage_descriptor_selection 'File(directory=" + self.model_files_path + ",prefix=\"model\"" + ",round=" + rond + ",cv_ind_id=" + self.ind_chunk_tmpl + ",cv_mon_id=" + self.mon_chunk_tmpl + ",cv_total=" + str(option_args.chunks) + ")' "

    self.log_file_path_tmpl = self.log_files_path + os.sep + self.job_file_tmpl + '.txt'
    self.variable_args += ' -logger File ' + self.log_file_path_tmpl
    self.merge_script = ''

    self.base_submission_cmd = ''
    self.variable_submission_cmd = ''
    local_cwd = os.getcwd()
    self.job_file_path = local_cwd + self.log_files_path[1:] + os.sep + ( self.job_file_tmpl if self.options.job_bundle_size == 1 else "bundle$$$BUNDLE$$$" ) + ".job.txt"
    ppn = 1
    if self.options.pthreads != None:
      ppn = self.options.pthreads
    #ppn *= self.options.job_bundle_size
    if self.options.pbs_jobs:
      extra_pbs_options = ''
      if self.options.pbs_resources:
        extra_pbs_options = ':' + self.options.pbs_resources
      self.base_submission_cmd += ' qsub -b30 -m n -d ' + local_cwd + ' -N "' + self.options.name + '" ' + \
        self.options.pbs_jobs + ' -j oe ' + ' -l nodes=' + str(self.options.job_bundle_size) + ':ppn=' + str(ppn) + \
        extra_pbs_options + ',mem=' + str(self.training_mb) + 'mb,walltime=' + self.writeWalltime()
      self.variable_submission_cmd = " -o '" + self.job_file_path + "' "
    if self.options.slurm_jobs:
      extra_slurm_options = ''
      if self.options.slurm_resources:
        extra_slurm_options = '--constraint="' + self.options.slurm_resources + '" ' 
      partition_flag = '--partition=' + self.options.slurm_partition if len(self.options.slurm_partition) else ''
      self.base_submission_cmd = ' srun ' + partition_flag + ' -vv -n ' + str(self.options.job_bundle_size) + ' -c ' + str(ppn) + \
                                 ' --mem-per-cpu=' + str(int(self.training_mb/ppn)) + ' -t ' + self.writeWalltime()
      self.variable_submission_cmd = " -o '" + self.job_file_path + "' "              

    # determine the # of keys and write out the status
    self.status_bar = self.CreateStatus(self.options)
    if self.options.no_write_scripts:
      return
    
    self.job_keys_to_run = []
    key = 0
    realkey = 0
    for ind in xrange(self.options.ind_id_range[0], self.options.ind_id_range[1] + 1):
      for mon in xrange(self.options.mon_id_range[0], self.options.mon_id_range[1] + 1):
        if (mon == ind) != self.options.mon_ind_same:
          continue
        for repeat in xrange(self.options.cv_repeats):
          if not self.options.complete:
            self.job_keys_to_run.append((ind,mon,key,repeat,realkey))
            key += 1
          else:
            if not self.checkLogFileFinished(self.replaceDynamicVariables(self.log_file_path_tmpl, ind, mon, key, repeat)):
              self.job_keys_to_run.append((ind,mon,key,repeat,realkey))
              key += 1
          realkey += 1
          
    self.number_jobs = len(self.job_keys_to_run)
    self.original_number_jobs = realkey
    self.all_submission_commands = []
    if self.options.pbs_jobs:
      self.writePbsScripts()
    elif self.options.slurm_jobs:
      self.writeSlurmScripts()
    else:
      self.writeLocalScripts()
    self.writeMergeScript()
    self.writeCommonCommands()
    self.success = False
    if (self.options.pbs_jobs or self.options.slurm_jobs) and not self.options.no_fail_prop:
      self.status_bar.abort_command = self.host_access_cmd + self.abort_script

def main():

  parser = BclCommandCreator.getParser()
  (option_args, remaining_args) = parser.parse_known_args()
  commander = BclCommandCreator(option_args, remaining_args)
  outputf = open(commander.log_files_path + os.sep + 'command.txt', 'w')
  outputf.write(' '.join(sys.argv) + '\n')
  outputf.close()

  commander.run()

  if commander.success or commander.options.dry_run or commander.options.no_write_scripts:
    sys.exit(0)
  else:
    sys.exit(-1)

if __name__ == '__main__':
    main()
