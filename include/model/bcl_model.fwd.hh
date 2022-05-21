// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons. 
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c) 
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_MODEL_FWD_HH_
#define BCL_MODEL_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the model namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace model
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AlignCutoff;
    class ApproximatorAlignCutoff;
    class ApproximatorApplicabilityDomainKohonen;
    class ApproximatorBase;
    class ApproximatorDecisionTree;
    class ApproximatorKappaNearestNeighbor;
    class ApproximatorKohonenNetwork;
    class ApproximatorLeverageMatrix;
    class ApproximatorLinearRegression;
    class ApproximatorNeuralNetwork;
    class ApproximatorNeuralNetworkSelective;
    class ApproximatorRestrictedBoltzmannMachine;
    class ApproximatorSupportVectorMachine;
    class ApproximatorSupportVectorMachineMultiOutput;
    class CollectFeaturesAbove;
    class CollectFeaturesInterface;
    class CollectFeaturesTop;
    class CrossValidationInfo;
    class DataSetLog;
    class DataSetMultiplied;
    class DataSetReducedToClusterCenters;
    class DataSetReducedToKMeans;
    class DataSetReducedToPrincipalComponents;
    class DataSetScore;
    class DataSetSelectColumns;
    class DataSetSqrt;
    class DataSetStatistics;
    class DecisionTree;
    class DescriptorSelectionBackwardElimination;
    class DescriptorSelectionExhaustive;
    class DescriptorSelectionFeatureForward;
    class DescriptorSelectionInterface;
    class DtreeBinaryPartition;
    class DtreeDataPartitionFunctionInterface;
    class DtreeGiniIndexDataPartitionFunction;
    class DtreeInformationGainDataPartitionFunction;
    class DtreeRocDataPartitionFunction;
    class DtreeSequenceDataPartitionFunction;
    class FeatureLabelSet;
    class FeatureResultAndState;
    class HasLabelsBase;
    class Interface;
    class InterfaceRetrieveFromFile;
    class InterfaceStoreInFile;
    class KappaNearestNeighbor;
    class KohonenNetworkApplicabilityDomain;
    class KohonenNetworkAverage;
    class KohonenNode;
    class LeverageMatrix;
    class MetaDataStorageFile;
    class MetaDataStorageInterface;
    class Model;
    class MultipleLinearRegression;
    class NeuralNetwork;
    class NeuralNetworkPerturbAttenuate;
    class NeuralNetworkPerturbMaxNorm;
    class NeuralNetworkPerturbationInterface;
    class NeuralNetworkSelectiveBackpropagationAccuracy;
    class NeuralNetworkSelectiveBackpropagationAdaptiveTolerance;
    class NeuralNetworkSelectiveBackpropagationBalanced;
    class NeuralNetworkSelectiveBackpropagationDefault;
    class NeuralNetworkSelectiveBackpropagationHybrid;
    class NeuralNetworkSelectiveBackpropagationInterface;
    class NeuralNetworkSelectiveBackpropagationLeadingSequence;
    class NeuralNetworkSelectiveBackpropagationTolerance;
    class NeuralNetworkUpdateWeightsBoundedSimplePropagation;
    class NeuralNetworkUpdateWeightsInterface;
    class NeuralNetworkUpdateWeightsResilientPropagation;
    class NeuralNetworkUpdateWeightsSimplePropagation;
    class ObjectiveFunctionAccuracy;
    class ObjectiveFunctionAccuracyWithExcludedRange;
    class ObjectiveFunctionAucRocCurve;
    class ObjectiveFunctionBinaryOperation;
    class ObjectiveFunctionBootstrap;
    class ObjectiveFunctionCategoricalMax;
    class ObjectiveFunctionConstant;
    class ObjectiveFunctionContingencyMatrixMeasure;
    class ObjectiveFunctionCutoffFromPercentile;
    class ObjectiveFunctionEnrichment;
    class ObjectiveFunctionEnrichmentAverage;
    class ObjectiveFunctionInformationGainRatio;
    class ObjectiveFunctionIntegralPrecisionFractionPredicted;
    class ObjectiveFunctionIntegralTnrTpr;
    class ObjectiveFunctionInterface;
    class ObjectiveFunctionMae;
    class ObjectiveFunctionPartial;
    class ObjectiveFunctionRmsd;
    class ObjectiveFunctionSegmentOverlap;
    class ObjectiveFunctionSelective;
    class ObjectiveFunctionWrapper;
    class PretrainNeuralNetworkFromFile;
    class PretrainNeuralNetworkInterface;
    class PretrainStackedAutoEncoder;
    class RescaleFeatureDataSet;
    class RestrictedBoltzmannMachineLayer;
    class RetrieveDataSetBalanced;
    class RetrieveDataSetBase;
    class RetrieveDataSetBootstrap;
    class RetrieveDataSetByFeature;
    class RetrieveDataSetById;
    class RetrieveDataSetByResult;
    class RetrieveDataSetChunk;
    class RetrieveDataSetCombined;
    class RetrieveDataSetEncodedByModel;
    class RetrieveDataSetFromDelimitedFile;
    class RetrieveDataSetFromFile;
    class RetrieveDataSetJoin;
    class RetrieveDataSetRandomized;
    class RetrieveDataSetRescaled;
    class RetrieveDataSetRows;
    class RetrieveDataSetYscramble;
    class RetrieveDatasetSubset;
    class RetrieveInterface;
    class ScoreDatasetBinaryOperation;
    class ScoreDatasetFScore;
    class ScoreDatasetInputSensitivity;
    class ScoreDatasetInputSensitivityDiscrete;
    class ScoreDatasetInterface;
    class ScoreDatasetNeuralNetworkInputSensitivity;
    class ScoreDatasetNeuralNetworkWeights;
    class ScoreDatasetNonRedundant;
    class ScoreDatasetPartition;
    class ScoreDatasetPearsonCorrelation;
    class ScoreDerivativeEnsemble;
    class StoreInterface;
    class SupportVectorKernelBase;
    class SupportVectorKernelPolynomial;
    class SupportVectorKernelRBF;
    class SupportVectorMachine;
    class SupportVectorMachineMultiOutput;
    class TrainRestrictedBoltzmannMachineLayer;
    class TrainingSchedule;
    class TransferFunctionInterface;
    class TransferGaussian;
    class TransferLinear;
    class TransferRectifier;
    class TransferSigmoid;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_DataType>
    class FeatureDataReference;

    template< typename t_DataType>
    class FeatureDataSet;

    template< typename t_DataType>
    class FeatureDataSetInterface;

    template< typename t_DataType>
    class FeatureInterface;

    template< typename t_DataType>
    class FeatureReference;

    template< typename t_DataType>
    class FeatureSimilarityMeasuresInterface;

  //////////////
  // typedefs //
  //////////////

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FWD_HH_
