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

#ifndef BCL_OPTI_FWD_HH_
#define BCL_OPTI_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the opti namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace opti
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class EnsembleFilter;
    class TrackerBase;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class ApproximatorEvolution;

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorGoldenSection;

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorModularBase;

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorModularInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorNelderMead;

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorPowell;

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorRootBisect;

    template< typename t_ArgumentResultType>
    class ApproximatorRootNewton;

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorRootRegulaFalsi;

    template< typename t_ArgumentResultType>
    class ApproximatorRootSecant;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionAccepted;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionAll;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionCombine;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionConvergenceArgument;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionConvergenceResult;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionDivergenceArgument;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionElapsedTime;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionFunction;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionImprovementRatio;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionNStep;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionNumberIterations;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionPhase;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionRejected;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionResultChanged;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionResultThreshold;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionSkippedSteps;

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionUnimproved;

    template< typename t_MemberType, typename t_FitnessType>
    class DeterministicTournamentSelector;

    template< typename t_DataType>
    class EnsembleNode;

    template< typename t_MemberType, typename t_FitnessType>
    struct EvolutionMemberCompareLess;

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionMemberUniqueInterface;

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionOperation;

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionOperationSelect;

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionPopulation;

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionPopulationBuilder;

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionPopulationMember;

    template< typename t_ArgumentType>
    class OptimizationIdentity;

    template< typename t_ArgumentType>
    class OptimizationInterface;

    template< typename t_ArgumentType>
    class Pipeline;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrintInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterArgumentToFile;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterDefault;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterWithCriterion;

    template< typename t_MemberType, typename t_FitnessType>
    class ProbabilisticSelection;

    template< typename t_ArgumentType>
    class ProcessorInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class Tracker;

    template< typename t_ArgumentType, typename t_ResultType>
    class TrackerWithHistory;

  //////////////
  // typedefs //
  //////////////

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_FWD_HH_
