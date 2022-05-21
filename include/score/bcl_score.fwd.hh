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

#ifndef BCL_SCORE_FWD_HH_
#define BCL_SCORE_FWD_HH_

// include the dependency file for this header
#include "bcl_score.depends.fwd.hh"

// This file contains forward declarations for the score namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace score
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AAAssignmentBLOSUM;
    class AAAssignmentBlastProfile;
    class AAAssignmentIdentity;
    class AAAssignmentMeanSimilarityMatrix;
    class AAAssignmentPAM;
    class AAAssignmentPHAT;
    class AAAssignmentProperty;
    class AAAssignmentSSPrediction;
    class AAAssignments;
    class AANeighborhoodDistances;
    class AANeighborhoodExposure;
    class AANeighborhoodExposurePrediction;
    class AANeighborhoodInterface;
    class AAPairAtomClash;
    class AAPairClash;
    class AAPairContact;
    class AAPairContactEnergy;
    class AAPairDistance;
    class AAPairDistanceFittedFunction;
    class AAPairDistanceInterface;
    class AAPairDistanceSmooth;
    class AAPairHiResClash;
    class AAPairSidechainInteraction;
    class AASequence;
    class AASequencePair;
    class Accessibility;
    class AccessibilityHydrophobicMoment;
    class AccessibilityHydrophobicMomentMagnitude;
    class AlignmentQuality;
    class AssignmentGapSimple;
    class BodyAssignment;
    class BodyConnectivityDensity;
    class BodyExtentAgreement;
    class BodyExtentPositionAgreement;
    class ConsensusEnrichment;
    class ContactOrder;
    class DataSetPairwiseBipolar;
    class DataSetPairwiseCoordinateExclusion;
    class DataSetPairwiseCoordinateTriangulation;
    class DataSetPairwiseDataDensity;
    class DataSetPairwiseDistanceChangeMagnitude;
    class DataSetPairwiseEuclidianDistance;
    class DataSetPairwiseResidueTypeExclusion;
    class DataSetPairwiseSSECenter;
    class DataSetPairwiseSSEConnection;
    class DataSetPairwiseSSESize;
    class DataSetPairwiseSSETerm;
    class DataSetPairwiseSequenceSeparation;
    class DataSetPairwiseSize;
    class DataSetPairwiseStructuralExposure;
    class DensityProfileSSEAgreement;
    class EPRAccessibility;
    class EnergyDistribution;
    class EnvironmentPredictions;
    class FuzzyLogicFilter;
    class LogNormalDistribution;
    class Loop;
    class LoopAngle;
    class LoopClosure;
    class PhiPsi;
    class PhiPsiWithSSPred;
    class PofR;
    class ProteinAtomDensity;
    class ProteinModel;
    class ProteinModelAANeighborhood;
    class ProteinModelAANeighborhoodDocking;
    class ProteinModelCompleteness;
    class ProteinModelDefinedLoops;
    class ProteinModelFragmentTopology;
    class ProteinModelGap;
    class ProteinModelInverted;
    class ProteinModelLessThan;
    class ProteinModelLoopDomainClosure;
    class ProteinModelMembraneTopology;
    class ProteinModelSSE;
    class ProteinModelSSEChirality;
    class ProteinModelSSECompleteness;
    class ProteinModelSSELinearLoopProximity;
    class ProteinModelSSENeighbors;
    class ProteinModelSSEPacking;
    class ProteinModelSSEPairs;
    class ProteinModelScoreSum;
    class ProteinModelTopology;
    class ProteinModelWrapper;
    class RadiusOfGyration;
    class ReadHistograms;
    class ResidualDipolarCouplingHistogram;
    class ResidualDipolarCouplingQValue;
    class RestraintAtomAttraction;
    class RestraintAtomDistance;
    class RestraintBodyProteinModel;
    class RestraintDistanceEPR;
    class RestraintDistanceSpinLabel;
    class RestraintEnergyWell;
    class RestraintNMRDistanceInterface;
    class RestraintNoeAttraction;
    class RestraintNoeKnowledgeBased;
    class RestraintPofr;
    class RestraintResidualDipolarCoupling;
    class RestraintSaxs;
    class RestraintXlink;
    class SSEMembraneAlignment;
    class SSEPackInterface;
    class SSEPairAngleDistance;
    class SSEPairClash;
    class SSEPairConnectivity;
    class SSEPairContact;
    class SSEPairGap;
    class SSEPairPacking;
    class SSEPairsFragments;
    class SSEPoolSSEs;
    class SSEPredictionInterface;
    class SSEPredictions;
    class SasType;
    class Score;
    class Scores;
    class StrandPairing;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_Member>
    class AlignmentAssignment;

    template< typename t_Member>
    class AlignmentGap;

    template< typename t_Member>
    class Assignment;

    template< typename t_Member>
    class AssignmentWithGap;

    template< typename t_ArgumentType>
    class Symmetry;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >, AAAssignments> AAAssignment;
    typedef math::FunctionInterfaceSerializable< restraint::AtomDistanceAssignment, double>                                       RestraintAtomDistanceAssignment;

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_FWD_HH_
