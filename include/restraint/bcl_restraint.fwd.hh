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

#ifndef BCL_RESTRAINT_FWD_HH_
#define BCL_RESTRAINT_FWD_HH_

// include the dependency file for this header
#include "bcl_restraint.depends.fwd.hh"

// This file contains forward declarations for the restraint namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace restraint
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AccessibilityAA;
    class AccessibilityAAAssignment;
    class AccessibilityProfile;
    class AccessibilityProfileAssignment;
    class AnalyzeAccessibilityChange;
    class AnalyzeAtomDistanceHeatmap;
    class AnalyzeAtomDistanceMeanSD;
    class AnalyzeAtomDistancePymol;
    class AnalyzeAtomDistanceScore;
    class AnalyzeAtomDistanceScoreHeatmap;
    class AnalyzeCoordinateDistanceDistribution;
    class AnalyzePerResidueRMSD;
    class AnalyzePerResidueRMSDBetweenEnsembles;
    class AnalyzeSas;
    class AtomDistance;
    class AtomDistanceAssignment;
    class Body;
    class ConeModel;
    class ContactData;
    class ContainsBodyOrigin;
    class DataPairwise;
    class DataSetPairwise;
    class Distance;
    class EPRAccessibilityData;
    class EPRDecay;
    class EPRDecaySimulation;
    class EPRDistanceAssignment;
    class EPRDistanceData;
    class HandlerAccessibilityAA;
    class HandlerAtomDistanceAssigned;
    class HandlerBody;
    class HandlerDataSetPairwiseIdentifiers;
    class HandlerDataSetPairwiseInterface;
    class HandlerEPRDecay;
    class HandlerInterface;
    class Interface;
    class LocatorCoordinatesFirstSideChainAtom;
    class LocatorCoordinatesHydrogen;
    class MutateDataSetPairwiseAdd;
    class MutateDataSetPairwiseFilterAAType;
    class MutateDataSetPairwiseFilterCoordinateExclusion;
    class MutateDataSetPairwiseFilterEuclidianDistance;
    class MutateDataSetPairwiseFilterExposure;
    class MutateDataSetPairwiseFilterSSESize;
    class MutateDataSetPairwiseFilterTriangulation;
    class MutateDataSetPairwiseRemove;
    class MutateTransformationMatrix3DNull;
    class MutateTransformationMatrix3DRotate;
    class NOEData;
    class PREData;
    class Piesa;
    class PofrData;
    class RDC;
    class RDCAssignment;
    class RDCData;
    class SasAnalysis;
    class SasDataParameters;
    class SasDebye;
    class SasDebyeInterface;
    class SasDensityData;
    class SasDistanceDensityPoint;
    class SasExperimentalAndCalculatedData;
    class SasExperimentalAndCalculatedDensity;
    class SasOptimization;
    class SasPofR;
    class SasPofRInterface;
    class SasScatteringData;
    class SasScatteringPoint;
    class SasTransformation;
    class SaxsData;
    class SaxsDataReduction;
    class SaxsOptiResult;
    class SaxsOptimization;
    class XlinkData;

  //////////////////////
  // template classes //
  //////////////////////

    template
    <
      typename t_Restraint,
      typename t_GroupIdentifier,
      typename t_GroupMember,
      typename t_GroupIdentifierCompare = std::less< t_GroupIdentifier>
    >
    class Assignment;

    template< typename t_DataType>
    class Group;

    template
    <
      typename t_GroupIdentifier,
      typename t_GroupMember,
      typename t_GroupIdentifierCompare = std::less< t_GroupIdentifier>
    >
    class GroupCollection;

    template< typename t_RestraintType>
    class HandlerBase;

    template
    <
      typename t_Argument,
      typename t_Restraint,
      typename t_GroupIdentifier,
      typename t_GroupMember,
      typename t_GroupIdentifierCompare = std::less< t_GroupIdentifier>
    >
    class RestraintInterface;

  //////////////
  // typedefs //
  //////////////

    typedef Assignment< util::ShPtrVector< assemble::SSEGeometryInterface>, size_t, assemble::SSE> SSEAssignment;
    typedef util::Implementation< Interface>                                                       Type;

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_FWD_HH_
