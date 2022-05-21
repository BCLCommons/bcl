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

#ifndef BCL_DESCRIPTOR_FWD_HH_
#define BCL_DESCRIPTOR_FWD_HH_

// include the dependency file for this header
#include "bcl_descriptor.depends.fwd.hh"

// This file contains forward declarations for the descriptor namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace descriptor
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AAAtomPosition;
    class AAAtomProperty;
    class AABlastProfile;
    class AABlastProfileEntropy;
    class AABlastWeightedProperty;
    class AADSSPInfo;
    class AADistanceFromAAType;
    class AADistanceToTMCenter;
    class AAHbondNeighbor;
    class AAIdentifier;
    class AAInfoFromSymMatrixFile;
    class AAInsideOutsideMembrane;
    class AALetterCode;
    class AAMoleculeProperty;
    class AAPairDistance;
    class AAPairProbability;
    class AAPairSSProbability;
    class AAPairSeparation;
    class AAPairType;
    class AAPhiPsi;
    class AAPoreOrientation;
    class AAPositionInTM;
    class AAProperty;
    class AASSEInfo;
    class AASSTMPrediction;
    class AASasa;
    class AASeqID;
    class AATMDirection;
    class AATripletHelixLocationType;
    class AATripletLookup;
    class AATripletType;
    class AtomAromaticityAxes;
    class AtomEffectivePolarizability;
    class AtomEstimatedSurfaceArea;
    class AtomFormalCharge;
    class AtomHBondInfo;
    class AtomIsSP;
    class AtomIsSP2;
    class AtomIsSP3;
    class AtomLonePairEN;
    class AtomMiscProperty;
    class AtomNeighborDirection;
    class AtomNumberValences;
    class AtomPiCharge;
    class AtomPolarizability;
    class AtomRelativePropertyScore;
    class AtomRingSize;
    class AtomSigmaCharge;
    class AtomStereocenters;
    class AtomSurfaceArea;
    class AtomTopologicalPolarSurfaceArea;
    class AtomTypeNumber;
    class AtomTypePropertyRetriever;
    class AtomVcharge;
    class AtomVolume;
    class AtomicNumber;
    class BondTypeCount;
    class BuserMetric;
    class CacheMap;
    class Central2DASign;
    class CheminfoProperties;
    class CoulombicForce;
    class Dataset;
    class ElementTypePropertyRetriever;
    class FourierAnalysisWindowCreator;
    class MACCS;
    class MolAlignPharmScore;
    class Molecule2DACode;
    class Molecule2DAMaxMin;
    class Molecule2DAMaxSign;
    class Molecule2DASmoothSignCode;
    class Molecule3DAClosestPairRealSpace;
    class Molecule3DAClosestPairRealSpaceAsymmetry;
    class Molecule3DACode;
    class Molecule3DAPairConvolution;
    class Molecule3DAPairConvolutionAsymmetry;
    class Molecule3DAPairRealSpace;
    class Molecule3DAPairRealSpaceAsymmetry;
    class Molecule3DAPairRealSpaceConvolution;
    class Molecule3DAPairRealSpaceConvolutionAsymmetry;
    class Molecule3DASmooth;
    class Molecule3DASmoothSignCode;
    class Molecule3DASmoothSignOcclusionCode;
    class Molecule3DASoftMax;
    class Molecule3DASoftMaxSign;
    class Molecule3DASoftMin;
    class Molecule3DDistribution;
    class Molecule3DSignDistribution;
    class MoleculeAsymmetry;
    class MoleculeAtomEnvironmentMap;
    class MoleculeCachedString;
    class MoleculeComplexity;
    class MoleculeDefault;
    class MoleculeDruglike;
    class MoleculeEntropyQHA;
    class MoleculeFragmentRescale;
    class MoleculeFragmentStatistics;
    class MoleculeGirth;
    class MoleculeHalogenatedAromaticRings;
    class MoleculeInterHBondCode;
    class MoleculeLipinskiViolations;
    class MoleculeLogP;
    class MoleculeLogP2008;
    class MoleculeMaximumFragmentStatistics;
    class MoleculeMiscProperty;
    class MoleculeName;
    class MoleculeOneFourClashScore;
    class MoleculeRDFCode;
    class MoleculeRDFGridCode;
    class MoleculeRDFMaxSignCode;
    class MoleculeRDFSignCode;
    class MoleculeRings;
    class MoleculeRotatableBonds;
    class MoleculeShapeMoments;
    class MoleculeSimilarity;
    class MoleculeStorageId;
    class MoleculeTotalBondEnergy;
    class MoleculeTriangulatorCode;
    class MoleculeVdwScore;
    class MoleculeXLogP;
    class MutationAADensity;
    class MutationAAProperty;
    class MutationDensity;
    class MutationId;
    class MutationNearestInSequence;
    class MutationNearestSpatially;
    class PairConvolutionCorrelationDNN;
    class ProteinId;
    class ProteinLigandCorrelationDNN;
    class ReactionStructureSearch;
    class SdfFileId;
    class SegmentFinder;
    class SegmentInfo;
    class SequenceSSECount;
    class StringSequence;
    class StructureCount;
    class StructureSearch;
    class Type;
    class UMol2D;
    class WindowWeightingInterface;
    class WindowWeights;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_DataType>
    class AtomPlanarity;

    template< typename t_DataType>
    class BandPassFilter;

    template< typename t_DataType, typename t_ReturnType>
    class Base;

    template< typename t_DataType, typename t_ReturnType>
    class BaseElement;

    template< typename t_DataType, typename t_ReturnType>
    class BaseElementOrSequence;

    template< typename t_DataType, typename t_ReturnType>
    class BasePair;

    template< typename t_DataType, typename t_ReturnType>
    class BaseSequence;

    template< typename t_DataType, typename t_ReturnType>
    class BaseTriplet;

    template< typename t_DataType>
    class BinaryOperation;

    template< typename t_DataType, typename t_ReturnType>
    class Combine;

    template< typename t_DataType, typename t_ReturnType>
    class Constants;

    template< typename t_DataType>
    class DatasetBuilder;

    template< typename t_DataType, typename t_ReturnType>
    class Define;

    template< typename t_DataType>
    class ElementHistogram1D;

    template< typename t_DataType>
    class ElementHistogram2D;

    template< typename t_DataType>
    class ElementPosition;

    template< typename t_DataType, typename t_ReturnType>
    class ForEach;

    template< typename t_BaseClass>
    class HasCache;

    template< typename t_DataType>
    class IterativePrediction;

    template< typename t_DataType>
    class Iterator;

    template< typename t_DataType>
    class KohonenMapInfo;

    template< typename t_DataType>
    class Limit;

    template< typename t_DataType, typename t_ReturnType>
    class Mapped;

    template< typename t_ReturnType>
    class MappedSequence;

    template< typename t_DataType>
    class MinMaxIndex;

    template< typename t_DataType, typename t_ReturnType>
    class Named;

    template< typename t_DataType, typename t_ReturnType>
    class NamedTemplate;

    template< typename t_DataType>
    class Numeric;

    template< typename t_DataType>
    class NumericString;

    template< typename t_DataType>
    class Offset;

    template< typename t_DataType>
    class OuterProduct;

    template< typename t_DataType, typename t_ReturnType>
    class Partial;

    template< typename t_DataType>
    class Periodogram;

    template< typename t_DataType, typename t_ReturnType>
    class Positional;

    template< typename t_DataType>
    class PowerSpectrum;

    template< typename t_DataType>
    class PowerSpectrumSequenceWidth;

    template< typename t_DataType>
    class Prediction;

    template< typename t_DataType>
    class PredictionInfo;

    template< typename t_DataType>
    class Rank;

    template< typename t_DataType>
    class ReplaceUndefinedValues;

    template< typename t_DataType>
    class Rescale;

    template< typename t_DataType>
    class SequenceHistogram1D;

    template< typename t_DataType>
    class SequenceInterface;

    template< typename t_DataType>
    class SequenceSegmentStatistics;

    template< typename t_DataType, typename t_ReturnType>
    class SequenceSize;

    template
    <
      typename t_DataType,
      typename t_Accumulator,
      const linal::Vector< float> &( t_Accumulator::*GetAccumulatedValues)() const
    >
    class SequenceStatistics;

    template
    <
      typename t_DataType,
      typename t_Accumulator,
      const linal::Vector< float> &( t_Accumulator::*GetAccumulatedValues)() const
    >
    class SequenceWeightedStatistics;

    template< typename t_DataType>
    class Sigmoid;

    template< typename t_DataType>
    class Sort;

    template
    <
      typename t_DataType,
      typename t_Accumulator,
      const float &( t_Accumulator::*GetAccumulatedValues)() const
    >
    class Statistics;

    template< typename t_DataType, typename t_ReturnType>
    class Template;

    template< typename t_DataType>
    class UnaryOperation;

    template< typename t_DataType>
    class UniformRandom;

    template< typename t_DataType, typename t_ReturnType>
    class Window;

    template< typename t_DataType>
    class WindowAverage;

    template< typename t_DataType>
    class WindowConditionalAverage;

    template< typename t_DataType>
    class WindowMinMax;

    template< typename t_DataType, typename t_ReturnType>
    class WindowPairSquare;

    template< typename t_DataType>
    class WindowSegmentStatistics;

    template< typename t_DataType>
    class WithinRange;

    template< typename t_DataType>
    class WithinRangeGaussian;

    template< typename t_DataType>
    class WithinRangeSmooth;

  //////////////
  // typedefs //
  //////////////

    typedef util::Implementation< Base< chemistry::AtomConformationalInterface, char> >  CheminfoID;
    typedef util::Implementation< Base< chemistry::AtomConformationalInterface, float> > CheminfoProperty;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_FWD_HH_
