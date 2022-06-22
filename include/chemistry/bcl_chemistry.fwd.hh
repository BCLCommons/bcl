// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
// (c) This file is part of the BCL software suite and is made available under license.
// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
// (c) For academic and non-profit users: 
// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
// (c) For commercial users: 
// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
// (c) For BCL developers at Vanderbilt University: 
// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
// (c)
// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the 
// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
// (c) Vanderbilt University. 
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org 
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)

#ifndef BCL_CHEMISTRY_FWD_HH_
#define BCL_CHEMISTRY_FWD_HH_

// include the dependency file for this header
#include "bcl_chemistry.depends.fwd.hh"

// This file contains forward declarations for the chemistry namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace chemistry
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AAFragmentComplete;
    class AtomClashScore;
    class AtomComplete;
    class AtomConfigurationalInterface;
    class AtomConfigurationalShared;
    class AtomConformationalInterface;
    class AtomConformationalShared;
    class AtomConstitutionalInterface;
    class AtomConstitutionalShared;
    class AtomEnvironment;
    class AtomEnvironmentBender;
    class AtomOneFourInteractionScore;
    class AtomTypeData;
    class AtomTypes;
    class AtomVdwScore;
    class AtomWithPositionInterface;
    class AtomsCompleteStandardizer;
    class BondAngleAssignment;
    class BondConfigurational;
    class BondConformational;
    class BondConstitutional;
    class BondDihedralAngles;
    class BondDihedralStericWeight;
    class BondIsometryHandler;
    class BondLengths;
    class CollectorValence;
    class ConfigurationGraphConverter;
    class ConfigurationInterface;
    class ConfigurationSet;
    class ConfigurationSetSameConstitution;
    class ConfigurationalBondTypeData;
    class ConfigurationalBondTypes;
    class ConformationComparisonByDihedralBins;
    class ConformationComparisonByDihedrals;
    class ConformationComparisonByFragments;
    class ConformationComparisonByProperty;
    class ConformationComparisonByRealSpaceRmsd;
    class ConformationComparisonByRmsd;
    class ConformationComparisonBySubstructure;
    class ConformationComparisonBySymmetryRmsd;
    class ConformationComparisonInterface;
    class ConformationComparisonMolAlignByParts;
    class ConformationComparisonMultiAlign;
    class ConformationComparisonPropertyFieldCorrelation;
    class ConformationComparisonPropertyRMSDX;
    class ConformationComparisonPsiField;
    class ConformationComparisonPsiFlexField;
    class ConformationGraphConverter;
    class ConformationInterface;
    class ConformationSet;
    class ConformationSetSameConfiguration;
    class ConformationSetSameConstitution;
    class ConstitutionGraphConverter;
    class ConstitutionInterface;
    class ConstitutionSet;
    class ConstitutionalBondTypeData;
    class ConstitutionalBondTypes;
    class CoulombicScore;
    class CsCode;
    class DescriptorToScoreAdaptor;
    class ElectronConfiguration;
    class ElementStructureFactor;
    class ElementTypeData;
    class ElementTypes;
    class FragmentAddMedChem;
    class FragmentAlchemy;
    class FragmentAlignToScaffold;
    class FragmentComplete;
    class FragmentConfigurationShared;
    class FragmentConformationShared;
    class FragmentConnect;
    class FragmentConnector;
    class FragmentConstitutionShared;
    class FragmentCyclize;
    class FragmentDockEngine;
    class FragmentEnsemble;
    class FragmentEvolveBase;
    class FragmentEvolveImplementations;
    class FragmentExtendWithLinker;
    class FragmentFeed;
    class FragmentFeedFromFile;
    class FragmentFeedFromStream;
    class FragmentFeedInterface;
    class FragmentFluorinate;
    class FragmentGraphMarker;
    class FragmentGrow;
    class FragmentHalogenate;
    class FragmentMakeConformers;
    class FragmentMapConformer;
    class FragmentMolecule;
    class FragmentMutateInterface;
    class FragmentMutateMCM;
    class FragmentMutateReact;
    class FragmentProbabilityScore;
    class FragmentReact;
    class FragmentRemoveAtom;
    class FragmentRemoveBond;
    class FragmentRingSwap;
    class FragmentSplitCommonSubstructure;
    class FragmentSplitConformations;
    class FragmentSplitECFPFragments;
    class FragmentSplitGADDFragments;
    class FragmentSplitInterface;
    class FragmentSplitIsolate;
    class FragmentSplitLargestComponent;
    class FragmentSplitLinearFragments;
    class FragmentSplitRigid;
    class FragmentSplitRings;
    class FragmentSplitRingsWithUnsaturatedSubstituents;
    class FragmentSplitScaffolds;
    class FragmentSplitUnbridgedRings;
    class FragmentStochasticPoseOptimizer;
    class FragmentTrackMutableAtoms;
    class HybridOrbitalTypeData;
    class HybridOrbitalTypes;
    class HydrogensHandler;
    class LigandDesignHelper;
    class LigandPocketFitScore;
    class MergeFragmentComplete;
    class MinimizeLigandInPocket;
    class MolAlignByParts;
    class MolecularConfigurationShared;
    class MolecularConformationShared;
    class MolecularConstitutionShared;
    class MoleculeComplete;
    class MoleculeEnsemble;
    class MoleculeEnvironment;
    class MoleculeFeatureMapper;
    class MoleculeFragmentRecombination;
    class MoleculeStorageFile;
    class MutateBondAngles;
    class MutateBondLengths;
    class MutateChirality;
    class MutateClashResolver;
    class MutateDihedralBond;
    class MutateDihedralsInterface;
    class MutateFragment;
    class MutateMoleculeGeneric;
    class MutateMultiFragment;
    class NCAAFragmentComplete;
    class PerturbMoleculePose;
    class PharmacophoreMapper;
    class PickAtomByElement;
    class PickAtomRandom;
    class PickFragmentPropertyWeighted;
    class PickFragmentRandom;
    class PossibleAtomTypesForAtom;
    class PriorityDihedralAngles;
    class ReactionComplete;
    class ReactionEnsemble;
    class ReactionSearch;
    class ReactionStructure;
    class ReactionWorker;
    class RingFragmentMap;
    class RotamerClusterCenter;
    class RotamerDihedralBondData;
    class RotamerEnsemble;
    class RotamerLibraryFile;
    class RotamerLibraryInterface;
    class SampleConformations;
    class ScoreFunctionGeneric;
    class SearchFragmentLibraryFromTree;
    class SmallMoleculeFragmentIsomorphism;
    class SmallMoleculeFragmentMapping;
    class SmallMoleculeMiscProperties;
    class SmallMoleculeQsarStorageFile;
    class SmallMoleculeStringPropertiesCached;
    class SmallMoleculeStringPropertiesMapped;
    class SmallMoleculeStringPropertiesName;
    class SmallMoleculeStringPropertiesNumeric;
    class SmallMoleculeStringPropertiesTypes;
    class StereocentersHandler;
    class StringPropertyInterface;
    class SubFragment;
    class SubstituentConformational;
    class ValenceHandler;
    class VoxelGridAtom;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_Atom>
    class AtomVector;

    template< typename t_BaseClass>
    class HasProperties;

    template< typename t_BaseClass>
    class HasPropertiesInterface;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< AtomTypeData, AtomTypes>                               AtomType;
    typedef util::Enum< ConfigurationalBondTypeData, ConfigurationalBondTypes> ConfigurationalBondType;
    typedef util::Enum< ConstitutionalBondTypeData, ConstitutionalBondTypes>   ConstitutionalBondType;
    typedef util::Enum< ElementTypeData, ElementTypes>                         ElementType;
    typedef util::Enum< HybridOrbitalTypeData, HybridOrbitalTypes>             HybridOrbitalType;

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FWD_HH_
