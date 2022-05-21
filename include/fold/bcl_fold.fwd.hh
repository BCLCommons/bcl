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

#ifndef BCL_FOLD_FWD_HH_
#define BCL_FOLD_FWD_HH_

// include the dependency file for this header
#include "bcl_fold.depends.fwd.hh"

// This file contains forward declarations for the fold namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace fold
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AddParabolicLoops;
    class CollectorLoopDomain;
    class CollectorLoopDomainAllNonRigid;
    class CollectorLoopDomainRandom;
    class CollectorUnconnectedSSE;
    class DefaultFlags;
    class DefaultMutates;
    class DefaultScores;
    class HandlerLocatorLoopDomain;
    class HandlerLocatorLoopSegment;
    class LocatorLoop;
    class LocatorLoopDomain;
    class LocatorLoopSegment;
    class LocatorMissingCoordinates;
    class LocatorUnconnectedSegments;
    class LoopDomain;
    class LoopDomainCToN;
    class LoopLibrary;
    class LoopParameters;
    class LoopSegment;
    class LoopSegmentSequenceOrder;
    class MutateAAPhi;
    class MutateAAPsi;
    class MutateAARotate;
    class MutateAASequenceGrow;
    class MutateAASetPhi;
    class MutateAASetPsi;
    class MutateDomainFlip;
    class MutateDomainMergeConsecutiveSSTypes;
    class MutateDomainSSEPairTrim;
    class MutateDomainSSESplit;
    class MutateDomainShuffle;
    class MutateDomainTransformation;
    class MutateLoopDomainDihedral;
    class MutateMembrane;
    class MutateMembraneChainMove;
    class MutateMultimer;
    class MutateProteinEnsembleAdd;
    class MutateProteinEnsembleRemove;
    class MutateProteinModel;
    class MutateProteinModelAddSheetFromTemplate;
    class MutateProteinModelChainMove;
    class MutateProteinModelCompress;
    class MutateProteinModelDomain;
    class MutateProteinModelDomainAdd;
    class MutateProteinModelFilterConformations;
    class MutateProteinModelFixLoopClosureWrapper;
    class MutateProteinModelGrowSSE;
    class MutateProteinModelLoopDomain;
    class MutateProteinModelLoopDomainCCD;
    class MutateProteinModelLoopDomainGrow;
    class MutateProteinModelLoopResize;
    class MutateProteinModelMoveAA;
    class MutateProteinModelMultipleGeometries;
    class MutateProteinModelPairStrands;
    class MutateProteinModelReplicateConformation;
    class MutateProteinModelSSE;
    class MutateProteinModelSSEAdd;
    class MutateProteinModelSSEAddMultiple;
    class MutateProteinModelSSEMove;
    class MutateProteinModelSSEPair;
    class MutateProteinModelSSEPairAlignAndPull;
    class MutateProteinModelSSEPairClash;
    class MutateProteinModelSSEPairFixLoopClosure;
    class MutateProteinModelSSEPairHinge;
    class MutateProteinModelSSERemove;
    class MutateProteinModelSSEResize;
    class MutateProteinModelSSESeed;
    class MutateProteinModelSSESplit;
    class MutateProteinModelSSESwap;
    class MutateProteinModelSSESwapBody;
    class MutateProteinModelSSESwapMultimer;
    class MutateProteinModelSSESwapWithPool;
    class MutateProteinModelSSESwapWithPoolOverlap;
    class MutateProteinModelStrandSwitchSheet;
    class MutateProteinModelSwitchConformation;
    class MutateProteinModelThreadSequence;
    class MutateSSEBendRamachandran;
    class MutateSSEBendRandom;
    class MutateSSEBendTemplate;
    class MutateSSEType;
    class MutateSheetCycle;
    class MutateSheetDivide;
    class MutateSheetFitToTemplate;
    class MutateSheetOrder;
    class MutateSheetRegisterFix;
    class MutateSheetRegisterShift;
    class MutateSheetSort;
    class MutateSheetTwist;
    class MutateTree;
    class Mutates;
    class MutationResidue;
    class PhiPsiGeneratorCCD;
    class PhiPsiGeneratorRamachandran;
    class PlacementDomain;
    class PlacementDomainInterface;
    class PlacementDomainUsingFoldTemplate;
    class PlacementSSEDistanceRestraint;
    class PlacementSSEIntoBody;
    class PlacementSSENextToSSE;
    class PlacementSSEShortLoop;
    class PlacementStrandNextToSheet;
    class ProteinGeometry;
    class ProtocolAssembly;
    class ProtocolCreate;
    class ProtocolDefault;
    class ProtocolDock;
    class ProtocolEM;
    class ProtocolEnsemble;
    class ProtocolEnsembleFilter;
    class ProtocolEnsembleReplicateConformation;
    class ProtocolEnsembleSwitchConformation;
    class ProtocolInterface;
    class ProtocolLoopClose;
    class ProtocolLoopCoordinateAdd;
    class ProtocolMembrane;
    class ProtocolMultimer;
    class ProtocolRefinement;
    class ProtocolRestraint;
    class ProtocolTemplate;
    class Protocols;
    class ScoreWeightSet;
    class Scores;
    class Setup;
    class StageFactory;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionLoopClosure;

    template< typename t_ObjectType, typename t_ArgumentType>
    class PlacementInterface;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< util::ShPtr< math::MutateInterface< assemble::ProteinModel> >, Mutates> Mutate;
    typedef util::Enum< util::SiPtr< ProtocolInterface>, Protocols>                             Protocol;
    typedef util::Enum< util::ShPtr< score::ProteinModel>, Scores>                              Score;

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_FWD_HH_
