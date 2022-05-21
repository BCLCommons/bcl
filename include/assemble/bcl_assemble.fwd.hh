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

#ifndef BCL_ASSEMBLE_FWD_HH_
#define BCL_ASSEMBLE_FWD_HH_

// include the dependency file for this header
#include "bcl_assemble.depends.fwd.hh"

// This file contains forward declarations for the assemble namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace assemble
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class  AAExposureInterface;
    class  AAExposures;
    class  AANeighborCount;
    class  AANeighborList;
    class  AANeighborListContainer;
    class  AANeighborListContainerGeneratorProteinModel;
    class  AANeighborListContainerGeneratorSSE;
    class  AANeighborListContainerGeneratorSSEPair;
    class  AANeighborListContainerPruner;
    class  AANeighborVector;
    class  AASasaOLS;
    class  AnalyzeChiAnglePairDistribution;
    class  AnalyzeChiAngleRecovery;
    class  AnalyzeProteinEnsembleAANeighborhood;
    class  AnalyzeProteinEnsembleInterface;
    class  Biomolecule;
    class  Chain;
    struct ChainLessThan;
    class  ChainMultiplier;
    struct ChainMultiplierLessThan;
    class  CollectorAASpecified;
    class  CollectorAAType;
    class  CollectorAllPossibleDomains;
    class  CollectorCommonAA;
    class  CollectorProteinModelConformationByScore;
    class  CollectorProteinModelConformations;
    class  CollectorSSE;
    class  CollectorSSEPaired;
    class  CollectorSSESize;
    class  CollectorSSEUnpaired;
    class  CollectorSheet;
    class  CollectorTopologyCombined;
    class  CollectorTopologyInterface;
    class  CollectorTopologySheet;
    class  Domain;
    class  DomainInterface;
    class  FoldTemplate;
    class  FoldTemplateHandler;
    class  LocatorAA;
    class  LocatorAtom;
    class  LocatorAtomCoordinatesInterface;
    class  LocatorChain;
    class  LocatorDomainRandom;
    class  LocatorDomainSSEPoolOverlapping;
    class  LocatorDomainSpecified;
    class  LocatorSSE;
    class  LocatorSSEFromProteinModelData;
    class  LocatorSSERandom;
    class  LocatorSSETerminusResidue;
    class  LocatorSSEUnpaired;
    class  LocatorSSEsRandom;
    class  LocatorSubDomainRandom;
    class  PickProteinModelConformationRandom;
    class  PickSSEFurthestEuclidean;
    class  PickSSERandom;
    class  PickSSEShortLoops;
    class  PickSSEsRandom;
    class  PrinterProteinModel;
    class  PrinterProteinModelEnsemble;
    class  PrinterProteinModelMovie;
    class  PrinterProteinModelMultimer;
    class  PrinterTrackerHistory;
    class  ProteinEnsemble;
    class  ProteinModel;
    class  ProteinModelData;
    class  ProteinModelInverter;
    class  ProteinModelMomentOfInertia;
    class  ProteinModelMultiplier;
    class  ProteinModelWithCache;
    class  ProteinModelWithMutations;
    class  ProteinStorageFile;
    class  ProteinWithCacheDatasetFromFile;
    class  ProteinWithCacheStorageFile;
    class  ProteinWithMutationsDatasetFromFile;
    class  Quality;
    class  QualityBatch;
    class  SSE;
    class  SSECompare;
    class  SSECompareByIdentity;
    class  SSECompareExtent;
    class  SSECompareOverlap;
    class  SSECompareType;
    class  SSEFactories;
    class  SSEFactoryConformation;
    class  SSEFactoryInterface;
    class  SSEFactoryMC;
    class  SSEGeometry;
    class  SSEGeometryInterface;
    class  SSEGeometryInterfaceLessThan;
    class  SSEGeometryPackerAllFragmentPairs;
    class  SSEGeometryPackerBestFragmentPair;
    class  SSEGeometryPackerBestFragmentPairs;
    class  SSEGeometryPackers;
    class  SSEGeometryPacking;
    class  SSEGeometryPackingCompareDistance;
    class  SSEGeometryPackingCompareInteractionWeight;
    class  SSEGeometryPackingCriteriaCombine;
    class  SSEGeometryPackingCriteriaContactType;
    class  SSEGeometryPackingCriteriaDistance;
    class  SSEGeometryPackingCriteriaDistancePerType;
    class  SSEGeometryPackingCriteriaInteractionWeight;
    class  SSEGeometryPackingCriteriaStrandWeight;
    class  SSEGeometryPackingListPickers;
    class  SSEGeometryPackingPickers;
    class  SSEGeometryPhiPsi;
    class  SSEGeometryPhiPsiLessThan;
    class  SSEGeometryWithinSizeTolerance;
    class  SSELessThan;
    class  SSELessThanBySize;
    class  SSELessThanNoOverlap;
    class  SSEPairTemplate;
    class  SSEPool;
    class  SSEPoolAgreement;
    class  SSEPoolInsertCoilIntoSSE;
    class  SSEPoolJoinSSEs;
    class  SSEPoolMoveAA;
    class  SSEPoolMutateSSE;
    class  SSEPoolSplitSSE;
    class  SSETransformer;
    class  SheetTemplateHandler;
    class  Topology;
    class  TopologyDistance;
    class  VoxelGridAA;
    class  VoxelGridAtom;
    class  VoxelGridMutation;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_DataType>
    class Ensemble;

    template< typename t_CriteriaType>
    class LocatorSSEFurthest;

    template< typename t_CriteriaType>
    class PickSSEFurthestEuclideanCenter;

    template< typename t_ReturnType>
    class SSEGeometryPackerCacheWrapper;

    template< typename t_ReturnType>
    class SSEGeometryPackerInterface;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< util::ShPtr< AAExposureInterface>, AAExposures>                                                                         AAExposure;
    typedef io::RetrieveInterface< util::ShPtr< ProteinModelWithCache>, util::ShPtrVector< ProteinModelWithCache> >                             RetrieveProteinModelWithCache;
    typedef util::Enum< util::ShPtr< SSEFactoryInterface>, SSEFactories>                                                                        SSEFactory;
    typedef util::Enum< util::ShPtr< SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > > >, SSEGeometryPackers> SSEGeometryPacker;
    typedef util::Enum< util::ShPtr< SSEGeometryPackerInterface< storage::List< SSEGeometryPacking> > >, SSEGeometryPackingListPickers>         SSEGeometryPackingListPicker;
    typedef util::Enum< util::ShPtr< SSEGeometryPackerInterface< SSEGeometryPacking> >, SSEGeometryPackingPickers>                              SSEGeometryPackingPicker;

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_FWD_HH_
