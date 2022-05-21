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

#ifndef BCL_BIOL_FWD_HH_
#define BCL_BIOL_FWD_HH_

// include the dependency file for this header
#include "bcl_biol.depends.fwd.hh"

// This file contains forward declarations for the biol namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace biol
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AA;
    class AABackBone;
    class AABackBoneCompleter;
    class AABase;
    class AACaCb;
    class AAClasses;
    class AACompareBySeqID;
    class AACompareData;
    class AACompareDataPtr;
    class AAComplete;
    class AAData;
    class AALessThanSeqID;
    class AASequence;
    class AASequenceFactory;
    class AASequenceFlexibility;
    class AASequencePhiPsi;
    class AASideChainFactory;
    class AATypeData;
    class AATypes;
    class AlignByAAData;
    class AlignByPdbID;
    class Atom;
    class AtomGroupTypeData;
    class AtomGroupTypes;
    class AtomTypeData;
    class AtomTypes;
    class BlastProfile;
    class BlastProfileHandler;
    class ChiAngle;
    class DSSP;
    class EnvironmentTypeData;
    class EnvironmentTypes;
    class ExposurePrediction;
    class Membrane;
    class Mutation;
    class ProteinCharge;
    class ProteinMutationSet;
    class ProteinParams;
    class Ramachandran;
    class Rotamer;
    class RotamerLibrary;
    class SSTypeData;
    class SSTypes;
    class SasaData;
    class SasaPoint;

  //////////////////////
  // template classes //
  //////////////////////

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< util::ShPtr< AABase>, AAClasses>       AAClass;
    typedef util::Enum< AATypeData, AATypes>                   AAType;
    typedef util::Enum< AtomGroupTypeData, AtomGroupTypes>     AtomGroupType;
    typedef util::Enum< AtomTypeData, AtomTypes>               AtomType;
    typedef util::Enum< EnvironmentTypeData, EnvironmentTypes> EnvironmentType;
    typedef util::Enum< SSTypeData, SSTypes>                   SSType;

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_FWD_HH_
