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

#ifndef BCL_SSPRED_FWD_HH_
#define BCL_SSPRED_FWD_HH_

// include the dependency file for this header
#include "bcl_sspred.depends.fwd.hh"

// This file contains forward declarations for the sspred namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace sspred
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class B2TMPRED;
    class BOCTOPUS;
    class CIPhiPsi;
    class CONPRED;
    class Dssp;
    class DsspStride;
    class JUFO;
    class JUFO9D;
    class JUFOANN;
    class Kaksi;
    class MASP;
    class Mahssmi;
    class MethodHandler;
    class MethodInterface;
    class Methods;
    class OCTOPUS;
    class PARTIFOLD;
    class PDB;
    class PROFTMB;
    class PROFphd;
    class PSIPRED;
    class Palsse;
    class SAM;
    class SSEFactoryHighest;
    class SSEFactoryThreshold;
    class Stride;
    class TALOS;
    class TMBETANET;
    class TMHMM;
    class TMMOD;

  //////////////////////
  // template classes //
  //////////////////////

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< util::ShPtr< MethodInterface>, Methods> Method;

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_FWD_HH_
