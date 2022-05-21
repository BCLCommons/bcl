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

#ifndef BCL_SDF_FWD_HH_
#define BCL_SDF_FWD_HH_

// include the dependency file for this header
#include "bcl_sdf.depends.fwd.hh"

// This file contains forward declarations for the sdf namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace sdf
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AtomInfo;
    class BondInfo;
    class CTabHandler;
    class Factory;
    class FragmentFactory;
    class MdlEntryTypeData;
    class MdlEntryTypes;
    class MdlHandler;
    class MdlHeader;
    class MdlProperty;
    class MolfileHandler;
    class RXNFactory;
    class RXNHandler;

  //////////////////////
  // template classes //
  //////////////////////

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< MdlEntryTypeData, MdlEntryTypes> MdlEntryType;

  } // namespace sdf
} // namespace bcl

#endif // BCL_SDF_FWD_HH_
