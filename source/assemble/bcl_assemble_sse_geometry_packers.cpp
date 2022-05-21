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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_sse_geometry_packers.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_packer_all_fragment_pairs.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief private constructor for constructing all SSEGeometryPackers
    SSEGeometryPackers::SSEGeometryPackers() :
      e_Default
      (
        AddEnum
        (
          "Default",
          SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > >::WrapPacker
          (
            util::ShPtr< SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > > >
            (
              new SSEGeometryPackerAllFragmentPairs
              (
                SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength(), SSEGeometryPackerAllFragmentPairs::GetDistanceCutoffMap()

              )
            ), true, false
          )
        )
      ),
      e_Default_NoCache
      (
        AddEnum
        (
          "Default_NoCache",
          SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > >::WrapPacker
          (
            util::ShPtr< SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > > >
            (
              new SSEGeometryPackerAllFragmentPairs
              (
                SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength(), SSEGeometryPackerAllFragmentPairs::GetDistanceCutoffMap()
              )
            ), false, false
          )
        )
      ),
      e_SSEClash
      (
        AddEnum
        (
          "SSEClash",
          SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > >::WrapPacker
          (
            util::ShPtr< SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > > >
            (
              new SSEGeometryPackerAllFragmentPairs
              (
                0.0, SSEGeometryPackerAllFragmentPairs::GetClashDistanceCutoffMap()
              )
            ), true, false
          )
        )
      ),
      e_SSEClash_NoCache
      (
        AddEnum
        (
          "SSEClash_NoCache",
          SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > >::WrapPacker
          (
            util::ShPtr< SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > > >
            (
              new SSEGeometryPackerAllFragmentPairs
              (
                0.0, SSEGeometryPackerAllFragmentPairs::GetClashDistanceCutoffMap()
              )
            ), false, false
          )
        )
      )
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEGeometryPackers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access static instance of SSEGeometryPackers enums
    //! @return static instance of SSEGeometryPackers enums
    const SSEGeometryPackers &GetSSEGeometryPackers()
    {
      return SSEGeometryPackers::GetEnums();
    }

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< assemble::SSEGeometryPackerInterface< storage::Vector< storage::List< assemble::SSEGeometryPacking> > > >, assemble::SSEGeometryPackers>;

  } // namespace util
} // namespace bcl
