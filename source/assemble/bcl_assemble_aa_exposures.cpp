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
#include "assemble/bcl_assemble_aa_exposures.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_aa_sasa_ols.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    //! @brief default constructor
    AAExposures::AAExposures() :
      e_NeighborCount ( AddEnum( "NeighborVector" , util::ShPtr< AAExposureInterface>( new AANeighborCount()))),
      e_NeighborVector( AddEnum( "NeighborCount"  , util::ShPtr< AAExposureInterface>( new AANeighborVector()))),
      e_AASasaOls     ( AddEnum( "SasaOls"        , util::ShPtr< AAExposureInterface>( new AASasaOLS())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAExposures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief construct on access function for all AAExposures
    //! @return reference to only instances of AAExposures
    AAExposures &GetAAExposures()
    {
      return AAExposures::GetEnums();
    }

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< assemble::AAExposureInterface>, assemble::AAExposures>;

  } // namespace util
} // namespace bcl
