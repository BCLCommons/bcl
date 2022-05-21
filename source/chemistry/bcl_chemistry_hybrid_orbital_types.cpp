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
#include "chemistry/bcl_chemistry_hybrid_orbital_types.h"
#include "util/bcl_util_enumerate.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    HybridOrbitalTypes::HybridOrbitalTypes() :
      util::Enumerate< HybridOrbitalTypeData, HybridOrbitalTypes>( false),
      e_Unhybridized( AddEnum( "Unhybridized", HybridOrbitalTypeData())),
      e_SP( AddEnum( "SP", HybridOrbitalTypeData( storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px)))),
      e_SP2( AddEnum( "SP2", HybridOrbitalTypeData( storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py)))),
      e_SP3( AddEnum( "SP3", HybridOrbitalTypeData( storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py, e_Pz))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &HybridOrbitalTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    const HybridOrbitalTypes &GetHybridOrbitalTypes()
    {
      return HybridOrbitalTypes::GetEnums();
    }

  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< chemistry::HybridOrbitalTypeData, chemistry::HybridOrbitalTypes>;

  } // namespace util
} // namespace bcl
