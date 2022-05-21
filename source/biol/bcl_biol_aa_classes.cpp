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
#include "biol/bcl_biol_aa_classes.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_back_bone.h"
#include "biol/bcl_biol_aa_complete.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct all AAClasses
    AAClasses::AAClasses() :
      e_AA(         AddEnum( "AA",         util::ShPtr< AABase>( new AA()))),
      e_AACaCb(     AddEnum( "AACaCb",     util::ShPtr< AABase>( new AACaCb()))),
      e_AABackBone( AddEnum( "AABackBone", util::ShPtr< AABase>( new AABackBone()))),
      e_AAComplete( AddEnum( "AAComplete", util::ShPtr< AABase>( new AAComplete())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAClasses::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief construct on access function for all AAClasses
    //! @return reference to only instances of AAClasses
    const AAClasses &GetAAClasses()
    {
      return AAClasses::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< biol::AABase>, biol::AAClasses>;

  } // namespace util
} // namespace bcl
