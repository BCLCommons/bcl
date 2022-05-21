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

#ifndef BCL_CHEMISTRY_HYBRID_ORBITAL_TYPES_H_
#define BCL_CHEMISTRY_HYBRID_ORBITAL_TYPES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_hybrid_orbital_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HybridOrbitalTypes
    //! @brief enumerates common hybrid orbital types
    //! @details Collection of common orbital hybrids ( sp, sp2, ...)
    //!
    //! @see @link example_chemistry_hybrid_orbital_types.cpp @endlink
    //! @author mueller, woetzen
    //! @date Aug 23, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HybridOrbitalTypes :
      public util::Enumerate< HybridOrbitalTypeData, HybridOrbitalTypes>
    {
      friend class util::Enumerate< HybridOrbitalTypeData, HybridOrbitalTypes>;

    public:

    //////////
    // data //
    //////////

      const HybridOrbitalType e_Unhybridized; //!< unhybridized - no di, tr, or te
      const HybridOrbitalType e_SP;  //!< sp-hybrid  - aka di : digonal
      const HybridOrbitalType e_SP2; //!< sp2-hybrid - aka tr : trigonal planar
      const HybridOrbitalType e_SP3; //!< sp3-hybrid - aka te : tetrahedral

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      HybridOrbitalTypes();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    }; // class HybridOrbitalTypes

    BCL_API const HybridOrbitalTypes &GetHybridOrbitalTypes();

  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< chemistry::HybridOrbitalTypeData, chemistry::HybridOrbitalTypes>;

  } // namespace util
} // namespace bcl

#endif // BCL_CHEMISTRY_HYBRID_ORBITAL_TYPES_H_
