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

#ifndef BCL_CHEMISTRY_ATOMIC_ORBITAL_TYPES_H_
#define BCL_CHEMISTRY_ATOMIC_ORBITAL_TYPES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_chemistry_atomic_orbital_types.h
  //! @brief enumerates basic atom orbital types without PrincipalQuantumNumber
  //!
  //! @remarks example unnecessary
  //! @author mueller
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace chemistry
  {

    //! @enum AtomicOrbitalTypes
    //! @brief enumerates basic atom orbital types without PrincipalQuantumNumber
    enum AtomicOrbitalTypes
    {
      e_S,
      e_Px,
      e_Py,
      e_Pz,
      e_Dxy,
      e_Dxz,
      e_Dyz,
      e_Dz2,
      e_Dx2y2,
      g_NumberOfAtomicOrbitalTypes
    };

    //! @brief GetEnumDescriptor provides the name of ENUM
    //! @param ENUM - the enum for which a name is desired
    BCL_API const std::string &GetAtomicOrbitalTypeString( const AtomicOrbitalTypes &ENUM);

    // Typedef of AtomicOrbitalTypesEnum to insulate the user from details of the enum itself
    typedef util::WrapperEnum< AtomicOrbitalTypes, &GetAtomicOrbitalTypeString, g_NumberOfAtomicOrbitalTypes>
      AtomicOrbitalTypesEnum;

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOMIC_ORBITAL_TYPES_H_
