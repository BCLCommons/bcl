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

#ifndef BCL_CHEMISTRY_BOND_ISOMETRY_H_
#define BCL_CHEMISTRY_BOND_ISOMETRY_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_chemistry_bond_isometry.h
  //! @brief enumerates isometry states for bonds
  //!
  //! @remarks example unnecessary
  //! @author mendenjl
  //! @date Feb 28, 2012
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace chemistry
  {

    //! Isometry of the bond type
    enum BondIsometry
    {
      e_NonIsometric,      //!< A bond that does not give rise to isometry
      e_EIsometry,         //!< A non-rotable bond with higher priority substituents on opposite sides
      e_ZIsometry,         //!< A non-rotable bond with higher priority substituents on same sides
      e_UnknownIsometry,   //!< Isometry that is unknown; may be any of the other types
      s_NumberOfIsometries
    };

    //! @brief Isometry as string
    //! @param ISOMETRY the Isometry whose name is desired
    //! @return the name as string
    BCL_API const std::string &GetIsometryName( const BondIsometry &ISOMETRY);

    //! IsometryEnum simplifies the usage of the Isometry enum of this class
    typedef util::WrapperEnum< BondIsometry, &GetIsometryName, s_NumberOfIsometries> BondIsometryEnum;

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_BOND_ISOMETRY_H_ 
