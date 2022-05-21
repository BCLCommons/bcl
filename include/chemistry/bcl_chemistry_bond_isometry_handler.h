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

#ifndef BCL_CHEMISTRY_BOND_ISOMETRY_HANDLER_H_
#define BCL_CHEMISTRY_BOND_ISOMETRY_HANDLER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BondIsometryHandler
    //! @brief assigns isometry E/Z for each isometric bond
    //! Known limitation:
    //!  Emergent isometric bonds are not considered by operator <
    //!    Emergent isometric bonds are formed when two otherwise identical substituents differ only in the chirality of
    //!    their substituent atoms or stereogenicity of their bonds
    //!
    //! @see @link example_chemistry_bond_isometry_handler.cpp @endlink
    //! @author mendenjl
    //! @date Mar 18, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BondIsometryHandler
    {
    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief Add E/Z isometry information to the bonds, given a conformation
      //! @param FRAGMENT Conformation upon which to add chirality information
      //! @param FORCE_RECALCULATION whether to relculate even if the bond isometry is already given
      static void AddIsometryInformation( AtomVector< AtomComplete> &FRAGMENT, const bool &FORCE_RECALCULATION = false);

    }; // class BondIsometryHandler

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_BOND_ISOMETRY_HANDLER_H_

