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

#ifndef BCL_CHEMISTRY_SOLVENTS_H_
#define BCL_CHEMISTRY_SOLVENTS_H_

// include the namespace header
#include "bcl_chemistry.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_chemistry_solvents.h
  //! @brief enumerates common solvents used for nmr on small molecules
  //!
  //! @see @link example_chemistry_solvents.cpp @endlink
  //! @author mueller
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace chemistry
  {
    //! solvent enumeration
    enum Solvents
    {
      e_UNDEFINED_SOLVENT, e_CDCl3, e_DMSO, e_METHANOL_D4, e_CCl4,
      e_ACETONE_D6, e_PYRIDIN_D5, e_BENZENE_D6, e_D2O, e_DMF_D7,
      s_NumberSolvents
    };

    //! check whether a solvent is in a defined state
    BCL_API
    bool IsDefined( const Solvents &SOLVENT);

    //! create solvent from string
    BCL_API Solvents SolventFromString( const std::string &SOLVENT);

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SOLVENTS_H_
