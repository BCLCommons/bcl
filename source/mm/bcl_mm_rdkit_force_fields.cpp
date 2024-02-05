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
#include "mm/bcl_mm_rdkit_force_fields.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace mm
  {
    //! @brief RdkitForceFields as string
    //! @param FF the force field desired
    //! @return the force field as string
    const std::string &GetRdkitForceFieldsName( const RdkitForceFields &FF)
    {
      static const std::string s_names[] =
      {
        "UFF",
        "MMFF94",
        "MMFF94s",
        GetStaticClassName< RdkitForceFields>()
      };
      return s_names[ FF];
    }

  } // namespace mm
} // namespace bcl
