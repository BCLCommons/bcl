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

#ifndef BCL_MM_RDKIT_FORCE_FIELDS_H_
#define BCL_MM_RDKIT_FORCE_FIELDS_H_

// include the namespace header
#include "bcl_mm.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_mm_rdkit_force_fields.h
  //! @brief enumerates force fields available through the RDKit library
  //!
  //! @remarks example unnecessary
  //! @author brownbp1
  //! @date Feb 04, 2024
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace mm
  {
    //! Force fields available in RDKit library
    enum RdkitForceFields
    {
      e_UFF,                //!< Universal Force Field
      e_MMFF94,             //!< Merck Molecular Force Field 94
      e_MMFF94s,            //!< Merck Molecular Force Field 94 static
      s_NumberForceFields
    };

    //! @brief RdkitForceFields as string
    //! @param FF the force field desired
    //! @return the force field as string
    BCL_API const std::string &GetRdkitForceFieldsName( const RdkitForceFields &FF);

    //! RdkitForceFieldsEnum simplifies the usage of the RdkitForceFields enum of this class
    typedef util::WrapperEnum< RdkitForceFields, &GetRdkitForceFieldsName, s_NumberForceFields> RdkitForceFieldsEnum;

  } // namespace mm
} // namespace bcl

#endif // BCL_MM_RDKIT_FORCE_FIELDS_H_
