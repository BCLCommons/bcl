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
#include "descriptor/bcl_descriptor_within_range.hpp"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> WithinRange< t_DataType>::s_Instance
    (
      util::Enumerated< Base< t_DataType, float> >::AddInstance
      (
        new WithinRange< t_DataType>()
      )
    );

    template class BCL_API WithinRange< chemistry::AtomConformationalInterface>;
    template class BCL_API WithinRange< biol::AABase>;
    template class BCL_API WithinRange< biol::Mutation>;
    template class BCL_API WithinRange< char>;

  } // namespace descriptor
} // namespace bcl
