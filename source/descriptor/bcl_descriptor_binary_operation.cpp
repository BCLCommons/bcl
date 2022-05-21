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
#include "descriptor/bcl_descriptor_binary_operation.hpp"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "math/bcl_math_assignments.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////
  // data //
  //////////

    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> BinaryOperation< t_DataType>::s_Instance
    (
      math::Assignments< float>::AddInstances< Base< t_DataType, float>, BinaryOperation< t_DataType> >()
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API BinaryOperation< chemistry::AtomConformationalInterface>;
    template class BCL_API BinaryOperation< biol::AABase>;
    template class BCL_API BinaryOperation< biol::Mutation>;
    template class BCL_API BinaryOperation< char>;

  } // namespace descriptor
} // namespace bcl
