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
#include "descriptor/bcl_descriptor_base.hpp"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_mutation.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    bool Base< t_DataType, t_ReturnType>::s_HaveUpdatedHelp
    (
      util::Implementation< Base< t_DataType, t_ReturnType> >::SetHelpWriter
      (
        &Base< t_DataType, t_ReturnType>::WriteInstancesHelp
      )
      &&
      util::Implementation< Base< t_DataType, t_ReturnType> >::SetUndefinedInstanceName
      (
        Base< t_DataType, t_ReturnType>::GetObjectName() + "/" + Base< t_DataType, t_ReturnType>::GetElementName()
        + ( type::Compare< t_ReturnType, float>::e_Same ? " Numeric " : " String ")
        + "Descriptor"
      )
      &&
      util::Enumerated< util::ImplementationInterface>::AddInstance
      (
        new util::Implementation< Base< t_DataType, t_ReturnType> >(),
        GetStaticClassName< Base< t_DataType, t_ReturnType> >()
      )
      &&
      util::Enumerated< util::ImplementationInterface>::AddInstance
      (
        new util::Implementation< Base< t_DataType, t_ReturnType> >(),
        util::Implementation< Base< t_DataType, t_ReturnType> >::GetUndefinedInstanceName()
      )
    );

    template class BCL_API Base< chemistry::AtomConformationalInterface, char>;
    template class BCL_API Base< chemistry::AtomConformationalInterface, float>;
    template class BCL_API Base< biol::AABase, char>;
    template class BCL_API Base< biol::AABase, float>;
    template class BCL_API Base< biol::Mutation, char>;
    template class BCL_API Base< biol::Mutation, float>;
    template class BCL_API Base< char, char>;
    template class BCL_API Base< char, float>;

  } // namespace descriptor
} // namespace bcl
