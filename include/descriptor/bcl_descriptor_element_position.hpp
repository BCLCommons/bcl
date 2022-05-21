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

// include header of this class
#include "bcl_descriptor_element_position.h"
// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "io/bcl_io_serializer.h"
#include "linal/bcl_linal_vector_3d.h"
#include "linal/bcl_linal_vector_reference.h"
#include "type/bcl_type_compare.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief Clone function
    //! @return pointer to new ElementPosition
    template< typename t_DataType>
    ElementPosition< t_DataType> *ElementPosition< t_DataType>::Clone() const
    {
      return new ElementPosition( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &ElementPosition< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &ElementPosition< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "Position");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    //! @note only one value is needed if this is a numeric descriptor, for char descriptors, assume that 99999 is the
    //! @note max, so 5 characters is sufficient
    template< typename t_DataType>
    size_t ElementPosition< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return 3;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer ElementPosition< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Returns the X,Y,Z coordinates of the " + this->GetElementName());

      return parameters;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void ElementPosition< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      const linal::Vector3D &position( ELEMENT->GetCenter());
      STORAGE( 0) = position.X();
      STORAGE( 1) = position.Y();
      STORAGE( 2) = position.Z();
    }

  } // namespace descriptor
} // namespace bcl
