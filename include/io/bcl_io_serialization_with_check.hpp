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

#ifndef BCL_IO_SERIALIZATION_WITH_CHECK_HPP_
#define BCL_IO_SERIALIZATION_WITH_CHECK_HPP_

// include the header of this class
#include "bcl_io_serialization_with_check.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from base class
    //! @param PARAMETER_CHECK any util::ParameterCheckInterface derived class
    //! @param DATA reference to the associated member variable
    template< typename t_DataType>
    SerializationWithCheck< t_DataType>::SerializationWithCheck
    (
      const command::ParameterCheckInterface &PARAMETER_CHECK,
      const t_DataType *DATA
    ) :
      SerializationBuiltin< t_DataType>( DATA),
      m_Check( PARAMETER_CHECK.Clone())
    {
    }

    //! @brief Clone function
    //! @return pointer to new SerializationWithCheck
    template< typename t_DataType>
    SerializationWithCheck< t_DataType> *SerializationWithCheck< t_DataType>::Clone() const
    {
      return new SerializationWithCheck( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &SerializationWithCheck< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief set the given object using the label
    //! @param OBJECT object to read
    //! @param LABEL label that is used to set the string
    //! @param ERR_STREAM stream to write errors to
    //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
    template< typename t_DataType>
    bool SerializationWithCheck< t_DataType>::TryReadObject
    (
      t_DataType &OBJECT,
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    ) const
    {
      if( !m_Check->IsAllowedParameter( LABEL.ToString(), GetStaticClassName< t_DataType>(), ERR_STREAM))
      {
        return false;
      }
      return SerializationBuiltin< t_DataType>::TryReadObject( OBJECT, LABEL, ERR_STREAM);
    }

    //! @brief writes the help for the label
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT the amount of indent
    //! @return the given stream to which the help was written to
    template< typename t_DataType>
    std::ostream &SerializationWithCheck< t_DataType>::WriteHelp
    (
      std::ostream &OSTREAM,
      const size_t INDENT
    ) const
    {
      return m_Check->WriteHelp( OSTREAM, INDENT);
    }

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_WITH_CHECK_HPP_

