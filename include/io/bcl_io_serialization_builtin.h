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

#ifndef BCL_IO_SERIALIZATION_BUILTIN_H_
#define BCL_IO_SERIALIZATION_BUILTIN_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializationBuiltin
    //! @brief serialization for builtin types (double, float, int, etc.) or character-based data types ( std::string, util::ObjectDataLabel, etc.)
    //!
    //! @see @link example_io_serialization_builtin.cpp @endlink
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class SerializationBuiltin :
      public SerializationBase< t_DataType>
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      explicit SerializationBuiltin( const t_DataType *DATA = NULL);

      //! @brief Clone function
      //! @return pointer to new SerializationBuiltin
      SerializationBuiltin *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief Get the label for an object of the given type
      //! @param OBJECT the object to get the label for
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      util::ObjectDataLabel GetLabelForObject( const t_DataType &OBJECT, const bool &WITH_DATA) const;

      //! @brief determine the type of value that the handler expects to parse
      util::DataType::Type GetSerializedType() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief set the given object using the label
      //! @param OBJECT object to read
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryReadObject( t_DataType &OBJECT, const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM) const;

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

    }; // class SerializationBuiltin

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< short>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< unsigned short>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< bool>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< signed char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< unsigned char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< std::string>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationBuiltin< util::ObjectDataLabel>;

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_BUILTIN_H_

