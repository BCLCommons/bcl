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

#ifndef BCL_IO_SERIALIZATION_WITH_MIN_MAX_H_
#define BCL_IO_SERIALIZATION_WITH_MIN_MAX_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_builtin.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializationWithMinMax
    //! @brief data labels for plain old data types (double, float, etc.) for which the number is desired within a
    //!        a specific range
    //!
    //! @see @link example_io_serialization_with_min_max.cpp @endlink
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class SerializationWithMinMax :
      public SerializationBuiltin< t_DataType>
    {
    private:

      t_DataType m_MinValue;
      t_DataType m_MaxValue;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from members
      //! @param MIN the minimum allowed character
      //! @param MAX the maximum allowed character
      SerializationWithMinMax( const t_DataType &MIN, const t_DataType &MAX);

      //! @brief constructor from members
      //! @param DATA reference to the associated member variable
      //! @param MIN the minimum allowed character
      //! @param MAX the maximum allowed character
      SerializationWithMinMax
      (
        const t_DataType *DATA,
        const t_DataType &MIN,
        const t_DataType &MAX
      );

      //! @brief Clone function
      //! @return pointer to new SerializationWithMinMax
      SerializationWithMinMax *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

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

    }; // class SerializationWithMinMax

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< short>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< unsigned short>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< signed char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithMinMax< unsigned char>;

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_WITH_MIN_MAX_H_

