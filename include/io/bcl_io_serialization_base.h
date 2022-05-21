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

#ifndef BCL_IO_SERIALIZATION_BASE_H_
#define BCL_IO_SERIALIZATION_BASE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_interface.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_assert.h"
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_object_data_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializationBase
    //! @brief adds a pointer to a templated data type onto the label handler interface, with corresponding get/set functions
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class SerializationBase :
      public SerializationInterface
    {

    private:

    //////////
    // data //
    //////////

      t_DataType *m_Data; //<! Reference to the associated data member

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from base class
      //! @param DATA reference to the associated member variable
      explicit SerializationBase( const t_DataType *DATA = NULL) :
        m_Data( DATA ? const_cast< t_DataType *>( DATA) : NULL)
      {
      }

      //! @brief Clone function
      //! @return pointer to new SerializationBase
      virtual SerializationBase *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the label containing only initialization parameters
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      util::ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const
      {
        BCL_Assert( m_Data, "Called GetLabel on " + GetStaticClassName< t_DataType>() + " without object!");
        return GetLabelForObject( *m_Data, WITH_DATA);
      }

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        BCL_Assert( m_Data, "Called TryRead on " + GetStaticClassName< t_DataType>() + " without object");
        return TryReadObject( *m_Data, LABEL, ERROR_STREAM);
      }

      //! @brief set the given object using the label
      //! @param OBJECT object to read
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      virtual bool TryReadObject
      (
        t_DataType &OBJECT,
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      ) const = 0;

      //! @brief Get the label for an object of the given type
      //! @param OBJECT the object to get the label for
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      virtual util::ObjectDataLabel GetLabelForObject( const t_DataType &OBJECT, const bool &WITH_DATA) const = 0;

      //! @brief Get a set of all class names used by the serializer. Useful for introspection
      //! @param TYPES set to insert type names into
      //! @param INCLUDE_OPTIONAL true to also count optional members
      //! @param INCLUDE_DATA true to also include data-containing members
      virtual void InsertDataTypes
      (
        storage::Map< std::string, size_t> &TYPES,
        const bool &INCLUDE_OPTIONAL = true,
        const bool &INCLUDE_DATA = false,
        const size_t &MAX_DEPTH = size_t( -1)
      ) const
      {
        BCL_Assert( m_Data, "Called InsertDataTypes on " + GetStaticClassName< t_DataType>() + " without object");
        InsertDataTypesForObject( *m_Data, TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA, MAX_DEPTH);
      }

      //! @brief Get a set of all class names used by the serializer. Useful for introspection
      //! @param TYPES set to insert type names into
      //! @param INCLUDE_OPTIONAL true to also count optional members
      //! @param INCLUDE_DATA true to also include data-containing members
      virtual void InsertDataTypesForObject
      (
        const t_DataType &OBJECT,
        storage::Map< std::string, size_t> &TYPES,
        const bool &INCLUDE_OPTIONAL,
        const bool &INCLUDE_DATA,
        const size_t &MAX_DEPTH
      ) const
      {
        ++TYPES[ GetStaticClassName< t_DataType>()];
      }

    protected:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief Get addresses of all objects serialized as part of this
      //! @notes used to ensure uniqueness of serialized objects
      std::vector< size_t> GetSerializationAddresses() const
      {
        return std::vector< size_t>( size_t( 1), size_t( m_Data));
      }

    }; // class SerializationBase

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_BASE_H_

