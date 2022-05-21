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

#ifndef BCL_IO_SERIALIZATION_INTERFACE_H_
#define BCL_IO_SERIALIZATION_INTERFACE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_io_validation_result.h"
#include "util/bcl_util_data_type.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <set>

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializationInterface
    //! @brief interface for ways of reading, writing, and getting help for something
    //!
    //! @note Derivation from object interface to be removed once all bcl objects only need
    //!       the core functionality given herein
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SerializationInterface :
      public virtual util::ObjectInterface
    {
    public:

      //! This is a typedef used to allow compile-time error detection for object serialization
      //! Generally, derived classes can override this class if they provide their own serialization methods
      typedef SerializationInterface t_SerializedType;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SerializationInterface
      virtual SerializationInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the label containing only initialization parameters
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      virtual util::ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const = 0;

      //! @brief get the basic data type that this handlers expects in TryRead and that is returned from GetLabel
      //! @return the basic data type that this handlers expects in TryRead and that is returned from GetLabel
      virtual util::DataType::Type GetSerializedType() const
      {
        // by default, derived objects can be serialized to / from a simple string
        return util::DataType::e_String;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      virtual bool TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM) = 0;

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return result of validation
      virtual ValidationResult ValidateRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      virtual std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const = 0;

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
      ) const;

      //! @brief Get addresses of all objects serialized as part of this
      //! @notes used to ensure uniqueness of serialized objects
      virtual std::vector< size_t> GetSerializationAddresses() const;

    }; // class SerializationInterface

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_INTERFACE_H_

