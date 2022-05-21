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
#include "io/bcl_io_serialization_interface.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    //! @brief set the value of the corresponding member based on the label
    //! @param LABEL label that is used to set the string
    //! @param ERROR_STREAM stream to write errors to
    //! @return result of validation
    ValidationResult SerializationInterface::ValidateRead
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      // check if the name or value is == the help string
      if( LABEL.IsScalar())
      {
        if
        (
          LABEL.GetName() == ValidationResult::GetHelpString()
          || ( LABEL.GetName().empty() && LABEL.GetValue() == ValidationResult::GetHelpString())
        )
        {
          WriteHelp( ERROR_STREAM);
          return ValidationResult( e_Help);
        }
      }
      // otherwise, just return the result of try read
      return ValidationResult( TryRead( LABEL, ERROR_STREAM));
    }

    //! @brief Get a set of all class names used by the serializer. Useful for introspection
    //! @param TYPES set to insert type names into
    //! @param INCLUDE_OPTIONAL true to also count optional members
    //! @param INCLUDE_DATA true to also include data-containing members
    void SerializationInterface::InsertDataTypes
    (
      storage::Map< std::string, size_t> &TYPES,
      const bool &INCLUDE_OPTIONAL,
      const bool &INCLUDE_DATA,
      const size_t &MAX_DEPTH
    ) const
    {
      ++TYPES[ this->GetClassIdentifier()];
    }

    //! @brief Get addresses of all objects serialized as part of this
    //! @notes used to ensure uniqueness of serialized objects
    std::vector< size_t> SerializationInterface::GetSerializationAddresses() const
    {
      return std::vector< size_t>( size_t( 1), size_t( this));
    }

  } // namespace io
} // namespace bcl
