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

#ifndef BCL_UTIL_IMPLEMENTATION_INTERFACE_H_
#define BCL_UTIL_IMPLEMENTATION_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_serializable_interface.h"
#include "io/bcl_io_fixed_line_width_writer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ImplementationInterface
    //! @brief Interface derived from by all implementation objects
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Apr 09, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ImplementationInterface :
      virtual public SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone the implementation
      //! @return pointer to new ImplementationInterface
      virtual ImplementationInterface *Clone() const = 0;

      //! @brief return an empty copy of the implementation
      //! @return pointer to an empty copy of the implementation
      virtual ImplementationInterface *Empty() const = 0;

      //! @brief assign to a different implementation via a string
      virtual ImplementationInterface &operator =( const std::string &DATA_LABEL)
      {
        AssertRead( ObjectDataLabel( DATA_LABEL));
        return *this;
      }

      //! @brief assign to a different implementation via a data label
      virtual ImplementationInterface &operator =( const ObjectDataLabel &DATA_LABEL)
      {
        AssertRead( DATA_LABEL);
        return *this;
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief determine the type of value that can be parsed
      DataType::Type GetSerializedType() const
      {
        return DataType::e_DynamicObject;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief remove any existing implementation
      virtual void Reset() = 0;

      //! @brief test whether the implementation is defined
      //! @return true if implementation is defined
      virtual bool IsDefined() const = 0;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

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
      ) const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // template class ImplementationInterface

    //! comparison of two implementations '<'
    //! used for ordering implementations in sets/maps
    inline bool operator <( const ImplementationInterface &NAME_A, const ImplementationInterface &NAME_B)
    {
      return NAME_A.GetLabel() < NAME_B.GetLabel();
    }

    //! comparison of two property interfaces '>'
    //! meaning the property interfaces on the right side are > property interfaces on the left side
    inline bool operator >( const ImplementationInterface &NAME_A, const ImplementationInterface &NAME_B)
    {
      return NAME_A.GetLabel() > NAME_B.GetLabel();
    }

    //! comparison of two property interfaces
    //! meaning the property interfaces on the right side are a subgroup of the property interfaces on the left side
    inline bool operator !=( const ImplementationInterface &NAME_A, const ImplementationInterface &NAME_B)
    {
      return NAME_B.GetLabel() != NAME_A.GetLabel();
    }

    //! equality of two property interfaces '=='
    inline bool operator ==( const ImplementationInterface &NAME_A, const ImplementationInterface &NAME_B)
    {
      return NAME_B.GetLabel() == NAME_A.GetLabel();
    }

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_IMPLEMENTATION_INTERFACE_H_

