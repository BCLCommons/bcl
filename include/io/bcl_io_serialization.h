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

#ifndef BCL_IO_SERIALIZATION_H_
#define BCL_IO_SERIALIZATION_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_base.h"
#include "bcl_io_serialization_builtin.h"
#include "bcl_io_serialization_container.h"
#include "bcl_io_serialization_interface.h"
#include "bcl_io_serialization_map.h"
#include "bcl_io_serialization_with_check.h"
#include "bcl_io_serialization_with_min_max.h"
#include "bcl_io_serialization_wrapper.h"
#include "bcl_io_serializer.h"
#include "math/bcl_math_limits.h"
#include "type/bcl_type_is_map.h"
#include "type/bcl_type_is_sequence.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Serialization
    //! @brief Helper methods to retrieve objects that handle serialization
    //!
    //! @see @link example_io_serialization.cpp @endlink
    //! @author mendenjl
    //! @date Oct 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Serialization
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      Serialization();

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class HasSerializedType
      //! @brief a simple helper class that determines whether a type has a typedef t_SerializedType
      //! @remarks Example unnecessary
      //!
      //! @author mendenjl
      //! @date Oct 29, 2012
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      template< typename t_Type>
      class HasSerializedType
      {
        //!< Helper typedef used for template meta programming
        typedef char   t_True;
        typedef double t_False;

        //! @brief a function declaration; used only to determine whether t_PseudoType has a typedef t_SerializedType
        //! @note this is an example of Substitution Failure Is Not An Error (SFINAE)
        template< typename t_PseudoType>
        static t_True Test( typename t_PseudoType::t_SerializedType *);

        //! @brief a function declaration; used only to determine whether t_PseudoType has a typedef t_SerializedType
        template< typename t_PseudoType>
        static t_False Test(...);

      public:

        enum { value = ( sizeof( Test< t_Type>( 0)) == sizeof( t_True))};
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Type
      //! @brief a simple helper class that determines the basic type of the templated object
      //! @remarks Example unnecessary
      //!
      //! @author mendenjl
      //! @date Oct 29, 2012
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      template< typename t_Type>
      struct Type
      {
        enum
        {
          e_Sequence =                                      //!< e_Sequence == true for basic sequences
            type::IsSequence< t_Type>::value
            && !type::IsMap< t_Type>::value
            && !HasSerializedType< t_Type>::value,          //!< Serialization-Interface derived object
          e_Serialized = HasSerializedType< t_Type>::value,
          e_Map =                                           //!< Map type objects
            type::IsMap< t_Type>::value && !HasSerializedType< t_Type>::value
        };
      };

    public:

    //////////////
    // GetAgent //
    //////////////

      //! @brief handler for string, no checking performed
      //! @param VALUE the data member to obtain a handler for
      //! @return a serialization handler for strings, no error checking or help
      static util::OwnPtr< SerializationBase< std::string> >
        GetAgent( const std::string *VALUE)
      {
        return
          util::OwnPtr< SerializationBase< std::string> >
          (
            new SerializationBuiltin< std::string>( VALUE)
          );
      }

      //! @brief handler for util::ObjectDataLabel, no checking performed
      //! @param VALUE the data member to obtain a handler for
      //! @return a serialization handler for util::ObjectDataLabel, no error checking or help
      static util::OwnPtr< SerializationBase< util::ObjectDataLabel> >
        GetAgent( const util::ObjectDataLabel *VALUE)
      {
        return
          util::OwnPtr< SerializationBase< util::ObjectDataLabel> >
          (
            new SerializationBuiltin< util::ObjectDataLabel>( VALUE)
          );
      }

      //! @brief handler for serializable agent interface; parameter checks are unnecessary due to TryRead in the interface
      //! @param OBJECT the data member to obtain a handler for
      //! @return a serialization handler for the given object
      static util::OwnPtr< SerializationInterface>
        GetAgent( const SerializationInterface *OBJECT)
      {
        return
         util::OwnPtr< SerializationInterface>
         (
           const_cast< SerializationInterface *>( OBJECT),
           false
         );
      }

      //! @brief handler for numeric types check that the values are valid for numeric values in the given range automatically
      //! @param VALUE a pointer to the data member to obtain a handler for
      //! @return a serialization handler for the given object
      template< typename t_DataType>
      static
      typename type::EnableIf
      <
        std::numeric_limits< t_DataType>::is_specialized,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type GetAgent( const t_DataType *VALUE)
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationBuiltin< t_DataType>( VALUE)
          );
      }

      //! @brief this overload matches simple sequence containers
      template< typename t_DataType>
      static
      typename type::EnableIf
      <
        Type< t_DataType>::e_Sequence,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type GetAgent( const t_DataType *CONTAINER)
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationContainer< t_DataType>
            (
              *GetAgentForContainedObject< typename SerializationContainer< t_DataType>::t_TypeContained>(),
              CONTAINER
            )
          );
      }

      //! @brief this overload matches maps
      template< typename t_DataType>
      static
      typename type::EnableIf
      <
        Type< t_DataType>::e_Map,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type GetAgent( const t_DataType *CONTAINER)
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationMap< t_DataType>
            (
              *GetAgentForContainedObject< typename t_DataType::key_type>(),
              *GetAgentForContainedObject< typename t_DataType::mapped_type>(),
              CONTAINER
            )
          );
      }

    /////////////////////////////////////////
    // handlers with additional validation //
    /////////////////////////////////////////

      //! @brief handler for a type that will use a particular check
      //! @param VALUE the data member to obtain a handler for
      //! @param CHECK the parameter check to use
      //! @return a serialization handler for the given object using the provided check
      template< typename t_DataType>
      static util::OwnPtr< SerializationBase< t_DataType> >
        GetAgentWithCheck( const t_DataType *VALUE, const command::ParameterCheckInterface &CHECK)
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationWithCheck< t_DataType>( CHECK, VALUE)
          );
      }

      //! @brief this overload matches simple sequence containers
      template< typename t_DataType>
      static
      typename type::EnableIf
      <
        Type< t_DataType>::e_Sequence,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type GetAgentWithSizeLimits
      (
        const t_DataType *CONTAINER,
        const size_t &MIN_SIZE,
        const size_t &MAX_SIZE = std::numeric_limits< size_t>::max()
      )
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationContainer< t_DataType>
            (
              *GetAgentForContainedObject< typename SerializationContainer< t_DataType>::t_TypeContained>(),
              CONTAINER,
              MIN_SIZE,
              MAX_SIZE
            )
          );
      }

      //! @brief this overload matches simple sequence containers
      template< typename t_DataType>
      static
      typename type::EnableIf
      <
        Type< t_DataType>::e_Map,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type GetAgentWithSizeLimits
      (
        const t_DataType *CONTAINER,
        const size_t &MIN_SIZE,
        const size_t &MAX_SIZE = std::numeric_limits< size_t>::max()
      )
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationMap< t_DataType>
            (
              *GetAgentForContainedObject< typename t_DataType::key_type>(),
              *GetAgentForContainedObject< typename t_DataType::value_type>(),
              CONTAINER,
              MIN_SIZE,
              MAX_SIZE
            )
          );
      }

      //! @brief handler for input files; checks that the file exists before setting
      //! @param VALUE the data member to obtain a handler for
      //! @return a serialization handler for the given object
      static util::OwnPtr< SerializationBase< std::string> >
      GetAgentInputFilename( const std::string *VALUE);

      //! @brief handler for numeric types check that the values are valid and within a certain range
      //! @param VALUE the data member to obtain a handler for
      //! @param MIN the minimum value to allow
      //! @param MAX the maximum value to allow
      //! @return a serialization handler for the given object
      template< typename t_DataType, typename t_OtherDataType>
      static
      typename type::EnableIf
      <
        std::numeric_limits< t_DataType>::is_specialized,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type
        GetAgentWithRange( const t_DataType *VALUE, const t_OtherDataType &MIN, const t_OtherDataType &MAX)
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationWithMinMax< t_DataType>( VALUE, t_DataType( MIN), t_DataType( MAX))
          );
      }

      //! @brief handler for numeric types check that the values are valid and within a certain range
      //! @param MIN the minimum value to allow
      //! @param MAX the maximum value to allow
      //! @return a serialization handler for a numeric value with the given validation
      template< typename t_DataType>
      static
      typename type::EnableIf
      <
        std::numeric_limits< t_DataType>::is_specialized,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type
        GetAgentWithRange( const t_DataType &MIN, const t_DataType &MAX)
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationWithMinMax< t_DataType>( MIN, MAX)
          );
      }

      //! @brief Convenience function for specifying a handler with only a min, unlimited max (inf for floating point types)
      //! @param VALUE a pointer to the data member to obtain a handler for
      //! @param MIN the minimum value to allow
      //! @return a serialization handler for the given object
      //! @note specialization for types with no infinity value
      template< typename t_DataType, typename t_OtherDataType>
      static util::OwnPtr< SerializationBase< t_DataType> >
      GetAgentWithMin( const t_DataType *VALUE, const t_OtherDataType &MIN)
      {
        // if the type has an infinity value, use that for the max, otherwise, just call max directly
        return GetAgentWithRange( VALUE, t_DataType( MIN), math::GetHighestUnboundedValue< t_DataType>());
      }

      //! @brief Convenience function for specifying a handler with only a min, unlimited max (inf for floating point types)
      //! @param MIN the minimum value to allow
      //! @return a serialization handler for the given object
      //! @note specialization for types with no infinity value
      template< typename t_DataType>
      static util::OwnPtr< SerializationBase< t_DataType> > GetAgentWithMin( const t_DataType &MIN)
      {
        // if the type has an infinity value, use that for the max, otherwise, just call max directly
        return GetAgentWithRange( MIN, math::GetHighestUnboundedValue< t_DataType>());
      }

      //! @brief Convenience function for specifying a handler with only a max, unlimited min (-inf for floating point types)
      //! @param VALUE a pointer to the data member to obtain a handler for
      //! @param MAX the maximum value to allow
      //! @return a serialization handler for the given object
      template< typename t_DataType, typename t_OtherDataType>
      static util::OwnPtr< SerializationBase< t_DataType> >
      GetAgentWithMax( const t_DataType *VALUE, const t_OtherDataType &MAX)
      {
        // if the type has an infinity value, use that for the max, otherwise, just call max directly
        return GetAgentWithRange( VALUE, math::GetLowestUnboundedValue< t_DataType>(), t_DataType( MAX));
      }

      //! @brief Convenience function for specifying a handler with only a max, unlimited min (-inf for floating point types)
      //! @param MAX the maximum value to allow
      //! @return a serialization handler for the given object
      template< typename t_DataType>
      static util::OwnPtr< SerializationBase< t_DataType> > GetAgentWithMax( const t_DataType &MAX)
      {
        // if the type has an infinity value, use that for the max, otherwise, just call max directly
        return GetAgentWithRange( math::GetLowestUnboundedValue< t_DataType>(), t_DataType( MAX));
      }

      template< typename t_DataType, typename t_OtherType>
      static
      typename type::EnableIf
      <
        Type< t_DataType>::e_Sequence,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type
      GetAgentContainerWithCheck
      (
        const t_DataType *CONTAINER,
        const util::OwnPtr< t_OtherType> &CHECK,
        const size_t &MIN_SIZE = size_t( 1),
        const size_t &MAX_SIZE = std::numeric_limits< size_t>::max()
      )
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationContainer< t_DataType>( *CHECK, CONTAINER, MIN_SIZE, MAX_SIZE)
          );
      }

    private:

      //! @brief method to get a handler for an object inside a container
      //! @return the handler for the contained object
      //! @note overload for serializable types
      template< typename t_DataType>
      static
      typename type::EnableIf
      <
        HasSerializedType< t_DataType>::value,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type
      GetAgentForContainedObject()
      {
        return
          util::OwnPtr< SerializationBase< t_DataType> >
          (
            new SerializationWrapper< t_DataType>()
          );
      }

      //! @brief method to get a handler for an object inside a container
      //! @return the handler for the contained object
      //! @note overload for types not inherited from Serialized object interface=
      template< typename t_DataType>
      static
      typename type::EnableIf
      <
        !HasSerializedType< t_DataType>::value,
        util::OwnPtr< SerializationBase< t_DataType> >
      >::Type
      GetAgentForContainedObject()
      {
        return GetAgent( ( const t_DataType *)( NULL));
      }

    }; // class Serialization

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_H_
