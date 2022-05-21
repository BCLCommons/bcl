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

#ifndef BCL_IO_SERIALIZATION_MAP_H_
#define BCL_IO_SERIALIZATION_MAP_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_base.h"
#include "bcl_io_serialize.h"
#include "iterate/bcl_iterate.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializationMap
    //! @brief data labels for any associative container (usually either a map or a hash map)
    //! Containers must have the following properties:
    //!   - A key_type and value_type typedef
    //!   - Functions called Begin, End that return const_iterators
    //!   - a constructor from input iterators
    //!
    //! @see @link example_io_serialization_map.cpp @endlink
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_ContainerType>
    class SerializationMap :
      public SerializationBase< t_ContainerType>
    {

    private:

    //////////
    // data //
    //////////

      //! import typedefs from the container type for the key and value for the type
      typedef typename t_ContainerType::key_type     key_type;
      typedef typename t_ContainerType::mapped_type  mapped_type;

      //! create a typedef for the type held internally by the object
      //! Because we only require a const_iterator typedef, it is most general
      //! to use the const_iterator type and just remove the const from it
      typedef std::pair< key_type, mapped_type>      t_DataType;

      util::OwnPtr< SerializationBase< key_type> >    m_KeyParameter;    //!< Controls how keys are set
      util::OwnPtr< SerializationBase< mapped_type> > m_MappedParameter; //!< Controls how mapped types are handled
      size_t                                          m_MinSize;         //!< Minimum # of objects in the container
      size_t                                          m_MaxSize;         //!< Maximum # of objects in the container

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from base class
      //! @param MEMBER_SETTER how to set up the objects that will be contained in this object
      //! @param CONTAINER reference to the associated container variable
      //! @param PARAMETER parameter to copy
      //! @param MIN_SIZE minimum number of objects the container should contain
      SerializationMap
      (
        const SerializationBase< key_type> &KEY_SETTER,
        const SerializationBase< mapped_type> &MAPPED_TYPE_SETTER,
        const t_ContainerType *CONTAINER = NULL,
        const size_t &MIN_SIZE = size_t( 1),
        const size_t &MAX_SIZE = std::numeric_limits< size_t>::max()
      ) :
        SerializationBase< t_ContainerType>( CONTAINER),
        m_KeyParameter( KEY_SETTER.Clone()),
        m_MappedParameter( MAPPED_TYPE_SETTER.Clone()),
        m_MinSize( MIN_SIZE),
        m_MaxSize( MAX_SIZE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new SerializationMap
      virtual SerializationMap *Clone() const
      {
        return new SerializationMap( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief Get the label for an object of the given type
      //! @param OBJ the object to get the label for
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      util::ObjectDataLabel GetLabelForObject( const t_ContainerType &OBJ, const bool &WITH_DATA) const
      {
        // make a vector to hold the labels of all objects that this container contains
        std::vector< util::ObjectDataLabel> labels;

        for( typename t_ContainerType::const_iterator itr( OBJ.Begin()), itr_end( OBJ.End()); itr != itr_end; ++itr)
        {
          labels.push_back
          (
            util::ObjectDataLabel
            (
              m_KeyParameter->GetLabelForObject( itr->first, WITH_DATA).ToString(),
              m_MappedParameter->GetLabelForObject( itr->second, WITH_DATA)
            )
          );
        }

        return util::ObjectDataLabel( labels);
      }

      //! @brief determine the type of value that the handler expects to parse
      util::DataType::Type GetSerializedType() const
      {
        return util::DataType::e_Map;
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
      bool TryReadObject( t_ContainerType &OBJECT, const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM) const
      {
        // check the number of arguments
        if( LABEL.GetNumberArguments() < m_MinSize || LABEL.GetNumberArguments() > m_MaxSize)
        {
          ERR_STREAM << "Expected a ";
          WriteDescription( ERR_STREAM);
          ERR_STREAM << " key-value pairs, but was given: " << LABEL.ToString();
          return false;
        }

        // common case, arguments were given in the label, these will form the container

        // create a vector to hold all values contained in the arguments of this label
        std::vector< t_DataType> holder( LABEL.GetNumberArguments());
        typename std::vector< t_DataType>::iterator itr_holder( holder.begin());

        // add every element to the vector
        for
        (
          util::ObjectDataLabel::const_iterator itr( LABEL.Begin()), itr_end( LABEL.End());
          itr != itr_end;
          ++itr, ++itr_holder
        )
        {
          if
          (
            !m_KeyParameter->TryReadObject( itr_holder->first, util::ObjectDataLabel( "", itr->GetName().empty() && itr->GetNumberArguments() ? itr->GetValue() : itr->GetName()), ERR_STREAM)
            || !m_MappedParameter->TryReadObject( itr_holder->second, *itr, ERR_STREAM)
          )
          {
            return false;
          }
        }

        OBJECT = t_ContainerType( holder.begin(), holder.end());

        return true;
      }

      //! @brief hook for label handlers to describe themselves
      //! @param OSTREAM stream to write the output to
      void WriteDescription( std::ostream &OSTREAM) const
      {
        OSTREAM << "Map";
        util::DataType::WriteSizeRequirements( OSTREAM, m_MinSize, m_MaxSize);
      }

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const
      {
        WriteDescription( OSTREAM);
        OSTREAM << '\n';
        // write any help for the contained object
        Serialize::InsertIndent( OSTREAM, INDENT + 1) << "Key info (LHS of = sign): ";
        m_KeyParameter->WriteHelp( OSTREAM, INDENT + 1);
        Serialize::InsertIndent( OSTREAM, INDENT + 1) << "Mapped value info (RHS of = sign): ";
        m_MappedParameter->WriteHelp( OSTREAM, INDENT + 1);

        // end
        return OSTREAM;
      }

      //! @brief Get a set of all class names used by the serializer. Useful for introspection
      //! @param TYPES set to insert type names into
      //! @param INCLUDE_OPTIONAL true to also count optional members
      //! @param INCLUDE_DATA true to also include data-containing members
      virtual void InsertDataTypesForObject
      (
        const t_ContainerType &OBJECT,
        storage::Map< std::string, size_t> &TYPES,
        const bool &INCLUDE_OPTIONAL,
        const bool &INCLUDE_DATA,
        const size_t &MAX_DEPTH
      ) const
      {
        if( m_MinSize != m_MaxSize)
        {
          ++TYPES[ GetStaticClassName< t_ContainerType>()];
        }
        if( MAX_DEPTH)
        {
          if( OBJECT.Begin() != OBJECT.End())
          {
            for( typename t_ContainerType::const_iterator itr( OBJECT.Begin()), itr_end( OBJECT.End()); itr != itr_end; ++itr)
            {
              m_KeyParameter->InsertDataTypesForObject( itr->first, TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA, MAX_DEPTH - 1);
              m_MappedParameter->InsertDataTypesForObject( itr->second, TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA, MAX_DEPTH - 1);
            }
          }
          else
          {
            for( size_t i( 0), nominal_sz( std::max( std::min( m_MaxSize, size_t( 3)), m_MinSize)); i < nominal_sz; ++i)
            {
              m_KeyParameter->InsertDataTypesForObject( key_type(), TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA, MAX_DEPTH - 1);
              m_MappedParameter->InsertDataTypesForObject( mapped_type(), TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA, MAX_DEPTH - 1);
            }
          }
        }
      }

    }; // class SerializationMap

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_MAP_H_

