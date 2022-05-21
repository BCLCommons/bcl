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

#ifndef BCL_UTIL_WRAPPER_ENUM_H_
#define BCL_UTIL_WRAPPER_ENUM_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_logger_interface.h"
#include "io/bcl_io_serialization_interface.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WrapperEnum
    //! @brief WrapperEnum wraps a normal enum class into a bcl object
    //!
    //! @see @link example_util_wrapper_enum.cpp @endlink
    //! @author mendenjl
    //! @date Oct 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template
    <
      typename t_Enum, // the enumerated type
      const std::string &( *GetEnumDescriptor)( const t_Enum &), // the function that returns a string given a t_Enum
      size_t s_NumberOfEnums
    >
    class WrapperEnum :
      public io::SerializationInterface
    {
    private:
    //////////
    // data //
    //////////

      //! actual enum value
      t_Enum m_Value;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor (initializes to first value in the enum)
      WrapperEnum() :
        m_Value( t_Enum( 0))
      {
      }

      //! @brief constructor from string
      //! @param VALUE the string to try to construct this enum from
      WrapperEnum( const std::string &VALUE) :
        m_Value( t_Enum( 0))
      {
        BCL_Assert( TryRead( VALUE, GetLogger()), "Bad enum value");
      }

      //! @brief constructor from enum
      //! @param ENUM the enum to use
      WrapperEnum( const t_Enum &VALUE) :
        m_Value( VALUE)
      {
      }

      //! @brief copy constructor
      WrapperEnum( const WrapperEnum &WRAPPER) :
        m_Value( WRAPPER.m_Value)
      {
      }

      //! @brief Clone function
      //! @return pointer to new SerializationBase
      WrapperEnum *Clone() const
      {
        return new WrapperEnum( *this);
      }

    /////////////////
    // Data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName< t_Enum>();
      }

      //! @brief get the value of the enum
      //! @return the raw enum value
      const t_Enum &GetEnum() const
      {
        return m_Value;
      }

      //! @brief get the value of the enum
      //! @return the raw enum value
      const std::string &GetString() const
      {
        return ( *GetEnumDescriptor)( m_Value);
      }

      //! @brief get the label containing only initialization parameters
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const
      {
        return ObjectDataLabel( "", GetString());
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief conversion to basic enum type
      operator const t_Enum &() const
      {
        return m_Value;
      }

      //! @brief conversion to string
      operator const std::string &() const
      {
        return GetString();
      }

      //! @brief assignment operator
      //! @param ENUM the enum to set this wrapper up using
      //! @return reference to this
      WrapperEnum &operator =( const WrapperEnum &ENUM)
      {
        m_Value = ENUM.m_Value;
        return *this;
      }

      //! @brief assignment operator
      //! @param ENUM the enum to set this wrapper up using
      //! @return reference to this
      WrapperEnum &operator =( const t_Enum &ENUM)
      {
        m_Value = ENUM;
        return *this;
      }

      //! @brief assignment operator
      //! @param STRING the string to set this wrapper up using
      //! @return reference to this
      WrapperEnum &operator =( const std::string &STRING)
      {
        BCL_Assert( TryRead( STRING, GetLogger()), "Bad enum value");
        return *this;
      }

      //! @brief equality comparison
      //! @param WRAPPER the other enum
      //! @return true if the enums are equal
      bool operator ==( const WrapperEnum &WRAPPER)
      {
        return m_Value == WRAPPER.m_Value;
      }

      //! @brief inequality comparison
      //! @param WRAPPER the other enum
      //! @return true if the enums are not equal
      bool operator !=( const WrapperEnum &WRAPPER)
      {
        return m_Value != WRAPPER.m_Value;
      }

      //! @brief equality comparison
      //! @param ENUM the other enum
      //! @return true if the enums are equal
      bool operator ==( const t_Enum &ENUM)
      {
        return m_Value == ENUM;
      }

      //! @brief inequality comparison
      //! @param ENUM the other enum
      //! @return true if the enums are not equal
      bool operator !=( const t_Enum &ENUM)
      {
        return m_Value != ENUM;
      }

      //! @brief inequality comparison
      //! @param ENUM the other enum
      //! @return true if the enums are not equal
      WrapperEnum &operator ++()
      {
        m_Value = t_Enum( long( m_Value) + 1);
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        return this->TryRead( LABEL.GetValue(), ERROR_STREAM);
      }

      //! @brief writes the help for the class
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const
      {
        OSTREAM << "Allowed values: {";
        // write list
        if( s_NumberOfEnums >= 1)
        {
          OSTREAM << ( *GetEnumDescriptor)( t_Enum( 0));
          //iterate over all enums from the GetFirst() to GetLast()
          for( size_t current_enum( 1); current_enum < s_NumberOfEnums; ++current_enum)
          {
            OSTREAM << ", " << ( *GetEnumDescriptor)( t_Enum( current_enum));
          }
        }
        OSTREAM << "}";
        if( s_NumberOfEnums > size_t( 10))
        {
          OSTREAM << '\n';
        }
        return OSTREAM;
      }

      //! @brief generate a vector of std::Strings for that enum
      //! @return Vector with strings representing each of the enums
      static const storage::Vector< std::string> &GetStringVector()
      {
        // generate vector of enum strings
        static const storage::Vector< std::string> s_string_vector
        (
          s_NumberOfEnums, &( *GetEnumDescriptor)( t_Enum( 0))
        );

        // return
        return s_string_vector;
      }

      //! @brief generate a vector of std::Strings for that enum
      //! @return Vector with strings representing each of the enums
      static const storage::Vector< WrapperEnum> &GetEnumVector()
      {
        // generate vector of enum strings
        static const storage::Vector< WrapperEnum> s_vector( GetEnumVectorImpl());

        // return
        return s_vector;
      }

    protected:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief set the value of the corresponding member based on the label
      //! @param VALUE value to set the enum up as
      //! @param ERR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryRead( const std::string &VALUE, std::ostream &ERR_STREAM)
      {
        //convert string to enum and assign to ENUMERATOR
        //iterate over all enums from the GetFirst() to GetLast()
        for( size_t current_enum( 0); current_enum < s_NumberOfEnums; ++current_enum)
        {
          //check that the STRING agrees with the description from the current enum
          if( ( *GetEnumDescriptor)( t_Enum( current_enum)) == VALUE)
          {
            //in case of agreement return the current enum
            m_Value = t_Enum( current_enum);
            return true;
          }
        }

        ERR_STREAM << "Tried to set enum of type " << GetStaticClassName< t_Enum>() << " with " << VALUE << ", ";
        WriteHelp( ERR_STREAM, 0);
        return false;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        io::Serialize::InsertIndent( OSTREAM, INDENT) << ( *GetEnumDescriptor)( m_Value);
        return OSTREAM;
      }

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        std::string value;
        ISTREAM >> value;
        BCL_Assert( TryRead( value, GetLogger()), "Bad enum value");
        return ISTREAM;
      }

      //! @brief read from std::istream
      //! This is the main function to be used for reading an object. It will read the first line assuming there is the
      //! class name as identifier - will be asserted, and after that it will call the overwritten Read function for the
      //! object.
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadObject( std::istream &ISTREAM)
      {
        // to maintain backwards-compatibility, do not write out the identifier
        return Read( ISTREAM);
      }

      //! @brief write to std::ostream
      //! This is the main function to be used for writing an object. It will write the classname as identifier,
      //! and after that it will call the overwritten Write function for the object
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &WriteObject( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // to maintain backwards-compatibility, do not write out the identifier
        return Write( OSTREAM, INDENT);
      }

      //! @brief create the enum vector
      //! @return the enum vector
      static storage::Vector< WrapperEnum> GetEnumVectorImpl()
      {
        storage::Vector< WrapperEnum> vector( s_NumberOfEnums);
        for( size_t i( 0); i < s_NumberOfEnums; ++i)
        {
          vector( i) = WrapperEnum( t_Enum( i));
        }
        return vector;
      }

    }; // template class WrapperEnum

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_WRAPPER_ENUM_H_

