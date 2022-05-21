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

#ifndef BCL_UTIL_ENUMERATE_HPP_
#define BCL_UTIL_ENUMERATE_HPP_

// includes from bcl - sorted alphabetically
#include "bcl_util.h"
#include "bcl_util_enums_instances.h"
#include "bcl_util_message.h"
#include "command/bcl_command_guesser.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <iostream>
#include <iterator>

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct Enums - only used by derived classes
    //! constructs a default constructed Enum as very first and undefined Enum
    //! @param IS_READ_WRITABLE declare that Enum derived instance as read and writable for the user over the
    //!        EnumsInstances and ObjectInstances mechanisms
    template< typename t_DataType, typename t_Derived>
    Enumerate< t_DataType, t_Derived>::Enumerate( const bool IS_READ_WRITABLE) :
      m_DataList(),
      m_EnumVector(),
      e_Undefined( const_cast< EnumDataType &>( GetUndefinedData()))
    {
      // if object is supposed to be read and writable, add it to the EnumsInstances (for reading) and in the Insert
      // function also to the ObjectInstances (for writing)
      if( IS_READ_WRITABLE)
      {
        // insert that instance
        EnumsInstances::GetEnumsInstances().Insert( *this, GetStaticClassName< t_Derived>());
      }
    }

    //! virtual copy constructor
    template< typename t_DataType, typename t_Derived>
    Enumerate< t_DataType, t_Derived> *Enumerate< t_DataType, t_Derived>::Clone() const
    {
      BCL_MessageCrt( "cloning of Enumerate is not supported")
      return NULL;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType, typename t_Derived>
    const std::string &Enumerate< t_DataType, t_Derived>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the only instance of the derived enumerators
    //! @return reference to singleton instance of t_Derived
    template< typename t_DataType, typename t_Derived>
    t_Derived &Enumerate< t_DataType, t_Derived>::GetEnums()
    {
      static t_Derived *s_Enum( new t_Derived());
      return *s_Enum;
    }

    //! @brief get strings of all enums in a vector
    //! @return std::vector filled with strings of all enums, sorted by enum index
    template< typename t_DataType, typename t_Derived>
    std::vector< std::string> Enumerate< t_DataType, t_Derived>::GetEnumStrings() const
    {
      return std::vector< std::string>( m_EnumVector.begin(), m_EnumVector.end());
    }

  ////////////////
  // operations //
  ////////////////

    //! protected constructor, no creation of objects outside of the class
    //! @param NAME   name of the current enum
    //! @param OBJECT object to be enumerated
    template< typename t_DataType, typename t_Derived>
    Enum< t_DataType, t_Derived> &Enumerate< t_DataType, t_Derived>::AddEnum
    (
      const std::string &NAME,
      const t_DataType &OBJECT
    )
    {
      // check that nobody uses name "Undefined"
      BCL_Assert
      (
        NAME != EnumDataType::GetUndefinedEnumName(),
        "name: \"" + EnumDataType::GetUndefinedEnumName() + "\" is reserved for default undefined Enum!"
      );

      // check that enum with this name does not exist yet
      if( HaveEnumWithName( NAME))
      {
        // Use std::cout since this loop could be reached while adding the logger enum,
        // in which case bcl_message will recurse in an infinite loop
        std::cerr << "Enumerator named: \"" + NAME + "\" already exists" << '\n';
        return *GetEnumIteratorFromName( NAME);
      }

      // add the enum to the list
      m_DataList.push_back( EnumDataType( m_DataList.size(), NAME, OBJECT));

      // add the enum to the vector
      m_EnumVector.push_back( Enum< t_DataType, t_Derived>( m_DataList.back()));

      // return refernce to the last inserted enum
      return m_EnumVector.back();
    }

    //! begin iterator for iterating over enums
    //! @return const_iterator to begin of all enums
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::const_iterator Enumerate< t_DataType, t_Derived>::Begin() const
    {
      return m_EnumVector.begin();
    }

    //! begin iterator for iterating over enums
    //! @return iterator to begin of all enums
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::iterator Enumerate< t_DataType, t_Derived>::Begin()
    {
      return m_EnumVector.begin();
    }

    //! end iterator for iterating over enums
    //! @return const_iterator to end of all enums
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::const_iterator Enumerate< t_DataType, t_Derived>::End() const
    {
      return m_EnumVector.end();
    }

    //! end iterator for iterating over enums
    //! @return iterator to end of all enums
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::iterator Enumerate< t_DataType, t_Derived>::End()
    {
      return m_EnumVector.end();
    }

    //! get number of enum objects
    //! @return number of enums
    template< typename t_DataType, typename t_Derived>
    size_t Enumerate< t_DataType, t_Derived>::GetEnumCount() const
    {
      return m_EnumVector.size();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get iterator to enum object from index of enum
    //! @param INDEX index of the enum to be searched for
    //! @return iterator to the enum with the INDEX - will point to end if there is no such enum with the INDEX
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::iterator Enumerate< t_DataType, t_Derived>::GetEnumIteratorFromIndex( const size_t INDEX)
    {
      // check that the INDEX is valid
      if( !HaveEnumWithIndex( INDEX))
      {
        if( INDEX != e_Undefined.GetIndex())
        {
          BCL_MessageCrt
          (
            "invalid enum INDEX for " + GetClassIdentifier() + ": "
            + Format()( INDEX) + " > " + Format()( GetEnumCount())
          );
        }
        // return to the end of the vector
        return m_EnumVector.end();
      }

      // return
      return m_EnumVector.begin() + INDEX;
    }

    //! @brief get const_iterator to enum object from index of enum
    //! @param INDEX index of the enum to be searched for
    //! @return iterator to the enum with the INDEX - will point to end if there is no such enum with the INDEX
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::const_iterator Enumerate< t_DataType, t_Derived>::GetEnumIteratorFromIndex( const size_t INDEX) const
    {
      // check that the INDEX is valid
      if( !HaveEnumWithIndex( INDEX))
      {
        if( INDEX != e_Undefined.GetIndex())
        {
          BCL_MessageCrt
          (
            "invalid enum INDEX for " + GetClassIdentifier() + ": "
            + Format()( INDEX) + " > " + Format()( GetEnumCount())
          );
        }
        // return the end of the vector
        return m_EnumVector.end();
      }

      // return
      return m_EnumVector.begin() + INDEX;
    }

    //! @brief get iterator to enum from name
    //! @param NAME name of the enum to be searched for
    //! @return iterator to the enum with the NAME - will point to end if there is no such enum with the NAME
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::iterator Enumerate< t_DataType, t_Derived>::GetEnumIteratorFromName( const std::string &NAME)
    {
      // find an enum with the name
      iterator itr
      (
        std::find_if( m_EnumVector.begin(), m_EnumVector.end(), std::bind2nd( std::equal_to< std::string>(), NAME))
      );

      // no iterator with the given NAME found - issue user warning
      if( itr == m_EnumVector.end() && NAME != EnumDataType::GetUndefinedEnumName())
      {
        BCL_MessageDbg( "there is no such enum with the name \"" + NAME + "\"");
      }

      // return iterator to end
      return itr;
    }

    //! @brief get const_iterator to enum from name
    //! @param NAME name of the enum to be searched for
    //! @return iterator to the enum with the NAME - will point to end if there is no such enum with the NAME
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::const_iterator Enumerate< t_DataType, t_Derived>::GetEnumIteratorFromName( const std::string &NAME) const
    {
      // find an enum with the name
      const_iterator itr
      (
        std::find_if( m_EnumVector.begin(), m_EnumVector.end(), std::bind2nd( std::equal_to< std::string>(), NAME))
      );

      // no iterator with the given NAME found - issue user warning
      if( itr == m_EnumVector.end() && NAME != EnumDataType::GetUndefinedEnumName())
      {
        BCL_MessageDbg( "there is no such enum with the name \"" + NAME + "\"");
      }

      // return iterator to end
      return itr;
    }

    //! @brief set an enum to the enum with a given name
    //! @param ENUM enum to set
    //! @param NAME name of the enum to be searched for
    //! @param ERR_STREAM error output stream, for if name is not an enum
    //! @return true on success
    template< typename t_DataType, typename t_Derived>
    bool Enumerate< t_DataType, t_Derived>::SetEnumFromName( EnumType &ENUM, const std::string &NAME, std::ostream &ERR_STREAM) const
    {
      // find an enum with the name
      const_iterator itr
      (
        std::find_if( m_EnumVector.begin(), m_EnumVector.end(), std::bind2nd( std::equal_to< std::string>(), NAME))
      );

      // no iterator with the given NAME found - issue user warning
      if( itr == m_EnumVector.end())
      {
        if( NAME != EnumDataType::GetUndefinedEnumName())
        {
          command::Guesser::GetDefaultGuesser().WriteGuesses
          (
            NAME,
            GetEnumStrings(),
            ERR_STREAM,
            GetStaticClassName< EnumType>()
          );

          return false;
        }
        else
        {
          ENUM = e_Undefined;
        }
      }
      else
      {
        ENUM = *itr;
      }

      // return true for successful set
      return true;
    }

    //! @brief Get an undefined EnumData of this type
    //! @return const reference to undefined Data
    template< typename t_DataType, typename t_Derived>
    const EnumData< t_DataType> &Enumerate< t_DataType, t_Derived>::GetUndefinedData()
    {
      return EnumData< t_DataType>::GetUndefinedData();
    }

    //! @brief Get an undefined Enum of this type
    //! @brief return reference to undefined enum, which is e_Undefined
    template< typename t_DataType, typename t_Derived>
    const Enum< t_DataType, t_Derived> &Enumerate< t_DataType, t_Derived>::GetUndefined() const
    {
      return e_Undefined;
    }

    //! @brief get iterator to enum from enum instance
    //! @param ENUM enum to be searched for
    //! @return iterator to the provided enum
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::iterator Enumerate< t_DataType, t_Derived>::GetEnumIteratorFromEnum( const Enum< t_DataType, t_Derived> &ENUM)
    {
      // use the index of that enum
      return GetEnumIteratorFromIndex( ENUM->m_Index);
    }

    //! @brief get iterator to enum from enum instance
    //! @param ENUM enum to be searched for
    //! @return const_iterator to the provided enum
    template< typename t_DataType, typename t_Derived>
    typename std::vector< Enum< t_DataType, t_Derived> >::const_iterator Enumerate< t_DataType, t_Derived>::GetEnumIteratorFromEnum( const Enum< t_DataType, t_Derived> &ENUM) const
    {
      // use the index of that enum
      return GetEnumIteratorFromIndex( ENUM->m_Index);
    }

    //! @brief get enum from index
    //! @param INDEX index of the enum
    //! @return const enum reference with the index, undefined if non such enum exists
    template< typename t_DataType, typename t_Derived>
    const Enum< t_DataType, t_Derived> &Enumerate< t_DataType, t_Derived>::GetEnumFromIndex( const size_t INDEX) const
    {
      // get iterator for that enum
      const_iterator itr( GetEnumIteratorFromIndex( INDEX));

      // enum with that index was not found
      if( itr == End())
      {
        return e_Undefined;
      }

      // if itr points to valid entry
      return *itr;
    }

    //! @brief get enum from index
    //! @param INDEX index of the enum
    //! @return enum reference with the index, undefined if non such enum exists
    template< typename t_DataType, typename t_Derived>
    Enum< t_DataType, t_Derived> &Enumerate< t_DataType, t_Derived>::GetEnumFromIndex( const size_t INDEX)
    {
      // get iterator for that enum
      iterator itr( GetEnumIteratorFromIndex( INDEX));

      // enum with that index was not found
      if( itr == End())
      {
        return e_Undefined;
      }

      // if itr points to valid entry
      return *itr;
    }

    //! @brief get enum from name
    //! @param NAME name of the enum to be searched for
    //! @return the const reference to enum with the NAME - will give undefined enum if no enum with NAME was found
    template< typename t_DataType, typename t_Derived>
    const Enum< t_DataType, t_Derived> &Enumerate< t_DataType, t_Derived>::GetEnumFromName( const std::string &NAME) const
    {
      // t_Derived will be statically initialized here if it isn't already
      const_iterator itr( GetEnumIteratorFromName( NAME));

      // enum with that name was not found
      if( itr == End())
      {
        return e_Undefined;
      }

      // if itr points to valid entry
      return *itr;
    }

    //! @brief get enum from name
    //! @param NAME name of the enum to be searched for
    //! @return the reference to enum with the NAME - will give undefined enum if no enum with NAME was found
    template< typename t_DataType, typename t_Derived>
    Enum< t_DataType, t_Derived> &Enumerate< t_DataType, t_Derived>::GetEnumFromName( const std::string &NAME)
    {
      // t_Derived will be statically initialized here if it isn't already
      iterator itr( GetEnumIteratorFromName( NAME));

      // enum with that name was not found
      if( itr == End())
      {
        return e_Undefined;
      }

      // if itr points to valid entry
      return *itr;
    }

    //! @brief check if an enum with that INDEX exists
    //! @param INDEX index of suspected enum
    //! @return true if an enum with INDEX exists, false if there is no such enum or if INDEX is undefined
    template< typename t_DataType, typename t_Derived>
    bool Enumerate< t_DataType, t_Derived>::HaveEnumWithIndex( const size_t INDEX) const
    {
      return INDEX < GetEnumCount();
    }

    //! @brief check if an enum with the NAME exists
    //! @param NAME name of suspected enum
    //! @return true if an enum with NAME exists alos, if it is the name of the undefined enum, false otherwise
    template< typename t_DataType, typename t_Derived>
    bool Enumerate< t_DataType, t_Derived>::HaveEnumWithName( const std::string &NAME) const
    {
      // check if the NAME is that of the undefined enum
      // otherwise, if there is an enum with that name in the enum vector
      return
      (
           NAME == EnumDataType::GetUndefinedEnumName()
        || (
             std::find_if
             (
               m_EnumVector.begin(), m_EnumVector.end(), std::bind2nd( std::equal_to< std::string>(), NAME)
             ) != m_EnumVector.end()
           )
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the list of enums
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the list was written to
    template< typename t_DataType, typename t_Derived>
    std::ostream &Enumerate< t_DataType, t_Derived>::WriteList( std::ostream &OSTREAM) const
    {
      OSTREAM << '{';
      // copy to output stream
      const_iterator itr( Begin()), itr_end( End());
      if( itr != itr_end)
      {
        OSTREAM << ' ' << itr->GetName();
        //iterate over all enums from the GetFirst() to GetLast()
        for( ++itr; itr != itr_end; ++itr)
        {
          OSTREAM << ", " << itr->GetName();
        }
      }
      // write closing parentheses
      OSTREAM << '}';
      return OSTREAM;
    }

    //! @brief write Enums to std::ostream
    //! @param OSTREAM ostream to write to
    //! @param INDENT the indent to be used for each line
    //! @return ref to the stream written to
    template< typename t_DataType, typename t_Derived>
    std::ostream &Enumerate< t_DataType, t_Derived>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member.  t_Derived will be statically initailized here if it is not already
      io::Serialize::Write( GetEnumCount(), OSTREAM, INDENT) << '\n';

      // write all enum data
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        io::Serialize::Write( ( *itr)->m_Name   , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( t_DataType( **itr), OSTREAM, INDENT) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief read Enums from io::IFStream
    //! @param ISTREAM istream to read from
    //! @return ref to the stream written to
    template< typename t_DataType, typename t_Derived>
    std::istream &Enumerate< t_DataType, t_Derived>::Read( std::istream &ISTREAM)
    {
      // read the number of enums to be read
      size_t enum_count;
      io::Serialize::Read( enum_count, ISTREAM);

      // count how many have been read so far
      for( size_t count( 0); ISTREAM.good() && count < enum_count; ++count)
      {
        std::string name;
        io::Serialize::Read( name, ISTREAM);

        t_DataType current_data;
        io::Serialize::Read( current_data, ISTREAM);

        // search for that enum with that name
        iterator itr( GetEnumIteratorFromName( name));

        // found already existing enum, overwrite the exiting enum with the current data
        if( itr != m_EnumVector.end())
        {
          **itr = current_data;
        }
        // no such EnumData was found - insert it as an additional Enum
        else
        {
          AddEnum( name, current_data);
        }
      }

      // end
      return ISTREAM;
    }

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_ENUMERATE_HPP_
