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
#include "descriptor/bcl_descriptor_cache_map.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CacheMap::s_Instance
    (
      GetObjectInstances().AddInstance( new CacheMap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from CacheMap and map of strings to strings, only takes properties that are numeric
    //! @param CacheMap the parent CacheMap to use; preferentially takes properties from here
    //! @param MAP map of keys to values
    CacheMap::CacheMap( const CacheMap &PARENT, const storage::Map< std::string, std::string> &MAP) :
      m_PropertyToValueMap( PARENT.m_PropertyToValueMap)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CacheMap
    CacheMap *CacheMap::Clone() const
    {
      return new CacheMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CacheMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the number of properties in the CacheMap
    //! @return the number of properties in the CacheMap
    size_t CacheMap::GetSize() const
    {
      return m_PropertyToValueMap.GetSize();
    }

    //! @brief get a specific property
    //! @param PROPERTY name of the desired property
    //! @return the vector of doubles for a specific property
    const CacheMap::value_type &CacheMap::Get( const CacheMap::key_type &NAME) const
    {
      // create an empty static vector for the case that the property is not found
      static const value_type s_empty_property;

      // find the property
      const_iterator itr_prop( m_PropertyToValueMap.Find( NAME));

      // if the property was not found, return the empty value
      if( itr_prop == m_PropertyToValueMap.End())
      {
        return s_empty_property;
      }

      // return the found property
      return itr_prop->second;
    }

    //! @brief find the property for a given name
    //! @param NAME name of property
    //! @return pointer to the value, NULL if NAME was not in CacheMap
    util::SiPtr< const CacheMap::value_type> CacheMap::Find( const key_type &NAME) const
    {
      // find the property
      const_iterator itr_prop( m_PropertyToValueMap.Find( NAME));

      // if the property was not found, return the empty value
      if( itr_prop == m_PropertyToValueMap.End())
      {
        return util::SiPtr< const CacheMap::value_type>();
      }

      return util::ToSiPtr( itr_prop->second);
    }

    //! @brief check if there a property in the CacheMap with the given name
    //! @param NAME name of property
    //! @return true if property with that name exists, false otherwise
    bool CacheMap::Has( const key_type &NAME) const
    {
      return m_PropertyToValueMap.Has( NAME);
    }

    //! @brief try to CacheMap all keys/values in the given map
    //! @param MAP map of keys to values
    void CacheMap::InsertNumericalMapItems( const storage::Map< std::string, std::string> &MAP)
    {
      for
      (
        storage::Map< std::string, std::string>::const_iterator itr( MAP.Begin()), itr_end( MAP.End());
        itr != itr_end;
        ++itr
      )
      {
        InsertIfNumerical( itr->first, itr->second);
      }
    }

    //! @brief try to CacheMap the given name and value string
    //! @param NAME the name string to attempt to CacheMap
    //! @param VALUE the value string to attempt to CacheMap
    //! @return true if the value could be CacheMapd
    bool CacheMap::InsertIfNumerical( const std::string &KEY, const std::string &VALUE)
    {
      // test whether the given values are numeric
      const storage::Vector< std::string> result
      (
        util::SplitString( util::TrimString( VALUE), " \t\n\r,")
      );
      value_type converted( result.GetSize());
      value_type::iterator itr_descriptor( converted.Begin());
      std::stringstream err_stream;
      // iterate over all strings and convert them to double
      for
      (
        storage::Vector< std::string>::const_iterator itr_str( result.Begin()), itr_str_end( result.End());
        itr_str != itr_str_end;
        ++itr_str, ++itr_descriptor
      )
      {
        // read the input value; if double, could be a nan so we need to use io::Serialize
        if( !util::TryConvertFromString( *itr_descriptor, *itr_str, err_stream))
        {
          return false;
        }
      }

      key_type label;
      if( !label.TryAssign( KEY, err_stream))
      {
        // if the label cannot be assigned, it means the name contains some special delimiters,
        // so set the value directly
        label = key_type();
        label.SetValue( KEY);
      }

      return m_PropertyToValueMap.Insert( std::make_pair( label, converted)).second;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set a property
    //! @brief NAME name of property
    //! @param PROPERTY vector of properties
    //! @return reference to the value inserted
    const CacheMap::value_type &CacheMap::Insert( const key_type &NAME, const value_type &PROPERTY)
    {
      return m_PropertyToValueMap[ NAME] = PROPERTY;
    }

    //! @brief remove any property with NAME
    //! @param NAME the name of the property to remove
    //! @param REMOVE_DEPENDENT also remove any descriptors that contained NAME anywhere within their label
    void CacheMap::Remove( const key_type &NAME, const bool &REMOVE_DEPENDENT)
    {
      storage::Map< key_type, value_type>::iterator itr( m_PropertyToValueMap.Find( NAME));
      if( itr != m_PropertyToValueMap.End())
      {
        m_PropertyToValueMap.RemoveElement( itr);
      }

      if( REMOVE_DEPENDENT)
      {
        for
        (
          storage::Map< key_type, value_type>::iterator
            itr( m_PropertyToValueMap.Begin()), itr_end( m_PropertyToValueMap.End());
          itr != itr_end;
          // iteration in loop
        )
        {
          if( itr->first.Find( NAME, true, false) != itr->first.End())
          {
            // property contained NAME internally; remove it
            storage::Map< key_type, value_type>::iterator itr_old( itr);
            ++itr;
            m_PropertyToValueMap.RemoveElement( itr_old);
          }
          else
          {
            // property is not a component of this descriptor; move on to the next
            ++itr;
          }
        }
      }
    }

    //! @brief extract any property with NAME into a new cache map
    //! @param NAME the name of the property to extract
    CacheMap CacheMap::ExtractRelatedProperties( const key_type &NAME)
    {
      CacheMap new_map;
      storage::Map< key_type, value_type>::iterator itr( m_PropertyToValueMap.Find( NAME));
      if( itr != m_PropertyToValueMap.End())
      {
        new_map.m_PropertyToValueMap.InsertElement( *itr);
        m_PropertyToValueMap.RemoveElement( itr);
      }

      for
      (
        storage::Map< key_type, value_type>::iterator
          itr( m_PropertyToValueMap.Begin()), itr_end( m_PropertyToValueMap.End());
        itr != itr_end;
        // iteration in loop
      )
      {
        if( itr->first.Find( NAME, true, false) != itr->first.End())
        {
          // property contained NAME internally; extract it
          new_map.m_PropertyToValueMap.InsertElement( *itr);
          storage::Map< key_type, value_type>::iterator itr_old( itr);
          ++itr;
          m_PropertyToValueMap.RemoveElement( itr_old);
        }
        else
        {
          // property is not a component of this descriptor; move on to the next
          ++itr;
        }
      }
      return new_map;
    }

    //! @brief merge a different cache map with this cache map.  If an entry is in both maps, its value is retained
    //! @param CACHE_MAP the old cache map
    void CacheMap::Merge( const CacheMap &MAP)
    {
      m_PropertyToValueMap.InsertElements( MAP.m_PropertyToValueMap.Begin(), MAP.m_PropertyToValueMap.End());
    }

    //! @brief Reset the CacheMap
    void CacheMap::Reset()
    {
      m_PropertyToValueMap.Reset();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CacheMap::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_PropertyToValueMap, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CacheMap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_PropertyToValueMap, OSTREAM, INDENT);
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace descriptor
} // namespace bcl
