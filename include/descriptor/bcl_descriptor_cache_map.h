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

#ifndef BCL_DESCRIPTOR_CACHE_MAP_H_
#define BCL_DESCRIPTOR_CACHE_MAP_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_object_data_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CacheMap
    //! @brief CacheMap for descriptors that are slower than O(N) to calculate
    //!
    //! @see @link example_descriptor_cache_map.cpp @endlink
    //! @author mendenjl
    //! @date Nov 06, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CacheMap :
      public util::ObjectInterface
    {
    public:

      //! typedef for the key type of the CacheMap
      typedef util::ObjectDataLabel key_type;

      //! typedef for the value type of the CacheMap
      typedef linal::Vector< float> value_type;

    private:

    //////////
    // data //
    //////////

      //! map of properties, where key is the property identifier/name and the value are the actual values
      storage::Map< key_type, value_type> m_PropertyToValueMap;

    public:

      //! typedef for the iterator into the map
      typedef storage::Map< key_type, value_type>::const_iterator const_iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CacheMap()
      {
      }

      //! @brief constructor from CacheMap and map of strings to strings, only takes properties that are numeric
      //! @param CacheMap the parent CacheMap to use; preferentially takes properties from here
      //! @param MAP map of keys to values
      CacheMap( const CacheMap &PARENT, const storage::Map< std::string, std::string> &MAP);

      //! @brief Clone function
      //! @return pointer to new CacheMap
      CacheMap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the number of properties in the CacheMap
      //! @return the number of properties in the CacheMap
      size_t GetSize() const;

      //! @brief get a specific property
      //! @param PROPERTY name of the desired property
      //! @return the vector of doubles for a specific property
      const value_type &Get( const key_type &NAME) const;

      //! @brief get the whole map
      const storage::Map< key_type, value_type> &GetMap() const
      {
        return m_PropertyToValueMap;
      }

      //! @brief find the property for a given name
      //! @param NAME name of property
      //! @return pointer to the value, NULL if NAME was not in CacheMap
      util::SiPtr< const value_type> Find( const key_type &NAME) const;

      //! @brief check if there a property in the CacheMap with the given name
      //! @param NAME name of property
      //! @return true if property with that name exists, false otherwise
      bool Has( const key_type &NAME) const;

      //! @brief try to CacheMap all keys/values in the given map
      //! @param MAP map of keys to values
      void InsertNumericalMapItems( const storage::Map< std::string, std::string> &MAP);

      //! @brief try to CacheMap the given name and value string if the property is not already in the CacheMap
      //! @param NAME the name string to attempt to CacheMap
      //! @param VALUE the value string to attempt to CacheMap
      //! @return true if the value could be CacheMapd
      bool InsertIfNumerical( const std::string &KEY, const std::string &VALUE);

    ////////////////
    // operations //
    ////////////////

      //! @brief set a property
      //! @brief NAME name of property
      //! @param PROPERTY vector of properties
      //! @return reference to the value inserted
      const value_type &Insert( const key_type &NAME, const value_type &PROPERTY);

      //! @brief remove any property with NAME
      //! @param NAME the name of the property to remove
      //! @param REMOVE_DEPENDENT also remove any descriptors that contained NAME anywhere within their label
      void Remove( const key_type &NAME, const bool &REMOVE_DEPENDENT = false);

      //! @brief extract any property with NAME into a new cache map
      //! @param NAME the name of the property to extract
      CacheMap ExtractRelatedProperties( const key_type &NAME);

      //! @brief merge a different cache map with this cache map.  If an entry is in both maps, its value is retained
      //! @param CACHE_MAP the old cache map
      void Merge( const CacheMap &MAP);

      //! @brief Reset the CacheMap
      void Reset();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class CacheMap

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_CACHE_MAP_H_

