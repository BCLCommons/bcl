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

#ifndef BCL_DESCRIPTOR_HAS_CACHE_H_
#define BCL_DESCRIPTOR_HAS_CACHE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from the bcl
#include "bcl_descriptor_cache_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HasCache
    //! @brief Mix-in interface for objects that have cached properties
    //!
    //! @see @link example_descriptor_has_cache.cpp @endlink
    //! @author mendenjl
    //! @date Nov 06, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_BaseClass>
    class HasCache :
      virtual public t_BaseClass
    {

    private:

    //////////
    // data //
    //////////

      mutable util::ShPtr< CacheMap> m_Cache; //!< Cached properties calculated for this object

    public:

    /////////////////
    // data access //
    /////////////////

      HasCache() :
        m_Cache( new CacheMap)
      {
      }
//
//      HasCache( const HasCache &A) :
//        m_Cache( new CacheMap)
//      {
//      }

      //! @brief cache a property
      //! @param NAME name of property
      //! @param PROPERTY vector of properties
      const CacheMap::value_type &Cache( const CacheMap::key_type &NAME, const CacheMap::value_type &PROPERTY) const
      {
        return m_Cache->Insert( NAME, PROPERTY);
      }

      //! @brief remove any property with NAME
      //! @param NAME the name of the property to remove
      void Uncache( const util::ObjectDataLabel &NAME) const
      {
        m_Cache->Remove( NAME);
      }

      //! @brief share the cache from another HasCache object
      void ShareCache( const HasCache &CACHE) const
      {
        m_Cache = CACHE.m_Cache;
      }

      //! @brief Reset the cache
      virtual void ResetCache() const
      {
        m_Cache->Reset();
      }

      //! @brief find the property for a given name
      //! @param NAME name of property
      //! @return pointer to the value, NULL if NAME was not in cache
      util::SiPtr< const CacheMap::value_type> FindInCache( const CacheMap::key_type &NAME) const
      {
        return m_Cache->Find( NAME);
      }

      //! @brief determine whether the given property exists either in the cache or store
      //! @param NAME name of property
      //! @return true iff the property exists for this molecule in the cache
      bool IsCached( const CacheMap::key_type &NAME) const
      {
        return m_Cache->Has( NAME);
      }

      //! @brief get a cached property
      //! @param NAME name of property
      //! @return values for that property
      const CacheMap::value_type &GetFromCache( const CacheMap::key_type &NAME) const
      {
        return m_Cache->Get( NAME);
      }

      //! @brief remove any property with NAME
      //! @param NAME the name of the property to remove
      //! @param REMOVE_DEPENDENT also remove any descriptors that contained NAME anywhere within their label
      void RemoveFromCache( const CacheMap::key_type &NAME, const bool &REMOVE_DEPENDENT = false) const
      {
        m_Cache->Remove( NAME, REMOVE_DEPENDENT);
      }

      //! @brief extract a given cache entry and any related descriptors from the cache
      //! @param NAME the name of the property to extract
      CacheMap ExtractRelatedCacheEntries( const CacheMap::key_type &NAME) const
      {
        return m_Cache->ExtractRelatedProperties( NAME);
      }

      //! @brief merge a different cache map with this cache map.  If an entry is in both maps, its value is retained
      //! @param CACHE_MAP the old cache map
      void MergeCache( const CacheMap &MAP) const
      {
        return m_Cache->Merge( MAP);
      }

      //! @brief try to cache all keys/values in the given map
      //! @param MAP map of keys to values
      void CacheNumeric( const storage::Map< std::string, std::string> &MAP) const
      {
        m_Cache->InsertNumericalMapItems( MAP);
      }

      //! @brief get the complete cache map
      util::ShPtr< CacheMap> &GetCacheMap() const
      {
        return m_Cache;
      }
//
//      //! @brief assignment operator
//      HasCache &operator=( const HasCache &A)
//      {
//        m_Cache = util::ShPtr< CacheMap>( new CacheMap);
//        return *this;
//      }

    }; // class HasCache

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_HAS_CACHE_H_
