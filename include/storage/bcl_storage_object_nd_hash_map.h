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

#ifndef BCL_STORAGE_OBJECT_ND_HASH_MAP_H_
#define BCL_STORAGE_OBJECT_ND_HASH_MAP_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_storage_hash_map.h"
#include "bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectNDHashMap
    //! @brief This is the template class for storing data for any number of objects in a Hashmap.
    //! @details attributes:
    //! a.) this is a specialized type of hashmap
    //! b.) key is a size_t which is generated as a bit shift and xor operation on the addresses of the objects desired
    //!     to be a part of the key
    //! c.) allows multiple objects to be hashed into a single key
    //!
    //! TODO: fix iterators so that they are the correct name type
    //!
    //! @see @link example_storage_object_nd_hash_map.cpp @endlink
    //! @author karakam
    //! @date 04/09/06
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! @tparam t_NumberOfObjects size_t which indicates the number of objects
    //! @tparam t_ObjectType which indicates the type of objects that correspond to the data
    //! @tparam t_DataType which indicates the type of data that will be held
    //! @tparam t_SymmetricOrNot bool which indicates whether the objects should be symmetric or not
    template< size_t t_NumberOfObjects, typename t_ObjectType, typename t_DataType>
    class ObjectNDHashMap :
      public HashMap< size_t, t_DataType>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      ObjectNDHashMap()
      {
      }

      //! copy constructor
      ObjectNDHashMap *Clone() const
      {
        return new ObjectNDHashMap( *this);
      }

      //! @brief destrcutor
      virtual ~ObjectNDHashMap()
      {
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief GetClassIdentifier returns class name
      //! @return returns the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Insert inserts data into HashMap
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > to be inserted
      //! @param DATA t_DataType to be inserted with OBJECTS
      //! @return returns a pair of iterator( position of insertion) and bool( success or not)
      std::pair< typename HashMap< size_t, t_DataType>::iterator, bool>
      Insert( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS, const t_DataType &DATA)
      {
        return Insert( HashKey( OBJECTS), DATA);
      }

      //! @brief Insert inserts data into HashMap with additional information used for hashkey
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > to be inserted
      //! @param INFO size_t which is the additional information used for the hashkey
      //! @param DATA t_DataType to be inserted with OBJECTS
      //! @return returns a pair of iterator( position of insertion) and bool( success or not)
      std::pair< typename HashMap< size_t, t_DataType>::iterator, bool>
      Insert
      (
        const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS,
        const size_t &INFO,
        const t_DataType &DATA
      )
      {
        return Insert( HashKey( OBJECTS, INFO), DATA);
      }

      //! Orders Objects and Inserts the data into HashMap
      //! @brief InsertAndOrderObjects orders Objects and Inserts the data into HashMap
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > to be inserted
      //! @param DATA t_DataType to be inserted with OBJECTS
      //! @return returns a pair of iterator( position of insertion) and bool( success or not)
      std::pair< typename HashMap< size_t, t_DataType>::iterator, bool>
      InsertAndOrderObjects
      (
        VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS,
        const t_DataType &DATA
      )
      {
        return Insert( HashKeyAndOrderObjects( OBJECTS), DATA);
      }

      //! Orders Objects and Inserts the data into HashMap with additional information used for hashkey
      //! @brief InsertAndOrderObjects orders Objects and Inserts the data into HashMap with additional information
      //! used for hashkey
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > to be inserted
      //! @param INFO size_t which is the additional information used for the hashkey
      //! @param DATA t_DataType to be inserted with OBJECTS
      //! @return returns a pair of iterator( position of insertion) and bool( success or not)
      std::pair< typename HashMap< size_t, t_DataType>::iterator, bool>
      InsertAndOrderObjects
      (
        VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS,
        const size_t &INFO,
        const t_DataType &DATA
      )
      {
        return Insert( HashKeyAndOrderObjects( OBJECTS, INFO), DATA);
      }

      //! @brief Find finds the data stored for objects
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > for which the data is desired
      //! @return returns an iterator to the desired object and data
      typename HashMap< size_t, t_DataType>::iterator
      Find( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS)
      {
        return HashMap< size_t, t_DataType>::m_Data.find( HashKey( OBJECTS));
      }

      //! @brief Find finds the data stored for objects and additional info
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > for which the data is desired
      //! @param INFO size_t which is the additional information desired
      //! @return returns an iterator to the desired object, info, and data
      typename HashMap< size_t, t_DataType>::iterator
      Find( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS, const size_t &INFO)
      {
        return HashMap< size_t, t_DataType>::m_Data.find( HashKey( OBJECTS, INFO));
      }

      //! @brief Find finds the data stored for objects
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > for which the data is desired
      //! @return returns a const_iterator to the desired object and data
      typename HashMap< size_t, t_DataType>::const_iterator
      Find( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS) const
      {
        return HashMap< size_t, t_DataType>::m_Data.find( HashKey( OBJECTS));
      }

      //! @brief Find finds the data stored for objects and additional info
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > for which the data is desired
      //! @param INFO size_t which is the additional information desired
      //! @return returns a const_iterator to the desired object, info, and data
      typename HashMap< size_t, t_DataType>::const_iterator
      Find( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS, const size_t &INFO) const
      {
        return HashMap< size_t, t_DataType>::m_Data.find( HashKey( OBJECTS, INFO));
      }

      // TODO: those functions should actually return a ref to the data, but if it cant find the data, it would need to
      // return a reference to an Undefined< t_DataType>::s_Undefined - but it is not implemented for most of datatypes yet
      //    //! Returns the data inserted for OBJECTS and INFO
      //    t_DataType GetData( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS) const
      //    {
      //      typename HashMap< size_t, t_DataType>::const_iterator itr = HashMap< size_t, t_DataType>::Find
      //      (
      //        HashKey( OBJECTS)
      //      );
      //      if( itr == HashMap< size_t, t_DataType>::End())
      //      {
      //        return t_DataType();
      //      }
      //      return itr->second;
      //    }
      //
      //    //! Returns the key inserted for OBJECTS and INFO
      //    t_DataType
      //    GetData
      //    (
      //      const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS, const size_t &INFO
      //    ) const
      //    {
      //      typename HashMap< size_t, t_DataType>::const_iterator itr = Find( HashKey( OBJECTS, INFO));
      //      if( itr == HashMap< size_t, t_DataType>::End())
      //      {
      //        return t_DataType();
      //      }
      //      return *itr;
      //    }

      //! @brief operator [] allows access to the data stored for a object analogous to operator of Hashmap
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > for which access is desired
      //! @return returns t_DataType which is the data accessed by OBJECTS
      virtual t_DataType &operator []( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS)
      {
        return HashMap< size_t, t_DataType>::operator[]( HashKey( OBJECTS));
      }

      //! @brief HashKey returns the Hashkey which made of address of t_NumberOfObjects objects
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > of the desired  HashKey
      //! @return returns size_t which is the address of OBJECTS
      static size_t HashKey( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS)
      {
        // create size_t "temp" which is an array to store addresses of t_NumberOfObjects objects
        size_t temp[ t_NumberOfObjects];

        // create VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > and order "ordered_objects" by address
        VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > ordered_objects( OBJECTS);

        // converts all the addresses to size_t
        for( size_t i( 0); i < t_NumberOfObjects; ++i)
        {
          temp[ i] = size_t( ordered_objects( i).GetPointer());
        }

        // return size_t which is address of OBJECTS
        return DoShiftForHashKey( temp);
      };

      //! @brief HashKey function for storing an additional information than just the addresses of objects in your key
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > for which access is desired
      //! @param INFO size_t which is the additional information
      //! @return returns size_t which is the address of OBJECTS
      static size_t
      HashKey( const VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS, const size_t &INFO)
      {
        // create size_t "temp" which is an array to store addresses of t_NumberOfObjects objects
        size_t temp[ t_NumberOfObjects + 1];

        // create VectorND "ordered_objects" and initialize with OBJECTS
        VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > ordered_objects( OBJECTS);

        // converts all the addresses to size_t
        for( size_t i( 0); i < t_NumberOfObjects; ++i)
        {
          temp[ i] = size_t( ordered_objects( i).GetPointer());
        }
        temp[ t_NumberOfObjects] = INFO;

        // return size_t array which is addresses of OBEJCTS
        return DoShiftForHashKey( temp, t_NumberOfObjects + 1);
      };

      //! @brief HashKeyAndOrderObjects orders objects and returns Hashkey made of address of t_NumberOfObjects objects
      //! @param OBJECTS VectorND to be ordered and hashkey given
      //! @return return size_t array which contains the addresses of the ordered objects in OBJECTS
      static size_t HashKeyAndOrderObjects( VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS)
      {
        // create size_t array "temp" to store addresses of t_NumberOfObjects objects
        size_t temp[ t_NumberOfObjects];

        // converts all the addresses to size_t
        for( size_t i( 0); i < t_NumberOfObjects; ++i)
        {
          temp[ i] = size_t( OBJECTS( i).GetPointer());
        }

        // do bit shift on addresses of OBJECTS and obtain a hashkey then return this hashkey size_t
        return DoShiftForHashKey( temp);
      };

      //! @brief HashKeyAndOrderObjects orders objects and overloaded HashKey function for cases where you want to store
      //! additional information beyond just the addresses of objects in your key
      //! @param OBJECTS VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > to be ordered and hashkey given
      //! @param INFO size_t which is the extra information
      //! @return returns size_t which is the array of addresses of the ordered OBJECTS
      static size_t
      HashKeyAndOrderObjects
      (
        VectorND< t_NumberOfObjects, util::SiPtr< const t_ObjectType> > &OBJECTS,
        const size_t &INFO
      )
      {
        // create size_t array "temp" to store addresses of t_NumberOfObjects + 1 INFO
        size_t temp[ t_NumberOfObjects + 1];

        // converts all the addresses to size_t
        for( size_t i( 0); i < t_NumberOfObjects; ++i)
        {
          temp[ i] = size_t( OBJECTS( i).GetPointer());
        }

        // set single value of array equal to INFO
        temp[ t_NumberOfObjects] = INFO;

        // do bit shift on addresses of OBJECTS and obtain a hashkey then return this hashkey size_t
        return DoShiftForHashKey( temp, t_NumberOfObjects + 1);
      };

    private:

      //! @brief DoShiftForHashKey shifts every address by a certain amount and returns their accumulated xor
      //! @param ARR contains the size_t addresses that need to be shifted
      //! @param NUMBER the number of addresses that needs to be shifted (defaulted to t_NumberOfObjects)
      //! @return returns size_t which is the accumulated xor combinations of ARR
      static size_t DoShiftForHashKey( size_t *ARR, const size_t NUMBER = t_NumberOfObjects)
      {
        // use bit rotations to store addresses in one number
        for( size_t i( 1); i < NUMBER; ++i)
        {
          // create const shorts "shift1" and "shift2" which define the bit shifts below
          const short shift1( ( ( i * 3) % 8) * g_SizeOfAddress), shift2( g_SizeOfSystem - shift1);

          // do bit shifts with either "shift1" or "shift2"
          ARR[ i] = ( ARR[ i] << shift1) | ( ARR[ i] >> shift2);

          // define first position of array as xor combination of ARR[i] and ARR[0]
          ARR[ 0] ^= ARR[ i];
        }

        // return first position of ARR which is a size_t of the accumulated xor combinations of ARR
        return ARR[ 0];
      };

      //! @brief operator [] hides this operator so that users are forced to use Insert function defined in this class
      //! @param KEY the key that should be inserted or located in the Hashmap
      //! @return returns reference to t_DataType
      virtual t_DataType &operator []( const size_t KEY)
      {
        return HashMap< size_t, t_DataType>::operator[]( KEY);
      }

      //! @brief Insert hides the virtual insert function defined in base class HashMap and
      //! forces user to use the insert function defined in this class
      //! @param KEY the key to be inserted into the Hashmap
      //! @param DATA the data that will be associated with KEY
      //! @return returns pair which is an iterator to off-the-end or KEY/DATA and a bool indicating the success or not
      virtual std::pair< typename HashMap< size_t, t_DataType>::iterator, bool>
      Insert( const size_t &KEY, const t_DataType &DATA)
      {
        return HashMap< size_t, t_DataType>::Insert( std::pair< size_t, t_DataType>( KEY, DATA));
      }

      //! @brief Insert hides the virtual insert function defined in base class HashMap and
      //! forces user to use the insert function defined in this class
      //! @param KEY_DATA_PAIR the pair of a key and a data to be inserted into the Hashmap
      //! @return returns pair which is an iterator to off-the-end or KEY/DATA and a bool indicating the success or not
      virtual std::pair< typename HashMap< size_t, t_DataType>::iterator, bool>
      Insert( const std::pair< size_t, t_DataType> &KEY_DATA_PAIR)
      {
        return HashMap< size_t, t_DataType>::Insert( KEY_DATA_PAIR);
      }

    };

    // instantiate s_Instance
    template< size_t t_NumberOfObjects, typename t_ObjectType, typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> ObjectNDHashMap< t_NumberOfObjects, t_ObjectType, t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new ObjectNDHashMap< t_NumberOfObjects, t_ObjectType, t_DataType>())
    );

    //! @brief operator << this is a dummy implementation for storage::Vector to compile. This operator cannot be used.
    //! @param OSTREAM the stream which is passed to the operator
    //! @param OBJECT_ND_HASHMAP is the Hashmap which would be written
    //! @return returns the stream, but gives BCL_Exit before this can happen
    template< size_t t_NumberOfObjects, typename t_ObjectType, typename t_DataType>
    inline std::ostream &operator <<
    (
      std::ostream &OSTREAM,
      const ObjectNDHashMap< t_NumberOfObjects, t_ObjectType, t_DataType> &OBJECT_ND_HASHMAP
    )
    {
      // ends program before STREAM is returned
      BCL_Exit( "this function cannot be used", -1);
      return OSTREAM;
    }

    //! @brief operator >> this is a dummy implementation for storage::Vector to compile. This operator cannot be used.
    //! @param ISTREAM the stream which is passed to the operator
    //! @param OBJECTNDHASHMAP is the Hashmap which would be input
    //! @return returns the stream, but gives BCL_Exit before this can happen
    template< size_t t_NumberOfObjects, typename t_ObjectType, typename t_DataType>
    inline std::istream &operator >>
    (
      std::istream &ISTREAM,
      const ObjectNDHashMap< t_NumberOfObjects, t_ObjectType, t_DataType> &OBJECTNDHASHMAP
    )
    {
      // ends program before STREAM is returned
      BCL_Exit( "this function cannot be used", -1);
      return ISTREAM;
    }

  } // namespace storage
} // namespace bcl

#endif // BCL_STORAGE_OBJECT_ND_HASH_MAP_H_
