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

#ifndef BCL_STORAGE_MAP_H_
#define BCL_STORAGE_MAP_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_storage_pair.h"
#include "bcl_storage_set.h"
#include "bcl_storage_vector.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically
#include <map>

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Map
    //! @brief This is a class for the Map container.
    //! @details Elements in this container are made up of data stored in association with a key.
    //!
    //! Attributes: (copied from http://www.sgi.com/tech/stl/Map.html)
    //! a.) stores pairs of key and value
    //! b.) values are access by key
    //! c.) no two elements have the same key
    //! d.) iterators are not invalidated unless the iterator was pointing to an element that was erased
    //! e.) always sorted in ascending order by key
    //!
    //! @see @link example_storage_map.cpp @endlink
    //! @author alexanns, jsmith
    //! @date 01/21/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! @tparam t_ContainerType indicates the type of standard container that will be used
    //! @tparam t_KeyType indicates the type of keys that will be used
    //! @tparam t_DataType indicates the type of data that the container will hold
    //! @tparam t_KeyCompare is the comparison that will be used to sort the elements; default is less than
    template< typename t_KeyType, typename t_DataType, typename t_KeyCompare>
    class Map :
      public util::ObjectInterface
    {

    protected:

    //////////
    // data //
    //////////

      //! standard container member of AssociativeContainer
      std::map< t_KeyType, t_DataType, t_KeyCompare> m_Data;

    public:

      //! typedef for iterator
      typedef typename std::map< t_KeyType, t_DataType, t_KeyCompare>::iterator               iterator;

      //! typedef for const_iterator
      typedef typename std::map< t_KeyType, t_DataType, t_KeyCompare>::const_iterator         const_iterator;

      //! typedef for reverse_iterator
      typedef typename std::map< t_KeyType, t_DataType, t_KeyCompare>::reverse_iterator       reverse_iterator;

      //! typedef for const_reverse_iterator
      typedef typename std::map< t_KeyType, t_DataType, t_KeyCompare>::const_reverse_iterator const_reverse_iterator;

      //! typedef for key type
      typedef typename std::map< t_KeyType, t_DataType, t_KeyCompare>::key_type               key_type;

      //! typedef for mapped type
      typedef typename std::map< t_KeyType, t_DataType, t_KeyCompare>::mapped_type            mapped_type;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Map() :
        m_Data()
      {
      }

      //! @brief constructor that constructs a map from a range of iterators [INCLUSIVE, NOT_INCLUSIVE)
      //! @param FIRST denotes the first element to be copied into the new map
      //! @param LAST denotes the end (won't be copied) of the range of elements to be copied into the new map
      template< typename t_Iterator>
      Map
      (
        const t_Iterator &FIRST,
        const t_Iterator &LAST
      ) :
        m_Data( FIRST, LAST)
      {
      }

      //! @brief copy constructor
      Map( const Map &MAP) :
        m_Data( MAP.m_Data)
      {
      }

      //! @brief move constructor
      Map( Map && MAP) :
        m_Data( std::move( MAP.m_Data))
      {
      }

      //! @brief clone virtual copy constructor
      //! @return pointer to a copy of the actual object
      Map< t_KeyType, t_DataType, t_KeyCompare> *Clone() const
      {
        return new Map< t_KeyType, t_DataType, t_KeyCompare>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object
      //! @return string with the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return iterator to reverse begin
      //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
      reverse_iterator ReverseBegin()
      {
        return m_Data.rbegin();
      }

      //! @brief return const_iterator to reverse begin
      //! @return const_reverse_iterator pointing to the beginning of the reversed container
      const_reverse_iterator ReverseBegin() const
      {
        return m_Data.rbegin();
      }

      //! @brief return iterator to reverse end
      //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
      reverse_iterator ReverseEnd()
      {
        return m_Data.rend();
      }

      //! @brief return const_iterator to reverse end
      //! @return const_reverse_iterator pointing to the end of the reversed container
      const_reverse_iterator ReverseEnd() const
      {
        return m_Data.rend();
      }

      //! @brief LowerBound gives an iterator to the place where a key would fit into the map.
      //! If the key already exists in the map then the iterator indicates the first occurrence.
      //! @param KEY for which the lower bound is desired
      //! @return gives an iterator to the lower bound denoted by KEY or off-the-end iterator if not possible
      iterator LowerBound( const t_KeyType &KEY)
      {
        return m_Data.lower_bound( KEY);
      }

      //! @brief LowerBound gives a const_iterator to the place where a key would fit into the map.
      //! If the key already exists in the map then the iterator indicates the first occurrence.
      //! @param KEY for which the lower bound is desired
      //! @return gives an iterator to the lower bound denoted by KEY or off-the-end iterator if not possible
      const_iterator LowerBound( const t_KeyType &KEY) const
      {
        return m_Data.lower_bound( KEY);
      }

      //! @brief UpperBound gives a const_iterator to the place where a key would fit into the map
      //! If the key already exists in the map then the iterator denotes one past the last occurrence.
      //! @param KEY for which the upper bound is desired
      //! @return gives an iterator to the upper bound denoted by KEY or an off-the-end iterator if not possible
      iterator UpperBound( const t_KeyType &KEY)
      {
        return m_Data.upper_bound( KEY);
      }

      //! @brief UpperBound gives a const_iterator to the place where a key would fit into the map
      //! If the key already exists in the map then the iterator denotes one past the last occurrence.
      //! @param KEY for which the upper bound is desired
      //! @return gives an iterator to the upper bound denoted by KEY or an off-the-end iterator if not possible
      const_iterator UpperBound( const t_KeyType &KEY) const
      {
        return m_Data.upper_bound( KEY);
      }

      //! @brief Erase removes all the elements which have the specified key
      //! @param KEY denotes the key which will be removed
      //! @return returns the number of keys which were removed
      size_t Erase( const t_KeyType &KEY)
      {
        return m_Data.erase( KEY);
      }

      //! @brief RemoveElement removes the element which is pointed to by the given iterator
      //! @param ITR iterator which denotes the element to be deleted
      void RemoveElement( iterator ITR)
      {
        m_Data.erase( ITR);
      }

      //! @brief Erase removes ell elements in range denoted by two iterators
      //! @param ITR_A iterator denoting the start of elements to be deleted
      //! @param ITR_B iterator denoting the end position of elements to be deleted (non-inclusive)
      void Erase( iterator ITR_A, iterator ITR_B)
      {
        m_Data.erase( ITR_A, ITR_B);
      }

      //! @brief Clear erases all elements of the associative container
      void Reset()
      {
        m_Data.clear();
      }

      //! @brief Find takes a key and returns an iterator to the element denoted by key, or an off the end iterator
      //! @param KEY the key for which an iterator is desired
      //! @return returns an iterator to the element
      iterator Find( const t_KeyType &KEY)
      {
        return m_Data.find( KEY);
      }

      //! @brief Find takes a key and returns a const_iterator to the element denoted by key, or off the end iterator
      //! @param KEY the key for which a const_iterator is desired
      //! @return returns an iterator to the element
      const_iterator Find( const t_KeyType &KEY) const
      {
        return m_Data.find( KEY);
      }

      //! @brief Count takes a key and gives the number of times it appears in the associative container
      //! @param KEY the key for which a count is desired
      //! @return returns the number of instances of KEY
      size_t Count( const t_KeyType &KEY) const
      {
        return m_Data.count( KEY);
      }

      //! @brief EqualRange gives a pair of iterators that denotes the range of the first and last instance of a key
      //! @param KEY for which the range is desired
      //! @return returns a pair of iterators denoting the first and last instance of a key
      std::pair< iterator, iterator>
      EqualRange( const t_KeyType &KEY)
      {
        return m_Data.equal_range( KEY);
      }

      //! @brief returns size of the container
      //! @return size, i.e. number of elements stored
      size_t GetSize() const
      {
        return m_Data.size();
      }

      //! @brief return the maximal size of the container
      //! @return maximum size, i.e. maximal number of elements to store
      size_t MaxSize() const
      {
        return m_Data.max_size();
      }

      //! @brief return all keys
      //! @return storage vector of all keys
      Set< t_KeyType, t_KeyCompare> GetKeys() const
      {
        // allocate vector for keys
        Set< t_KeyType, t_KeyCompare> keys;

        typename Set< t_KeyType, t_KeyCompare>::iterator last_insert( keys.Begin());

        // iterate over all elements
        for( const_iterator itr( m_Data.begin()), itr_end( m_Data.end()); itr != itr_end; ++itr)
        {
          last_insert = keys.Insert( last_insert, itr->first);
        }

        // end
        return keys;
      }

      //! @brief return all keys
      //! @return storage vector of all keys
      Vector< t_KeyType> GetKeysAsVector() const
      {
        // allocate vector for keys
        Vector< t_KeyType> keys;
        keys.AllocateMemory( m_Data.size());

        // iterate over all elements
        for( const_iterator itr( m_Data.begin()), itr_end( m_Data.end()); itr != itr_end; ++itr)
        {
          keys.PushBack( itr->first);
        }

        // end
        return keys;
      }

      //! @brief return all values
      //! @return storage vector of all values
      Vector< t_DataType> GetMappedValues() const
      {
        // allocate vector for keys
        Vector< t_DataType> values;
        values.AllocateMemory( m_Data.size());

        // iterate over all elements
        for( const_iterator itr( m_Data.begin()), itr_end( m_Data.end()); itr != itr_end; ++itr)
        {
          values.PushBack( itr->second);
        }

        // end
        return values;
      }

      //! @brief returns a const reference to the used stl container
      //! @return const reference to the internal stl container
      std::map< t_KeyType, t_DataType, t_KeyCompare> const &InternalData() const
      {
        return m_Data;
      }

      //! @brief returns a changeable reference to the used stl container
      //! @return changeable reference to the internal stl container
      std::map< t_KeyType, t_DataType, t_KeyCompare> &InternalData()
      {
        return m_Data;
      }

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin()
      {
        return m_Data.begin();
      }

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const
      {
        return m_Data.begin();
      }

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End()
      {
        return m_Data.end();
      }

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const
      {
        return m_Data.end();
      }

      //! @brief Empty tests whether the container is empty or not
      //! @return boolean value indicating if the container is empty or not
      bool IsEmpty() const
      {
        return m_Data.empty();
      }

      //! @brief Swap exchanges the contents of two associative containers
      //! @param SWAPPER is the container with which the contents of m_Data will be swapped
      void Swap( Map< t_KeyType, t_DataType, t_KeyCompare> &SWAPPER)
      {
        m_Data.swap( SWAPPER.m_Data);
      }

      //! @brief Insert takes a standard pair and inserts it into the map
      //! @param STD_PAIR is the pair that is desired to be inserted into the map
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      std::pair< iterator, bool>
      Insert( const std::pair< t_KeyType, t_DataType> &STD_PAIR)
      {
        return m_Data.insert( STD_PAIR);
      }

      //! @brief Insert takes a standard pair and inserts it into the map
      //! @param BCL_PAIR is the pair that is desired to be inserted into the map
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      std::pair< iterator, bool>
      Insert( const Pair< t_KeyType, t_DataType> &BCL_PAIR)
      {
        return m_Data.insert
        (
          std::pair< t_KeyType, t_DataType>( BCL_PAIR.First(), BCL_PAIR.Second())
        );
      }

      //! @brief InsertElement takes a standard pair and inserts it into the map
      //! @param STD_PAIR is the pair that is desired to be inserted into the map
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      void
      InsertElement( const std::pair< t_KeyType, t_DataType> &STD_PAIR)
      {
        m_Data.insert( STD_PAIR);
      }

      //! @brief Insert takes a BCL pair and inserts it into the map
      //! @param BCL_PAIR is the pair that is desired to be inserted into the map
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      void
      InsertElement( const Pair< t_KeyType, t_DataType> &BCL_PAIR)
      {
        m_Data.insert
        (
          std::pair< t_KeyType, t_DataType>( BCL_PAIR.First(), BCL_PAIR.Second())
        );
      }

      //! @brief InsertElements inserts a number of elements into the container [INCLUSIVE, NOT_INCLUSIVE)
      //! @param ITR_A iterator denoting the first element to be inserted
      //! @param ITR_B iterator denoting the end of the elements to be inserted (not included in the insertions)
      template< typename t_Iterator>
      void InsertElements( const t_Iterator ITR_A, const t_Iterator ITR_B)
      {
        m_Data.insert( ITR_A, ITR_B);
      }

      //! @brief InsertElements insert all elements of argument MAP
      //! @param MAP container containing objects of t_ValueType that are inserted
      void
      InsertElements( const Map< t_KeyType, t_DataType, t_KeyCompare> &MAP)
      {
        m_Data.insert
        (
          MAP.Begin(), MAP.End()
        );
      }

      //! @brief gets the value for a given KEY
      //! @param KEY key to be searched for
      //! @return the value for the given KEY
      const t_DataType &GetValue( const t_KeyType &KEY) const
      {
        // search for the key
        const const_iterator itr( m_Data.find( KEY));

        BCL_Assert
        (
          itr != m_Data.end(),
          "Error, could not locate key: " + util::Format()( KEY)
          + " in constant map: " + util::Format()( *this)
        );

        return itr->second;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator [] gives a reference to the data associated with a specified key
      //! @param KEY is the key for which the data is desired
      //! @return a reference to the data associated with KEY
      t_DataType &operator []( const t_KeyType &KEY)
      {
        return m_Data[ KEY];
      }

      //! @brief assignment operator
      Map &operator =( const Map &MAP)
      {
        m_Data = MAP.m_Data;
        return *this;
      }

      //! @brief move assignment operator
      Map &operator =( Map && MAP)
      {
        m_Data = std::move( MAP.m_Data);
        return *this;
      }

      //! @brief Has checks whether the map has the given KEY
      //! @param KEY is the key to check for
      //! @return true if KEY is found in the map
      bool Has( const t_KeyType &KEY) const
      {
        return m_Data.find( KEY) != m_Data.end();
      }

      //! @brief InsertElement takes a standard pair and inserts it into the UniqueSortedAssociativeContainer
      //! with the passed iterator as a suggestion on where to start looking for the insertion
      //! @tparam t_Iterator is the type of iterator that suggests where the insert should start looking
      //! @param ITR is the iterator denoting the position to start looking for the correct insertion point
      //! @param STD_PAIR is the pair that is desired to be inserted into the UniqueSortedAssociativeContainer
      //! @return an iterator pointing to the place of insertion
      iterator
      InsertElement( iterator ITR, const std::pair< t_KeyType, t_DataType> &STD_PAIR)
      {
        return m_Data.insert( ITR, STD_PAIR);
      }

      //! @brief InsertElement takes a BCL pair and inserts it into the UniqueAssociativeContainer
      //! with the passed iterator as a suggestion on where to start looking for the insertion
      //! @tparam t_Iterator is the type of iterator that suggests where the insert should start looking
      //! @param ITR is the iterator denoting the position to start looking for the correct insertion point
      //! @param BCL_PAIR is the pair that is desired to be inserted into the UniqueSortedAssociativeContainer
      //! @return an iterator pointing to the place of insertion
      iterator
      InsertElement( iterator ITR, const Pair< t_KeyType, t_DataType> &BCL_PAIR)
      {
        return m_Data.insert
        (
          ITR,
          std::pair< t_KeyType, t_DataType>( BCL_PAIR.First(), BCL_PAIR.Second())
        );
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief creates the storage map from one element
      //! @param DATA_A first element
      static Map< t_KeyType, t_DataType, t_KeyCompare> Create
      (
        const std::pair< t_KeyType, t_DataType> &DATA_A
      )
      {
        // create map from two elements
        Map< t_KeyType, t_DataType, t_KeyCompare> map;
        map.InsertElement( DATA_A);

        // return map
        return map;
      }

      //! @brief creates the storage map from two elements
      //! @param DATA_A first element
      //! @param DATA_B second element
      static Map< t_KeyType, t_DataType, t_KeyCompare> Create
      (
        const std::pair< t_KeyType, t_DataType> &DATA_A,
        const std::pair< t_KeyType, t_DataType> &DATA_B
      )
      {
        // create map from two elements
        Map< t_KeyType, t_DataType, t_KeyCompare> map;
        map.InsertElement( DATA_A);
        map.InsertElement( DATA_B);

        // return map
        return map;
      }

      //! @brief creates the storage map from three elements
      //! @param DATA_A first element
      //! @param DATA_B second element
      //! @param DATA_C third element
      static Map< t_KeyType, t_DataType, t_KeyCompare> Create
      (
        const std::pair< t_KeyType, t_DataType> &DATA_A,
        const std::pair< t_KeyType, t_DataType> &DATA_B,
        const std::pair< t_KeyType, t_DataType> &DATA_C
      )
      {
        // create map from three elements
        Map< t_KeyType, t_DataType, t_KeyCompare> map;
        map.InsertElement( DATA_A);
        map.InsertElement( DATA_B);
        map.InsertElement( DATA_C);

        // return map
        return map;
      }

      //! @brief creates the storage map from four elements
      //! @param DATA_A first element
      //! @param DATA_B second element
      //! @param DATA_C third element
      //! @param DATA_D fourth element
      static Map< t_KeyType, t_DataType, t_KeyCompare> Create
      (
        const std::pair< t_KeyType, t_DataType> &DATA_A,
        const std::pair< t_KeyType, t_DataType> &DATA_B,
        const std::pair< t_KeyType, t_DataType> &DATA_C,
        const std::pair< t_KeyType, t_DataType> &DATA_D
      )
      {
        // create map from four elements
        Map< t_KeyType, t_DataType, t_KeyCompare> map;
        map.InsertElement( DATA_A);
        map.InsertElement( DATA_B);
        map.InsertElement( DATA_C);
        map.InsertElement( DATA_D);

        // return map
        return map;
      }

      //! @brief creates the storage map from five elements
      //! @param DATA_A first element
      //! @param DATA_B second element
      //! @param DATA_C third element
      //! @param DATA_D fourth element
      //! @param DATA_E fifth element
      static Map< t_KeyType, t_DataType, t_KeyCompare> Create
      (
        const std::pair< t_KeyType, t_DataType> &DATA_A,
        const std::pair< t_KeyType, t_DataType> &DATA_B,
        const std::pair< t_KeyType, t_DataType> &DATA_C,
        const std::pair< t_KeyType, t_DataType> &DATA_D,
        const std::pair< t_KeyType, t_DataType> &DATA_E
      )
      {
        // create map from four elements
        Map< t_KeyType, t_DataType, t_KeyCompare> map;
        map.InsertElement( DATA_A);
        map.InsertElement( DATA_B);
        map.InsertElement( DATA_C);
        map.InsertElement( DATA_D);
        map.InsertElement( DATA_E);

        // return map
        return map;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write container to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write size of container
        io::Serialize::Write( GetSize(), OSTREAM, INDENT) << '\n';

        // write each element of container to new line
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          io::Serialize::Write( std::pair< t_KeyType, t_DataType>( *itr), OSTREAM, INDENT) << '\n';
        }

        // end
        return OSTREAM;
      }

      //! @brief read container from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // empty the container
        Reset();

        // create int "size" for holding the size of the container
        size_t size;

        // read in size of container
        io::Serialize::Read( size, ISTREAM);

        // read in each element from stream
        for( size_t count( 0); count < size; ++count)
        {
          std::pair< t_KeyType, t_DataType> current;

          BCL_Assert( io::Serialize::Read( current, ISTREAM), "Error reading element!");

          // check that insert was successful
          BCL_Assert( Insert( current).second, "error inserting");
        }

        // return
        return ISTREAM;
      }

    }; // class Map

    // instantiate s_Instance
    template< typename t_KeyType, typename t_DataType, typename t_KeyCompare>
    const util::SiPtr< const util::ObjectInterface> Map< t_KeyType, t_DataType, t_KeyCompare>::s_Instance
    (
      GetObjectInstances().AddInstance( new Map< t_KeyType, t_DataType, t_KeyCompare>())
    );

    //! @brief operator == checks if two Maps are the same
    //! @param MAP_A first container
    //! @param MAP_B second container
    //! @return true, if Maps are identical
    template< typename t_KeyType, typename t_DataType, typename t_CompareType>
    inline bool operator ==
    (
      const Map< t_KeyType, t_DataType, t_CompareType> &MAP_A,
      const Map< t_KeyType, t_DataType, t_CompareType> &MAP_B
    )
    {
      return MAP_A.InternalData() == MAP_B.InternalData();
    }

    //! @brief operator == checks if two Maps are the same
    //! @param MAP_A first container
    //! @param MAP_B second container
    //! @return true, if Maps are identical
    template< typename t_KeyType, typename t_DataType, typename t_CompareType>
    inline bool operator !=
    (
      const Map< t_KeyType, t_DataType, t_CompareType> &MAP_A,
      const Map< t_KeyType, t_DataType, t_CompareType> &MAP_B
    )
    {
      return MAP_A.InternalData() != MAP_B.InternalData();
    }

    //! @brief operator < allows for ordering of maps
    //! @param MAP_A first container
    //! @param MAP_B second container
    //! @return true, if Maps are identical
    template< typename t_KeyType, typename t_DataType, typename t_CompareType>
    inline bool operator <
    (
      const Map< t_KeyType, t_DataType, t_CompareType> &MAP_A,
      const Map< t_KeyType, t_DataType, t_CompareType> &MAP_B
    )
    {
      typename Map< t_KeyType, t_DataType, t_CompareType>::const_iterator
        itr_a( MAP_A.Begin()), itr_a_end( MAP_A.End()), itr_b( MAP_B.Begin()), itr_b_end( MAP_B.End());
      for( ; itr_a != itr_a_end && itr_b != itr_b_end; ++itr_a, ++itr_b)
      {
        if( itr_a->first != itr_b->first)
        {
          return itr_a->first < itr_b->first;
        }
        else if( itr_a->second != itr_b->second)
        {
          return itr_a->second < itr_b->second;
        }
      }
      return itr_a == itr_a_end && itr_b != itr_b_end;
    }

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_MAP_H_
