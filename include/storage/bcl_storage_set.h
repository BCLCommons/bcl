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

#ifndef BCL_STORAGE_SET_H_
#define BCL_STORAGE_SET_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <set>

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Set
    //! @brief This is a class for the Set container.
    //! @details Elements in this container are made up of data stored as a key.
    //!
    //! Attributes: (copied from http://www.sgi.com/tech/stl/Set.html)
    //! a.) stores key (as opposed to key value pairs as in Map and MultiMap)
    //! b.) key is also its value type
    //! c.) no two elements are the same
    //! d.) iterators are not invalidated unless the iterator was pointing to an element that was erased
    //! e.) always sorted in ascending order
    //!
    //! @see @link example_storage_set.cpp @endlink
    //! @author alexanns
    //! @date 01/21/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! @tparam t_ContainerType indicates the type of standard container that will be used
    //! @tparam t_KeyType indicates the type of keys that will be used
    //! @tparam t_KeyCompare is the comparison that will be used to sort the elements; default is less than
    template< typename t_KeyType, typename t_KeyCompare>
    class Set :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      std::set< t_KeyType, t_KeyCompare> m_Data;

    public:

      //! typedef for iterator
      typedef typename std::set< t_KeyType, t_KeyCompare>::iterator               iterator;

      //! typedef for const_iterator
      typedef typename std::set< t_KeyType, t_KeyCompare>::const_iterator         const_iterator;

      //! typedef for reverse_iterator
      typedef typename std::set< t_KeyType, t_KeyCompare>::reverse_iterator       reverse_iterator;

      //! typedef for const_reverse_iterator
      typedef typename std::set< t_KeyType, t_KeyCompare>::const_reverse_iterator const_reverse_iterator;

      //! typedef for const_reference
      typedef typename std::set< t_KeyType, t_KeyCompare>::const_reference        const_reference;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Set() :
        m_Data()
      {
      }

      //! @brief construct from single t_KeyType object
      //! @param OBJECT t_KeyType object which will be inserted in to the newly constructed Set
      Set( const t_KeyType &OBJECT) :
        m_Data()
      {
        m_Data.insert( OBJECT);
      }

      //! @brief construct from two t_KeyType objects
      //! @param OBJECT_A the first t_KeyType object which will be inserted into the newly constructed Set
      //! @param OBJECT_B the second t_KeyType object which will be inserted into the newly constructed Set
      Set( const t_KeyType &OBJECT_A, const t_KeyType &OBJECT_B) :
        m_Data()
      {
        m_Data.insert( OBJECT_A);
        m_Data.insert( OBJECT_B);
      }

      //! @brief construct from three t_KeyType objects
      //! @param OBJECT_A the first t_KeyType object which will be inserted into the newly constructed Set
      //! @param OBJECT_B the second t_KeyType object which will be inserted into the newly constructed Set
      //! @param OBJECT_C the third t_KeyType object which will be inserted into the newly constructed Set
      Set( const t_KeyType &OBJECT_A, const t_KeyType &OBJECT_B, const t_KeyType &OBJECT_C) :
        m_Data()
      {
        m_Data.insert( OBJECT_A);
        m_Data.insert( OBJECT_B);
        m_Data.insert( OBJECT_C);
      }

      //! @brief construct from four t_KeyType objects
      //! @param OBJECT_A the first t_KeyType object which will be inserted into the newly constructed Set
      //! @param OBJECT_B the second t_KeyType object which will be inserted into the newly constructed Set
      //! @param OBJECT_C the third t_KeyType object which will be inserted into the newly constructed Set
      //! @param OBJECT_D the fourth t_KeyType object which will be inserted into the newly constructed Set
      Set
      (
        const t_KeyType &OBJECT_A, const t_KeyType &OBJECT_B, const t_KeyType &OBJECT_C, const t_KeyType &OBJECT_D
      ) :
        m_Data()
      {
        m_Data.insert( OBJECT_A);
        m_Data.insert( OBJECT_B);
        m_Data.insert( OBJECT_C);
        m_Data.insert( OBJECT_D);
      }

      //! @brief construct a List from iterator [FIRST, LAST) range
      //! @tparam t_Iterator which is the type of iterator that indicates the range
      //! @param FIRST t_Iterator to the first element to be copied
      //! @param LAST t_Iterator to the first element after the last copied element
      template< typename t_Iterator>
      Set( t_Iterator FIRST, t_Iterator LAST) :
        m_Data()
      {
        m_Data.insert( FIRST, LAST);
      }

      //! @brief copy constructor
      Set
      (
        const Set< t_KeyType, t_KeyCompare> &SET
      ) :
        m_Data( SET.m_Data)
      {
      }

      //! @brief move constructor
      Set( Set< t_KeyType, t_KeyCompare> && SET) :
        m_Data( std::move( SET.m_Data))
      {
      }

      //! @brief clone copy constructor
      //! @return pointer to a copy of the actual object
      Set< t_KeyType, t_KeyCompare> *Clone() const
      {
        return new Set< t_KeyType, t_KeyCompare>( *this);
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
      //! If the key already exists in the map then the iterator indicates the first occurance.
      //! @param KEY for which the lowerbound is desired
      //! @return gives an iterator to the lowerbound denoted by KEY or off-the-end iterator if not possible
      iterator LowerBound( const t_KeyType &KEY)
      {
        return m_Data.lower_bound( KEY);
      }

      //! @brief LowerBound gives a const_iterator to the place where a key would fit into the map.
      //! If the key already exists in the map then the iterator indicates the first occurance.
      //! @param KEY for which the lowerbound is desired
      //! @return gives an iterator to the lowerbound denoted by KEY or off-the-end iterator if not possible
      const_iterator LowerBound( const t_KeyType &KEY) const
      {
        return m_Data.lower_bound( KEY);
      }

      //! @brief UpperBound gives a const_iterator to the place where a key would fit into the map
      //! If the key already exists in the map then the iterator denotes one past the last occurance.
      //! @param KEY for which the upperbound is desired
      //! @return gives an iterator to the upperbound denoted by KEY or an off-the-end iterator if not possible
      iterator UpperBound( const t_KeyType &KEY)
      {
        return m_Data.upper_bound( KEY);
      }

      //! @brief UpperBound gives a const_iterator to the place where a key would fit into the map
      //! If the key already exists in the map then the iterator denotes one past the last occurance.
      //! @param KEY for which the upperbound is desired
      //! @return gives an iterator to the upperbound denoted by KEY or an off-the-end iterator if not possible
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

      //! @brief uses an iterator range of t_KeyType to remove elements from the set
      //!        The t_KeyType object behind each iterator is used to remove that element
      //! @param ITR the first t_KeyType object iterator
      //! @param ITR_END one past the last t_KeyType object iterator that will be removed
      template< typename t_Iterator>
      void EraseKeys( t_Iterator ITR, t_Iterator ITR_END)
      {
        for( ; ITR != ITR_END; ++ITR)
        {
          Erase( *ITR);
        }
      }

      //! @brief Clear erases all elements of the associative container
      void Reset()
      {
        m_Data.clear();
      }

      //! @brief InsertElement insert ELEMENT into the container
      //! @param ELEMENT an object of t_DataType that is inserted
      void InsertElement( const t_KeyType &ELEMENT)
      {
        m_Data.insert( ELEMENT);
      }

      //! @brief InsertElements insert all elements of argument SET
      //! @param SET container containing objects of t_DataType that are inserted
      void InsertElements( const Set< t_KeyType, t_KeyCompare> &SET)
      {
        m_Data.insert( SET.Begin(), SET.End());
      }

      //! @brief InsertElements inserts a number of elements into the container [INCLUSIVE, NOT_INCLUSIVE)
      //! @param ITR_A iterator denoting the first element to be inserted
      //! @param ITR_B iterator denoting the end of the elements to be inserted (not included in the insertions)
      template< typename t_Iterator>
      void InsertElements( t_Iterator ITR_A, t_Iterator ITR_B)
      {
        m_Data.insert( ITR_A, ITR_B);
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

      //! @brief Contains returns true if KEY is is in the set, and false otherwise
      //! @param KEY the key to test for membership in the set
      //! @return true if KEY is is in the set, and false otherwise
      bool Contains( const t_KeyType &KEY) const
      {
        return m_Data.find( KEY) != m_Data.end();
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

      //! @brief returns a const reference to the used stl container
      //! @return const reference to the internal stl container
      std::set< t_KeyType, t_KeyCompare> const &InternalData() const
      {
        return m_Data;
      }

      //! @brief returns a changeable reference to the used stl container
      //! @return changeable reference to the internal stl container
      std::set< t_KeyType, t_KeyCompare> &InternalData()
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

      //! @brief Swap exchanges the contents of two sets
      //! @param SWAP_SET is the set with which the contents of m_Data will be swapped
      void Swap( Set< t_KeyType, t_KeyCompare> &SWAP_SET)
      {
        m_Data.swap( SWAP_SET.m_Data);
      }

      //! @brief test whether this set is a subset of another set
      //! @param SUPERSET the putative superset
      //! @return true if every member of this is in SUPERSET
      bool IsSubsetOf( const Set< t_KeyType, t_KeyCompare> &SUPERSET) const
      {
        const size_t this_size( GetSize()), superset_size( SUPERSET.GetSize());
        if( this_size > superset_size)
        {
          return false;
        }
        else if( this_size == superset_size)
        {
          return InternalData() == SUPERSET.InternalData();
        }
        // *this is smaller than SUPERSET, so need to check by elements. Use std::includes since Set sorts the ranges
        return std::includes( SUPERSET.Begin(), SUPERSET.End(), Begin(), End(), t_KeyCompare());
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief Insert takes a standard pair and inserts it into the UniqueSortedAssociativeContainer
      //! with the passed iterator as a suggestion on where to start looking for the insertion
      //! @tparam t_Iterator is the type of iterator that suggests where the insert should start looking
      //! @param ITR is the iterator denoting the position to start looking for the correct insertion point
      //! @param KEY is the pair that is desired to be inserted into the UniqueSortedAssociativeContainer
      //! @return an iterator pointing to the place of insertion
      iterator
      Insert( iterator ITR, const t_KeyType &KEY)
      {
        return m_Data.insert( ITR, KEY);
      }

      //! @brief Insert takes a standard pair and inserts it into the UniqueSortedAssociativeContainer
      //! with the passed iterator as a suggestion on where to start looking for the insertion
      //! @tparam t_Iterator is the type of iterator that suggests where the insert should start looking
      //! @param KEY is the pair that is desired to be inserted into the UniqueSortedAssociativeContainer
      //! @return an iterator pointing to the place of insertion
      std::pair< iterator, bool>
      Insert( const t_KeyType &KEY)
      {
        return m_Data.insert( KEY);
      }

      //! @brief assignment operator
      Set &operator =( const Set &SET)
      {
        m_Data = SET.m_Data;
        return *this;
      }

      //! @brief move assignment operator
      Set &operator =( Set && SET)
      {
        m_Data = std::move( SET.m_Data);
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

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
          t_KeyType current;

          BCL_Assert( io::Serialize::Read( current, ISTREAM), "Error reading element!");

          // check that insert was successful
          BCL_Assert( Insert( current).second, "error inserting");
        }

        // return
        return ISTREAM;
      }

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
          io::Serialize::Write( *itr, OSTREAM, INDENT) << '\n';
        }

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////
    public:

      //! @brief creates the Set from one element (for consistency)
      static Set< t_KeyType, t_KeyCompare> Create( const t_KeyType &DATA)
      {
        // return set
        return Set< t_KeyType, t_KeyCompare>( DATA);
      }

      //! @brief creates the Set from two elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B
      )
      {
        // create set from two elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);

        // return set
        return set;
      }

      //! @brief creates the Set from three elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C
      )
      {
        // create set from three elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);

        // return set
        return set;
      }

      //! @brief creates the Set from four elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D
      )
      {
        // create set from four elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);

        // return set
        return set;
      }

      //! @brief creates the Set from five elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E
      )
      {
        // create set from five elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);

        // return set
        return set;
      }

      //! @brief creates the Set from six elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);

        // return set
        return set;
      }

      //! @brief creates the Set from seven elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G
      )
      {
        // create set from 7 elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);

        // return set
        return set;
      }

      //! @brief creates the Set from eight elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H
      )
      {
        // create set from 8 elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);

        // return set
        return set;
      }

      //! @brief creates the Set from nine elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I
      )
      {
        // create set from 9 elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);

        // return set
        return set;
      }

      //! @brief creates the Set from 10 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J
      )
      {
        // create set from 10 elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);

        // return set
        return set;
      }

      //! @brief creates the Set from 11 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K
      )
      {
        // create set from 11 elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);

        // return set
        return set;
      }

      //! @brief creates the Set from 12 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L
      )
      {
        // create set from 12 elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);

        // return set
        return set;
      }

      //! @brief creates the Set from 13 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M
      )
      {
        // create set from 13 elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);

        // return set
        return set;
      }

      //! @brief creates the Set from 14 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);

        // return set
        return set;
      }

      //! @brief creates the Set from 15 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);

        // return set
        return set;
      }

      //! @brief creates the Set from 16 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O,
        const t_KeyType &DATA_P
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);
        set.Insert( DATA_P);

        // return set
        return set;
      }

      //! @brief creates the Set from 17 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O,
        const t_KeyType &DATA_P,
        const t_KeyType &DATA_Q
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);
        set.Insert( DATA_P);
        set.Insert( DATA_Q);

        // return set
        return set;
      }

      //! @brief creates the Set from 18 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O,
        const t_KeyType &DATA_P,
        const t_KeyType &DATA_Q,
        const t_KeyType &DATA_R
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);
        set.Insert( DATA_P);
        set.Insert( DATA_Q);
        set.Insert( DATA_R);

        // return set
        return set;
      }

      //! @brief creates the Set from 19 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O,
        const t_KeyType &DATA_P,
        const t_KeyType &DATA_Q,
        const t_KeyType &DATA_R,
        const t_KeyType &DATA_S
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);
        set.Insert( DATA_P);
        set.Insert( DATA_Q);
        set.Insert( DATA_R);
        set.Insert( DATA_S);

        // return set
        return set;
      }

      //! @brief creates the Set from 20 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O,
        const t_KeyType &DATA_P,
        const t_KeyType &DATA_Q,
        const t_KeyType &DATA_R,
        const t_KeyType &DATA_S,
        const t_KeyType &DATA_T
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);
        set.Insert( DATA_P);
        set.Insert( DATA_Q);
        set.Insert( DATA_R);
        set.Insert( DATA_S);
        set.Insert( DATA_T);

        // return set
        return set;
      }

      //! @brief creates the Set from 21 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O,
        const t_KeyType &DATA_P,
        const t_KeyType &DATA_Q,
        const t_KeyType &DATA_R,
        const t_KeyType &DATA_S,
        const t_KeyType &DATA_T,
        const t_KeyType &DATA_U
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);
        set.Insert( DATA_P);
        set.Insert( DATA_Q);
        set.Insert( DATA_R);
        set.Insert( DATA_S);
        set.Insert( DATA_T);
        set.Insert( DATA_U);

        // return set
        return set;
      }

      //! @brief creates the Set from 22 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O,
        const t_KeyType &DATA_P,
        const t_KeyType &DATA_Q,
        const t_KeyType &DATA_R,
        const t_KeyType &DATA_S,
        const t_KeyType &DATA_T,
        const t_KeyType &DATA_U,
        const t_KeyType &DATA_V
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);
        set.Insert( DATA_P);
        set.Insert( DATA_Q);
        set.Insert( DATA_R);
        set.Insert( DATA_S);
        set.Insert( DATA_T);
        set.Insert( DATA_U);
        set.Insert( DATA_V);

        // return set
        return set;
      }

      //! @brief creates the Set from 24 elements
      static Set< t_KeyType, t_KeyCompare> Create
      (
        const t_KeyType &DATA_A,
        const t_KeyType &DATA_B,
        const t_KeyType &DATA_C,
        const t_KeyType &DATA_D,
        const t_KeyType &DATA_E,
        const t_KeyType &DATA_F,
        const t_KeyType &DATA_G,
        const t_KeyType &DATA_H,
        const t_KeyType &DATA_I,
        const t_KeyType &DATA_J,
        const t_KeyType &DATA_K,
        const t_KeyType &DATA_L,
        const t_KeyType &DATA_M,
        const t_KeyType &DATA_N,
        const t_KeyType &DATA_O,
        const t_KeyType &DATA_P,
        const t_KeyType &DATA_Q,
        const t_KeyType &DATA_R,
        const t_KeyType &DATA_S,
        const t_KeyType &DATA_T,
        const t_KeyType &DATA_U,
        const t_KeyType &DATA_V,
        const t_KeyType &DATA_W,
        const t_KeyType &DATA_X
      )
      {
        // create set from six elements
        Set< t_KeyType, t_KeyCompare> set;

        set.Insert( DATA_A);
        set.Insert( DATA_B);
        set.Insert( DATA_C);
        set.Insert( DATA_D);
        set.Insert( DATA_E);
        set.Insert( DATA_F);
        set.Insert( DATA_G);
        set.Insert( DATA_H);
        set.Insert( DATA_I);
        set.Insert( DATA_J);
        set.Insert( DATA_K);
        set.Insert( DATA_L);
        set.Insert( DATA_M);
        set.Insert( DATA_N);
        set.Insert( DATA_O);
        set.Insert( DATA_P);
        set.Insert( DATA_Q);
        set.Insert( DATA_R);
        set.Insert( DATA_S);
        set.Insert( DATA_T);
        set.Insert( DATA_U);
        set.Insert( DATA_V);
        set.Insert( DATA_W);
        set.Insert( DATA_X);

        // return set
        return set;
      }

    }; // class Set

    // instantiate s_Instance
    template< typename t_KeyType, typename t_KeyCompare>
    const util::SiPtr< const util::ObjectInterface> Set< t_KeyType, t_KeyCompare>::s_Instance
    (
      GetObjectInstances().AddInstance( new Set< t_KeyType, t_KeyCompare>())
    );

    //! @brief operator == checks if two containers are the same
    //! @param SET_A first container
    //! @param SET_B second container
    //! @return returns boolean value
    template< typename t_DataType>
    inline bool operator ==
    (
      const Set< t_DataType> &SET_A,
      const Set< t_DataType> &SET_B
    )
    {
      return SET_A.InternalData() == SET_B.InternalData();
    }

    //! @brief operator == checks if two containers are the same
    //! @param SET_A first container
    //! @param SET_B second container
    //! @return returns boolean value
    template< typename t_DataType>
    inline bool operator <
    (
      const Set< t_DataType> &SET_A,
      const Set< t_DataType> &SET_B
    )
    {
      return SET_A.InternalData() < SET_B.InternalData();
    }

    //! @brief operator + checks if two containers are the same
    //! @param SET_A first container
    //! @param SET_B second container
    //! @return returns boolean value
    template< typename t_DataType>
    Set< t_DataType> operator +
    (
      const Set< t_DataType> &A,
      const Set< t_DataType> &B
    )
    {
      Set< t_DataType> new_set;
      std::set_union( A.Begin(), A.End(), B.Begin(), B.End(), std::inserter( new_set.InternalData(), new_set.Begin()));
      return new_set;
    }

    //! @brief operator + checks if two containers are the same
    //! @param SET_A first container
    //! @param SET_B second container
    //! @return returns boolean value
    template< typename t_DataType>
    Set< t_DataType> operator -
    (
      const Set< t_DataType> &A,
      const Set< t_DataType> &B
    )
    {
      Set< t_DataType> new_set;
      std::set_difference( A.Begin(), A.End(), B.Begin(), B.End(), std::inserter( new_set.InternalData(), new_set.Begin()));
      return new_set;
    }

    //! @brief operator + checks if two containers are the same
    //! @param SET_A first container
    //! @param SET_B second container
    //! @return returns boolean value
    template< typename t_DataType>
    Set< t_DataType> operator ^
    (
      const Set< t_DataType> &A,
      const Set< t_DataType> &B
    )
    {
      Set< t_DataType> new_set;
      std::set_intersection( A.Begin(), A.End(), B.Begin(), B.End(), std::inserter( new_set.InternalData(), new_set.Begin()));
      return new_set;
    }

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_SET_H_
