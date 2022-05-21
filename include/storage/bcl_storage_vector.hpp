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

#ifndef BCL_STORAGE_VECTOR_HPP_
#define BCL_STORAGE_VECTOR_HPP_

// include the header of this class
#include "bcl_storage_vector.h"

// includes from bcl - sorted alphabetically
#include "bcl_storage.h"
#include "io/bcl_io_serialize.h"
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically
#include <algorithm>

namespace bcl
{
  namespace storage
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Vector< t_DataType>::Vector() :
      m_Data()
    {
    }

    //! @brief construct a Vector from optional size and single element
    //! @param SIZE optional size, default 0
    //! @param VALUE optional single element which is assigned to all stored elements, default is t_DataType()
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const size_t SIZE, const t_DataType &VALUE) :
      m_Data( SIZE, VALUE)
    {
    }

    //! @brief construct a Vector from number of size and pointer to data
    //! @param SIZE number of elements behind the pointer *DATA
    //! @param DATA pointer to the elements
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const size_t SIZE, const t_DataType *DATA) :
      m_Data( DATA, DATA + SIZE)
    {
    }

    //! @brief construct a Vector from another Vector
    //! @param VECTOR the Deque to copy
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const Vector< t_DataType> &VECTOR) :
      m_Data( VECTOR.m_Data)
    {
    }

    //! @brief construct a Vector from another Vector
    //! @param VECTOR the Deque to copy
    template< typename t_DataType>
    Vector< t_DataType>::Vector( Vector< t_DataType> && VECTOR) :
      m_Data( std::move( VECTOR.m_Data))
    {
    }

    //! @brief construct a Vector from another Vector
    //! @param VECTOR the Deque to copy
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const std::vector< t_DataType> &VECTOR) :
      m_Data( VECTOR)
    {
    }

    //! @brief constuct a Vector from a subsequence [ POS, POS + LENGTH) of another Vector
    //! @param VECTOR the Vector to copy
    //! @param POS the first element to copy
    //! @param LENGTH the (maximum) number of elements to use
    //! in no event will elements beyond the end of VECTOR be copied
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const Vector< t_DataType> &VECTOR, const size_t POS, const size_t LENGTH)
    {
      if( POS < VECTOR.GetSize()) // if POS was beyond the length of the vector, then this vector will be empty
      {
        m_Data =
          std::vector< t_DataType>
          (
            VECTOR.Begin() + POS,
            VECTOR.Begin() + POS + std::min( VECTOR.GetSize() - POS, LENGTH)
          );
      }
    }

    //! @brief copy constructor
    //! @return pointer to a copy of the actual object
    template< typename t_DataType>
    Vector< t_DataType> *Vector< t_DataType>::Clone() const
    {
      return new Vector< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object
    //! @return string with the class name
    template< typename t_DataType>
    const std::string &Vector< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return changeable reference to element at POS
    //! @param POS the position of the returned element
    //! @return changeable reference to element at POS
    template< typename t_DataType>
    t_DataType &Vector< t_DataType>::operator()( const size_t POS)
    {
      AssertValidPosition( POS);
      return m_Data[ POS];
    }

    //! @brief return const_reference to element at POS
    //! @param POS the position of the returned element
    //! @return const_reference to element at POS
    template< typename t_DataType>
    const t_DataType &Vector< t_DataType>::operator()( const size_t POS) const
    {
      AssertValidPosition( POS);
      return m_Data[ POS];
    }

    //! @brief return changeable iterator pointing to element at POS
    //! @param POS the position of the element the iterator is pointing to
    //! @return changeable iterator pointing to element at POS
    template< typename t_DataType>
    typename Vector< t_DataType>::iterator Vector< t_DataType>::operator []( const size_t POS)
    {
      AssertValidIteratorPosition( POS);
      return m_Data.begin() + POS;
    }

    //! @brief return const_iterator pointing to element at POS
    //! @param POS the position of the element the iterator is pointing to
    //! @return const_iterator pointing to element at POS
    template< typename t_DataType>
    typename Vector< t_DataType>::const_iterator Vector< t_DataType>::operator []( const size_t POS) const
    {
      AssertValidIteratorPosition( POS);
      return m_Data.begin() + POS;
    }

    //! @brief returns size of the container
    //! @return size, i.e. number of elements stored
    template< typename t_DataType>
    size_t Vector< t_DataType>::GetSize() const
    {
      return m_Data.size();
    }

    //! @brief return the maximal size of the container
    //! @return maximum size, i.e. maximal number of elements to store
    template< typename t_DataType>
    size_t Vector< t_DataType>::MaxSize() const
    {
      return m_Data.max_size();
    }

    //! @brief returns a const reference to the used stl container
    //! @return const reference to the internal stl container
    template< typename t_DataType>
    const std::vector< t_DataType> &Vector< t_DataType>::InternalData() const
    {
      return m_Data;
    }

    //! @brief returns a changeable reference to the used stl container
    //! @return changeable reference to the internal stl container
    template< typename t_DataType>
    std::vector< t_DataType> &Vector< t_DataType>::InternalData()
    {
      return m_Data;
    }

    //! @brief return iterator on begin
    //! @return iterator pointing to the beginning of the container, i.e. the first element
    template< typename t_DataType>
    typename Vector< t_DataType>::iterator Vector< t_DataType>::Begin()
    {
      return m_Data.begin();
    }

    //! @brief return const_iterator on begin
    //! @return const_iterator pointing to the beginning of the container, i.e. the first element
    template< typename t_DataType>
    typename Vector< t_DataType>::const_iterator Vector< t_DataType>::Begin() const
    {
      return m_Data.begin();
    }

    //! @brief return changeable iterator to last element
    //! @return changeable iterator to last element
    template< typename t_DataType>
    typename Vector< t_DataType>::iterator Vector< t_DataType>::Last()
    {
      BCL_Assert( !IsEmpty(), "Cannot return iterator to last element from empty sequence container");
      return --m_Data.end();
    }

    //! @brief return const iterator to last element
    //! @return const iterator to last element
    template< typename t_DataType>
    typename Vector< t_DataType>::const_iterator Vector< t_DataType>::Last() const
    {
      BCL_Assert( !IsEmpty(), "Cannot return iterator to last element from empty sequence container");
      return --m_Data.end();
    }

    //! @brief return iterator on end
    //! @return iterator pointing to the end of the container, i.e. behind the last element
    template< typename t_DataType>
    typename Vector< t_DataType>::iterator Vector< t_DataType>::End()
    {
      return m_Data.end();
    }

    //! @brief return const_iterator on end
    //! @return const_iterator pointing to the end of the container, i.e. behind the last element
    template< typename t_DataType>
    typename Vector< t_DataType>::const_iterator Vector< t_DataType>::End() const
    {
      return m_Data.end();
    }

    //! @brief return iterator to reverse begin
    //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
    template< typename t_DataType>
    typename Vector< t_DataType>::reverse_iterator Vector< t_DataType>::ReverseBegin()
    {
      return m_Data.rbegin();
    }

    //! @brief return const_iterator to reverse begin
    //! @return const_reverse_iterator pointing to the beginning of the reversed container
    template< typename t_DataType>
    typename Vector< t_DataType>::const_reverse_iterator Vector< t_DataType>::ReverseBegin() const
    {
      return m_Data.rbegin();
    }

    //! @brief return iterator to reverse end
    //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
    template< typename t_DataType>
    typename Vector< t_DataType>::reverse_iterator Vector< t_DataType>::ReverseEnd()
    {
      return m_Data.rend();
    }

    //! @brief return const_iterator to reverse end
    //! @return const_reverse_iterator pointing to the end of the reversed container
    template< typename t_DataType>
    typename Vector< t_DataType>::const_reverse_iterator Vector< t_DataType>::ReverseEnd() const
    {
      return m_Data.rend();
    }

    //! @brief return const reference to first element
    //! @return const reference to first element of type t_DataType
    template< typename t_DataType>
    const t_DataType &Vector< t_DataType>::FirstElement() const
    {
      BCL_Assert( !IsEmpty(), "no first element in empty vector");
      return m_Data.front();
    }

    //! @brief return a changeable reference to first element
    template< typename t_DataType>
    t_DataType &Vector< t_DataType>::FirstElement()
    {
      BCL_Assert( !IsEmpty(), "no first element in empty vector");
      return m_Data.front();
    }

    //! @brief return const reference to last element
    //! @return const reference to last element of type t_DataType
    template< typename t_DataType>
    const t_DataType &Vector< t_DataType>::LastElement() const
    {
      BCL_Assert( !IsEmpty(), "no last element in empty vector");
      return m_Data.back();
    }

    //! @brief return changeable reference to last element
    //! @return changeable reference to last element of type t_DataType
    template< typename t_DataType>
    t_DataType &Vector< t_DataType>::LastElement()
    {
      BCL_Assert( !IsEmpty(), "no last element in empty vector");
      return m_Data.back();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief move operator
    template< typename t_DataType>
    Vector< t_DataType> &Vector< t_DataType>::operator =( Vector< t_DataType> && VECTOR)
    {
      // check for different argument
      if( this != &VECTOR)
      {
        // set members
        m_Data = std::move( VECTOR.m_Data);
      }

      // return
      return *this;
    }

    //! @brief equal operator
    template< typename t_DataType>
    Vector< t_DataType> &Vector< t_DataType>::operator =( const Vector< t_DataType> &VECTOR)
    {
      // check for different argument
      if( this != &VECTOR)
      {
        // set members
        m_Data = VECTOR.m_Data;
      }

      // return
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief checks whether container is empty
    //! @return if the container is empty or not
    template< typename t_DataType>
    bool Vector< t_DataType>::IsEmpty() const
    {
      return m_Data.empty();
    }

    //! @brief delete all elements
    template< typename t_DataType>
    void Vector< t_DataType>::Reset()
    {
      m_Data.clear();
    }

    //! @brief insert ELEMENT into the container
    //! @param ELEMENT an object of t_DataType that is inserted
    template< typename t_DataType>
    void Vector< t_DataType>::InsertElement( const t_DataType &ELEMENT)
    {
      m_Data.push_back( ELEMENT);
    }

    //! @brief copies all elements of VECTOR into this object at POS
    //! @param POS the position at which VECTOR is inserted
    //! @param VECTOR the sequence with the elements for insertion
    template< typename t_DataType>
    void Vector< t_DataType>::InsertElements( const size_t POS, const Vector< t_DataType> &VECTOR)
    {
      m_Data.reserve( m_Data.size() + VECTOR.GetSize());
      m_Data.insert( operator[]( POS), VECTOR.Begin(), VECTOR.End());
    }

    //! @brief insert all elements of argument CONTAINER
    //! @param CONTAINER container containing objects of t_DataType that are inserted
    //! @param ITR_POS iterator
    template< typename t_DataType>
    void Vector< t_DataType>::InsertElements( iterator ITR_POS, const Vector< t_DataType> &CONTAINER)
    {
      size_t position( ITR_POS - m_Data.begin());
      m_Data.reserve( m_Data.size() + CONTAINER.GetSize());
      m_Data.insert( m_Data.begin() + position, CONTAINER.Begin(), CONTAINER.End());
    }

    //! @brief insert NUMBER_ELEMENTS copies of argument ELEMENT at POS
    //! @param POS the position for insertion
    //! @param ELEMENT the element to be inserted
    //! @param NUMBER_ELEMENTS number of copies of element to be inserted, 1 by default
    template< typename t_DataType>
    void Vector< t_DataType>::InsertElements( const size_t POS, const t_DataType &ELEMENT, const size_t NUMBER_ELEMENTS)
    {
      m_Data.insert( operator[]( POS), NUMBER_ELEMENTS, ELEMENT);
    }

    //! @brief delete NUMBER_ELEMENTS elements starting at position POS
    //! @param POS position of first element to be removed
    //! @param NUMBER_ELEMENTS number of elements to be removed, 1 by default
    template< typename t_DataType>
    void Vector< t_DataType>::RemoveElements( const size_t POS, const size_t NUMBER_ELEMENTS)
    {
      // check starting position
      AssertValidPosition( POS);

      // check end position if different
      const size_t endpos( POS + NUMBER_ELEMENTS);
      if( NUMBER_ELEMENTS != 1)
      {
        AssertValidPosition( endpos - 1);
      }

      // remove
      m_Data.erase( Begin() + POS, Begin() + endpos);
    }

    //! @brief delete and return element at random position
    //! @return the removed element
    template< typename t_DataType>
    t_DataType Vector< t_DataType>::RemoveRandomElement()
    {
      BCL_Assert( !( IsEmpty()), "Cannot remove element from empty container");
      iterator itr_rem( random::GetGlobalRandom().Iterator( Begin(), End(), GetSize()));
      t_DataType copy_of_element( *itr_rem);
      m_Data.erase( itr_rem);
      return copy_of_element;
    }

    //! @brief allocate memory for a needed length
    //! @param NUMBER_ELEMENTS number of elements for which memory is allocated
    template< typename t_DataType>
    void Vector< t_DataType>::AllocateMemory( const size_t NUMBER_ELEMENTS)
    {
      m_Data.reserve( NUMBER_ELEMENTS);
    }

    //! @brief resizes the sequence to the length NUMBER_ELEMENTS, eventually fill with ELEMENTs
    //! @param NUMBER_ELEMENTS number of elements to resize to
    //! @param ELEMENT optional element of t_DataType to fill the Sequence Container with
    template< typename t_DataType>
    void Vector< t_DataType>::Resize( const size_t NUMBER_ELEMENTS, const t_DataType &ELEMENT)
    {
      m_Data.resize( NUMBER_ELEMENTS, ELEMENT);
    }

    //! @brief set all elements of the sequence to one given value X
    //! @param ELEMENT the element of type t_DataType which will be assigned to all elements, default is t_DataType()
    template< typename t_DataType>
    void Vector< t_DataType>::SetAllElements( const t_DataType &ELEMENT)
    {
      for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        ( *itr) = ELEMENT;
      }
    }

    //! @brief delete single element at ITR
    //! @return iterator pointing to the element immediatly following the deleted one
    //! @param ITR iterator pointing to the element that will be destroyed
    template< typename t_DataType>
    typename Vector< t_DataType>::iterator Vector< t_DataType>::Remove( iterator ITR)
    {
      return ITR == End() ? End() : m_Data.erase( ITR);
    }

    //! @brief delete single element at ITR
    //! @param ITR iterator pointing to the element that will be destroyed
    template< typename t_DataType>
    void Vector< t_DataType>::RemoveElement( iterator ITR)
    {
      ITR == End() ? End() : m_Data.erase( ITR);
    }

    //! @brief append a ELEMENT to the end of *this
    //! @param ELEMENT an object of t_DataType that is appended
    template< typename t_DataType>
    Vector< t_DataType> &Vector< t_DataType>::Append( const t_DataType &ELEMENT)
    {
      InsertElement( ELEMENT);
      return *this;
    }

    //! @brief append a CONTAINER to the end of *this
    //! @param CONTAINER container containing objects of t_DataType that are appended
    template< typename t_DataType>
    Vector< t_DataType> &Vector< t_DataType>::Append( const Vector< t_DataType> &CONTAINER)
    {
      InsertElements( End(), CONTAINER);
      return *this;
    }

    //! @brief inserts ELEMENT to the end
    //! @param ELEMENT element to be inserted, default element is t_DataType()
    template< typename t_DataType>
    void Vector< t_DataType>::PushBack( const t_DataType &ELEMENT)
    {
      m_Data.push_back( ELEMENT);
    }

    //! @brief inserts ELEMENT to the end
    //! @param ELEMENT element to be inserted, default element is t_DataType()
    template< typename t_DataType>
    void Vector< t_DataType>::PushBack( t_DataType && ELEMENT)
    {
      m_Data.push_back( std::move( ELEMENT));
    }

    //! @brief removes last element
    template< typename t_DataType>
    void Vector< t_DataType>::PopBack()
    {
      m_Data.pop_back();
    }

    //! @brief Randomly permutes the order of all elements in the vector.
    //! After a call to Shuffle, each element has an
    //! equal chance of being at any given location in the vector
    template< typename t_DataType>
    void Vector< t_DataType>::Shuffle()
    {
      std::random_shuffle( Begin(), End(), random::GetRandomSizeT);
    }

    //! @brief Puts the elements of a vector in a different order
    //! @param NEW_ORDER Indices to map the indices of this vector onto
    //! @note  this algorithm is suboptimal in that it goes to the trouble of copying
    //! @note  this vectors elements.  Use/Write a different algorithm if memory is an
    //! @note  issue or copy construction is slow
    template< typename t_DataType>
    void Vector< t_DataType>::Reorder( const Vector< size_t> &NEW_ORDER)
    {
      // copy this vector
      Vector< t_DataType> newbie( *this);

      // resize the vector, if necessary
      if( NEW_ORDER.GetSize() != GetSize())
      {
        Resize( NEW_ORDER.GetSize());
      }

      // copy the elements back according to the NEW_ORDER
      for( size_t i( 0), sz( m_Data.size()); i < sz; ++i)
      {
        m_Data[ i] = newbie( NEW_ORDER( i));
      }
    }

    //! @brief creates the storagevector from two elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A
    )
    {
      // return vector
      return Vector< t_DataType>( size_t( 1), DATA_A);
    }

    //! @brief creates the storagevector from two elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B
    )
    {
      // create vector from two elements
      Vector< t_DataType> vector( 2, DATA_A);
      vector( 1) = DATA_B;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from three elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C
    )
    {
      // create vector from three elements
      Vector< t_DataType> vector( 3, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from four elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D
    )
    {
      // create vector from four elements
      Vector< t_DataType> vector( 4, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from five elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E
    )
    {
      // create vector from five elements
      Vector< t_DataType> vector( 5, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from six elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F
    )
    {
      // create vector from six elements
      Vector< t_DataType> vector( 6, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;
      vector( 5) = DATA_F;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from seven elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G
    )
    {
      // create vector from seven elements
      Vector< t_DataType> vector( 7, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;
      vector( 5) = DATA_F;
      vector( 6) = DATA_G;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from eight elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G,
      const t_DataType &DATA_H
    )
    {
      // create vector from eight elements
      Vector< t_DataType> vector( 8, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;
      vector( 5) = DATA_F;
      vector( 6) = DATA_G;
      vector( 7) = DATA_H;

      // return vector
      return vector;

    }

    //! @brief creates the storagevector from nine elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G,
      const t_DataType &DATA_H,
      const t_DataType &DATA_I
    )
    {
      // create vector from nine elements
      Vector< t_DataType> vector( 9, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;
      vector( 5) = DATA_F;
      vector( 6) = DATA_G;
      vector( 7) = DATA_H;
      vector( 8) = DATA_I;

      // return vector
      return vector;

    }

    //! @brief creates the storagevector from ten elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G,
      const t_DataType &DATA_H,
      const t_DataType &DATA_I,
      const t_DataType &DATA_J
    )
    {
      // create vector from ten elements
      Vector< t_DataType> vector( 10, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;
      vector( 5) = DATA_F;
      vector( 6) = DATA_G;
      vector( 7) = DATA_H;
      vector( 8) = DATA_I;
      vector( 9) = DATA_J;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from 11 elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G,
      const t_DataType &DATA_H,
      const t_DataType &DATA_I,
      const t_DataType &DATA_J,
      const t_DataType &DATA_K
    )
    {
      // create vector from ten elements
      Vector< t_DataType> vector( 11, DATA_A);
      vector(  1) = DATA_B;
      vector(  2) = DATA_C;
      vector(  3) = DATA_D;
      vector(  4) = DATA_E;
      vector(  5) = DATA_F;
      vector(  6) = DATA_G;
      vector(  7) = DATA_H;
      vector(  8) = DATA_I;
      vector(  9) = DATA_J;
      vector( 10) = DATA_K;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from 11 elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G,
      const t_DataType &DATA_H,
      const t_DataType &DATA_I,
      const t_DataType &DATA_J,
      const t_DataType &DATA_K,
      const t_DataType &DATA_L
    )
    {
      // create vector from ten elements
      Vector< t_DataType> vector( 12, DATA_A);
      vector(  1) = DATA_B;
      vector(  2) = DATA_C;
      vector(  3) = DATA_D;
      vector(  4) = DATA_E;
      vector(  5) = DATA_F;
      vector(  6) = DATA_G;
      vector(  7) = DATA_H;
      vector(  8) = DATA_I;
      vector(  9) = DATA_J;
      vector( 10) = DATA_K;
      vector( 11) = DATA_L;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from 15 elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G,
      const t_DataType &DATA_H,
      const t_DataType &DATA_I,
      const t_DataType &DATA_J,
      const t_DataType &DATA_K,
      const t_DataType &DATA_L,
      const t_DataType &DATA_M,
      const t_DataType &DATA_N,
      const t_DataType &DATA_O
    )
    {
      // create vector from 15 elements
      Vector< t_DataType> vector( 15, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;
      vector( 5) = DATA_F;
      vector( 6) = DATA_G;
      vector( 7) = DATA_H;
      vector( 8) = DATA_I;
      vector( 9) = DATA_J;
      vector( 10) = DATA_K;
      vector( 11) = DATA_L;
      vector( 12) = DATA_M;
      vector( 13) = DATA_N;
      vector( 14) = DATA_O;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from 16 elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G,
      const t_DataType &DATA_H,
      const t_DataType &DATA_I,
      const t_DataType &DATA_J,
      const t_DataType &DATA_K,
      const t_DataType &DATA_L,
      const t_DataType &DATA_M,
      const t_DataType &DATA_N,
      const t_DataType &DATA_O,
      const t_DataType &DATA_P
    )
    {
      // create vector from 16 elements
      Vector< t_DataType> vector( 16, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;
      vector( 5) = DATA_F;
      vector( 6) = DATA_G;
      vector( 7) = DATA_H;
      vector( 8) = DATA_I;
      vector( 9) = DATA_J;
      vector( 10) = DATA_K;
      vector( 11) = DATA_L;
      vector( 12) = DATA_M;
      vector( 13) = DATA_N;
      vector( 14) = DATA_O;
      vector( 15) = DATA_P;

      // return vector
      return vector;
    }

    //! @brief creates the storagevector from 22 elements
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::Create
    (
      const t_DataType &DATA_A,
      const t_DataType &DATA_B,
      const t_DataType &DATA_C,
      const t_DataType &DATA_D,
      const t_DataType &DATA_E,
      const t_DataType &DATA_F,
      const t_DataType &DATA_G,
      const t_DataType &DATA_H,
      const t_DataType &DATA_I,
      const t_DataType &DATA_J,
      const t_DataType &DATA_K,
      const t_DataType &DATA_L,
      const t_DataType &DATA_M,
      const t_DataType &DATA_N,
      const t_DataType &DATA_O,
      const t_DataType &DATA_P,
      const t_DataType &DATA_Q,
      const t_DataType &DATA_R,
      const t_DataType &DATA_S,
      const t_DataType &DATA_T,
      const t_DataType &DATA_U,
      const t_DataType &DATA_V
    )
    {
      // create vector from 22 elements
      Vector< t_DataType> vector( 22, DATA_A);
      vector( 1) = DATA_B;
      vector( 2) = DATA_C;
      vector( 3) = DATA_D;
      vector( 4) = DATA_E;
      vector( 5) = DATA_F;
      vector( 6) = DATA_G;
      vector( 7) = DATA_H;
      vector( 8) = DATA_I;
      vector( 9) = DATA_J;
      vector( 10) = DATA_K;
      vector( 11) = DATA_L;
      vector( 12) = DATA_M;
      vector( 13) = DATA_N;
      vector( 14) = DATA_O;
      vector( 15) = DATA_P;
      vector( 16) = DATA_Q;
      vector( 17) = DATA_R;
      vector( 18) = DATA_S;
      vector( 19) = DATA_T;
      vector( 20) = DATA_U;
      vector( 21) = DATA_V;

      // return vector
      return vector;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief exits if given an invalid position
    //! @param POS the position to check
    template< typename t_DataType>
    void Vector< t_DataType>::AssertValidPosition( const size_t POS) const
    {
      BCL_Assert
      (
        POS < GetSize(),
        util::Format()( POS) + " is outside range: [ 0 .. " + util::Format()( GetSize()) + " ) for type: " + GetStaticClassName< t_DataType>()
      );
    }

    //! @brief exits if given an invalid position
    //! @param POS the position to check
    template< typename t_DataType>
    void Vector< t_DataType>::AssertValidIteratorPosition( const size_t POS) const
    {
      BCL_Assert
      (
        POS <= GetSize(),
        util::Format()( POS) + " is outside range: [ 0 .. " + util::Format()( GetSize()) + " ]"
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write container to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &Vector< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write size of container
      io::Serialize::Write( GetSize(), OSTREAM, INDENT) << '\n';

      // write each element of container to new line
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        io::Serialize::Write( *itr, OSTREAM, INDENT) << '\n';
      }

      OSTREAM << std::flush;

      // end
      return OSTREAM;
    }

    //! @brief read container from io::IFStream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Vector< t_DataType>::Read( std::istream &ISTREAM)
    {
      // create int "size" for holding the size of the container
      int size;

      // read in size of container
      io::Serialize::Read( size, ISTREAM);

      // match container size to size of incoming data
      Resize( size);

      // read in each element from stream
      for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        io::Serialize::Read( *itr, ISTREAM);
        BCL_Assert
        (
          ISTREAM.good() || ISTREAM.eof(),
          "Error reading element! This was read:\n" + util::Format()( ( *itr))
          + " at position " + util::Format()( size_t( itr - Begin()))
        );
      }

      // return
      return ISTREAM;
    }

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Vector< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Vector< t_DataType>())
    );

  } // namespace storage
} // namespace bcl

#endif // BCL_STORAGE_VECTOR_HPP_
