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

#ifndef BCL_STORAGE_VECTOR_H_
#define BCL_STORAGE_VECTOR_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <limits>
#include <vector>

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Vector
    //! @brief This is a template class for Vector Containers.
    //! @details Attributes: (copied from http://www.sgi.com/tech/stl/Vector.html)
    //! - random access to elements
    //! - constant time insertion and removal of elements at end of sequence
    //! - linear time insertion and removal of elements in the middle or in the beginning
    //! - iterators are invalidated when memory is reallocated
    //! - inserting or deleting element invalidates all iterators after that point
    //!
    //! @tparam t_DataType the type of data to be held by the container
    //!
    //! @see @link example_storage_vector.cpp @endlink
    //! @author heinzes1
    //! @date 9/5/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Vector :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      std::vector< t_DataType> m_Data;

    public:

      //! typedef for iterator
      typedef typename std::vector< t_DataType>::iterator               iterator;

      //! typedef for const_iterator
      typedef typename std::vector< t_DataType>::const_iterator         const_iterator;

      //! typedef for reverse_iterator
      typedef typename std::vector< t_DataType>::reverse_iterator       reverse_iterator;

      //! typedef for const_reverse_iterator
      typedef typename std::vector< t_DataType>::const_reverse_iterator const_reverse_iterator;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Vector();

      //! @brief construct a Vector from optional size and single element
      //! @param SIZE optional size, default 0
      //! @param VALUE optional single element which is assigned to all stored elements, default is t_DataType()
      Vector( const size_t SIZE, const t_DataType &VALUE = t_DataType());

      //! @brief construct a Vector from number of size and pointer to data
      //! @param SIZE number of elements behind the pointer *DATA
      //! @param DATA pointer to the elements
      Vector( const size_t SIZE, const t_DataType *DATA);

      //! @brief construct a Vector from another Vector
      //! @param VECTOR the Deque to copy
      Vector( const Vector< t_DataType> &VECTOR);

      //! @brief move operator
      //! @param VECTOR the Deque to copy
      Vector( Vector< t_DataType> && VECTOR);

      //! @brief construct a Vector from another Vector
      //! @param VECTOR the Deque to copy
      Vector( const std::vector< t_DataType> &VECTOR);

      //! @brief constuct a Vector from a subsequence [ POS, POS + LENGTH) of another Vector
      //! @param VECTOR the Vector to copy
      //! @param POS the first element to copy
      //! @param LENGTH the (maximum) number of elements to use
      //! in no event will elements beyond the end of VECTOR be copied
      Vector
      (
        const Vector< t_DataType> &VECTOR,
        const size_t POS,
        const size_t LENGTH = std::numeric_limits< size_t>::max()
      );

      //! @brief construct a Vector from iterator [FIRST, LAST) range
      //! @param FIRST iterator to the first element to be copied
      //! @param LAST iterator to the first element after the last copied element
      template< typename t_OtherDataType>
      explicit Vector( const Vector< t_OtherDataType> &VECTOR) :
        m_Data( VECTOR.Begin(), VECTOR.End())
      {
      }

      //! @brief construct a Vector from iterator [FIRST, LAST) range
      //! @param FIRST iterator to the first element to be copied
      //! @param LAST iterator to the first element after the last copied element
      template< typename t_InputIterator>
      Vector( const t_InputIterator &FIRST, const t_InputIterator &LAST) :
        m_Data( FIRST, LAST)
      {
      }

      //! @brief copy constructor
      //! @return pointer to a copy of the actual object
      Vector< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object
      //! @return string with the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return changeable reference to element at POS
      //! @param POS the position of the returned element
      //! @return changeable reference to element at POS
      t_DataType &operator()( const size_t POS);

      //! @brief return const_reference to element at POS
      //! @param POS the position of the returned element
      //! @return const_reference to element at POS
      const t_DataType &operator()( const size_t POS) const;

      //! @brief return changeable iterator pointing to element at POS
      //! @param POS the position of the element the iterator is pointing to
      //! @return changeable iterator pointing to element at POS
      iterator operator []( const size_t POS);

      //! @brief return const_iterator pointing to element at POS
      //! @param POS the position of the element the iterator is pointing to
      //! @return const_iterator pointing to element at POS
      const_iterator operator []( const size_t POS) const;

      //! @brief returns size of the container
      //! @return size, i.e. number of elements stored
      size_t GetSize() const;

      //! @brief return the maximal size of the container
      //! @return maximum size, i.e. maximal number of elements to store
      size_t MaxSize() const;

      //! @brief returns a const reference to the used stl container
      //! @return const reference to the internal stl container
      const std::vector< t_DataType> &InternalData() const;

      //! @brief returns a changeable reference to the used stl container
      //! @return changeable reference to the internal stl container
      std::vector< t_DataType> &InternalData();

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin();

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const;

      //! @brief return changeable iterator to last element
      //! @return changeable iterator to last element
      iterator Last();

      //! @brief return const iterator to last element
      //! @return const iterator to last element
      const_iterator Last() const;

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End();

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const;

      //! @brief return iterator to reverse begin
      //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
      reverse_iterator ReverseBegin();

      //! @brief return const_iterator to reverse begin
      //! @return const_reverse_iterator pointing to the beginning of the reversed container
      const_reverse_iterator ReverseBegin() const;

      //! @brief return iterator to reverse end
      //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
      reverse_iterator ReverseEnd();

      //! @brief return const_iterator to reverse end
      //! @return const_reverse_iterator pointing to the end of the reversed container
      const_reverse_iterator ReverseEnd() const;

      //! @brief return const reference to first element
      //! @return const reference to first element of type t_DataType
      t_DataType const &FirstElement() const;

      //! @brief return a changeable reference to first element
      t_DataType &FirstElement();

      //! @brief return const reference to last element
      //! @return const reference to last element of type t_DataType
      t_DataType const &LastElement() const;

      //! @brief return changeable reference to last element
      //! @return changeable reference to last element of type t_DataType
      t_DataType &LastElement();

    ///////////////
    // operators //
    ///////////////

      //! @brief equal operator
      Vector< t_DataType> &operator =( const Vector< t_DataType> &VECTOR);

      //! @brief move assignment operator
      Vector< t_DataType> &operator =( Vector< t_DataType> && VECTOR);

    ////////////////
    // operations //
    ////////////////

      //! @brief checks whether container is empty
      //! @return if the container is empty or not
      bool IsEmpty() const;

      //! @brief delete all elements
      void Reset();

      //! @brief insert ELEMENT into the container
      //! @param ELEMENT an object of t_DataType that is inserted
      void InsertElement( const t_DataType &ELEMENT);

      //! @brief insert ELEMENT after ITR into the container
      //! @param ITR iterator on the Sequence Container whereafter ELEMENT is inserted
      //! @param ELEMENT an object of t_DataType that is inserted
      template< typename t_InputIterator>
      void InsertElement( const t_InputIterator ITR, const t_DataType &ELEMENT)
      {
        m_Data.insert( ITR, ELEMENT);
      }

      //! @brief copies all elements of VECTOR into this object at POS
      //! @param POS the position at which VECTOR is inserted
      //! @param VECTOR the sequence with the elements for insertion
      void InsertElements( const size_t POS, const Vector< t_DataType> &VECTOR);

      //! @brief insert all elements of argument CONTAINER
      //! @param CONTAINER container containing objects of t_DataType that are inserted
      //! @param ITR_POS iterator
      void InsertElements( iterator ITR_POS, const Vector< t_DataType> &CONTAINER);
      //! @brief insert all elements of argument CONTAINER
      //! @param ITR_POS iterator
      //! @param BEGIN first element
      //! @param END last element
      template< typename t_InputIterator>
      void InsertElements( iterator ITR_POS, t_InputIterator BEGIN, t_InputIterator END)
      {
        m_Data.insert( ITR_POS, BEGIN, END);
      }

      //! @brief insert NUMBER_ELEMENTS copies of argument ELEMENT at POS
      //! @param POS the position for insertion
      //! @param ELEMENT the element to be inserted
      //! @param NUMBER_ELEMENTS number of copies of element to be inserted, 1 by default
      void InsertElements( const size_t POS, const t_DataType &ELEMENT, const size_t NUMBER_ELEMENTS = 1);

      //! @brief delete NUMBER_ELEMENTS elements starting at position POS
      //! @param POS position of first element to be removed
      //! @param NUMBER_ELEMENTS number of elements to be removed, 1 by default
      void RemoveElements( const size_t POS, const size_t NUMBER_ELEMENTS = 1);

      //! @brief delete and return element at random position
      //! @return the removed element
      t_DataType RemoveRandomElement();

      //! @brief allocate memory for a needed length
      //! @param NUMBER_ELEMENTS number of elements for which memory is allocated
      void AllocateMemory( const size_t NUMBER_ELEMENTS);

      //! @brief resizes the sequence to the length NUMBER_ELEMENTS, eventually fill with ELEMENTs
      //! @param NUMBER_ELEMENTS number of elements to resize to
      //! @param ELEMENT optional element of t_DataType to fill the Sequence Container with
      void Resize( const size_t NUMBER_ELEMENTS, const t_DataType &ELEMENT = t_DataType());

      //! @brief set all elements of the sequence to one given value X
      //! @param ELEMENT the element of type t_DataType which will be assigned to all elements, default is t_DataType()
      void SetAllElements( const t_DataType &ELEMENT = t_DataType());

      //! @brief delete single element at ITR
      //! @return iterator pointing to the element immediatly following the deleted one
      //! @param ITR iterator pointing to the element that will be destroyed
      iterator Remove( iterator ITR);

      //! @brief delete single element at ITR
      //! @param ITR iterator pointing to the element that will be destroyed
      void RemoveElement( iterator ITR);

      //! @brief append a ELEMENT to the end of *this
      //! @param ELEMENT an object of t_DataType that is appended
      Vector< t_DataType> &Append( const t_DataType &ELEMENT);

      //! @brief append a CONTAINER to the end of *this
      //! @param CONTAINER container containing objects of t_DataType that are appended
      Vector< t_DataType> &Append( const Vector< t_DataType> &CONTAINER);

      //! @brief inserts ELEMENT to the end
      //! @param ELEMENT element to be inserted, default element is t_DataType()
      void PushBack( const t_DataType &ELEMENT = t_DataType());

      //! @brief inserts ELEMENT to the end
      //! @param ELEMENT element to be inserted, default element is t_DataType()
      void PushBack( t_DataType && ELEMENT);

      //! @brief removes last element
      void PopBack();

      //! @brief Randomly permutes the order of all elements in the vector.
      //! After a call to Shuffle, each element has an
      //! equal chance of being at any given location in the vector
      void Shuffle();

      //! @brief Sort sorts the vector according to a binary predicate
      //! in approx N log N, where N is the length of the list
      //! @tparam t_BinaryPredicate type of binary predicate that will be used for comparison
      //! @param COMPARABLE t_BinaryPredicate used to compare two elements of the vector
      template< typename t_BinaryPredicate>
      void Sort( t_BinaryPredicate COMPARABLE)
      {
        std::sort( Begin(), End(), COMPARABLE);
      }

      //! @brief find an element in the vector
      //! @param ELEMENT what to search for.  t_OtherDataType must be comparable (via ==) to t_DataType
      //! @param START_POSITION the index to begin searching for ELEMENT at
      //! @return the index of the vector where ELEMENT is located.  If ELEMENT is not in the vector, this will be the
      //! @return size of the vector
      //! @note  test whether Find was successful on a vector (V) with V.Find(...) < V.GetSize()
      template< typename t_OtherDataType>
      size_t Find( const t_OtherDataType &ELEMENT, const size_t START_POSITION = 0) const
      {
        return size_t( std::find( operator[]( START_POSITION), End(), ELEMENT) - Begin());
      }

      //! @brief Puts the elements of a vector in a different order
      //! @param NEW_ORDER Indices to map the indices of this vector onto
      //! @note  this algorithm is suboptimal in that it goes to the trouble of copying
      //! @note  this vectors elements.  Use/Write a different algorithm if memory is an
      //! @note  issue or copy construction is slow
      void Reorder( const Vector< size_t> &NEW_ORDER);

      //! @brief creates the storagevector from one element
      static Vector< t_DataType> Create( const t_DataType &DATA_A);

      //! @brief creates the storagevector from two elements
      static Vector< t_DataType> Create
      (
        const t_DataType &DATA_A,
        const t_DataType &DATA_B
      );

      //! @brief creates the storagevector from three elements
      static Vector< t_DataType> Create
      (
        const t_DataType &DATA_A,
        const t_DataType &DATA_B,
        const t_DataType &DATA_C
      );

      //! @brief creates the storagevector from four elements
      static Vector< t_DataType> Create
      (
        const t_DataType &DATA_A,
        const t_DataType &DATA_B,
        const t_DataType &DATA_C,
        const t_DataType &DATA_D
      );

      //! @brief creates the storagevector from five elements
      static Vector< t_DataType> Create
      (
        const t_DataType &DATA_A,
        const t_DataType &DATA_B,
        const t_DataType &DATA_C,
        const t_DataType &DATA_D,
        const t_DataType &DATA_E
      );

      //! @brief creates the storagevector from six elements
      static Vector< t_DataType> Create
      (
        const t_DataType &DATA_A,
        const t_DataType &DATA_B,
        const t_DataType &DATA_C,
        const t_DataType &DATA_D,
        const t_DataType &DATA_E,
        const t_DataType &DATA_F
      );

      //! @brief creates the storagevector from seven elements
      static Vector< t_DataType> Create
      (
        const t_DataType &DATA_A,
        const t_DataType &DATA_B,
        const t_DataType &DATA_C,
        const t_DataType &DATA_D,
        const t_DataType &DATA_E,
        const t_DataType &DATA_F,
        const t_DataType &DATA_G
      );

      //! @brief creates the storagevector from eight elements
      static Vector< t_DataType> Create
      (
        const t_DataType &DATA_A,
        const t_DataType &DATA_B,
        const t_DataType &DATA_C,
        const t_DataType &DATA_D,
        const t_DataType &DATA_E,
        const t_DataType &DATA_F,
        const t_DataType &DATA_G,
        const t_DataType &DATA_H
      );

      //! @brief creates the storagevector from nine elements
      static Vector< t_DataType> Create
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
      );

      //! @brief creates the storagevector from ten elements
      static Vector< t_DataType> Create
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
      );

      //! @brief creates the storagevector from 11 elements
      static Vector< t_DataType> Create
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
      );

      //! @brief creates the storagevector from 11 elements
      static Vector< t_DataType> Create
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
      );

      //! @brief creates the storagevector from 15 elements
      static Vector< t_DataType> Create
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
      );

      //! @brief creates the storagevector from 16 elements
      static Vector< t_DataType> Create
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
      );

      //! @brief creates the storagevector from 22 elements
      static Vector< t_DataType> Create
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
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief exits if given an invalid position for dereferencing (POS >= GetSize())
      //! @param POS the position to check
      void AssertValidPosition( const size_t POS) const;

      //! @brief exits if given an invalid position for an iterator (POS > GetSize())
      //! @param POS the position to check
      void AssertValidIteratorPosition( const size_t POS) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief write container to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read container from io::IFStream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

    }; // template class Vector

    //! @brief operator == checks if two Vectors are the same
    //! @param VECTOR_A first container
    //! @param VECTOR_B second container
    //! @return true, if vectors are identical
    template< typename t_DataType>
    bool operator ==
    (
      const Vector< t_DataType> &VECTOR_A,
      const Vector< t_DataType> &VECTOR_B
    )
    {
      return VECTOR_A.InternalData() == VECTOR_B.InternalData();
    }

    //! @brief operator == checks if two Vectors are the same
    //! @param VECTOR_A first container
    //! @param VECTOR_B second container
    //! @return true, if vectors are identical
    template< typename t_DataType>
    bool operator !=
    (
      const Vector< t_DataType> &VECTOR_A,
      const Vector< t_DataType> &VECTOR_B
    )
    {
      return VECTOR_A.InternalData() != VECTOR_B.InternalData();
    }

    //! @brief operator < to compare to vectors
    //! @param VECTOR_LHS left hand side operand
    //! @param VECTOR_RHS right hand side operand
    //! @return true if the first position in VECTOR_LHS LHS that differs from VECTOR_RHS is less
    template< typename t_DataType>
    bool operator <( const Vector< t_DataType> &VECTOR_LHS, const Vector< t_DataType> &VECTOR_RHS)
    {
      // compare the first non-equal elements according to <
      for( size_t i( 0), vec_size( std::min( VECTOR_LHS.GetSize(), VECTOR_RHS.GetSize())); i < vec_size; ++i)
      {
        if( !( VECTOR_LHS( i) == VECTOR_RHS( i)))
        {
          return VECTOR_LHS( i) < VECTOR_RHS( i);
        }
      }
      return VECTOR_LHS.GetSize() < VECTOR_RHS.GetSize();
    }

    //! @brief create an index vector consisting of N integers [0, N), representing indices
    //! @param N the upper end of the index vector range
    //! @return vector containing 0,...,N-1
    BCL_API Vector< size_t> CreateIndexVector( const size_t &N);

    //! @brief create an index vector consisting of consecutive integers in range [MIN, MAX), representing indices
    //! @param MIN the lower end of the range
    //! @param MAX the upper end of the range (exclusive)
    //! @return vector containing MIN,...,MAX-1
    BCL_API Vector< size_t> CreateIndexVector( const size_t &MIN, const size_t &MAX);

    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< std::string>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< size_t>;
  } // namespace storage
} // namespace bcl

// include template implementation file
#include "bcl_storage_vector.hpp"

#endif // BCL_STORAGE_VECTOR_H_
