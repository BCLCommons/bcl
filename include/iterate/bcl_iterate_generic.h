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

#ifndef BCL_ITERATE_GENERIC_H_
#define BCL_ITERATE_GENERIC_H_

// include the namespace header
#include "bcl_iterate.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_iterate_with_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace iterate
  {
    /////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Generic
    //! @tparam t_Return type that the iterator returns when dereferenced ( e.g. qualified return type of *iterator)
    //!
    //! This class encapsulates common iterator functionality without requiring the user to know the container type of
    //! the iterator.  Several convenience functions extend the classes functionality
    //!
    //! @see @link example_iterate_generic.cpp @endlink
    //! @author mendenjl
    //! @date Jul 14, 2012
    //!
    /////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Return>
    class Generic :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      Interface< t_Return> *m_Iterator; //!< Pointer to a Generic iterator WithRange

    public:

      typedef t_Return&    reference;  //!< type returned by operator*
      typedef t_Return     value_type; //!< type needed to hold the result of an operator* (never a reference)
      typedef t_Return*    pointer;    //!< type returned by operator->
      typedef long         difference_type; //!< type returned by std::distance for this iterator
      typedef typename std::bidirectional_iterator_tag iterator_category; //!< Type of iterator

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Generic() :
        m_Iterator( new WithRange< pointer, t_Return>( NULL, NULL))
      {
      }

      //! @brief copy constructor
      Generic( const Generic &ORIGINAL) :
        m_Iterator( ORIGINAL.m_Iterator ? ORIGINAL.m_Iterator->Clone() : NULL)
      {
      }

      //! @brief constructor from interface
      template< typename t_OtherIteratorType>
      explicit Generic( const t_OtherIteratorType &ORIGINAL) :
        m_Iterator( ORIGINAL.CloneIterator())
      {
      }

      //! @brief constructor from GenericBase constructed with new
      explicit Generic( Interface< t_Return> *const ORIGINAL) :
        m_Iterator( ORIGINAL)
      {
      }

      //! @brief construct from begin and end of iterated-upon range
      //! @param BEGIN beginning of iterated range
      //! @param END end of the iterated range
      template< typename t_Iterator>
      Generic( t_Iterator BEGIN, t_Iterator END) :
        m_Iterator( new WithRange< t_Iterator, t_Return>( BEGIN, END))
      {
      }

      //! @brief construct from a non-const container
      //! @param BEGIN beginning of iterated range
      //! @param END end of the iterated range
      //! @param ITR the position at which the iterator should start
      template< typename t_Iterator>
      Generic( t_Iterator BEGIN, t_Iterator END, t_Iterator ITR) :
        m_Iterator( new WithRange< t_Iterator, t_Return>( BEGIN, END, ITR))
      {
      }

      //! @brief destructor
      ~Generic()
      {
        delete m_Iterator;
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief virtual copy constructor
      Generic *Clone() const
      {
        return new Generic( *this);
      }

      //! @brief determine whether the iterator is at the beginning
      //! @return true iff the iterator is at the beginning
      bool AtBegin() const
      {
        return m_Iterator->AtBegin();
      }

      //! @return true iff the iterator is at the end of the array
      bool NotAtEnd() const
      {
        return m_Iterator->NotAtEnd();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset the underlying iterator to the start of its range (also resets displacement)
      //! @param TO_BEGINNING true -> start at beginning, false -> start at reverse beginning (last element)
      void Restart()
      {
        m_Iterator->Restart();
      }

      //! @brief reset the underlying iterator to the start of its range (also resets displacement)
      Generic &GotoBegin()
      {
        m_Iterator->Restart();
        return *this;
      }

      //! @brief reset the underlying iterator to the end of its range (also determines displacement)
      Generic &GotoEnd()
      {
        m_Iterator->GotoEnd();
        return *this;
      }

      //! @brief return an iterator to the beginning of the range
      Generic Begin() const
      {
        return Generic( m_Iterator->Begin());
      }

      //! @brief return an iterator to the end of the range
      Generic End() const
      {
        return Generic( m_Iterator->End());
      }

      //! @brief determine how many elements are between the iterator and the first element of the range
      //! @return a long indicating the iterators distance from Begin()
      //! @note return type has to be a long because iterator could be before Begin()
      //! @note O(1) for vectors, O( size of range) for everything else
      long GetPosition() const
      {
        return m_Iterator->GetPosition();
      }

      //! @brief determine how many elements are between the iterator and the last element of the range
      //! @return a long indicating the iterators distance from the last element ()
      //! @note return type has to be a long because iterator could be before End()
      //! @note O(1) for vectors, O( size of range) for everything else
      long GetReversePosition() const
      {
        return m_Iterator->GetReversePosition();
      }

      //! @brief get a hard copy of the interface object
      Interface< t_Return> *CloneIterator() const
      {
        return m_Iterator ? m_Iterator->Clone() : NULL;
      }

      //! @return the size of the range
      //! @note O(1) for vectors, O( size of range) for everything else
      size_t GetSize() const
      {
        return m_Iterator->GetSize();
      }

      //! @return true if both iterators are pointing to the same chunk in memory
      //! @note yields compile time error if the two iterators have different pointer types
      template< typename t_OtherIterator>
      bool operator ==( t_OtherIterator const &ITR) const
      {
        return m_Iterator->operator ==( ITR);
      }

      //! @return true if both iterators are pointing to the same chunk in memory
      //! @note yields compile time error if the two iterators have different pointer types
      template< typename t_OtherIterator>
      bool operator !=( t_OtherIterator const &ITR) const
      {
        return m_Iterator->operator !=( ITR);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Asserts that the iterator is valid and prints a message if it is not
      void AssertNotAtEnd() const
      {
        m_Iterator->AssertNotAtEnd();
      }

      //! @brief moves the iterator to the position provided or the end of the range, whichever comes first
      //! @param POSITION the position to go to
      //! @return whether the iterator is still valid
      bool GotoPosition( const size_t &POSITION)
      {
        return m_Iterator->GotoPosition( POSITION);
      }

      //! @brief iterator forward until an element that satisfies CONDITION is reached
      //! @return true if a later element satisfied the condition
      template< typename t_Condition>
      bool GotoNextElementThatSatisfies( t_Condition CONDITION)
      {
        return m_Iterator->GotoNextElementThatSatisfies( CONDITION);
      }

      //! @brief iterator forward until an element that satisfies CONDITION is reached
      //! @return true if an earlier element satisfied the condition
      template< typename t_Condition>
      bool GotoPrevElementThatSatisfies( t_Condition CONDITION)
      {
        return m_Iterator->GotoPrevElementThatSatisfies( CONDITION);
      }

      //! @brief go to a random position in within the range
      //! @return the position that the iterator went to
      size_t GotoRandomPosition()
      {
        return m_Iterator->GotoRandomPosition();
      }

      //! @param CONDITION the condition to use to decide how whether to count a given element
      //! @return the number of elements in the range that satisfy the condition
      template< typename t_Condition>
      size_t CountElementsThatSatisfy( t_Condition CONDITION) const
      {
        return m_Iterator->CountElementsThatSatisfy( CONDITION);
      }

      //! @param CONDITION the condition to use to decide how whether to count a given element
      //! @return the number of elements in the range that satisfy the condition
      template< typename t_Condition, typename t_Container>
      void InsertElementsThatSatisfyConditionInto( t_Condition CONDITION, t_Container &CONTAINER) const
      {
        m_Iterator->CountElementsThatSatisfy( CONDITION, CONTAINER);
      }

      //! @param CONDITION the condition to use to decide how whether to count a given element
      //! @return the number of elements in the range that satisfy the condition
      //! Not called GetSize because the size is not known by the iterator
      template< typename t_Condition>
      size_t RandomlyPickElementThatSatisfies( t_Condition CONDITION)
      {
        return m_Iterator->RandomlyPickElementThatSatisfies( CONDITION);
      }

      //! @brief add all elements following the current iterator to a certain container
      //! @param BASE a Generic iterator base that contains things to iterator over
      //! @param CONTAINER the container to store the objects in
      template< class t_Container>
      void InsertElementsFromIteratorInto( t_Container CONTAINER) const
      {
        m_Iterator->InsertElementsFromIteratorInto( CONTAINER);
      }

    ///////////////
    // operators //
    ///////////////

      //! @return the data member
      reference operator *() const
      {
        return **m_Iterator;
      }

      //! @brief get the pointer to the data member
      //! @return a pointer to the data member
      //! If the return value of operator* is a temporary value,
      //! then the return value will be converted to void * to prevent dereferencing
      pointer operator ->() const
      {
        return m_Iterator->operator->();
      }

      //! @brief assignment operator
      //! @param ORIGINAL the original iterator
      //! @return reference to this iterator
      Generic &operator =( const Generic< t_Return> &ORIGINAL)
      {
        if( m_Iterator != ORIGINAL.m_Iterator)
        {
          delete m_Iterator;
          m_Iterator = ORIGINAL.m_Iterator ? ORIGINAL.m_Iterator->Clone() : NULL;
        }
        return *this;
      }

      //! @brief operator ++ (prefix, e.g. ++a)
      //! @return a reference to the iterator after incrementing
      Generic &operator ++()
      {
        ++*m_Iterator;
        return *this;
      }

      //! @brief operator -- (prefix, e.g. --a)
      //! @return a reference to the iterator after decrementing
      Generic &operator --()
      {
        --*m_Iterator;
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! Like native iterators, Generic iterators cannot be written in any meaningful way

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        BCL_Exit( "Cannot read a or write a" + this->GetClassIdentifier(), -1);
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        BCL_Exit( "Cannot read or write a " + this->GetClassIdentifier(), -1);
        return OSTREAM;
      }
    };

  } // namespace iterate
} // namespace bcl

#endif // BCL_ITERATE_GENERIC_H_
