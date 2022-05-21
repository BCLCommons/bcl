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

#ifndef BCL_ITERATE_INTERFACE_H_
#define BCL_ITERATE_INTERFACE_H_

// include the namespace header
#include "bcl_iterate.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace iterate
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Interface
    //! @brief Abstract base class for iterators.  Defines the Interface and basic functions used by the Interface
    //!
    //! @tparam t_Return the type of object returned by operator*(), skipping the reference
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Jul 14, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Return>
    class Interface :
      public util::ObjectInterface
    {
    public:

    /////////////
    // friends //
    /////////////

      friend class Generic< t_Return>; //!< generic iterator needs to be a friend to call Begin()/End()

    //////////////
    // typedefs //
    //////////////

      // These typedefs are common for all iterators
      typedef long      difference_type; //!< type returned by std::distance for this iterator
      typedef t_Return& reference;       //!< type returned by operator*
      typedef t_Return  value_type;      //!< type needed to hold the result of an operator*
      typedef t_Return* pointer;         //!< type returned by operator->
      typedef typename std::bidirectional_iterator_tag iterator_category; //!< Type of iterator

    /////////////////
    // data access //
    /////////////////

      //! Virtual clone
      virtual Interface *Clone() const = 0;

      //! @brief get a hard copy of the interface object
      virtual Interface *CloneIterator() const = 0;

      //! @brief determine how many elements are between the iterator and the first element of the range
      //! @return a long indicating the iterators distance from Begin()
      //! @note return type has to be a long because iterator could be before Begin()
      //! @note O(1) for vectors, O( size of range) for everything else
      virtual long GetPosition() const = 0;

      //! @brief determine how many elements are between the iterator and the last element of the range
      //! @return a long indicating the iterators distance from the last element ()
      //! @note return type has to be a long because iterator could be before End()
      //! @note O(1) for vectors, O( size of range) for everything else
      virtual long GetReversePosition() const = 0;

      //! @return the size of the range
      //! @note O(1) for vectors, O( size of range) for everything else
      virtual size_t GetSize() const = 0;

      //! @brief test whether this iterator is at the end of its range
      //! @return true iff the iterator is not at the end
      virtual bool NotAtEnd() const = 0;

      //! @brief determine whether the iterator is at the beginning
      //! @return true iff the iterator is at the beginning
      virtual bool AtBegin() const = 0;

    /////////////////////////
    // abstract operations //
    /////////////////////////

      //! @brief reset the underlying iterator to the start of its range (also resets displacement)
      virtual void Restart() = 0;

      //! @brief reset the underlying iterator to the start of its range (also resets displacement)
      virtual Interface &GotoBegin() = 0;

      //! @brief reset the underlying iterator to the start of its range (also resets displacement)
      virtual Interface &GotoEnd() = 0;

    ////////////////////////
    // abstract operators //
    ////////////////////////

      //! @return the data member
      virtual reference operator *() const = 0;

      //! @return a pointer to the data member
      //! If the return value of operator* is a temporary value,
      //! then the return value will be converted to void * to prevent dereferencing
      virtual pointer operator ->() const = 0;

      //! @brief operator ++ (prefix, e.g. ++a)
      //! @return a reference to the iterator after incrementing
      virtual Interface &operator ++() = 0;

      //! @brief operator -- (prefix, e.g. --a)
      //! @return a reference to the iterator after decrementing
      virtual Interface &operator --() = 0;

      //! @return true if both iterators are pointing to the same chunk in memory
      //! @note yields compile time error if the two iterators have different pointer types
      template< typename t_OtherIterator>
      bool operator ==( t_OtherIterator const &ITR) const
      {
        return operator->() == &*ITR;
      }

      //! @return true if both iterators are pointing to the same chunk in memory
      //! @note yields compile time error if the two iterators have different pointer types
      template< typename t_OtherIterator>
      bool operator !=( t_OtherIterator const &ITR) const
      {
        return operator->() != &*ITR;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief moves the iterator to the position provided or the end of the range, whichever comes first
      //! @param POSITION the position to go to
      //! @return whether the iterator is still valid
      virtual bool GotoPosition( const size_t &POSITION) = 0;

      //! @brief Asserts that the iterator is valid and prints a message if it is not
      void AssertNotAtEnd() const
      {
        BCL_Assert
        (
          NotAtEnd(),
          "Tried to perform operation requiring a valid iterator with an iterator at the end of its range"
        );
      }

      //! @brief iterator forward until an element that satisfies CONDITION is reached
      //! @return true if a later element satisfied the condition
      template< typename t_Condition>
      bool GotoNextElementThatSatisfies( t_Condition CONDITION)
      {
        while( NotAtEnd() && !CONDITION( operator*())) // so long as the condition is false and we haven't reached end
        {
          // increment this
          operator++();
        }

        return NotAtEnd();
      }

      //! @brief iterator forward until an element that satisfies CONDITION is reached
      //! @return true if an earlier element satisfied the condition
      template< typename t_Condition>
      bool GotoPrevElementThatSatisfies( t_Condition CONDITION)
      {
        while( NotAtEnd() && !CONDITION( operator*())) // so long as the condition is false and we haven't reached end
        {
          // decrement this
          operator--();
        }

        return NotAtEnd();
      }

      //! @brief go to a random position in within the range
      //! @return the position that the iterator went to
      size_t GotoRandomPosition()
      {
        const size_t position( random::GetGlobalRandom().Random( GetSize() - 1));
        GotoPosition( position);
        return position;
      }

      //! @param CONDITION the condition to use to decide how whether to count a given element
      //! @return the number of elements in the range that satisfy the condition
      template< typename t_Condition>
      size_t CountElementsThatSatisfy( t_Condition CONDITION) const
      {
        size_t count( 0);

        // make a clone of this object
        util::OwnPtr< Interface< t_Return> > itr( Clone());
        itr->Restart();

        while( itr->GotoNextElementThatSatisfies( CONDITION)) // continue until no more elements satisfy the condition
        {
          ++*itr;
          ++count;
        }

        return count;
      }

      //! @param CONDITION the condition to use to decide how whether to count a given element
      //! @return the number of elements in the range that satisfy the condition
      template< typename t_Condition, typename t_Container>
      void InsertElementsThatSatisfyConditionInto( t_Condition CONDITION, t_Container &CONTAINER) const
      {
        // make a clone of this object
        util::OwnPtr< Interface< t_Return> > itr( Clone());
        itr->Restart();

        while( itr->GotoNextElementThatSatisfies( CONDITION)) // continue until we reach the end
        {
          CONTAINER.Insert( **itr);
          ++*itr;
        }
      }

      //! @param CONDITION the condition to use to decide how whether to count a given element
      //! @return the number of elements in the range that satisfy the condition
      //! Not called GetSize because the size is not known by the iterator
      template< typename t_Condition>
      size_t RandomlyPickElementThatSatisfies( t_Condition CONDITION)
      {
        Restart(); // start off at the start

        //! find out how many objects satisfy the given condition
        const size_t total( CountElementsThatSatisfy( CONDITION));

        //! number of objects satisfying the condition that will preceed the chosen object in the range
        const size_t random_object_number( random::GetGlobalRandom().Random( total));

        size_t position( 0);

        for( size_t count( 0); position < random_object_number; operator++(), ++position, ++count)
        {
          while( !CONDITION( operator*())) // continue the condition is satisfeid
          {
            operator++();
            ++position;
          }
        }

        return position;
      }

      //! @brief add all elements following the current iterator to a certain container
      //! @param BASE a generic iterator base that contains things to iterator over
      //! @param CONTAINER the container to store the objects in
      template< class t_Container>
      void InsertElementsFromIteratorInto( t_Container CONTAINER) const
      {
        // make a clone of this object
        util::OwnPtr< Interface< t_Return> > itr( *this);
        itr->Reset();

        while( itr->NotAtEnd()) // continue until we reach the end
        {
          CONTAINER.Insert( **itr);
          ++*itr;
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! Iterators cannot be serialized, since the bounds of the containers they iterate on are addresses that
      //! are not written

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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return a pointer to a new iterator base
      //! @return a generic iterator base to be used for constructing a new generic iterator
      virtual Interface *Begin() const = 0;

      //! @brief return an iterator to the end of the range
      //! @return a generic iterator base to be used for constructing a new generic iterator
      virtual Interface *End() const = 0;

    };

  } // namespace iterate
} // namespace bcl

#endif // BCL_ITERATE_INTERFACE_H_
