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

#ifndef BCL_ITERATE_WITH_RANGE_H_
#define BCL_ITERATE_WITH_RANGE_H_

// include the namespace header
#include "bcl_iterate.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_iterate_interface.h"
#include "type/bcl_type_enable_if.h"
#include "type/bcl_type_is_a.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace iterate
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WithRange
    //!
    //! @brief an iterator with a built-in range and convenience functions, such as moving to a random position,
    //!        and getting the current position
    //!
    //! uses a normal iterator, which can be retrieved by calling GetNativeIterator() for use in e.g. Remove() functions
    //! in most containers
    //!
    //! @tparam t_Iterator the type of iterator to use
    //! @tparam t_Return const-qualified type of object returned by operator * on this iterator;
    //!         t_Return defaults to the type returned by *t_Iterator
    //!
    //! @author mendenjl
    //! @see @link example_iterate_with_range.cpp @endlink
    //! @date Jul 14, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Iterator, typename t_Return>
    class WithRange :
      public Interface< t_Return>
    {
    private:

    //////////
    // data //
    //////////

      t_Iterator m_Begin; //!< the beginning of the range iterated upon
      t_Iterator m_Itr;   //!< current position of the iterator
      t_Iterator m_End;   //!< the end of the range being iterated upon

    public:

      //! Value, pointer, and reference types of this iterator
      //! These typedefs enable one to use generic iterators with std:: algorithms, functions, and other
      //! classes in this namespace
      //! All std iterators have these typedefs

      //! non-const, non-reference type (appropriate for storing the result of dereferencing the iterator)
      typedef typename type::RemoveConst< t_Return>::Type value_type;

      //! type returned by operator->
      typedef t_Return *pointer;

      //! type returned by operator*
      typedef t_Return &reference;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from an iterator range
      WithRange( t_Iterator const &BEGIN, t_Iterator const &END) :
        m_Begin( BEGIN),
        m_Itr( BEGIN),
        m_End( END)
      {
      }

      //! @brief constructor from an iterator range and iterator, usually within that range
      WithRange
      (
        t_Iterator const &BEGIN,
        t_Iterator const &END,
        t_Iterator const &ITR
      ) :
        m_Begin( BEGIN),
        m_Itr( ITR),
        m_End( END)
      {
      }

      //! @brief Clone function
      //! @return pointer to new WithRange
      WithRange *Clone() const
      {
        return new WithRange( *this);
      }

      //! @brief get a hard copy of the interface object
      Interface< t_Return> *CloneIterator() const
      {
        return Clone();
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

      //! @brief GetBegin returns the beginning
      //! @return a constant iterator to the position after the last element of the range
      t_Iterator const &GetBegin() const
      {
        return m_Begin;
      }

      //! @brief GetEnd returns the end
      //! @return a constant iterator to the position after the last element of the range
      t_Iterator const &GetEnd() const
      {
        return m_End;
      }

      //! @return the primary iterator
      t_Iterator const &GetNativeIterator() const
      {
        return m_Itr;
      }

      //! @brief SetRange set the range of this iterator
      //! @param BEGIN iterator to the first element of the range
      //! @param END iterator to the last element of the range
      void SetRange( const t_Iterator &BEGIN, const t_Iterator &END)
      {
        m_Itr = m_Begin = BEGIN;
        m_End = END;
      }

      //! @brief determine whether the iterator is at the end
      //! @return true iff the iterator is not at the end
      bool NotAtEnd() const
      {
        return m_Itr != m_End;
      }

      //! @brief determine whether the iterator is at the beginning
      //! @return true iff the iterator is at the beginning
      bool AtBegin() const
      {
        return m_Itr == m_Begin;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset the underlying iterator to the start of its range
      void Restart()
      {
        m_Itr = m_Begin;
      }

      //! @brief reset the underlying iterator to the start of its range
      WithRange &GotoBegin()
      {
        m_Itr = m_Begin;
        return *this;
      }

      //! @brief moves the iterator to the position provided or the end of the range, whichever comes first
      //! @param POSITION the position to go to
      //! @return whether the iterator is still valid
      bool GotoPosition( const size_t &POSITION)
      {
        long difference( long( POSITION) - GetPosition());
        if( difference > 0)
        {
          // advance, but not past the end of the container
          std::advance( m_Itr, std::min( long( GetReversePosition()), difference));
        }
        else if( difference < 0)
        {
          // go backwards; no danger of going before the beginning since POSITION is unsigned
          std::advance( m_Itr, difference);
        }
        return NotAtEnd();
      }

      //! @brief reset the underlying iterator to the end of its range
      WithRange &GotoEnd()
      {
        m_Itr = m_End;
        return *this;
      }

      //! @brief determine how many elements are between the iterator and the first element of the range
      //! @return a long indicating the iterators distance from Begin()
      //! @note return type has to be a long because iterator could be before Begin()
      //! @note O(1) for vectors, O( size of range) for everything else
      long GetPosition() const
      {
        return std::distance( m_Begin, m_Itr); // determine how far iterator is from the beginning
      }

      //! @brief determine how many elements are between the iterator and the last element of the range
      //! @return a long indicating the iterators distance from the last element ()
      //! @note return type has to be a long because iterator could be before End()
      //! @note O(1) for vectors, O( size of range) for everything else
      long GetReversePosition() const
      {
        return std::distance( m_Itr, m_End);
      }

      //! @return the size of the range
      //! @note O(1) for vectors, O( size of range) for everything else
      size_t GetSize() const
      {
        return std::distance( m_Begin, m_End);
      }

    ///////////////
    // operators //
    ///////////////

      //! @return the result of the view applied to the dereferenced iterator
      reference operator *() const
      {
        return Dereference< t_Iterator>();
      }

      //! @return the result of the view applied to the dereferenced iterator
      pointer operator ->() const
      {
        return &Dereference< t_Iterator>();
      }

      //! @brief Const ref access to the iterator
      //! @return the internal iterator
      operator t_Iterator const &() const
      {
        return m_Itr;
      }

      //! @brief assignment operator
      WithRange &operator =( const WithRange &ORIGINAL)
      {
        m_Begin = ORIGINAL.m_Begin;
        m_Itr = ORIGINAL.m_Itr;
        m_End = ORIGINAL.m_End;
        return *this;
      }

      //! @brief assignment operator
      WithRange &operator =( const t_Iterator &ITR)
      {
        m_Itr = ITR;
        return *this;
      }

      //! @brief operator ++ (prefix, e.g. ++a)
      //! @return a reference to the iterator after incrementing
      WithRange &operator ++()
      {
        ++m_Itr;
        return *this;
      }

      //! @brief operator ++ (postfix, e.g. a++)
      //! @param dummy parameter to make this the postfix ++
      //! @return a copy of the iterator before operator++ was applied
      WithRange operator ++( int)
      {
        WithRange old( *this);
        ++m_Itr;
        return old;
      }

      //! @brief operator -- (prefix, e.g. --a)
      //! @return a reference to the iterator after decrementing
      WithRange &operator --()
      {
        --m_Itr;
        return *this;
      }

      //! @brief operator -- (postfix, e.g. a--)
      //! @param dummy parameter to make this the postfix --
      //! @return a copy of the interator before operator-- was applied
      WithRange operator --( int)
      {
        WithRange old( *this);
        --m_Itr;
        return old;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return an iterator to the beginning of the range
      //! @return a generic iterator base to be used for constructing a new generic iterator
      Interface< t_Return> *Begin() const
      {
        return new WithRange( m_Begin, m_End);
      }

      //! @brief return an iterator to the end of the range
      //! @return a generic iterator base to be used for constructing a new generic iterator
      Interface< t_Return> *End() const
      {
        return new WithRange( m_Begin, m_End, m_End);
      }

      //! @brief a special function to dereference the iterator
      //! @return the reference type
      //! @note specialization for when dereferencing the iterator directly returns the correct type
      template< typename t_IteratorOther>
      typename type::EnableIf
      <
        type::IsA< typename std::iterator_traits< t_IteratorOther>::reference, value_type>::value,
        reference
      >::Type
      Dereference() const
      {
        return *m_Itr;
      }

      //! @brief a special function to dereference the iterator
      //! @return the reference type
      //! @note specialization for when dereferencing the iterator and the underlying object returns the correct type
      template< typename t_IteratorOther>
      typename type::EnableIf
      <
        !type::IsA< typename std::iterator_traits< t_IteratorOther>::reference, value_type>::value,
        reference
      >::Type
      Dereference() const
      {
        return **m_Itr;
      }
    };

  } // namespace iterate
} // namespace bcl

#endif // BCL_ITERATE_WITH_RANGE_H_

