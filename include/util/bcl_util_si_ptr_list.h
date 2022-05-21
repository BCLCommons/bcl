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

#ifndef BCL_UTIL_SI_PTR_LIST_H_
#define BCL_UTIL_SI_PTR_LIST_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_si_ptr.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SiPtrList
    //! @brief The SiPtrList class is a list of simple pointers.
    //!
    //! @see @link example_util_si_ptr_list.cpp @endlink
    //! @author alexanns
    //! @date 01/21/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class SiPtrList :
      public storage::List< SiPtr< t_DataType> >
    {

    private:

    //////////
    // data //
    //////////

    public:

      typedef typename storage::List< SiPtr< t_DataType> >::iterator iterator;
      typedef typename storage::List< SiPtr< t_DataType> >::const_iterator const_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct a SiPtrList from optional size and single element
      //! @param SIZE optional size, default 0
      //! @param DATA optional single element which is assigned to all stored elements, default is empty SiPtr
      SiPtrList( const size_t SIZE = 0, const SiPtr< t_DataType> &DATA = SiPtr< t_DataType>()) :
        storage::List< SiPtr< t_DataType> >( SIZE, DATA)
      {
      }

      //! @brief construct a SiPtrList from number of size and pointer to data
      //! @param SIZE number of elements behind the pointer *DATA
      //! @param DATA pointer to the elements
      SiPtrList( const size_t SIZE, t_DataType *const DATA) :
        storage::List< SiPtr< t_DataType> >( SIZE)
      {
        // create a a pointer to the DATA
        t_DataType *ptr_data = DATA;

        // iterate over the DATA and insert it into the SiPtrList
        for
        (
          typename storage::List< SiPtr< t_DataType> >::iterator
            itr_begin( storage::List< SiPtr< t_DataType> >::Begin()),
            itr_end( storage::List< SiPtr< t_DataType> >::End());
          itr_begin != itr_end;
          ++itr_begin
        )
        {
          ( *itr_begin) = SiPtr< t_DataType>( ptr_data++);
        }
      }

      //! @brief copy constructor construct a SiPtrList from another SiPtrList
      //! @param SI_PTR_LIST the SiPtrList to copy
      SiPtrList( const SiPtrList< t_DataType> &SI_PTR_LIST) :
        storage::List< SiPtr< t_DataType> >( SI_PTR_LIST)
      {
      }

      //! @brief move constructor construct a SiPtrList from another SiPtrList
      //! @param SI_PTR_LIST the SiPtrList to move
      SiPtrList( SiPtrList< t_DataType> && SI_PTR_LIST) :
        storage::List< SiPtr< t_DataType> >( std::move( SI_PTR_LIST))
      {
      }

      //! @brief construct a SiPtrList from iterator [FIRST, LAST) range
      //! @tparam t_Iterator which is the type of iterator that indicates the range
      //! @param FIRST t_Iterator to the first element to be copied
      //! @param LAST t_Iterator to the first element after the last copied element
      template< typename t_Iterator>
      SiPtrList( const t_Iterator &FIRST, const t_Iterator &LAST) :
        storage::List< SiPtr< t_DataType> >( FIRST, LAST)
      {
      }

      //! @brief virtual copy constructor
      //! @return pointer to a copy of the actual object
      SiPtrList< t_DataType> *Clone() const
      {
        return new SiPtrList< t_DataType>( *this);
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

    ////////////////
    // operations //
    ////////////////

      using storage::List< SiPtr< t_DataType> >::PushFront;

      //! @brief PushFront adds an element to the beginning of the SiPtrList
      //! @param POINTER is a pointer to t_DataType to be inserted into the SiPtrList
      void PushFront( t_DataType *const POINTER)
      {
        storage::List< SiPtr< t_DataType> >::PushFront( SiPtr< t_DataType>( POINTER));
      }

      using storage::List< SiPtr< t_DataType> >::PushBack;

      //! @brief PushBack adds an element to the end of the SiPtrList
      //! @param POINTER is a pointer to t_DataType which will be added to the end of the list
      void PushBack( t_DataType *const POINTER)
      {
        storage::List< SiPtr< t_DataType> >::PushBack( SiPtr< t_DataType>( POINTER));
      }

      //! checks whether at at least one position ShPtr is undefined
      bool IsDefined() const
      {
        //check each individual pointer if it is defined
        for
        (
          typename storage::List< SiPtr< t_DataType> >::const_iterator
            itr( storage::List< SiPtr< t_DataType> >::Begin()),
            itr_end( storage::List< SiPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          //return false if one individual pointer is not defined
          if( !itr->IsDefined())
          {
            return false;
          }
        }

        //return true if every pointer is defined
        return true;
      }

      //! @brief assign operator for assigning given SiPtrVector to this one
      //! @param SI_PTR_VECTOR SiPtrVector to be copied
      //! @return pointer to this SiPtrVector after being assigned to given SiPtrVector
      SiPtrList &operator =( const SiPtrList &SI_PTR_LIST)
      {
        // class base class equal
        storage::List< SiPtr< t_DataType> >::operator =( SI_PTR_LIST);

        // return
        return *this;
      }

      //! @brief move assign operator for assigning given SiPtrVector to this one
      //! @param SI_PTR_VECTOR SiPtrVector to be copied
      //! @return pointer to this SiPtrVector after being assigned to given SiPtrVector
      SiPtrList< t_DataType> &operator =( SiPtrList && SI_PTR_LIST)
      {
        // class base class equal
        storage::List< SiPtr< t_DataType> >::operator =( std::move( SI_PTR_LIST));

        // return
        return *this;
      }

    }; // template class SiPtrList

    //! @brief convert data of t_DataType denoted by two iterators [FIRST, LAST) into a SiPtrVector
    //! @tparam t_Iterator which is the type of iterator that indicates the range
    //! @param FIRST t_Iterator to the first element to be copied
    //! @param LAST t_Iterator to the first element after the last copied element
    template< typename t_DataType, typename t_Iterator>
    SiPtrList< t_DataType> ConvertToSiPtrList( t_Iterator FIRST, t_Iterator LAST)
    {
      SiPtrList< t_DataType> si_ptr_list;
      while( FIRST != LAST)
      {
        si_ptr_list.PushBack( &( *FIRST));
        ++FIRST;
      }

      return si_ptr_list;
    }

    //! @brief return a SiPtrList< const t_DataType> constructed from an iterator rage od t_DataType
    //! @param FIRST t_Iterator to the first element to be pointed to
    //! @param LAST t_Iterator to the first element after the last copied element
    //! @return SiPtrList< const t_DataType> constructed from range
    template< typename t_DataType, typename t_Iterator>
    SiPtrList< const t_DataType>
    ConvertToConstSiPtrList
    (
      t_Iterator FIRST,
      t_Iterator LAST
    )
    {
      // initialize new SiPtrList
      SiPtrList< const t_DataType> new_siptrlist;

      // iterate over all elements
      for( ; FIRST != LAST; ++FIRST)
      {
        // construct a SiPtr to this element and insert it into the new SiPtrList
        new_siptrlist.PushBack( SiPtr< const t_DataType>( &*FIRST));
      }

      // end
      return new_siptrlist;
    }

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_SI_PTR_LIST_H_
