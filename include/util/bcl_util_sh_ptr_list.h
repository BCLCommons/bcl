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

#ifndef BCL_UTIL_SH_PTR_LIST_H_
#define BCL_UTIL_SH_PTR_LIST_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_sh_ptr.h"
#include "bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ShPtrList
    //! @brief The ShPtrList class is a list of shared pointers.
    //!
    //! @see @link example_util_sh_ptr_list.cpp @endlink
    //! @author alexanns
    //! @date 01/28/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class ShPtrList :
      public storage::List< ShPtr< t_DataType> >
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct a ShPtrList from optional size and single element
      //! @param SIZE optional size, default 0
      //! @param DATA optional single element which is assigned to all stored elements, default is empty ShPtr
      ShPtrList( const size_t SIZE = 0, const ShPtr< t_DataType> &DATA = ShPtr< t_DataType>()) :
        storage::List< ShPtr< t_DataType> >( SIZE, DATA)
      {
      }

      //! @brief construct a ShPtrList from number of size and pointer to data
      //! @param SIZE number of elements behind the pointer *DATA
      //! @param DATA pointer to the elements
      ShPtrList( const size_t SIZE, t_DataType *const DATA) :
        storage::List< ShPtr< t_DataType> >( SIZE)
      {
        // create a a pointer to the DATA
        t_DataType *ptr_data = DATA;

        // iterate over the DATA and insert it into the ShPtrList
        for
        (
          typename storage::List< ShPtr< t_DataType> >::iterator
            itr_begin( storage::List< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::List< ShPtr< t_DataType> >::End());
          itr_begin != itr_end;
          ++itr_begin
        )
        {
          ( *itr_begin) = ShPtr< t_DataType>( PtrHardCopy( ptr_data++));
        }
      }

      //! @brief construct a ShPtrList from iterator [FIRST, LAST) range
      //! @tparam t_Iterator which is the type of iterator that indicates the range
      //! @param FIRST t_Iterator to the first element to be copied
      //! @param LAST t_Iterator to the first element after the last copied element
      template< typename t_Iterator>
      ShPtrList( const t_Iterator &FIRST, const t_Iterator &LAST) :
        storage::List< ShPtr< t_DataType> >( FIRST, LAST)
      {
      }

      //! copy constructor
      ShPtrList( const ShPtrList &A) :
        storage::List< ShPtr< t_DataType> >( A)
      {
      }

      //! copy constructor
      ShPtrList( ShPtrList && A) :
        storage::List< ShPtr< t_DataType> >( std::move( A))
      {
      }

      //! @brief virtual copy constructor
      //! @return pointer to a copy of the actual object
      ShPtrList< t_DataType> *Clone() const
      {
        return new ShPtrList< t_DataType>( *this);
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

      //! @brief HardCopy returns a new ShPtrList with new object behind the shared pointers
      //! @return returns a ShPtrList to of new objects
      ShPtrList< t_DataType> HardCopy() const
      {
        // instantiate ShPtrList "hard_copy" to hold new copies of objects
        ShPtrList< t_DataType> hard_copy;

        // iterate over all elements of "hard_copy" and assign a hardcopy to the sharedpointers
        for
        (
          typename storage::List< ShPtr< t_DataType> >::const_iterator itr_source( storage::List< ShPtr< t_DataType> >::Begin()),
            itr_end_source( storage::List< ShPtr< t_DataType> >::End());
          itr_source != itr_end_source;
          ++itr_source
        )
        {
          // make a HardCopy of the ShPtr behind the iterator and pass it to "hard_copy"
          hard_copy.PushBack( itr_source->HardCopy());
        }

        // return hard copied ShPtrList
        return hard_copy;
      }

      //! @brief GetSiPtrList creates a SiPtrList out of the ShPtrList
      //! @return SiPtrList< t_DataType> created from the objects of the ShPtrList
      operator SiPtrList< t_DataType>()
      {
        // create "si_ptr_list" to hold SiPtrs to the t_DataTypes of the ShPtrList
        SiPtrList< t_DataType> si_ptr_list;

        // iterate over the ShPtrList in order to add elements to "si_ptr_list"
        for
        (
          typename storage::List< ShPtr< t_DataType> >::iterator itr_source( storage::List< ShPtr< t_DataType> >::Begin()),
            itr_end_source( storage::List< ShPtr< t_DataType> >::End());
          itr_source != itr_end_source;
          ++itr_source
        )
        {
          // add current element of ShPtr to "si_ptr_list"
          si_ptr_list.PushBack( *itr_source);
        }

        // return "si_ptr_list"
        return si_ptr_list;
      }

      //! @brief equal operator
      ShPtrList &operator =( const ShPtrList &LIST)
      {
        storage::List< ShPtr< t_DataType> >::operator =( LIST);
        return *this;
      }

      //! @brief move assignment operator
      ShPtrList &operator =( ShPtrList && LIST)
      {
        storage::List< ShPtr< t_DataType> >::operator =( std::move( LIST));
        return *this;
      }

      //! checks whether at at least one position ShPtr is undefined
      bool IsDefined() const
      {
        // check each individual pointer if it is defined
        for
        (
          typename storage::List< ShPtr< t_DataType> >::const_iterator itr( storage::List< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::List< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          // return false if one individual pointer is not defined
          if( !itr->IsDefined())
          {
            return false;
          }
        }

        // return true if every pointer is defined
        return true;
      }

    }; // class ShPtrList

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_SH_PTR_LIST_H_
