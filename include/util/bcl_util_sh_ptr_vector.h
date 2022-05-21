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

#ifndef BCL_UTIL_SH_PTR_VECTOR_H_
#define BCL_UTIL_SH_PTR_VECTOR_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_sh_ptr.h"
#include "bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ShPtrVector
    //! @brief This is the ShPtrVector class
    //!
    //! @see @link example_util_sh_ptr_vector.cpp @endlink
    //! @author woetzen, meilerj
    //! @date 03.10.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class ShPtrVector :
      public storage::Vector< ShPtr< t_DataType> >
    {

    public:

      typedef typename storage::Vector< ShPtr< t_DataType> >::iterator iterator;
      typedef typename storage::Vector< ShPtr< t_DataType> >::const_iterator const_iterator;

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ShPtrVector() :
        storage::Vector< ShPtr< t_DataType> >()
      {
      }

      //! construct ShPtrVector from optional LENGTH
      ShPtrVector( const size_t LENGTH) :
        storage::Vector< ShPtr< t_DataType> >( LENGTH)
      {
      }

      //! construct ShPtrVector of LENGTH and ELEMENT
      ShPtrVector( const size_t LENGTH, const t_DataType &ELEMENT) :
        storage::Vector< ShPtr< t_DataType> >( LENGTH)
      {
        Initialize( ELEMENT);
      }

      //! construct ShPtrVector of LENGTH and ShPtr to ELEMENT
      ShPtrVector( const size_t LENGTH, const ShPtr< t_DataType> &SP_ELEMENT) :
        storage::Vector< ShPtr< t_DataType> >( LENGTH, SP_ELEMENT)
      {
      }

      //! construct ShPtrVector of existing Objects from pointer on datafield of N objects
      ShPtrVector( const size_t SIZE, t_DataType *const DATA) :
        storage::Vector< ShPtr< t_DataType> >( SIZE)
      {
        t_DataType *dat = DATA;
        for
        (
          iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          ( *itr) = ShPtr< t_DataType>( ( dat++)->Clone());
        }
      }

      //! @brief construct a ShPtrVector from iterator [FIRST, LAST) range
      //! @tparam t_Iterator which is the type of iterator that indicates the range
      //! @param FIRST t_Iterator to the first element to be copied
      //! @param LAST t_Iterator to the first element after the last copied element
      template< typename t_Iterator>
      ShPtrVector( const t_Iterator &FIRST, const t_Iterator &LAST) :
        storage::Vector< ShPtr< t_DataType> >( FIRST, LAST)
      {
      }

      //! copy constructor
      ShPtrVector( const ShPtrVector< t_DataType> &SHAREDPOINTERVECTOR) :
        storage::Vector< ShPtr< t_DataType> >( SHAREDPOINTERVECTOR)
      {
      }

      //! move constructor
      ShPtrVector( ShPtrVector< t_DataType> && SHAREDPOINTERVECTOR) :
        storage::Vector< ShPtr< t_DataType> >( std::move( SHAREDPOINTERVECTOR))
      {
      }

      //! construct ShPtrVector from storage::Vector of ShPtr but only makes a soft copy
      ShPtrVector( const storage::Vector< ShPtr< t_DataType> > &STORAGEVECTORSHAREDPOINTER) :
        storage::Vector< ShPtr< t_DataType> >( STORAGEVECTORSHAREDPOINTER)
      {
      };

      //! construct ShPtrVector from SharedpointerVector of derived classes from t_DataType
      template< typename t_DerivedDataType>
      ShPtrVector( const ShPtrVector< t_DerivedDataType> &SHAREDPOINTERVECTOR) :
        storage::Vector< ShPtr< t_DataType> >( SHAREDPOINTERVECTOR.GetSize())
      {
        iterator itr1 = storage::Vector< ShPtr< t_DataType> >::Begin();
        typename std::vector< ShPtr< t_DerivedDataType> >::const_iterator itr2 = SHAREDPOINTERVECTOR.Begin();
        for( size_t i( 0); i < storage::Vector< ShPtr< t_DataType> >::GetSize(); ++i, ++itr1, ++itr2)
        {
          *itr1 = ShPtr< t_DataType>( *itr2);
        }
      }

      //! copy constructor
      ShPtrVector< t_DataType> *Clone() const
      {
        return new ShPtrVector< t_DataType>( *this);
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

      //! return sub ShPtrVector from POS of LENGTH
      ShPtrVector< t_DataType> SubShPtrVector( const size_t POS, const size_t LENGTH) const
      {
        return storage::Vector< ShPtr< t_DataType> >( *this, POS, LENGTH);
      }

      //! return sub ShPtrVector from ( POS1) to ( POS2)
      ShPtrVector< t_DataType> SubShPtrVector( const std::pair< size_t, size_t> POS) const
      {
        return SubShPtrVector( POS.first, ( POS.second - POS.first + 1));
      }

      //! Initialize new objects the smart pointers are pointing to and return reference to this
      ShPtrVector< t_DataType> &Initialize( const t_DataType &ELEMENT = t_DataType())
      {
        for
        (
          iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          ( *itr) = ShPtr< t_DataType>( ELEMENT.Clone());
        }

        return *this;
      }

      using storage::Vector< ShPtr< t_DataType> >::InsertElement;

      //! copies ShPtr of element into ShPtrVector at position ( POS) (SoftCopy)
      void InsertElement( const size_t POS, const ShPtr< t_DataType> &SHAREDPOINTER)
      {
        storage::Vector< ShPtr< t_DataType> >::InsertElements( POS, SHAREDPOINTER);
      }

      //! insert ShPtr to copy of the element to ShPtrVector at position ( POS) by default empty Object
      void InsertElement( const size_t POS, const t_DataType &ELEMENT = t_DataType())
      {
        storage::Vector< ShPtr< t_DataType> >::InsertElements( POS, ShPtr< t_DataType>( ELEMENT.Clone()));
      }

      //! @brief Given a vector and an element to look up, return index of the element
      //! @param ELEMENT The element we're looking up
      //! @return the position of the element, or an invalid position if nonexistent
      size_t GetIndex( const ShPtr< t_DataType> &ELEMENT) const
      {
        if( !storage::Vector< ShPtr< t_DataType> >::GetSize())
        {
          return ( size_t)( storage::Vector< ShPtr< t_DataType> >::GetSize());
        }

        // first, see if ELEMENT occupies a memory location in the array
        // If so, take the difference with Begin, and return the index
        if
        (
          &ELEMENT >= &*storage::Vector< ShPtr< t_DataType> >::Begin()
          &&
          &ELEMENT < &storage::Vector< ShPtr< t_DataType> >::LastElement()
        )
        {
          return &ELEMENT - &*storage::Vector< ShPtr< t_DataType> >::Begin();
        }

        // otherwise, loop through the vector and compare addresses
        for
        (
          const_iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          // need to use addresses, since some types may not have operator == defined
          if( &( **itr) == &( *ELEMENT))
          {
            return ( size_t)( itr - storage::Vector< ShPtr< t_DataType> >::Begin());
          }
        }

        // didn't find it, so return the size which is an invalid index
        return ( size_t)( storage::Vector< ShPtr< t_DataType> >::GetSize());
      }

      //! @brief Given a vector and an element to look up, return index of the element
      //! @param ELEMENT The element we're looking up
      //! @return the position of the element, or an invalid position if nonexistent
      //! @note this is an overload of GetIndex( ShPtr) function for the case when we have a pointer
      //! @note rather than a shared object.  This happens when we use SiPtrs to point to objects held by
      //! @note ShPtrs.
      size_t GetIndex( const t_DataType *ELEMENT) const
      {
        // loop through the vector and compare addresses
        for
        (
          const_iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          // need to use addresses, since some types may not have operator == defined
          if( &( **itr) == ELEMENT)
          {
            return ( size_t)( itr - storage::Vector< ShPtr< t_DataType> >::Begin());
          }
        }

        // didn't find it, so return the size which is an invalid index
        return ( size_t)( storage::Vector< ShPtr< t_DataType> >::GetSize());
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief equal operator
      ShPtrVector< t_DataType> &operator =( const ShPtrVector< t_DataType> &SH_PTR_VECTOR)
      {
        // call base class operator =
        storage::Vector< ShPtr< t_DataType> >::operator =( SH_PTR_VECTOR);

        // return
        return *this;
      }

      //! @brief move assignment operator
      ShPtrVector &operator =( ShPtrVector && VECTOR)
      {
        storage::Vector< ShPtr< t_DataType> >::operator =( std::move( VECTOR));
        return *this;
      }

      //! @brief conversion operator to simple-pointer vector
      operator SiPtrVector< t_DataType>()
      {
        SiPtrVector< t_DataType> new_simplepointervector;
        new_simplepointervector.AllocateMemory( storage::Vector< ShPtr< t_DataType> >::GetSize());
        for
        (
          iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          new_simplepointervector.PushBack( SiPtr< t_DataType>( itr->operator->()));
        }

        return new_simplepointervector;
      }

      //! @brief conversion operator to const simple-pointer vector
      operator SiPtrVector< const t_DataType>() const
      {
        SiPtrVector< const t_DataType> new_simplepointervector;
        new_simplepointervector.AllocateMemory( storage::Vector< ShPtr< t_DataType> >::GetSize());
        for
        (
          const_iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          new_simplepointervector.PushBack( SiPtr< const t_DataType>( itr->operator->()));
        }

        return new_simplepointervector;
      }

      //! @brief conversion operator to const simple-pointer vector of alternate type
      template< typename t_OtherDataType>
      operator SiPtrVector< t_OtherDataType>()
      {
        SiPtrVector< t_OtherDataType> new_simplepointervector;
        new_simplepointervector.AllocateMemory( storage::Vector< ShPtr< t_DataType> >::GetSize());
        for
        (
          iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          new_simplepointervector.PushBack( SiPtr< t_OtherDataType>( itr->operator->()));
        }

        return new_simplepointervector;
      }

      //! @brief conversion operator to cost simple-pointer vector of alternate type
      template< typename t_OtherDataType>
      operator SiPtrVector< const t_OtherDataType>() const
      {
        SiPtrVector< const t_OtherDataType> new_simplepointervector;
        new_simplepointervector.AllocateMemory( storage::Vector< ShPtr< t_DataType> >::GetSize());
        for
        (
          const_iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
          itr != itr_end;
          ++itr
        )
        {
          new_simplepointervector.PushBack( SiPtr< const t_OtherDataType>( itr->operator->()));
        }

        return new_simplepointervector;
      }

      //! returns hard copy of object
      ShPtrVector< t_DataType> HardCopy() const
      {
        //instantiate sharedPointervector with hardcopies of each object
        ShPtrVector< t_DataType> hardcopy( storage::Vector< ShPtr< t_DataType> >::GetSize());

        //iterator to the first element of this vector
        const_iterator           itr_source = storage::Vector< ShPtr< t_DataType> >::Begin();
        //iterate over all elements of the hardcopyvecotor and assign a hardcopy to the sharedpointers
        for
        (
          iterator itr_hardcopy( hardcopy.Begin()), itr_hardcopy_end( hardcopy.End());
          itr_hardcopy != itr_hardcopy_end;
          ++itr_hardcopy, ++itr_source
        )
        {
          ( *itr_hardcopy) = ShPtr< t_DataType>( ( *itr_source)->Clone()); //make a HardCopy of the ShPtr behind the Pointer and pass it to the hardcopy
        }

        //return hardcopied ShPtrVector
        return hardcopy;
      }

      //! checks whether at at least one position ShPtr is undefined
      bool IsDefined() const
      {
        //check each individual pointer if it is defined
        for
        (
          const_iterator
            itr( storage::Vector< ShPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< ShPtr< t_DataType> >::End());
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

    }; // templat class ShPtrVector

    // instantiate s_Instance
    template< typename t_DataType>
    const SiPtr< const ObjectInterface> ShPtrVector< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new ShPtrVector< t_DataType>())
    );

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_SH_PTR_VECTOR_H_

