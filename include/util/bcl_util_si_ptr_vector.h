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

#ifndef BCL_UTIL_SI_PTR_VECTOR_H_
#define BCL_UTIL_SI_PTR_VECTOR_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_si_ptr.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SiPtrVector
    //! @brief This is the SiPtrVector class
    //!
    //! @see @link example_util_si_ptr_vector.cpp @endlink
    //! @author woetzen, meilerj, karakam
    //! @date 23.09.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class SiPtrVector :
      public storage::Vector< SiPtr< t_DataType> >
    {

    public:

    //////////
    // data //
    //////////

      //! typedef for iterator
      typedef typename storage::Vector< SiPtr< t_DataType> >::iterator iterator;
      //! typedef for const_iterator
      typedef typename storage::Vector< SiPtr< t_DataType> >::const_iterator const_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct SiPtrVector from optional SIZE will be filled with default t_DataTypes
      //! @param SIZE requested size of vector
      SiPtrVector( const size_t SIZE = 0) :
        storage::Vector< SiPtr< t_DataType> >( SIZE)
      {
      }

      //! @brief construct a SiPtrVector of existing objects of length SIZE from a pointer on a data field
      //! @param SIZE requested size of vector
      //! @param DATA pointer to beginning of data field that contains objects to be copied
      SiPtrVector( const size_t SIZE, t_DataType *const DATA) :
        storage::Vector< SiPtr< t_DataType> >( SIZE)
      {
        // initialize pointer to first element of DATA
        t_DataType *data = DATA;

        // iterate over given DATA and this SiPtrVector
        for
        (
          typename std::vector< SiPtr< t_DataType> >::iterator
            itr( storage::Vector< SiPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< SiPtr< t_DataType> >::End());
          itr != itr_end; ++itr
        )
        {
          // copy the pointer to element
          ( *itr) = SiPtr< t_DataType>( data++);
        }
      }

      //! @brief copy constructor
      //! @param SI_PTR_VECTOR_RHS SiPtrVector to be copied
      SiPtrVector( const SiPtrVector< t_DataType> &SI_PTR_VECTOR_RHS) :
        storage::Vector< SiPtr< t_DataType> >( SI_PTR_VECTOR_RHS)
      {
      }

      //! @brief move constructor
      //! @param SI_PTR_VECTOR_RHS SiPtrVector to be moved
      SiPtrVector( SiPtrVector< t_DataType> && SI_PTR_VECTOR_RHS) :
        storage::Vector< SiPtr< t_DataType> >( std::move( SI_PTR_VECTOR_RHS))
      {
      }

      //! @brief construct SiPtrVector from iterator FIRST and LAST range
      //! @param FIRST iterator to beginning of the range to be copied
      //! @param LAST iterator to the end of the range to be copied
      template< typename II>
      SiPtrVector( const II &FIRST, const II &LAST) :
        storage::Vector< SiPtr< t_DataType> >( FIRST, LAST)
      {
      }

      //! @brief construct SiPtrVector from SiPtrVector of another data type
      //! @param SI_PTR_VECTOR_RHS SiPtrVector of type t_OtherDataType to be copied
      template< typename t_OtherDataType>
      SiPtrVector( const SiPtrVector< t_OtherDataType> &SI_PTR_VECTOR_RHS) :
        storage::Vector< SiPtr< t_DataType> >( SI_PTR_VECTOR_RHS.GetSize())
      {
        // initialize iterators to the beginning and end of this SiPtrVector
        typename std::vector< SiPtr< t_DataType> >::iterator
          itr1( storage::Vector< SiPtr< t_DataType> >::Begin()),
          itr1_end( storage::Vector< SiPtr< t_DataType> >::End());

        // initialize iterators to the beginning and end of given SI_PTR_VECTOR_RHS
        typename std::vector< SiPtr< t_OtherDataType> >::const_iterator
          itr2( SI_PTR_VECTOR_RHS.Begin()),
          itr2_end( SI_PTR_VECTOR_RHS.End());

        // iterate over this and the given SiPtrVector
        for( ; itr1 != itr1_end && itr2 != itr2_end; ++itr1, ++itr2)
        {
          // copy elements over
          ( *( itr1)) = ( *( itr2));
        }
      }

      //! @brief construct SiPtrVector from SiPtrVector of another data type
      //! @param SI_PTR_VECTOR_RHS SiPtrVector of type t_OtherDataType to be copied
      template< typename t_OtherDataType>
      SiPtrVector( SiPtrVector< t_OtherDataType> &SI_PTR_VECTOR_RHS) :
        storage::Vector< SiPtr< t_DataType> >( SI_PTR_VECTOR_RHS.GetSize())
      {
        // initialize iterators to the beginning and end of this SiPtrVector
        typename std::vector< SiPtr< t_DataType> >::iterator
          itr1( storage::Vector< SiPtr< t_DataType> >::Begin()),
          itr1_end( storage::Vector< SiPtr< t_DataType> >::End());

        // initialize iterators to the beginning and end of given SI_PTR_VECTOR_RHS
        typename std::vector< SiPtr< t_OtherDataType> >::iterator
          itr2( SI_PTR_VECTOR_RHS.Begin()),
          itr2_end( SI_PTR_VECTOR_RHS.End());

        // iterate over this and the given SiPtrVector
        for( ; itr1 != itr1_end && itr2 != itr2_end; ++itr1, ++itr2)
        {
          // copy elements over
          ( *( itr1)) = ( *( itr2));
        }
      }

      //! @brief construct SiPtrVector from a single SiPtr
      //! @param SIMPLE_POINTER SiPtr to be copied
      SiPtrVector( const SiPtr< t_DataType> &SIMPLE_POINTER) :
         storage::Vector< SiPtr< t_DataType> >( 1, SIMPLE_POINTER)
      {
      }

      //! @brief construct SiPtrVector from two SiPtrs
      //! @brief SIMPLE_POINTER_A first SiPtr to be copied
      //! @brief SIMPLE_POINTER_B second SiPtr to be copied
      SiPtrVector( const SiPtr< t_DataType> &SIMPLE_POINTER_A, const SiPtr< t_DataType> &SIMPLE_POINTER_B) :
         storage::Vector< SiPtr< t_DataType> >( 2)
      {
        storage::Vector< SiPtr< t_DataType> >::operator()( 0) = SIMPLE_POINTER_A;
        storage::Vector< SiPtr< t_DataType> >::operator()( 1) = SIMPLE_POINTER_B;
      }

      //! @brief copy constructor
      //! @return pointer to a new SiPtrVector copied from this one
      SiPtrVector< t_DataType> *Clone() const
      {
        return new SiPtrVector< t_DataType>( *this);
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

      //! @brief inserts SiPtr to element to the end
      //! @param SIMPLE_POINTER SiPtr to element to be inserted
      void PushBack( const SiPtr< t_DataType> &SIMPLE_POINTER)
      {
        storage::Vector< SiPtr< t_DataType> >::PushBack( SIMPLE_POINTER);
      }

      //! @brief inserts element to the end
      //! @param POINTER pointer to elements to be inserted
      void PushBack( t_DataType *const POINTER)
      {
        storage::Vector< SiPtr< t_DataType> >::PushBack( SiPtr< t_DataType>( POINTER));
      }

      //! @brief return sub SiPtrVector from given position POS of length LENGTH
      //! @param POS starting position of the subvector
      //! @param LENGTH length of the subvector
      //! @return requested sub SiPtrVector
      SiPtrVector< t_DataType> SubSiPtrVector( const size_t POS, const size_t LENGTH) const
      {
        return SiPtrVector< t_DataType>
        (
          storage::Vector< SiPtr< t_DataType> >::Begin() + POS,
          storage::Vector< SiPtr< t_DataType> >::Begin() + POS + LENGTH
        );
      }

      //! @brief return sub SiPtrVector specified by a pair of beginning and end position
      //! @param POS pair of beginning and end positions
      //! @return requested sub SiPtrVector
      SiPtrVector< t_DataType> SubSiPtrVector( const std::pair< size_t, size_t> &POS) const
      {
        return SubSiPtrVector( POS.first, ( POS.second - POS.first + 1));
      }

      //! @brief copies SiPtr to object and inserts it at position POS
      //! @param POS position at which the new elements is going to be inserted
      //! @param SIMPLE_POINTER pointer to element to be inserted
      void InsertElement( const size_t POS, const SiPtr< t_DataType> SIMPLE_POINTER)
      {
        storage::Vector< SiPtr< t_DataType> >::InsertElements( POS, SIMPLE_POINTER);
      }

      //! @brief copies given element and inserts it at position POS
      //! @param POS position at which the new elements is going to be inserted
      //! @param POINTER pointer to element to be inserted
      void InsertElement( const size_t POS, t_DataType *const POINTER)
      {
        InsertElement( POS, SiPtr< t_DataType>( POINTER));
      }

      //! @brief Given a vector and an element to look up, return index of the element
      //! @param ELEMENT The element we're looking up
      //! @return the position of the element, or an invalid position if nonexistent
      size_t GetIndex
      (
        const SiPtr< t_DataType> &ELEMENT
      )
      {
        // loop through the vector and compare addresses
        for
        (
          typename std::vector< SiPtr< t_DataType> >::iterator
            itr( storage::Vector< SiPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< SiPtr< t_DataType> >::End());
          itr != itr_end; ++itr
        )
        {
          // need to use addresses, since some types may not have operator == defined
          if( &( **itr) == &( *ELEMENT))
            return ( size_t)( itr - storage::Vector< SiPtr< t_DataType> >::Begin());
        }
        // didn't find it, so return the size which is an invalid index
        return ( size_t)( storage::Vector< SiPtr< t_DataType> >::GetSize());
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief assign operator for assigning given SiPtrVector to this one
      //! @param SI_PTR_VECTOR SiPtrVector to be copied
      //! @return pointer to this SiPtrVector after being assigned to given SiPtrVector
      SiPtrVector< t_DataType> &operator =( const SiPtrVector< t_DataType> &SI_PTR_VECTOR)
      {
        // class base class equal
        storage::Vector< SiPtr< t_DataType> >::operator =( SI_PTR_VECTOR);

        // return
        return *this;
      }

      //! @brief move assign operator for assigning given SiPtrVector to this one
      //! @param SI_PTR_VECTOR SiPtrVector to be copied
      //! @return pointer to this SiPtrVector after being assigned to given SiPtrVector
      SiPtrVector< t_DataType> &operator =( SiPtrVector< t_DataType> && SI_PTR_VECTOR)
      {
        // class base class equal
        storage::Vector< SiPtr< t_DataType> >::operator =( std::move( SI_PTR_VECTOR));

        // return
        return *this;
      }

      //! @brief check that all SiPtrs in this SiPtrVector are defined (!= NULL)
      //! @return whether all SiPtrs in this SiPtrVector are defined (!= NULL)
      bool IsDefined() const
      {
        // iterate over elements
        for
        (
          typename std::vector< SiPtr< t_DataType> >::const_iterator
            itr( storage::Vector< SiPtr< t_DataType> >::Begin()),
            itr_end( storage::Vector< SiPtr< t_DataType> >::End());
          itr != itr_end; ++itr
        )
        {
          // if this pointer is not defined return false
          if( !( itr->IsDefined()))
          {
            return false;
          }
        }
        // otherwise all SiPtrs are valid, therefore return true
        return true;
       };

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
             itr( storage::Vector< SiPtr< t_DataType> >::Begin()),
             itr_end( storage::Vector< SiPtr< t_DataType> >::End());
           itr != itr_end;
           ++itr
         )
         {
           // need to use addresses, since some types may not have operator == defined
           if( &( **itr) == ELEMENT)
           {
             return ( size_t)( itr - storage::Vector< SiPtr< t_DataType> >::Begin());
           }
         }

         // didn't find it, so return the size which is an invalid index
         return ( size_t)( storage::Vector< SiPtr< t_DataType> >::GetSize());
       }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read SiPtrVector from std::istream
      //! @param ISTREAM input stream
      //! @return istream which SiPtrVector was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write SiPtrVector to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which SiPtrVector was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    ////////////////
    // operations //
    ////////////////
    public:

      //! @brief creates the simple pointer vector from two elements
      //! @param DATA_A first data
      //! @param DATA_B second data
      //! @return SiPtrVector constructed
      static SiPtrVector< t_DataType> Create
      (
        t_DataType &DATA_A,
        t_DataType &DATA_B
      )
      {
        // create vector from two elements
        SiPtrVector< t_DataType> vector( 2, &DATA_A);
        vector( 1) = &DATA_B;

        // return vector
        return vector;
      }

      //! @brief creates the simple pointer vector from three elements
      //! @param DATA_A first data
      //! @param DATA_B second data
      //! @param DATA_C third data
      //! @return SiPtrVector constructed
      static SiPtrVector< t_DataType> Create
      (
        t_DataType &DATA_A,
        t_DataType &DATA_B,
        t_DataType &DATA_C
      )
      {
        // create vector from three elements
        SiPtrVector< t_DataType> vector( 3, &DATA_A);
        vector( 1) = &DATA_B;
        vector( 2) = &DATA_C;

        // return vector
        return vector;
      }

      //! @brief creates the simple pointer vector from four elements
      //! @param DATA_A first data
      //! @param DATA_B second data
      //! @param DATA_C third data
      //! @param DATA_D fourth data
      //! @return SiPtrVector constructed
      static SiPtrVector< t_DataType> Create
      (
        t_DataType &DATA_A,
        t_DataType &DATA_B,
        t_DataType &DATA_C,
        t_DataType &DATA_D
      )
      {
        // create vector from four elements
        SiPtrVector< t_DataType> vector( 4, &DATA_A);
        vector( 1) = &DATA_B;
        vector( 2) = &DATA_C;
        vector( 3) = &DATA_D;

        // return vector
        return vector;
      }

      //! @brief creates the simple pointer vector from five elements
      //! @param DATA_A first data
      //! @param DATA_B second data
      //! @param DATA_C third data
      //! @param DATA_D fourth data
      //! @param DATA_E fifth data
      //! @return SiPtrVector constructed
      static SiPtrVector< t_DataType> Create
      (
        t_DataType &DATA_A,
        t_DataType &DATA_B,
        t_DataType &DATA_C,
        t_DataType &DATA_D,
        t_DataType &DATA_E
      )
      {
        // create vector from five elements
        SiPtrVector< t_DataType> vector( 5, &DATA_A);
        vector( 1) = &DATA_B;
        vector( 2) = &DATA_C;
        vector( 3) = &DATA_D;
        vector( 4) = &DATA_E;

        // return vector
        return vector;
      }

      //! @brief creates the simple pointer vector from six elements
      //! @param DATA_A first data
      //! @param DATA_B second data
      //! @param DATA_C third data
      //! @param DATA_D fourth data
      //! @param DATA_E fifth data
      //! @param DATA_F sixth data
      //! @return SiPtrVector constructed
      static SiPtrVector< t_DataType> Create
      (
        t_DataType &DATA_A,
        t_DataType &DATA_B,
        t_DataType &DATA_C,
        t_DataType &DATA_D,
        t_DataType &DATA_E,
        t_DataType &DATA_F
      )
      {
        // create vector from six elements
        SiPtrVector< t_DataType> vector( 6, &DATA_A);
        vector( 1) = &DATA_B;
        vector( 2) = &DATA_C;
        vector( 3) = &DATA_D;
        vector( 4) = &DATA_E;
        vector( 5) = &DATA_F;

        // return vector
        return vector;
      }

    }; // SiPtrVector Class

    //! @brief read SiPtrVector from std::istream
    //! @param ISTREAM input stream
    //! @return istream which SiPtrVector was read from
    template< typename t_DataType>
    std::istream &
    SiPtrVector< t_DataType>::Read( std::istream &ISTREAM)
    {
      return storage::Vector< SiPtr< t_DataType> >::Read( ISTREAM);
    }

    //! @brief write SiPtrVector to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which SiPtrVector was written to
    template< typename t_DataType>
    std::ostream &
    SiPtrVector< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return storage::Vector< SiPtr< t_DataType> >::Write( OSTREAM, INDENT);
    }

    //! this is a dummy *  operator to allow usage of SumFunction class with std::pair< AA, AA>
    template< typename t_DataType, class T2> inline SiPtrVector< t_DataType> operator *( const T2 &VALUE, const SiPtrVector< t_DataType> &SIMPLEPOINTERVECTOR)
    {
      BCL_Exit( "This dummy function should never be called!", -1);
      return SIMPLEPOINTERVECTOR;
    }

    //! this is a dummy += operator to allow usage of SumFunction class with std::pair< AA, AA>
    template< typename t_DataType> inline SiPtrVector< t_DataType> operator +=( const SiPtrVector< t_DataType> &SIMPLEPOINTERVECTOR_A, const SiPtrVector< t_DataType> &SIMPLEPOINTERVECTOR_B)
    {
      BCL_Exit( "This dummy function should never be called!", -1);
      return SIMPLEPOINTERVECTOR_A;
    }

    //! @brief construct SiPtrVector from storage::Vector
    //! @param STORAGE_VECTOR storage::Vector from which SiPtrVector is going to be constructed
    //! @return SiPtrVector constructor from given storage::VECTOR
    template< typename t_DataType>
    SiPtrVector< t_DataType>
    ConvertToSiPtrVector
    (
      storage::Vector< t_DataType> &STORAGE_VECTOR
    )
    {
      // initialize a new SiPtrVector and allocate the memory
      SiPtrVector< t_DataType> new_siptrvector;
      new_siptrvector.AllocateMemory( STORAGE_VECTOR.GetSize());

      // iterate over the STORAGE_VECTOR
      for
      (
        typename std::vector< t_DataType>::iterator
          arg_itr( STORAGE_VECTOR.Begin()), arg_itr_end( STORAGE_VECTOR.End());
        arg_itr != arg_itr_end;
        ++arg_itr
      )
      {
        // insert a SiPtr to the current element into the new SiPtrVector constructed
        new_siptrvector.PushBack( SiPtr< t_DataType>( *arg_itr));
      }

      // end
      return new_siptrvector;
    }

    //! @brief return a storage::Vector with copies of all objects in the SiPtrVector
    //! @param SIMPLE_POINTER_VECTOR SiPtrVector from which elements will be copied
    //! @return storage::Vector with copies of all objects in the SiPtrVector
    template< typename t_DataType, typename t_OtherDataType>
    storage::Vector< t_DataType>
    ConvertToStorageVector
    (
      const SiPtrVector< t_OtherDataType> &SIMPLE_POINTER_VECTOR
    )
    {
      // initialize new storage::Vector and allocate memory
      storage::Vector< t_DataType> new_storagevector;
      new_storagevector.AllocateMemory( SIMPLE_POINTER_VECTOR.GetSize());

      // iterate over all elements in the SiPtrVector
      for
      (
        typename std::vector< SiPtr< t_OtherDataType> >::const_iterator
          this_itr( SIMPLE_POINTER_VECTOR.Begin()), this_itr_end( SIMPLE_POINTER_VECTOR.End());
        this_itr != this_itr_end;
        ++this_itr
      )
      {
        // make a copy of the element behind this SiPtr and insert it into the new storage::Vector
        new_storagevector.PushBack( t_DataType( **this_itr));
      }

      // end
      return new_storagevector;
    }

    //! @brief return a SiPtrVector< const t_DataType> constructed from storage::Vector< t_DataType>
    //! @param STORAGE_VECTOR storage::Vector from which SiPtrVector< const t_DataType> is going to be created
    //! @return SiPtrVector< const t_DataType> constructed from storage::Vector< t_DataType>
    template< typename t_DataType>
    SiPtrVector< const t_DataType>
    ConvertToConstSiPtrVector
    (
      const storage::Vector< t_DataType> &STORAGE_VECTOR
    )
    {
      // initialize new SiPtrVector and allocate memory
      SiPtrVector< const t_DataType> new_siptrvector;
      new_siptrvector.AllocateMemory( STORAGE_VECTOR.GetSize());

      // iterate over all elements in STORAGE_VECTOR
      for
      (
        typename std::vector< t_DataType>::const_iterator
          stv_itr( STORAGE_VECTOR.Begin()), stv_itr_end( STORAGE_VECTOR.End());
        stv_itr != stv_itr_end;
        ++stv_itr
      )
      {
        // construct a SiPtr to this element and insert it into the new SiPtrVector
        new_siptrvector.PushBack( SiPtr< const t_DataType>( &*stv_itr));
      }

      // end
      return new_siptrvector;
    }

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_SI_PTR_VECTOR_H_
