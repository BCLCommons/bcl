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

#ifndef BCL_LINAL_VECTOR_INTERFACE_HPP_
#define BCL_LINAL_VECTOR_INTERFACE_HPP_

// include the header of this class
#include "bcl_linal_vector_interface.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_const_reference.h"
#include "bcl_linal_vector_reference.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  ////////////////
  // operations //
  ////////////////

    //! @brief normalize the vector (such that inner product == 1)
    template< typename t_DataType>
    VectorInterface< t_DataType> &VectorInterface< t_DataType>::Normalize()
    {
      math::Statistics::Normalize( Begin(), End());
      return *this;
    }

    //! set to sum SUM (sum of all elements = SUM)
    template< typename t_DataType>
    VectorInterface< t_DataType> &VectorInterface< t_DataType>::SetToSum( const t_DataType &SUM)
    {
      math::Statistics::SetToSum( Begin(), End(), SUM);
      return *this;
    }

    //! @brief set all elements in vector to random numbers
    //! @param MIN minimum value to set elements in the vector to
    //! @param MAX maximum value to set elements in the vector to
    template< typename t_DataType>
    VectorInterface< t_DataType> &VectorInterface< t_DataType>::SetRand( const t_DataType MIN, const t_DataType MAX)
    {
      math::Statistics::SetRand( Begin(), End(), MIN, MAX);
      return *this;
    }

    //! @brief copies elements of argument VECTOR into this object at position ( POS)
    //! @param POS beginning of range to replace elements
    //! @param VECTOR vector of elements to replace elements in this vector between POS and POS + VECTOR.GetSize()
    template< typename t_DataType>
    VectorInterface< t_DataType> &VectorInterface< t_DataType>::ReplaceElements
    (
      const size_t POS,
      const VectorConstInterface< t_DataType> &VECTOR
    )
    {
      if( VECTOR.GetSize())
      {
        // check whether the end position is valid; the start position will be validated by operator [] (below)
        BCL_Assert
        (
          POS + VECTOR.GetSize() <= GetSize(),
          "cannot replace elements beyond end of vector " + util::Format()( POS) + " > " + util::Format()( GetSize())
        );

        // copy elements
        std::copy( VECTOR.Begin(), VECTOR.End(), operator[]( POS));
      }
      return *this;
    }

    //! @brief return a vector reference to a given slice of the vector
    //! @param POS_START position of the first element to be referenced
    //! @param POS_END position of the last element to be referenced; will be truncated if beyond the end of the vector
    //! @return const reference to elements [ POS_START, min(POS_END,GetSize()-1]
    template< typename t_DataType>
    VectorReference< t_DataType> VectorInterface< t_DataType>::Slice( const size_t POS_START, const size_t POS_END)
    {
      BCL_Assert( POS_START <= POS_END, "Slice called with POS_START > POS_END!");
      return
        POS_END >= GetSize()
        ?  VectorReference< t_DataType>( GetSize() - POS_START, Begin() + POS_START)
        :  VectorReference< t_DataType>( POS_END + 1 - POS_START, Begin() + POS_START);
    }

    //! @brief return a vector reference to a given slice of the vector
    //! @param POS_START position of the first element to be referenced
    //! @param SIZE max number elements to be referenced; will be truncated if beyond the end of the vector
    //! @return const reference to elements [ POS_START, min(POS_START + SIZE,GetSize())
    template< typename t_DataType>
    VectorReference< t_DataType> VectorInterface< t_DataType>::Reference( const size_t POS_START, const size_t SIZE)
    {
      BCL_Assert( POS_START < GetSize(), "Slice called with POS_START >= GetSize()!");
      return
        SIZE >= GetSize() - POS_START
        ? VectorReference< t_DataType>( GetSize() - POS_START, Begin() + POS_START)
        : VectorReference< t_DataType>( SIZE, Begin() + POS_START);
    }

    //! @brief return a vector reference to a given slice of the vector
    //! @param POS_START position of the first element to be referenced
    //! @param POS_END position of the last element to be referenced; will be truncated if beyond the end of the vector
    //! @return const reference to elements [ POS_START, min(POS_END,GetSize()-1]
    template< typename t_DataType>
    VectorConstReference< t_DataType> VectorInterface< t_DataType>::Slice
    (
      const size_t POS_START,
      const size_t POS_END
    ) const
    {
      BCL_Assert( POS_START <= POS_END, "Slice called with POS_START > POS_END!");
      return
        POS_END >= GetSize()
        ?  VectorConstReference< t_DataType>( GetSize() - POS_START, Begin() + POS_START)
        :  VectorConstReference< t_DataType>( POS_END + 1 - POS_START, Begin() + POS_START);
    }

    //! @brief return a vector reference to a given slice of the vector
    //! @param POS_START position of the first element to be referenced
    //! @param SIZE max number elements to be referenced; will be truncated if beyond the end of the vector
    //! @return const reference to elements [ POS_START, min(POS_START + SIZE,GetSize())
    template< typename t_DataType>
    VectorConstReference< t_DataType> VectorInterface< t_DataType>::Reference
    (
      const size_t POS_START,
      const size_t SIZE
    ) const
    {
      BCL_Assert( POS_START < GetSize(), "Slice called with POS_START >= GetSize()!");
      return
        SIZE >= GetSize() - POS_START
        ? VectorConstReference< t_DataType>( GetSize() - POS_START, Begin() + POS_START)
        : VectorConstReference< t_DataType>( SIZE, Begin() + POS_START);
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_INTERFACE_HPP_
