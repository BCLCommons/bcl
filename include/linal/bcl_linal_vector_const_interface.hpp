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

#ifndef BCL_LINAL_VECTOR_CONST_INTERFACE_HPP_
#define BCL_LINAL_VECTOR_CONST_INTERFACE_HPP_

// include source of this class
#include "bcl_linal_vector_const_interface.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_const_reference.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual destructor
    template< typename t_DataType>
    VectorConstInterface< t_DataType>::~VectorConstInterface()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief test whether the container is empty
    //! @return true if the container is empty
    template< typename t_DataType>
    bool VectorConstInterface< t_DataType>::IsEmpty() const
    {
      return GetSize() == size_t( 0);
    }

    //! @brief test whether all elements in the container are defined
    //! @return true if all elements in the container are defined
    template< typename t_DataType>
    bool VectorConstInterface< t_DataType>::IsDefined() const
    {
      return math::Statistics::IsDefined( Begin(), End());
    }

    //! @brief reference to first element
    //! @return const reference to first element in range containing all elements of Vector
    template< typename t_DataType>
    const t_DataType &VectorConstInterface< t_DataType>::First() const
    {
      BCL_Assert( GetSize(), "Cannot get first element of empty " + GetClassIdentifier());
      return *Begin();
    }

    //! @brief reference to last element
    //! @return const reference to last element in range containing all elements of Vector
    template< typename t_DataType>
    const t_DataType &VectorConstInterface< t_DataType>::Last() const
    {
      BCL_Assert( GetSize(), "Cannot get last element of empty " + GetClassIdentifier());
      return *( End() - 1);
    }

    //! @brief return a vector reference to a given slice of the vector
    //! @param POS_START position of the first element to be referenced
    //! @param POS_END position of the last element to be referenced; will be truncated if beyond the end of the vector
    //! @return const reference to elements [ POS_START, min(POS_END,GetSize()-1]
    template< typename t_DataType>
    VectorConstReference< t_DataType> VectorConstInterface< t_DataType>::Slice
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
    VectorConstReference< t_DataType> VectorConstInterface< t_DataType>::Reference
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

  ////////////////
  // operations //
  ////////////////

    //! @brief norm = length of vector
    //! @return length of vector
    template< typename t_DataType>
    t_DataType VectorConstInterface< t_DataType>::Norm() const
    {
      return math::Statistics::Norm( Begin(), End());
    }

    //! @brief square norm = square length of vector
    //! @return square length of vector
    template< typename t_DataType>
    t_DataType VectorConstInterface< t_DataType>::SquareNorm() const
    {
      return math::Statistics::SquareNorm( Begin(), End());
    }

    //! @brief sum up all elements
    //! @return sum of all elements
    template< typename t_DataType>
    t_DataType VectorConstInterface< t_DataType>::Sum() const
    {
      return math::Statistics::Sum( Begin(), End());
    }

    //! @brief min of all elements
    //! @return min of all elements
    template< typename t_DataType>
    t_DataType VectorConstInterface< t_DataType>::Min() const
    {
      return math::Statistics::MinimumValue( Begin(), End());
    }

    //! @brief max of all elements
    //! @return max of all elements
    template< typename t_DataType>
    t_DataType VectorConstInterface< t_DataType>::Max() const
    {
      return math::Statistics::MaximumValue( Begin(), End());
    }

    //! @brief compare two vectors for equality
    //! @param RHS rhs vector to compare to
    //! @return true if size and content match
    template< typename t_DataType>
    bool VectorConstInterface< t_DataType>::operator ==( const VectorConstInterface< t_DataType> &RHS) const
    {
      if( GetSize() != RHS.GetSize())
      {
        return false;
      }

      return std::equal( Begin(), End(), RHS.Begin());
    }

    //! @brief compare two vectors for inequality
    //! @param RHS rhs vector to compare to
    //! @return true if size and/or content do not match
    template< typename t_DataType>
    bool VectorConstInterface< t_DataType>::operator !=( const VectorConstInterface< t_DataType> &RHS) const
    {
      return !operator ==( RHS);
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_CONST_INTERFACE_HPP_
