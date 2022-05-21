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

#ifndef BCL_LINAL_VECTOR_INTERFACE_H_
#define BCL_LINAL_VECTOR_INTERFACE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_const_interface.h"
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorInterface
    //! @brief TODO: add a brief comment to this class
    //!
    //! @see @link example_linal_vector_interface.cpp @endlink
    //! @author woetzen
    //! @date Aug 11, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class VectorInterface :
      public VectorConstInterface< t_DataType>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new VectorInterface< t_DataType, size_t N>
      virtual VectorInterface< t_DataType> *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      typedef t_DataType*       iterator;
      typedef t_DataType&       reference;

      //! @brief size of vector
      //! @return size of Vector
      virtual size_t GetSize() const = 0;

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of Vector
      virtual const t_DataType *Begin() const = 0;

      //! @brief pointer to First Element
      //! @return pointer to first element in range containing all elements of Vector
      virtual t_DataType *Begin() = 0;

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Vector
      virtual const t_DataType *End() const = 0;

      //! @brief pointer to end of range
      //! @return pointer to address one after last element in Vector
      virtual t_DataType *End() = 0;

      //! @brief return a vector reference to a given slice of the vector
      //! @param POS_START position of the first element to be referenced
      //! @param POS_END position of the last element to be referenced; will be truncated if beyond the end of the vector
      //! @return const reference to elements [ POS_START, min(POS_END,GetSize()-1]
      virtual VectorReference< t_DataType> Slice( const size_t POS_START, const size_t POS_END = size_t( -1));

      //! @brief return a vector reference to a given slice of the vector
      //! @param POS_START position of the first element to be referenced
      //! @param SIZE max number elements to be referenced; will be truncated if beyond the end of the vector
      //! @return const reference to elements [ POS_START, min(POS_START + SIZE,GetSize())
      virtual VectorReference< t_DataType> Reference( const size_t POS_START, const size_t SIZE = size_t( -1));

      //! @brief return a vector reference to a given slice of the vector
      //! @param POS_START position of the first element to be referenced
      //! @param POS_END position of the last element to be referenced; will be truncated if beyond the end of the vector
      //! @return const reference to elements [ POS_START, min(POS_END,GetSize()-1]
      virtual VectorConstReference< t_DataType> Slice( const size_t POS_START, const size_t POS_END = size_t( -1)) const;

      //! @brief return a vector reference to a given slice of the vector
      //! @param POS_START position of the first element to be referenced
      //! @param SIZE max number elements to be referenced; will be truncated if beyond the end of the vector
      //! @return const reference to elements [ POS_START, min(POS_START + SIZE,GetSize())
      virtual VectorConstReference< t_DataType> Reference( const size_t POS_START, const size_t SIZE = size_t( -1)) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief return non-const reference to element ( POS)
      //! @param POS position of the element requested
      //! @return non-const reference to element ( POS)
      virtual t_DataType &operator()( const size_t POS) = 0;

      //! @brief return const reference to element ( POS)
      //! @param POS position of the element requested
      //! @return const reference to element ( POS)
      virtual const t_DataType &operator()( const size_t POS) const = 0;

      //! @brief return pointer to const element ( POS)
      //! @param POS position of the element requested
      //! @return pointer to const element ( POS)
      virtual const t_DataType *operator[]( const size_t POS) const = 0;

      //! @brief return pointer to element ( POS)
      //! @param POS position of the element requested
      //! @return pointer to element ( POS)
      virtual t_DataType *operator[]( const size_t POS) = 0;

      //! @brief equal operator
      //! @param VALUE all elements are set to that value
      //! @return reference to this assigned Vector
      virtual VectorInterface< t_DataType> &operator =( const t_DataType &VALUE) = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief normalize the vector (such that inner product == 1)
      VectorInterface< t_DataType> &Normalize();

      //! @brief set to sum SUM (sum of all elements = SUM)
      VectorInterface< t_DataType> &SetToSum( const t_DataType &SUM = t_DataType( 1.0));

      //! @brief set all elements in vector to random numbers
      //! @param MIN minimum value to set elements in the vector to
      //! @param MAX maximum value to set elements in the vector to
      VectorInterface< t_DataType> &SetRand( const t_DataType MIN, const t_DataType MAX);

      //! @brief copies elements of argument VECTOR into this object at position ( POS)
      //! @param POS beginning of range to replace elements
      //! @param VECTOR vector of elements to replace elements in this vector between POS and POS + VECTOR.GetSize()
      VectorInterface< t_DataType> &ReplaceElements( const size_t POS, const VectorConstInterface< t_DataType> &VECTOR);

    }; // template class VectorInterface

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API VectorInterface< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorInterface< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorInterface< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorInterface< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorInterface< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorInterface< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorInterface< char>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_INTERFACE_H_
