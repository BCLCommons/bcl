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

#ifndef BCL_LINAL_VECTOR_CONST_INTERFACE_H_
#define BCL_LINAL_VECTOR_CONST_INTERFACE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorConstInterface
    //! @brief TODO: add a brief comment to this class
    //!
    //! @see @link example_linal_vector_const_interface.cpp @endlink
    //! @author loweew, mendenjl
    //! @date Feb 23, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class VectorConstInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new VectorConstInterface< t_DataType, size_t N>
      virtual VectorConstInterface< t_DataType> *Clone() const = 0;

      //! @brief virtual destructor
      virtual ~VectorConstInterface();

    /////////////////
    // data access //
    /////////////////

      typedef const t_DataType* const_iterator;
      typedef const t_DataType& const_reference;

      //! @brief size of vector
      //! @return size of Vector
      virtual size_t GetSize() const = 0;

      //! @brief test whether the container is empty
      //! @return true if the container is empty
      virtual bool IsEmpty() const;

      //! @brief test whether all elements in the container are defined
      //! @return true if all elements in the container are defined
      virtual bool IsDefined() const;

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of Vector
      virtual const t_DataType *Begin() const = 0;

      //! @brief reference to first element
      //! @return const reference to first element in range containing all elements of Vector
      virtual const t_DataType &First() const;

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Vector
      virtual const t_DataType *End() const = 0;

      //! @brief reference to last element
      //! @return const reference to last element in range containing all elements of Vector
      virtual const t_DataType &Last() const;

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

      //! @brief return const reference to element ( POS)
      //! @param POS position of the element requested
      //! @return const reference to element ( POS)
      virtual const t_DataType &operator()( const size_t POS) const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief norm = length of vector
      //! @return length of vector
      t_DataType Norm() const;

      //! @brief square norm = square length of vector
      //! @return square length of vector
      t_DataType SquareNorm() const;

      //! @brief sum up all elements
      //! @return sum of all elements
      t_DataType Sum() const;

      //! @brief get the min
      //! @return the min value
      t_DataType Min() const;

      //! @brief get the max
      //! @return the max value
      t_DataType Max() const;

      //! @brief compare two vectors for equality
      //! @param RHS rhs vector to compare to
      //! @return true if size and content match
      bool operator ==( const VectorConstInterface< t_DataType> &RHS) const;

      //! @brief compare two vectors for inequality
      //! @param RHS rhs vector to compare to
      //! @return true if size and/or content do not match
      bool operator !=( const VectorConstInterface< t_DataType> &RHS) const;

    }; // template class VectorConstInterface

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API VectorConstInterface< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorConstInterface< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorConstInterface< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorConstInterface< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorConstInterface< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorConstInterface< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API VectorConstInterface< char>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_CONST_INTERFACE_H_
