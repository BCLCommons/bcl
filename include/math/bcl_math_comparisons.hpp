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

#ifndef BCL_MATH_COMPARISONS_HPP_
#define BCL_MATH_COMPARISONS_HPP_

// include header of this class
#include "bcl_math_comparisons.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_binary_function_bind_second.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct all functors
    template< typename t_DataType>
    Comparisons< t_DataType>::Comparisons() :
      e_Less
      (
        AddEnum
        (
          "less",
          util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >
          ( new util::BinaryFunctionSTLWrapper< std::less< t_DataType> >())
        )
      ), //!< "<" operator
      e_LessEqual
      (
        AddEnum
        (
          "less_equal",
          util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >
          ( new util::BinaryFunctionSTLWrapper< std::less_equal< t_DataType> >())
        )
      ), //!< "<=" operator
      e_Greater
      (
        AddEnum
        (
          "greater",
          util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >
          ( new util::BinaryFunctionSTLWrapper< std::greater< t_DataType> >())
        )
      ), //!< ">" operator
      e_GreaterEqual
      (
        AddEnum
        (
          "greater_equal",
          util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >
          ( new util::BinaryFunctionSTLWrapper< std::greater_equal< t_DataType> >())
        )
      ), //!< ">=" operator
      e_Equal
      (
        AddEnum
        (
          "equal",
          util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >
          ( new util::BinaryFunctionSTLWrapper< std::equal_to< t_DataType> >())
        )
      ), //!< "==" operator
      e_NotEqual
      (
        AddEnum
        (
          "not_equal",
          util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >
          ( new util::BinaryFunctionSTLWrapper< std::not_equal_to< t_DataType> >())
        )
      ) //!< "!=" operator
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &Comparisons< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief create a unary predicate function for the comparison and a LHS VALUE
    //! it return a function of the form f(x) = x {comp} VALUE
    //! @param COMPARISON the comparison function to use
    //! @param VALUE the right hand side value
    //! @return unary function which can be used as a unary predicate
    template< typename t_DataType>
    BinaryFunctionBindSecond< t_DataType, t_DataType, bool> Comparisons< t_DataType>::CreateUnaryPredicate
    (
      const Comparison &COMPARISON,
      const t_DataType &VALUE
    ) const
    {
      return BinaryFunctionBindSecond< t_DataType, t_DataType, bool>( *COMPARISON, VALUE);
    }

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_COMPARISONS_HPP_ 
