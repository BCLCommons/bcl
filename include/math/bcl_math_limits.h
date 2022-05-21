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

#ifndef BCL_MATH_LIMITS_H_
#define BCL_MATH_LIMITS_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_math_limits.h
    //! @brief Set of templated functions to retrieve the highest and lowest possible values for different types
    //! @see @link example_math_limits.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Oct 18, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! @brief convenience function to determine the most negative value that is bounded (e.g. non -inf) for a type
    //! @note specialization for integral types
    template< typename t_DataType>
    typename type::EnableIf< std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type GetLowestBoundedValue()
    {
      return std::numeric_limits< t_DataType>::min();
    }

    //! @brief convenience function to determine the most negative value that is bounded (e.g. non -inf) for a type
    //! @note specialization for integral types
    template< typename t_DataType>
    typename type::EnableIf< !std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type GetLowestBoundedValue()
    {
      return -std::numeric_limits< t_DataType>::max();
    }

    //! @brief convenience function to determine the most negative / smallest  value (possibly -inf) for a type
    //! @note specialization for integral types
    template< typename t_DataType>
    typename type::EnableIf< std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type GetLowestUnboundedValue()
    {
      return std::numeric_limits< t_DataType>::min();
    }

    //! @brief convenience function to determine the most negative / smallest value (possibly -inf) for a type
    //! @note specialization for float types
    template< typename t_DataType>
    typename type::EnableIf< !std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type GetLowestUnboundedValue()
    {
      return -std::numeric_limits< t_DataType>::infinity();
    }

    //! @brief convenience function to determine the most positive value that is bounded (e.g. non +inf) for a type
    template< typename t_DataType>
    t_DataType GetHighestBoundedValue()
    {
      return std::numeric_limits< t_DataType>::max();
    }

    //! @brief convenience function to determine the most positive value (possibly +inf) for a type
    //! @note specialization for integral types
    template< typename t_DataType>
    typename type::EnableIf< std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type GetHighestUnboundedValue()
    {
      return std::numeric_limits< t_DataType>::max();
    }

    //! @brief convenience function to determine the most positive value (possibly +inf) for a type
    //! @note specialization for float types
    template< typename t_DataType>
    typename type::EnableIf< !std::numeric_limits< t_DataType>::is_integer, t_DataType>::Type GetHighestUnboundedValue()
    {
      return std::numeric_limits< t_DataType>::infinity();
    }

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_LIMITS_H_
