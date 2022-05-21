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

#ifndef BCL_MATH_RUNNING_STATISTICS_H_
#define BCL_MATH_RUNNING_STATISTICS_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_running_average.h"
#include "bcl_math_running_average_sd.h"
#include "bcl_math_running_min_max.h"
#include "bcl_math_running_sum.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RunningStatistics
    //! @brief Provides functionality to automatically add all math::Running... classes to an enumerated instance
    //!
    //! @see @link example_math_running_statistics.cpp @endlink
    //! @author mendenjl
    //! @date Feb 07, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class RunningStatistics
    {

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! undefined default constructor because no one should be constructing this class
      RunningStatistics();

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief function to add all instances to the enumerated interface for
      //! @tparam t_InterfaceType the enumerated type to add all newly-created t_ClassTypes to
      //! @tparam t_FirstParam is the first template parameter of the class
      //! @tparam t_ClassType the actual class instance; first parameter can be anything, 2nd is the accumulator type
      //!         and the third is the function to call on that accumulator to get the results
      template
      <
        typename t_InterfaceType,
        typename t_FirstParam,
        template
        <
          typename,
          typename t_Accumulator,
          const t_DataType &( t_Accumulator::*t_Function)() const
        > class t_ClassType
      >
      static bool AddInstances()
      {
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningMinMax< t_DataType>, &RunningMinMax< t_DataType>::GetMax>()
        );
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningMinMax< t_DataType>, &RunningMinMax< t_DataType>::GetMin>()
        );
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningMinMax< t_DataType>, &RunningMinMax< t_DataType>::GetRange>()
        );
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningAverageSD< t_DataType>, &RunningAverageSD< t_DataType>::GetStandardDeviation>()
        );
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningAverage< t_DataType>, &RunningAverage< t_DataType>::GetAverage>()
        );
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningSum< t_DataType>, &RunningSum< t_DataType>::GetSum>()
        );
        return true;
      }

      //! @brief function to add all instances that support weighting to the enumerated interface for
      //! @tparam t_InterfaceType the enumerated type to add all newly-created t_ClassTypes to
      //! @tparam t_FirstParam is the first template parameter of the class
      //! @tparam t_ClassType the actual class instance; first parameter can be anything, 2nd is the accumulator type
      //!         and the third is the function to call on that accumulator to get the results
      template
      <
        typename t_InterfaceType,
        typename t_FirstParam,
        template
        <
          typename,
          typename t_Accumulator,
          const t_DataType &( t_Accumulator::*t_Function)() const
        > class t_ClassType
      >
      static bool AddWeightedInstances()
      {
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningAverageSD< t_DataType>, &RunningAverageSD< t_DataType>::GetStandardDeviation>()
        );
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningAverage< t_DataType>, &RunningAverage< t_DataType>::GetAverage>()
        );
        util::Enumerated< t_InterfaceType>::AddInstance
        (
          new t_ClassType< t_FirstParam, RunningSum< t_DataType>, &RunningSum< t_DataType>::GetSum>()
        );
        return true;
      }

    }; // template class RunningStatistics

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RUNNING_STATISTICS_H_
