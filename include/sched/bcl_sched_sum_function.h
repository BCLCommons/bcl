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

#ifndef BCL_SCHED_SUM_FUNCTION_H_
#define BCL_SCHED_SUM_FUNCTION_H_

// include the namespace header
#include "bcl_sched.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_sum_function.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SumFunction
    //! @brief This class is a template for t_ResultType = Function( t_ArgumentType&) allowing addition
    //! of functions with coefficients.
    //! @details Assuming y being element of t_ResultType, x being element of t_ArgumentType and s being scalar
    //! this class can represent the following functions:
    //! y1 = s1 * g1(x)
    //! y2 = s2 * g2(x)
    //!
    //! f(x) = y0 + y1 + y2 = y0 + s1 * g1(x) + s2 * g2(x)
    //!
    //! operators provide the functionality to divide the sumfunction, add an absolute to y0, add coefficients
    //! s3 * g3(x), add another sumfunction and to create sumfunctions from function interfaces and coefficients.
    //!
    //! @see @link example_sched_sum_function.cpp @endlink
    //! @author woetzen
    //! @date Aug 21, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class SumFunction :
      public math::SumFunction< t_ArgumentType, t_ResultType>
    {
    public:

      //! @typedef single term with scalar
      typedef typename math::SumFunction< t_ArgumentType, t_ResultType> Base;
      typedef typename Base::Term Term;

      //! @typedef all terms
      typedef storage::Vector< Term> Functions;

      //! @typedef const_iterator for iterating over function
      typedef typename Functions::const_iterator const_iterator;
      //! @typedef iterator for iterating over function
      typedef typename Functions::iterator iterator;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SumFunction( const std::string &SCHEME = "");

      //! @brief construct from ShPtr on Function, optional coefficient, and absolute value
      //! y = y0 + s1*g1(x)
      //! @param ABSOLUTE    y0 as absolute in sumfunction - default 0
      //! @param COEFFICIENT s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param SP_FUNCTION g1(x) as part of sumfunction as ShPtr to FunctionInterface
      SumFunction
      (
        const t_ResultType &ABSOLUTE,
        const double &COEFFICIENT,
        const util::Implementation< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > &SP_FUNCTION,
        const std::string &SCHEME = ""
      );

      //! @brief copy constructor
      //! @return pointer to new SumFunction which is a copy of this
      SumFunction< t_ArgumentType, t_ResultType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief actual f(x) operator taking ARGUMENT and returning the sum of all function values
      //! @brief ARGUMENT x to f(x)
      //! @return the value of f(x) = sum
      t_ResultType
      operator()
      (
        const t_ArgumentType &ARGUMENT
      ) const;

    }; // class SumFunction

  } // namespace sched
} // namespace bcl

// include template implementation file
#include "bcl_sched_sum_function.hpp"

#endif // BCL_SCHED_SUM_FUNCTION_H_
