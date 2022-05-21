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

#ifndef BCL_MATH_COMPARISONS_H_
#define BCL_MATH_COMPARISONS_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Comparisons
    //! @brief collects different comparison operations, mainly for sorting purposes (<, >, <=, >=)
    //!
    //! @see @link example_math_comparisons.cpp @endlink
    //! @author woetzen
    //! @date Oct 7, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Comparisons :
      public util::Enumerate
                   <
                     util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >,
                     Comparisons< t_DataType>
                   >
    {

      friend class util::Enumerate
                    <
                      util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >,
                      Comparisons< t_DataType>
                    >;

    public:

      typedef typename util::Enumerate
                             <
                               util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >,
                               Comparisons< t_DataType>
                             >::EnumType Comparison;

      typedef typename util::Enumerate
                             <
                               util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >,
                               Comparisons< t_DataType>
                             >::EnumDataType ComparisonData;

      using util::Enumerate
            <
              util::ShPtr< util::BinaryFunctionInterface< t_DataType, t_DataType, bool> >,
              Comparisons< t_DataType>
            >::AddEnum;

    //////////
    // data //
    //////////

      const Comparison e_Less;         //!< "<" operator
      const Comparison e_LessEqual;    //!< "<=" operator
      const Comparison e_Greater;      //!< ">" operator
      const Comparison e_GreaterEqual; //!< ">=" operator
      const Comparison e_Equal;        //!< "==" operator
      const Comparison e_NotEqual;     //!< "!=" operator

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct all functors
      Comparisons();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief create a unary predicate function for the comparison and a LHS VALUE
      //! it return a function of the form f(x) = x {comp} VALUE
      //! @param COMPARISON the comparison function to use
      //! @param VALUE the right hand side value
      //! @return unary function which can be used as a unary predicate
      BinaryFunctionBindSecond< t_DataType, t_DataType, bool> CreateUnaryPredicate
      (
        const Comparison &COMPARISON,
        const t_DataType &VALUE
      ) const;

    }; // template class Comparisons

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Comparisons< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Comparisons< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Comparisons< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Comparisons< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Comparisons< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Comparisons< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Comparisons< bool>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Comparisons< char>;

  } // namespace math

  namespace util
  {

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<             double,             double, bool> >, math::Comparisons<             double> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<              float,              float, bool> >, math::Comparisons<              float> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<                int,                int, bool> >, math::Comparisons<                int> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<       unsigned int,       unsigned int, bool> >, math::Comparisons<       unsigned int> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<      unsigned long,      unsigned long, bool> >, math::Comparisons<      unsigned long> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface< unsigned long long, unsigned long long, bool> >, math::Comparisons< unsigned long long> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<               bool,               bool, bool> >, math::Comparisons<               bool> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<               char,               char, bool> >, math::Comparisons<               char> >;

  } // namespace util
} // namespace bcl

#endif // BCL_MATH_COMPARISONS_H_ 
