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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math_comparisons.hpp"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Comparisons< double>;
    template class BCL_API Comparisons< float>;
    template class BCL_API Comparisons< int>;
    template class BCL_API Comparisons< unsigned int>;
    template class BCL_API Comparisons< unsigned long>;
    template class BCL_API Comparisons< unsigned long long>;
    template class BCL_API Comparisons< bool>;
    template class BCL_API Comparisons< char>;

  } // namespace math

  namespace util
  {

    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<             double,             double, bool> >, math::Comparisons<             double> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<              float,              float, bool> >, math::Comparisons<              float> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<                int,                int, bool> >, math::Comparisons<                int> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<       unsigned int,       unsigned int, bool> >, math::Comparisons<       unsigned int> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<      unsigned long,      unsigned long, bool> >, math::Comparisons<      unsigned long> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface< unsigned long long, unsigned long long, bool> >, math::Comparisons< unsigned long long> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<               bool,               bool, bool> >, math::Comparisons<               bool> >;
    template class BCL_API Enumerate< ShPtr< BinaryFunctionInterface<               char,               char, bool> >, math::Comparisons<               char> >;

  } // namespace util
} // namespace bcl
