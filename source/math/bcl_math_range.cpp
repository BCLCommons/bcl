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
#include "math/bcl_math_range.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    const std::string &RangeBorders::GetConditionLeftChars()
    {
      static const std::string s_condition_left_char( "[(");
      return s_condition_left_char;
    }
    const std::string &RangeBorders::GetConditionRightChars()
    {
      static const std::string s_condition_right_char( "])");
      return s_condition_right_char;
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Range< double>;
    template class BCL_API Range< float>;
    template class BCL_API Range< int>;
    template class BCL_API Range< unsigned int>;
    template class BCL_API Range< unsigned long>;
    template class BCL_API Range< unsigned long long>;
    template class BCL_API Range< bool>;
    template class BCL_API Range< char>;

  } // namespace math
} // namespace bcl
