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
#include "linal/bcl_linal_vector.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Vector< double>;
    template class BCL_API Vector< float>;
    template class BCL_API Vector< int>;
    // these are necessary to ensure that we get both the cl_uint type, which is used in
    // some opencl classes, as well as size_t.  It is possible that cl_uint and size_t are the same
    // type on some platforms.
    template class BCL_API Vector< unsigned int>;
    template class BCL_API Vector< unsigned long>;
    template class BCL_API Vector< unsigned long long>;
    template class BCL_API Vector< char>;

  } // namespace linal
} // namespace bcl
