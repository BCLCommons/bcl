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
#include "io/bcl_io_serialization_builtin.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API SerializationBuiltin< double>;
    template class BCL_API SerializationBuiltin< float>;
    template class BCL_API SerializationBuiltin< short>;
    template class BCL_API SerializationBuiltin< int>;
    template class BCL_API SerializationBuiltin< long>;
    template class BCL_API SerializationBuiltin< long long>;
    template class BCL_API SerializationBuiltin< unsigned short>;
    template class BCL_API SerializationBuiltin< unsigned int>;
    template class BCL_API SerializationBuiltin< unsigned long>;
    template class BCL_API SerializationBuiltin< unsigned long long>;
    template class BCL_API SerializationBuiltin< bool>;
    template class BCL_API SerializationBuiltin< char>;
    template class BCL_API SerializationBuiltin< signed char>;
    template class BCL_API SerializationBuiltin< unsigned char>;
    template class BCL_API SerializationBuiltin< std::string>;
    template class BCL_API SerializationBuiltin< util::ObjectDataLabel>;

  } // namespace io
} // namespace bcl
