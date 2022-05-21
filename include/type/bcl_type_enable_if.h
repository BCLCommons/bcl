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

#ifndef BCL_TYPE_ENABLE_IF_H_
#define BCL_TYPE_ENABLE_IF_H_

// include the namespace header
#include "bcl_type.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace type
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnableIf
    //! @brief a simple helper class that can be used as a return type to enable a function if a parameter was true
    //! @see @link example_type_enable_if.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< bool ENABLE, typename t_Type>
    struct EnableIf :
      public Base
    {
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnableIf
    //! @brief specialization for ENABLE == true; then define t_Type
    //! @see @link example_type_enable_if.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Type>
    struct EnableIf< true, t_Type> :
      public Base
    {
      typedef t_Type Type;
    };

  } // namespace type
} // namespace bcl

#endif //BCL_TYPE_ENABLE_IF_H_
