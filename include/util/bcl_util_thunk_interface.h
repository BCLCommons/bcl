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

#ifndef BCL_UTIL_THUNK_INTERFACE_H_
#define BCL_UTIL_THUNK_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ThunkInterface
    //! @brief This class is an Interface class for functions with no arguments (thunks)
    //!
    //! @remarks example unnecessary
    //! @author riddeljs
    //! @date 10.23.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ResultType>
    class ThunkInterface :
      public ObjectInterface
    {
    public:

    ///////////////
    // operators //
    ///////////////

      //! virtual operator taking no argument and returning a t_ResultType object
      virtual t_ResultType operator()() = 0;

    }; //end template class ThunkInterface

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_THUNK_INTERFACE_H_
