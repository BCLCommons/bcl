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

#ifndef BCL_TYPE_IS_MAP_H_
#define BCL_TYPE_IS_MAP_H_

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
    //! @class IsMap
    //! @brief a simple helper class that determines whether a type has typedefs value_type and key_type
    //! @see @link example_type_is_map.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Jul 14, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Type>
    class IsMap :
      public Base
    {
      //! @brief a function declaration; used only to determine whether t_PseudoType has map-required typedefs
      //! @note this is an example of Substitution Failure Is Not An Error (SFINAE)
      template< typename t_PseudoType>
      static Yes Test( typename t_PseudoType::mapped_type *, typename t_PseudoType::key_type *);

      //! @brief a function declaration; used only to determine whether t_PseudoType has map-required typedefs
      template< typename t_PseudoType>
      static No Test(...);

    public:

      enum { value = ( sizeof( Test< t_Type>( 0, 0)) == sizeof( Yes))};
    };

  } // namespace type
} // namespace bcl

#endif //BCL_TYPE_IS_MAP_H_
