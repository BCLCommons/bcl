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

#ifndef BCL_TYPE_IS_SEQUENCE_H_
#define BCL_TYPE_IS_SEQUENCE_H_

// include the namespace header
#include "bcl_type.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_type_compare.h"
#include "bcl_type_is_derived_from.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace type
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class IsSequence
    //! @brief a helper class that determines whether a type has a typedef const_iterator, e.g. is a sequence
    //! @see @link example_type_is_sequence.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Jul 14, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Type>
    class IsSequence :
      public Base
    {
      //! @brief a function declaration; used only to determine whether t_PseudoType has a typedef const_iterate
      //! @note this is an example of Substitution Failure Is Not An Error (SFINAE)
      template< typename t_PseudoType>
      static Yes Test( typename t_PseudoType::const_iterator *);

      //! @brief a function declaration; used only to determine whether t_PseudoType has a typedef const_iterate
      template< typename t_PseudoType>
      static No Test(...);

    public:

      enum { value = ( sizeof( Test< t_Type>( 0)) == sizeof( Yes))};
    };

  } // namespace type
} // namespace bcl

#endif //BCL_TYPE_IS_SEQUENCE_H_
