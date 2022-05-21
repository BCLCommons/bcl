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

#ifndef BCL_TYPE_REMOVE_CONST_REF_H_
#define BCL_TYPE_REMOVE_CONST_REF_H_

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
    //! @class RemoveReference
    //! @brief a simple helper class that removes the & from a type
    //! @see @link example_type_remove_const_ref.cpp @endlink
    //!
    //! @details RemoveReference< char>::Type is char
    //!          RemoveReference< char &>::Type is char
    //!          This can be useful when the type passed to a template may or may not be a reference type, and we need
    //!          the non-referenced type
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Type>
    struct RemoveReference :
      public Base
    {
      typedef t_Type Type;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RemoveReference
    //! @brief specialization for reference types
    //! @see @link example_type_remove_const_ref.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Type>
    struct RemoveReference< t_Type &> :
      public Base
    {
      typedef t_Type Type;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RemoveConst
    //! @brief a simple helper class that removes the const from a type
    //! @see @link example_type_remove_const_ref.cpp @endlink
    //!
    //! @details RemoveConst< char>::Type is char
    //!          RemoveConst< const char>::Type is char
    //!          RemoveConst< char const *>::Type is const char * because the type (pointer to const char) is non-const
    //!          RemoveConst< char const * const>::Type is const char *
    //!          This can be useful when the type passed to a template may or may not be const, and we need the
    //!          non-const type
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Type>
    struct RemoveConst
    {
      typedef t_Type Type;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RemoveConst
    //! @brief specialization for const types
    //! @see @link example_type_remove_const_ref.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Type>
    struct RemoveConst< t_Type const>
    {
      typedef t_Type Type;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RemoveConstRef
    //! @brief combination of remove const with remove reference
    //! @see @link example_type_remove_const_ref.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Type>
    struct RemoveConstRef :
      public RemoveConst< typename RemoveReference< t_Type>::Type>
    {
    };

  } // namespace type
} // namespace bcl

#endif //BCL_TYPE_REMOVE_CONST_REF_H_
