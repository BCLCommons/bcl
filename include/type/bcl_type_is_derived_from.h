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

#ifndef BCL_TYPE_IS_DERIVED_FROM_H_
#define BCL_TYPE_IS_DERIVED_FROM_H_

// include the namespace header
#include "bcl_type.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_type_compare.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace type
  {
    // anonymous namespace, because the internal class should only be used by this file
    namespace
    {
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class IsDerived
      //! @brief this class is used by IsDerivedFrom to test whether t_Derived is derived from t_Base
      //!
      //! @see @link example_type_is_derived_from.cpp @endlink
      //!
      //! @author mendenjl
      //! @date Nov 29, 2012
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      template< typename t_Derived, typename t_Base>
      struct IsDerived :
        public Base
      {
      private:

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //!
        //! @class Tester
        //! @brief a trivial class that offers conversion functions to t_Base and t_Derived
        //! @remarks example unnecessary
        //!
        //! @author mendenjl
        //! @date Nov 29, 2012
        //!
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        struct Tester
        {
          operator t_Base*() const;
          operator t_Derived*();
        };

      public:

        //! @brief a function declaration; used only to determine whether t_Derived is derived from t_Base
        //! @note this is an example of Substitution Failure Is Not An Error (SFINAE)
        template< typename t_ArbitraryType>
        static Yes Test( t_Derived *, t_ArbitraryType);
        static No Test( t_Base *, int);

        // If you are really dying to know why this works
        //! @see @link http://stackoverflow.com/questions/2910979/how-is-base-of-works @endlink
        enum { value = sizeof( Test( Tester(), int())) == sizeof( Yes)};
      };
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class IsDerivedFrom
    //! @brief a simple helper class that can be used to determine at compile time whether one class is derived from the other
    //! @see @link example_type_is_derived_from.cpp @endlink
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Derived, typename t_Base>
    struct IsDerivedFrom :
      public Base
    {
      // calculate the actual value after removing const & from both t_Derived and t_Base
      enum
      {
        value = IsDerived< typename RemoveConstRef< t_Derived>::Type, typename RemoveConstRef< t_Base>::Type>::value
      };
    };

  } // namespace type
} // namespace bcl

#endif //BCL_TYPE_IS_DERIVED_FROM_H_
