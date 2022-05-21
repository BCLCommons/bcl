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

#ifndef BCL_TYPE_H_
#define BCL_TYPE_H_

// include the namespace forward header
#include "bcl_type.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically
#include <string>

namespace bcl
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_type.h
  //! @namespace namespace for template meta-programming classes that determine compile-time available type information
  //! Classes in this namespace share several features:
  //! 1. they have absolutely no defined functions
  //! 2. None of the classes can ever be constructed, copied, or changed in any way
  //! 3. All classes will have at least one instantiation with exactly one of the following:
  //!    i.  An externally visible typedef
  //!    ii. An anonymous enum with some value which can be used in other contexts
  //! As an implication of #2, no forward headers should be made for classes in this namespace
  //!
  //! @brief For more detailed explanation of template metaprogramming, look on wikipedia, or see
  //! @brief Modern C++ Design: Generic Programming and Design Patterns Applied. Addison-Wesley. ISBN 3-8266-1347-3.
  //!
  //! @see @link example_type.cpp @endlink
  //! @author mendenjl
  //! @date Nov 29, 2012
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace type
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API const std::string &GetNamespaceIdentifier();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Base
    //! @brief a base class for all objects in this namespace, renders them non-copyable and non-constructable
    //! @remarks example unnecessary
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    struct BCL_API Base
    {
    private:

      //! @brief undefined default constructor
      Base();

      //! @brief undefined copy constructor
      Base( const Base &);

      //! @brief undefined assignment operator
      Base &operator =( const Base &);
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Yes
    //! @brief represents a true statement
    //! @remarks example unnecessary
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    struct BCL_API Yes :
      public Base
    {
    private:

      long m_Weight[ 4]; //!< Data member, used solely to make sizeof(Yes) != sizeof(No)

    public:
      //! enum for logical value of the class; e.g. t_Type::e_Yes, where t_Type is either Yes or No,
      //! and get back the value as a bool
      enum
      {
        e_Yes = 1,
        e_No  = 0
      };
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class No
    //! @brief represents a false statement
    //! @remarks example unnecessary
    //!
    //! @author mendenjl
    //! @date Nov 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    struct BCL_API No :
      public Base
    {
      //! enum for logical value of the class; e.g. t_Type::e_Yes, where t_Type is either Yes or No,
      //! and get back the value as a bool
      enum
      {
        e_Yes = 0,
        e_No  = 1
      };
    };
  } // namespace type
} // namespace bcl

#endif //BCL_TYPE_H_
