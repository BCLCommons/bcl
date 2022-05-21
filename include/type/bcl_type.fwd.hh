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

#ifndef BCL_TYPE_FWD_HH_
#define BCL_TYPE_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the type namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace type
  {
  /////////////////////
  // regular classes //
  /////////////////////

    struct Base;
    struct No;
    struct Yes;

  //////////////////////
  // template classes //
  //////////////////////

    template< bool CHOICE, typename t_TypeA, typename t_TypeB>
    struct Chooser;

    template< typename t_TypeA, typename t_TypeB>
    struct Compare;

    template< bool ENABLE, typename t_Type>
    struct EnableIf;

    template< typename t_TypeA, typename t_TypeB>
    struct FuzzyCompare;

    template< typename t_Type, typename t_OtherType>
    struct IsA;

    template< typename t_Derived, typename t_Base>
    struct IsDerivedFrom;

    template< typename t_Type>
    class IsMap;

    template< typename t_Type>
    class IsSequence;

    template< typename t_Type>
    struct RemoveConst;

    template< typename t_Type>
    struct RemoveConstRef;

    template< typename t_Type>
    struct RemoveReference;

  //////////////
  // typedefs //
  //////////////

  } // namespace type
} // namespace bcl

#endif // BCL_TYPE_FWD_HH_
