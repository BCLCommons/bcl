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

#ifndef BCL_FUNCTION_FWD_HH_
#define BCL_FUNCTION_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the function namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace function
  {
  /////////////////////
  // regular classes //
  /////////////////////

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_IntermediateTypeLeft, typename t_IntermediateTypeRight, typename t_ResultType>
    class BinaryAdapter;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryCached;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryInterface;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinarySum;

    template< typename t_ClassType, typename t_ResultType>
    class Member;

    template< typename t_ClassType, typename t_ResultType>
    class MemberConst;

    template< typename t_ClassType, typename t_ArgumentType, typename t_ResultType>
    class MemberUnary;

    template< typename t_ClassType, typename t_ArgumentType, typename t_ResultType>
    class MemberUnaryConst;

    template< typename t_ArgumentType, typename t_IntermediateType, typename t_ResultType>
    class UnaryAdapter;

    template< typename t_ArgumentType, typename t_ResultType>
    class UnaryCached;

    template< typename t_ArgumentType, typename t_ResultType>
    class UnaryInterface;

  //////////////
  // typedefs //
  //////////////

  } // namespace function
} // namespace bcl

#endif // BCL_FUNCTION_FWD_HH_
