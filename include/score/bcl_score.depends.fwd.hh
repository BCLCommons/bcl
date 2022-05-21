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

#ifndef BCL_SCORE_DEPENDS_FWD_HH_
#define BCL_SCORE_DEPENDS_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// external includes - sorted alphabetically

// AUTOMATICALLY GENERATED FILE -- DO NOT EDIT!
// This file contains forward declarations needed by the forward header for the score namespace

namespace bcl
{
////////////////////////////////
// class forward declarations //
////////////////////////////////

  namespace biol
  {
    class AABase;
  } // namespace biol
  namespace function
  {
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryInterface;
  } // namespace function
  namespace math
  {
    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterfaceSerializable;
  } // namespace math
  namespace restraint
  {
    class AtomDistanceAssignment;
  } // namespace restraint
  namespace util
  {
    template< typename t_DataType>
    class ShPtr;

    template< typename t_DataType, typename t_Derived>
    class Enum;
  } // namespace util

///////////////////////
// external typedefs //
////////////////////////

} // namespace bcl

#endif // BCL_SCORE_DEPENDS_FWD_HH_
