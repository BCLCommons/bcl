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

#ifndef BCL_FIND_FWD_HH_
#define BCL_FIND_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the find namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace find
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class PickBodyExtent;
    class PickBodyRandom;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_ArgumentType>
    class CollectorCriteriaCombined;

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class CollectorCriteriaInterface;

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class CollectorCriteriaWrapper;

    template< typename t_ReturnType, typename t_ArgumentType>
    class CollectorInterface;

    template< typename t_ReturnType, typename t_ArgumentType, typename t_IntermediateResultType>
    class Locator;

    template< typename t_Argument>
    class LocatorCoordinatesAverage;

    template< typename t_Argument>
    class LocatorCoordinatesInterface;

    template< typename t_ArgumentType>
    class LocatorCoordinatesKnown;

    template< typename t_ArgumentType>
    class LocatorCoordinatesTetrahedral;

    template< typename t_ArgumentType>
    class LocatorCoordinatesTrigonal;

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType, typename t_IntermediateResultType>
    class LocatorCriteria;

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class LocatorCriteriaInterface;

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class LocatorCriteriaWrapper;

    template< typename t_ReturnType, typename t_ArgumentType>
    class LocatorInterface;

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class PickCriteriaInterface;

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class PickCriteriaWrapper;

    template< typename t_ReturnType, typename t_ArgumentType>
    class PickInterface;

  //////////////
  // typedefs //
  //////////////

  } // namespace find
} // namespace bcl

#endif // BCL_FIND_FWD_HH_
