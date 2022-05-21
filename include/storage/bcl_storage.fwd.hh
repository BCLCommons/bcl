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

#ifndef BCL_STORAGE_FWD_HH_
#define BCL_STORAGE_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// external includes - sorted alphabetically
#include <cstddef>
#include <functional>

// This file contains forward declarations for the storage namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace storage
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class GreaterThanThird;
    class LessThanFirst;
    class LessThanSecond;
    class LessThanThird;
    class TableHeader;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_KeyType, typename t_DataType>
    class HashMap;

    template< typename t_DataType>
    class List;

    template< typename t_KeyType, typename t_DataType, typename t_KeyCompare = std::less< t_KeyType> >
    class Map;

    template< typename t_KeyType, typename t_DataType, typename t_KeyCompare = std::less< t_KeyType> >
    class MultiMap;

    template< typename t_KeyType, typename t_KeyCompare = std::less< t_KeyType> >
    class MultiSet;

    template< size_t t_NumberOfObjects, typename t_ObjectType, typename t_DataType>
    class ObjectNDHashMap;

    template< typename t_First, typename t_Second>
    class Pair;

    template< typename t_First, typename t_Second>
    class PairBinaryPredicateFirst;

    template< typename t_First, typename t_Second>
    class PairBinaryPredicateSecond;

    template< typename t_First>
    class PairEqualFirst;

    template< typename t_DataType>
    class Row;

    template< typename t_DataType>
    struct RowComparison;

    template< typename t_KeyType, typename t_KeyCompare = std::less< t_KeyType> >
    class Set;

    template< typename t_DataType>
    class SymmetricMatrix;

    template< typename t_DataType>
    class Table;

    template< typename t_First, typename t_Second, typename t_Third>
    class Triplet;

    template< typename t_DataType>
    class Vector;

    template< unsigned int N, typename t_DataType>
    class VectorND;

  //////////////
  // typedefs //
  //////////////

  } // namespace storage
} // namespace bcl

#endif // BCL_STORAGE_FWD_HH_
