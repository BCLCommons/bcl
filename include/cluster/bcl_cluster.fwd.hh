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

#ifndef BCL_CLUSTER_FWD_HH_
#define BCL_CLUSTER_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the cluster namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace cluster
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class InputFeatures;
    class NodeDescriptionFromFile;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class Dendrogram;

    template< typename t_DataType, typename t_PrecisionType>
    class DendrogramUncached;

    template< typename t_PrecisionType>
    class DistancesEuclidean;

    template< typename t_DataType, typename t_PrecisionType>
    class DistancesStored;

    template< typename t_DataType, typename t_PrecisionType>
    class InputClasses;

    template< typename t_DataType, typename t_PrecisionType>
    class InputInterface;

    template< typename t_PrecisionType>
    class InputPairwiseList;

    template< typename t_PrecisionType>
    class InputTable;

    template< typename t_PrecisionType>
    class KeyCompare;

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageAverage;

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageClasses;

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageComplete;

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageInterface;

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageSingle;

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageTotal;

    template< typename t_DataType, typename t_PrecisionType>
    class Node;

    template< typename t_DataType, typename t_PrecisionType>
    class NodeColorer;

    template< typename t_DataType, typename t_PrecisionType>
    class NodeDescriptionAverage;

    template< typename t_DataType, typename t_PrecisionType>
    class OutputCenters;

    template< typename t_DataType, typename t_PrecisionType>
    class OutputClasses;

    template< typename t_DataType, typename t_PrecisionType>
    class OutputInterface;

    template< typename t_DataType, typename t_PrecisionType>
    class OutputPymol;

    template< typename t_PrecisionType>
    class OutputPymolLabelProteinModelFromString;

    template< typename t_PrecisionType>
    class OutputPymolLabelSmallMolecule;

    template< typename t_PrecisionType>
    class OutputPymolLabelString;

    template< typename t_DataType, typename t_PrecisionType>
    class OutputRows;

    template< typename t_DataType, typename t_PrecisionType>
    class OutputSortedMatrix;

  //////////////
  // typedefs //
  //////////////

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_FWD_HH_
