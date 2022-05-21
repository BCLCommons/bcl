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

#ifndef BCL_CLUSTER_OUTPUT_INTERFACE_H_
#define BCL_CLUSTER_OUTPUT_INTERFACE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically
#include "io/bcl_io.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputInterface
    //! @brief is the interface from which classes that are used to output the results of clustering
    //!        should be derived
    //!
    //! @remarks example unnecessary
    //! @author alexanns
    //! @date October 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class OutputInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief WriteOutput writes some data of a list of Nodes to some output stream
      //! @param FILENAME
      //! @param NODE_LIST the list of nodes which is going to be output
      //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
      virtual void WriteOutput
      (
        const std::string &FILENAME, const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST
      ) const = 0;

    }; // OutputInterface

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_OUTPUT_INTERFACE_H_ 
