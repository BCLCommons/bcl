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

#ifndef BCL_MODEL_DTREE_INFORMATION_GAIN_DATA_PARTITION_FUNCTION_H_
#define BCL_MODEL_DTREE_INFORMATION_GAIN_DATA_PARTITION_FUNCTION_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_dtree_data_partition_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DtreeInformationGainDataPartitionFunction
    //! @brief is derived from DtreeDataPartitionFunctionInterface and
    //!        implements a specific way to partition a dataset of a corresponding decision tree using
    //!        an InformationGain splitting approach.
    //!
    //! @see @link example_model_dtree_information_gain_data_partition_function.cpp @endlink
    //! @author lemmonwa
    //! @date 07/02/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DtreeInformationGainDataPartitionFunction :
      public DtreeDataPartitionFunctionInterface
    {

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! copy constructor
      DtreeInformationGainDataPartitionFunction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator determines the best partition of the data
      //! @param DATA vector of feature result references to consider
      //! @return function returns how to partition to the data
      DtreeBinaryPartition operator()( const storage::Vector< FeatureResultAndState> &DATA) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief calculates information gain from counts of states in each segment
      //! @param STATE_COUNTS counts of each state (e.g. true, false, or 0,1,2, etc.) on one segment of the partition
      //! @param SIZE total number of items in the partition
      //! @return the information gain metric (unnormalized)
      float RatePartition
      (
        const linal::Vector< size_t> &STATE_COUNTS,
        const size_t &SIZE
      ) const;

    }; // class DtreeInformationGainDataPartitionFunction

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_DTREE_INFORMATION_GAIN_DATA_PARTITION_FUNCTION_H_
