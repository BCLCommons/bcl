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

#ifndef BCL_MODEL_DTREE_ROC_DATA_PARTITION_FUNCTION_H_
#define BCL_MODEL_DTREE_ROC_DATA_PARTITION_FUNCTION_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_dtree_data_partition_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DtreeRocDataPartitionFunction
    //! @brief is derived from DtreeDataPartitionFunctionInterface and implements
    //!        a specific way to partition a dataset of a corresponding decision tree using
    //!        a Roc-tree splitting approach. It uses a ROC curve for every descriptor value range to rank them.
    //!        The best pick is used for partitioning. See ROC-tree paper (Hossain 2006).
    //!        AUC of Roc curve is determined by distance from diagonal since only destinguishing of two groups is
    //!        required instead is specifying actives and inactives. To prevent overfitting the AUC is limited to a
    //!        range of 0.5 to 0.95.
    //!
    //! @see @link example_model_dtree_roc_data_partition_function.cpp @endlink
    //! @author teixeipl, lemmonwa
    //! @date 07/17/2009; 07/02/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DtreeRocDataPartitionFunction :
      public DtreeDataPartitionFunctionInterface
    {

    //////////
    // data //
    //////////

    private:

      //! minimum area under the ROC curve (AuROC) before termination of iterations. default value is 0.5
      float m_MinimumAuROC;

      //! maximum area under the ROC curve (AuROC) before termination of iterations, based on Hossain paper to avoid
      //! overfitting. See ROC-tree paper (Hossain 2006). Default value is 0.95.
      float m_MaximumAuROC;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DtreeRocDataPartitionFunction();

      //! @brief constructor taking all parameters
      //! @param MAX maximum of allowed ROC AUC
      //! @param MIN minimum of allowed ROC AUC
      DtreeRocDataPartitionFunction( const float MAX, const float MIN);

      //! @brief copy constructor
      DtreeRocDataPartitionFunction *Clone() const;

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

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief rates by misclassification
      //! @param STATE_COUNTS counts of each state (e.g. true, false, or 0,1,2, etc.) on one segment of the partition
      //! @param SIZE total number of items in the partition
      //! @return the # of correct classifications on this side
      float RatePartition
      (
        const linal::Vector< size_t> &STATE_COUNTS,
        const size_t &SIZE
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class DtreeRocDataPartitionFunction

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_DTREE_ROC_DATA_PARTITION_FUNCTION_H_
