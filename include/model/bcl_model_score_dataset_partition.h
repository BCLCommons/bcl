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

#ifndef BCL_MODEL_SCORE_DATASET_PARTITION_H
#define BCL_MODEL_SCORE_DATASET_PARTITION_H

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_dtree_data_partition_function_interface.h"
#include "bcl_model_score_dataset_interface.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreDatasetPartition
    //! @brief calculate f-score values for every column in combined active and inactive dataset
    //!
    //!                                     (X_avg_pos - X_avg)^2 + (X_avg_neg - X_avg)^2
    //! f(i)= -------------------------------------------------------------------------------------------------------
    //!       1/(n_pos - 1) SUM k=1 to n_pos (X_k_pos - X_avg)^2 + 1/(n_neg - 1) SUM k=1 to n_neg (X_k_neg - X_avg)^2
    //!
    //!
    //! @author mendenjl
    //! @see @link example_model_score_dataset_partition.cpp @endlink
    //! @date Apr 06, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreDatasetPartition :
      public ScoreDatasetInterface
    {
    private:

    //////////
    // data //
    //////////

      util::Implementation< DtreeDataPartitionFunctionInterface> m_Partitioner;
      float m_Cutoff; //!< Cutoff between actives/inactives

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ScoreDatasetPartition
      ScoreDatasetPartition *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief score a given dataset
      //! @param DATASET dataset of interest
      //! @return scores of the dataset
      linal::Vector< float> Score( const descriptor::Dataset &DATASET) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ScoreDatasetPartition

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SCORE_DATASET_PARTITION_H

