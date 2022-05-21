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

#ifndef BCL_MODEL_SCORE_DATASET_NON_REDUNDANT_H
#define BCL_MODEL_SCORE_DATASET_NON_REDUNDANT_H

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_score_dataset_interface.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreDatasetNonRedundant
    //! @brief Returns 1 for each descriptor that is not a scaling/shifting of another previously seen descriptor, 0
    //!        otherwise
    //!
    //!
    //! @author mendenjl
    //! @see @link example_model_score_dataset_non_redundant.cpp @endlink
    //! @date Jan 08, 2015
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreDatasetNonRedundant :
      public ScoreDatasetInterface
    {

    //////////
    // data //
    //////////

      //! Max Z-Score distance for all pairs of two descriptor values can be and yet still be considered redundant
      float m_ZScoreTolerance;

      //! Maximum number of outliers to the z-score tolerance before declaring the two descriptor values unique
      size_t m_MaxOutliers;

      //! Whether to allow up to one constant, non-zero, descriptor in a dataset (useful to allow an offset for linear regression)
      bool   m_AllowOneConstant;

      //! Minimum range between the highest and lowest value seen for a descriptor to consider it noteworthy
      float m_MinSpan;

      //! Minimum (absolute) standard deviation for the descriptor to consider it relevant
      float m_MinStd;

      //! Minimum absolute cross-correlation for removal
      float m_MinAbsCrossCorrelation;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ScoreDatasetNonRedundant();

      //! @brief Clone function
      //! @return pointer to new ScoreDatasetNonRedundant
      ScoreDatasetNonRedundant *Clone() const;

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

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ScoreDatasetNonRedundant

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SCORE_DATASET_NON_REDUNDANT_H

