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

#ifndef BCL_MODEL_DATA_SET_SCORE_H_
#define BCL_MODEL_DATA_SET_SCORE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_label_set.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetScore
    //! @brief interface to retrieve data sets independent of source
    //!
    //! @see @link example_model_data_set_score.cpp @endlink
    //! @author mendenjl
    //! @date Apr 02, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetScore :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      FeatureLabelSet       m_FeatureCodeLabel; //!< FeatureLabels considered by score
      linal::Vector< float> m_Scores;           //!< Scores for each feature or eigenvector, if known

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new DataSetScore
      virtual DataSetScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief Get the features of the data set
      //! @return the features of the data set
      const FeatureLabelSet &GetFeatures() const;

      //! @brief Get the scores for the data set
      //! @return the scores for each feature of the data set
      const linal::Vector< float> &GetScores() const;

      //! @brief Set the features of the data set
      //! @param FEATURES the new features
      void SetFeatures( const FeatureLabelSet &FEATURES);

      //! @brief Set the scores for each feature
      //! @param SCORES the new scores
      void SetScores( const linal::Vector< float> &SCORES);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class DataSetScore

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_DATA_SET_SCORE_H_ 
