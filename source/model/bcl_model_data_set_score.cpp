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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "model/bcl_model_data_set_score.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new DataSetScore
    DataSetScore *DataSetScore::Clone() const
    {
      return new DataSetScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get the features of the data set
    //! @return the features of the data set
    const FeatureLabelSet &DataSetScore::GetFeatures() const
    {
      return m_FeatureCodeLabel;
    }

    //! @brief Get the scores for the data set
    //! @return the scores for each feature of the data set
    const linal::Vector< float> &DataSetScore::GetScores() const
    {
      return m_Scores;
    }

    //! @brief Set the features of the data set
    //! @param FEATURES the new features
    void DataSetScore::SetFeatures( const FeatureLabelSet &FEATURES)
    {
      m_FeatureCodeLabel = FEATURES;
    }

    //! @brief Set the scores for each feature
    //! @param SCORES the new scores
    void DataSetScore::SetScores( const linal::Vector< float> &SCORES)
    {
      m_Scores = SCORES;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetScore::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_FeatureCodeLabel, ISTREAM);
      io::Serialize::Read( m_Scores, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_FeatureCodeLabel, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scores, OSTREAM, INDENT);
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
