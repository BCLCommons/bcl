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
#include "model/bcl_model_align_cutoff.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> AlignCutoff::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignCutoff())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AlignCutoff::AlignCutoff() :
      m_InternalModel(),
      m_ExperimentalCutoff(),
      m_AdjustedCutoff(),
      m_SlopeAboveCutoff(),
      m_SlopeBelowCutoff(),
      m_RescaleInput(),
      m_RescaleOutput()
    {
    }

    //! @brief constructor with a model, cutoff as parameters
    AlignCutoff::AlignCutoff
    (
      const util::ShPtr< Interface> &MODEL,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_OUTPUT,
      const float EXPERIMENTAL_CUTOFF,
      const float ADJUSTED_CUTOFF,
      const float SLOPE_ABOVE_CUTOFF,
      const float SLOPE_BELOW_CUTOFF
    ) :
      m_InternalModel( MODEL),
      m_ExperimentalCutoff( EXPERIMENTAL_CUTOFF),
      m_AdjustedCutoff( ADJUSTED_CUTOFF),
      m_SlopeAboveCutoff( SLOPE_ABOVE_CUTOFF),
      m_SlopeBelowCutoff( SLOPE_BELOW_CUTOFF),
      m_RescaleInput( RESCALE_INPUT),
      m_RescaleOutput( RESCALE_OUTPUT)
    {
    }

    //! copy constructor
    AlignCutoff *AlignCutoff::Clone() const
    {
      return new AlignCutoff( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AlignCutoff::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &AlignCutoff::GetAlias() const
    {
      static const std::string s_Name( "AlignCutoff");
      return s_Name;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get number outputs
    //! @return number of prediction outputs equivalent to the number of result columns
    size_t AlignCutoff::GetNumberOutputs() const
    {
      return m_InternalModel->GetNumberOutputs();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void AlignCutoff::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      if( !m_RescaleInput.IsDefined())
      {
        FEATURE.DeScale();
        return;
      }
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FEATURE.DeScale();
        FEATURE.Rescale( *m_RescaleInput);
      }
    }

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> AlignCutoff::PredictWithoutRescaling
    (
      const FeatureDataSetInterface< float> &FEATURE
    ) const
    {
      FeatureDataSet< float> predictions( m_InternalModel->PredictWithoutRescaling( FEATURE));

      // number of rows in predictions
      const size_t num_rows( predictions.GetNumberFeatures());

      // number cols in predictions
      const size_t num_cols( predictions.GetFeatureSize());

//      BCL_Message
//      (
//        util::Message::e_Verbose,
//        "m_AdjustedCutoff: " + util::Format()( m_AdjustedCutoff)
//        + " --- Slope above: " + util::Format()( m_SlopeBelowCutoff)
//        + " Slope Below: " + util::Format()( m_SlopeAboveCutoff)
//      );

      // get the raw matrix
      linal::Matrix< float> &classifications( predictions.GetRawMatrix());
      // subtract the adjusted cutoff
      classifications -= m_AdjustedCutoff;
      for( size_t row_id( 0); row_id < num_rows; ++row_id)
      {
        // get the row out of the matrix
        float *classification_row( classifications[ row_id]);
        for( size_t col_id( 0); col_id < num_cols; ++col_id)
        {
          // multiply by the appropriate slope
          classification_row[ col_id] *= classification_row[ col_id] > 0.0f ? m_SlopeAboveCutoff : m_SlopeBelowCutoff;
        }
      }
      // add in the experimental cutoff
      classifications += m_ExperimentalCutoff;

      return
        m_RescaleOutput.IsDefined()
        ? FeatureDataSet< float>( classifications, *m_RescaleOutput)
        : FeatureDataSet< float>( classifications);
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> AlignCutoff::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( m_RescaleInput.IsDefined())
      {
        if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
        {
          FeatureDataSet< float> feature( FEATURE);
          feature.Rescale( *m_RescaleInput);
          return PredictWithoutRescaling( feature).DeScale();
        }
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read AlignCutoff from std::istream
    std::istream &AlignCutoff::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_InternalModel, ISTREAM);
      io::Serialize::Read( m_ExperimentalCutoff, ISTREAM);
      io::Serialize::Read( m_AdjustedCutoff,     ISTREAM);
      io::Serialize::Read( m_SlopeAboveCutoff,   ISTREAM);
      io::Serialize::Read( m_SlopeBelowCutoff,   ISTREAM);
      io::Serialize::Read( m_RescaleInput, ISTREAM);
      io::Serialize::Read( m_RescaleOutput, ISTREAM);

      // end
      return ISTREAM;
    }

    //! write AlignCutoff into std::ostream
    std::ostream &AlignCutoff::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_InternalModel, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExperimentalCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AdjustedCutoff,     OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SlopeAboveCutoff,   OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SlopeBelowCutoff,   OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleInput, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleOutput, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
