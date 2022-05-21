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

#ifndef BCL_MODEL_ALIGN_CUTOFF_H_
#define BCL_MODEL_ALIGN_CUTOFF_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignCutoff
    //! @brief Wrapper model that alters threshold and slope on each side of that threshold
    //! @details Often used for combining outputs, since outputs from regularly trained models have different thresholds
    //!          at which they are most precise
    //!
    //! @see @link example_model_align_cutoff.cpp @endlink
    //! @author butkiem1, mendenjl
    //! @date Jul 07, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AlignCutoff :
      public Interface
    {
    public:

    //////////
    // data //
    //////////

      //! internal model used for training
      util::ShPtr< Interface> m_InternalModel;

      //! original cutoff for experimental data
      float m_ExperimentalCutoff;

      //! model adjusted prediction cutoff
      float m_AdjustedCutoff;

      //! Slope of values above the cutoff
      float m_SlopeAboveCutoff;

      //! Slope of values below the cutoff
      float m_SlopeBelowCutoff;

      //! rescale input to transfer function input range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleInput;

      //! rescale output from transfer function output range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleOutput;

    private:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AlignCutoff();

      //! @brief constructor with a model and cutoff as parameters
      AlignCutoff
      (
        const util::ShPtr< Interface> &MODEL,
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT,
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE_OUTPUT,
        const float EXPERIMENTAL_CUTOFF,
        const float ADJUSTED_CUTOFF,
        const float SLOPE_ABOVE_CUTOFF,
        const float SLOPE_BELOW_CUTOFF
      );

      //! copy constructor
      AlignCutoff *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! get number output neurons
      size_t GetNumberOutputs() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( FeatureDataSet< float> &FEATURE) const;

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> PredictWithoutRescaling
      (
        const FeatureDataSetInterface< float> &FEATURE
      ) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read AlignCutoff from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write AlignCutoff into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class AlignCutoff

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_ALIGN_CUTOFF_H_
