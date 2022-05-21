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

#ifndef BCL_MODEL_LEVERAGE_MATRIX_H_
#define BCL_MODEL_LEVERAGE_MATRIX_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_interface.h"
#include "bcl_model_rescale_feature_data_set.h"
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LeverageMatrix
    //! @brief is for creating a linear regression model for a set of data and using the
    //!        determined weights to calculate values based on the model
    //!
    //! @see @link example_model_leverage_matrix.cpp @endlink
    //! @author mendenjl
    //! @date May 06, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LeverageMatrix :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      linal::Matrix< float> m_Weights; //!< Weights determined by iterate
      util::ShPtr< RescaleFeatureDataSet> m_Rescaler; //!< Rescale function

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LeverageMatrix();

      //! @brief parameter constructor
      //! @param WEIGHTS weights used in the optimization
      LeverageMatrix
      (
        const linal::MatrixConstInterface< float> &WEIGHTS,
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE = util::ShPtr< RescaleFeatureDataSet>()
      );

      //! @brief Clone is the virtual copy constructor
      LeverageMatrix *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get the output feature size for this model
      //! @return the output feature size for this model
      size_t GetNumberOutputs() const;

      //! @brief get the weights used for this model
      //! @return the weights used for this model
      const linal::Matrix< float> &GetWeights() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( FeatureDataSet< float> &FEATURE) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class LeverageMatrix

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_LEVERAGE_MATRIX_H_
