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
#include "model/bcl_model_multiple_linear_regression.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MultipleLinearRegression::s_Instance
    (
      GetObjectInstances().AddInstance( new MultipleLinearRegression())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MultipleLinearRegression::MultipleLinearRegression()
    {
    }

    //! @brief parameter constructor
    //! @param WEIGHTS weights used in the optimization
    MultipleLinearRegression::MultipleLinearRegression( const linal::MatrixConstInterface< float> &WEIGHTS) :
      m_Weights( WEIGHTS)
    {
    }

    //! @brief Clone is the virtual copy constructor
    MultipleLinearRegression *MultipleLinearRegression::Clone() const
    {
      return new MultipleLinearRegression( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MultipleLinearRegression::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MultipleLinearRegression::GetAlias() const
    {
      static const std::string s_alias( "Regression");
      return s_alias;
    }

    //! @brief get the output feature size for this model
    //! @return the output feature size for this model
    size_t MultipleLinearRegression::GetNumberOutputs() const
    {
      return m_Weights.GetNumberCols();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void MultipleLinearRegression::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      // no rescaling is ever performed by linear regression, so just remove any rescaling on the dataset
      FEATURE.DeScale();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> MultipleLinearRegression::PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURE) const
    {
      BCL_Assert
      (
        FEATURE.GetFeatureSize() == m_Weights.GetNumberRows(),
        "Weights are of wrong dimensions (" + util::Format()( m_Weights.GetNumberRows())
        + ") for the given feature set (" + util::Format()( FEATURE.GetFeatureSize()) + " rows)"
      );
      return FeatureDataSet< float>( FEATURE.GetMatrix() * m_Weights);
    }

    //! @brief get the weights used for this model
    //! @return the weights used for this model
    const linal::Matrix< float> &MultipleLinearRegression::GetWeights() const
    {
      return m_Weights;
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> MultipleLinearRegression::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      return PredictWithoutRescaling( FEATURE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MultipleLinearRegression::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Weights, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MultipleLinearRegression::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Weights, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MultipleLinearRegression::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Multiplies a given feature matrix by another matrix to determine output"
      );
      // TODO: Add data member
      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace model
} // namespace bcl

