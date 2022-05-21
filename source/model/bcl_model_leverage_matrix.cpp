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
#include "linal/bcl_linal_vector_operations.h"
#include "model/bcl_model_leverage_matrix.h"

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
    const util::SiPtr< const util::ObjectInterface> LeverageMatrix::s_Instance
    (
      GetObjectInstances().AddInstance( new LeverageMatrix())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LeverageMatrix::LeverageMatrix()
    {
    }

    //! @brief parameter constructor
    //! @param WEIGHTS weights used in the optimization
    LeverageMatrix::LeverageMatrix
    (
      const linal::MatrixConstInterface< float> &WEIGHTS,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE
    ) :
      m_Weights( WEIGHTS),
      m_Rescaler( RESCALE)
    {
    }

    //! @brief Clone is the virtual copy constructor
    LeverageMatrix *LeverageMatrix::Clone() const
    {
      return new LeverageMatrix( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LeverageMatrix::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &LeverageMatrix::GetAlias() const
    {
      static const std::string s_alias( "Leverage");
      return s_alias;
    }

    //! @brief get the output feature size for this model
    //! @return the output feature size for this model
    size_t LeverageMatrix::GetNumberOutputs() const
    {
      return 1;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void LeverageMatrix::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      // no rescaling is ever performed by linear regression, so just remove any rescaling on the dataset
      if( m_Rescaler.IsDefined())
      {
        FEATURE.Rescale( *m_Rescaler);
      }
      else
      {
        FEATURE.DeScale();
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> LeverageMatrix::PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURE) const
    {
      linal::Matrix< float> hat_product( FEATURE.GetMatrix() * m_Weights);
      linal::Matrix< float> result( FEATURE.GetMatrix().GetNumberRows(), size_t( 1), float( 0.0));
      for( size_t i( 0), n_rows( result.GetNumberRows()); i < n_rows; ++i)
      {
        result( i, 0) = linal::ScalarProduct( hat_product.GetRow( i), FEATURE.GetMatrix().GetRow( i));
      }
      return result;
    }

    //! @brief get the weights used for this model
    //! @return the weights used for this model
    const linal::Matrix< float> &LeverageMatrix::GetWeights() const
    {
      return m_Weights;
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> LeverageMatrix::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( m_Rescaler.IsDefined() && ( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_Rescaler))
      {
        FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( *m_Rescaler);
        return PredictWithoutRescaling( feature);
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LeverageMatrix::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Weights, ISTREAM);
      io::Serialize::Read( m_Rescaler, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &LeverageMatrix::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Weights, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Rescaler, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LeverageMatrix::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes the leverage of given features with respect to the original training dataset"
      );
      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace model
} // namespace bcl

