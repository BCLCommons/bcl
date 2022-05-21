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
#include "model/bcl_model_objective_function_cutoff_from_percentile.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average_sd.h"
#include "model/bcl_model_data_set_select_columns.h"
#include "model/bcl_model_feature_label_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionCutoffFromPercentile::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionCutoffFromPercentile()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! copy constructor
    ObjectiveFunctionCutoffFromPercentile *ObjectiveFunctionCutoffFromPercentile::Clone() const
    {
      return new ObjectiveFunctionCutoffFromPercentile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief default constructor
    ObjectiveFunctionCutoffFromPercentile::ObjectiveFunctionCutoffFromPercentile() :
      m_Percentile( 0.5)
    {
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ObjectiveFunctionCutoffFromPercentile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief determine what sign of the derivative of this objective function indicates improvement
    //! @return the sign of the derivative of this objective function indicates improvement
    opti::ImprovementType ObjectiveFunctionCutoffFromPercentile::GetImprovementType() const
    {
      if( m_Objective.IsDefined())
      {
        return m_Objective->GetImprovementType();
      }
      return opti::e_LargerIsBetter;
    }

    //! @brief get the overall goal of the objective function
    //! @return the goal of the objective function
    ObjectiveFunctionInterface::Goal ObjectiveFunctionCutoffFromPercentile::GetGoalType() const
    {
      if( m_Objective.IsDefined())
      {
        return m_Objective->GetGoalType();
      }
      return e_Other;
    }

    //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
    //!        case indicates whether the prediction was true (upper) or false (lower)
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
    linal::Matrix< char> ObjectiveFunctionCutoffFromPercentile::GetFeaturePredictionClassifications
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      return m_Objective->GetFeaturePredictionClassifications( EXPERIMENTAL, PREDICTED);
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionCutoffFromPercentile::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      storage::Vector< float> data_copy( DATA.GetMatrix().Begin(), DATA.GetMatrix().End());
      if( m_Objective->GetRankingParity())
      {
        data_copy.Sort( std::greater< float>());
      }
      else
      {
        data_copy.Sort( std::less< float>());
      }
      float threshold
      (
        data_copy( std::min( size_t( m_Percentile * data_copy.GetSize()), data_copy.GetSize() - size_t( 1)))
      );
      auto scaling_ptr( DATA.GetScaling());
      if( scaling_ptr.IsDefined())
      {
        threshold = scaling_ptr->DescaleValue( 0, threshold);
      }
      m_Objective->SetThreshold( threshold);
      m_Objective->SetData( DATA, IDS);
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionCutoffFromPercentile::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      return m_Objective->operator()( EXPERIMENTAL, PREDICTED);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionCutoffFromPercentile::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Sets the threshold/cutoff/division between classes of a user-specified objective function based on a percentile "
        "of the data"
      );

      parameters.AddInitializer
      (
        "percentile",
        "percentile of the data to consider in the positive class for setting the cutoff",
        io::Serialization::GetAgentWithRange( &m_Percentile, 0.0, 1.0),
        "0.5"
      );
      parameters.AddInitializer
      (
        "",
        "Objective function to set the threshold of",
        io::Serialization::GetAgent( &m_Objective)
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
