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
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_objective_function_constant.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionConstant::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance( new ObjectiveFunctionConstant())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionConstant::ObjectiveFunctionConstant() :
      m_ConstantValue( float( 0)),
      m_ImprovementType( opti::e_SmallerEqualIsBetter)
    {
    }

    //! @brief default constructor
    ObjectiveFunctionConstant::ObjectiveFunctionConstant
    (
      const float CONST_VALUE,
      opti::ImprovementType IMPROVE,
      GoalEnum GOAL,
      const float &CUTOFF
    ) :
      m_ConstantValue( CONST_VALUE),
      m_ImprovementType( IMPROVE),
      m_GoalType( GOAL),
      m_Threshold( CUTOFF)
    {
    }

    //! copy constructor
    ObjectiveFunctionConstant *ObjectiveFunctionConstant::Clone() const
    {
      return new ObjectiveFunctionConstant( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ObjectiveFunctionConstant::GetAlias() const
    {
      static const std::string s_Name( "Constant");
      return s_Name;
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionConstant::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // no data available, return constant result
      return m_ConstantValue;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionConstant::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Returns a constant value as objective function result; used to weight objective functions or provide an offset"
      );
      parameters.AddInitializer
      (
        "value",
        "constant value to return",
        io::Serialization::GetAgent( &m_ConstantValue),
        "3.40282e+38"
      );
      parameters.AddInitializer
      (
        "direction",
        "Determines whether an increasing or decreasing value indicates improvement",
        io::Serialization::GetAgent( &m_ImprovementType),
        "SmallerIsBetter"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
