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
#include "model/bcl_model_score_dataset_binary_operation.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_divide_equals.h"
#include "math/bcl_math_minus_equals.h"
#include "math/bcl_math_plus_equals.h"
#include "math/bcl_math_power_equals.h"
#include "math/bcl_math_times_equals.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
  //////////
  // data //
  //////////

    const util::SiPtr< const util::ObjectInterface> ScoreDatasetBinaryOperation::s_Add
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance
      (
        new ScoreDatasetBinaryOperation( math::PlusEquals< float>())
      )
    );

    const util::SiPtr< const util::ObjectInterface> ScoreDatasetBinaryOperation::s_Subtract
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance
      (
        new ScoreDatasetBinaryOperation( math::MinusEquals< float>())
      )
    );

    const util::SiPtr< const util::ObjectInterface> ScoreDatasetBinaryOperation::s_Divide
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance
      (
        new ScoreDatasetBinaryOperation( math::DivideEquals< float>())
      )
    );

    const util::SiPtr< const util::ObjectInterface> ScoreDatasetBinaryOperation::s_Multiply
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance
      (
        new ScoreDatasetBinaryOperation( math::TimesEquals< float>())
      )
    );

    const util::SiPtr< const util::ObjectInterface> ScoreDatasetBinaryOperation::s_Exponentiate
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance
      (
        new ScoreDatasetBinaryOperation( math::PowerEquals< float>())
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from an operation
    ScoreDatasetBinaryOperation::ScoreDatasetBinaryOperation
    (
      const math::AssignmentOperationInterface< float> &OPERATION
    ) :
      m_Scores( 2),
      m_Op( OPERATION)
    {
    }

    //! @brief virtual copy constructor
    ScoreDatasetBinaryOperation *ScoreDatasetBinaryOperation::Clone() const
    {
      return new ScoreDatasetBinaryOperation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ScoreDatasetBinaryOperation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &ScoreDatasetBinaryOperation::GetAlias() const
    {
      return m_Op->GetVerb();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief score a given dataset
    //! @param DATASET dataset of interest
    //! @return scores of the dataset
    linal::Vector< float> ScoreDatasetBinaryOperation::Score( const descriptor::Dataset &DATASET) const
    {
      // instantiate vector to store the result of the next property
      if( m_Scores.GetSize() < size_t( 2))
      {
        return linal::Vector< float>();
      }

      // ensure that the properties are defined
      if( !m_Scores( 0).IsDefined() || !m_Scores( 1).IsDefined())
      {
        return linal::Vector< float>();
      }

      // initialize the descriptor return vector with the first descriptor
      linal::Vector< float> result( m_Scores( 0)->Score( DATASET));

      // get a reference to the operation
      const math::AssignmentOperationInterface< float> &operation( *m_Op);

      for
      (
        storage::Vector< util::Implementation< ScoreDatasetInterface> >::const_iterator
          itr_scores( m_Scores.Begin() + 1), itr_scores_end( m_Scores.End());
        itr_scores != itr_scores_end;
        ++itr_scores
      )
      {
        if( !itr_scores->IsDefined())
        {
          return linal::Vector< float>();
        }

        // calculate the next property
        linal::Vector< float> rhs( ( *itr_scores)->Score( DATASET));

        // check for vectors of incorrect size -> indicates scores that were not calculated
        if( rhs.GetSize() != result.GetSize())
        {
          return linal::Vector< float>();
        }

        for // perform the operation on each element of the vector
        (
          linal::Vector< float>::iterator itr_result( result.Begin()), itr_rhs( rhs.Begin()), itr_rhs_end( rhs.End());
          itr_rhs != itr_rhs_end;
          ++itr_rhs, ++itr_result
        )
        {
          operation( *itr_result, *itr_rhs);
        }
      }

      // return vector of added properties.
      return result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreDatasetBinaryOperation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( m_Op->GetVerb() + " two scores");

      // determine whether the lhs and rhs are interchangeable
      if( m_Op.GetAlias() == "+" || m_Op.GetAlias() == "*")
      {
        parameters.AddInitializer
        (
          "",
          "Scores to " + m_Op->GetVerb(),
          io::Serialization::GetAgentWithSizeLimits( &m_Scores, size_t( 2)) // require at least two scores
        );
      }
      else
      {
        parameters.AddInitializer
        (
          "lhs",
          "argument for the left hand side of the operation",
          io::Serialization::GetAgent( &m_Scores( 0))
        );
        parameters.AddInitializer
        (
          "rhs",
          "argument for the right hand side of the operation",
          io::Serialization::GetAgent( &m_Scores( 1))
        );
      }

      return parameters;
    }

  } // namespace model
} // namespace bcl
