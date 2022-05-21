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
#include "score/bcl_score_sas_type.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_statistics.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    //! @brief ScoreFunction as string
    //! @param SCORE_FUNCTION the ScoreFunction
    //! @return the string for the ScoreFunction
    const std::string &SasType::GetFunctionDescriptor( const ScoreFunction &SCORE_FUNCTION)
    {
      static const std::string s_descriptors[] =
      {
        "chi",
        "cumulative",
        "stovgaard",
        GetStaticClassName< ScoreFunction>()
      };

      return s_descriptors[ SCORE_FUNCTION];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasType::s_Instance
    (
      GetObjectInstances().AddInstance( new SasType())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Constructor
    SasType::SasType() : m_UseErrors( false), m_ScoreType( e_chi)
    {
    }

    //! @brief Alternative Constructor
    SasType::SasType( const bool &USE_ERRORS, const ScoreFunctionEnum &SCORE_TYPE)
    :
        m_UseErrors( USE_ERRORS),
        m_ScoreType( SCORE_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasType
    SasType *SasType::Clone() const
    {
      return new SasType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to calculate the chi score for different SAXS intensity curves
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return derivative chi score
    double SasType::CalculateChiScore( const restraint::SasExperimentalAndCalculatedData &SAXS_DATA) const
    {
      return SAXS_DATA.ComputeScoringFunction( m_UseErrors);
    }

    //! @brief function to calculate cumulative integral score for different SAXS intensity curves
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return cumulative integral score
    double SasType::CalculateCumulativeIntegralScore
    (
      const restraint::SasExperimentalAndCalculatedData &SAXS_DATA
    ) const
    {
      // Convert the experimental intensity to list of doubles
      storage::List< double> experimental_intensity
      (
        restraint::SasAnalysis::ConvertIntensityDataToList( SAXS_DATA.GetExperimentalData())
      );

      // Convert the calculated intensity to list of doubles
      storage::List< double> calculated_intensity
      (
        restraint::SasAnalysis::ConvertIntensityDataToList( SAXS_DATA.GetCalculatedData())
      );

      // Calculate the Score
      double cumulative_score
      (
        math::Statistics::CumulativeEuclidian
        (
          experimental_intensity.Begin(),
          experimental_intensity.End(),
          calculated_intensity.Begin(),
          calculated_intensity.End()
        )
      );
      return cumulative_score;
    }

    //! @brief function to calculate Stovgaard score for different SAXS intensity curves
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return Stovgaard score
    double SasType::CalculateStovgaardScore( const restraint::SasExperimentalAndCalculatedData &SAXS_DATA) const
    {
      // initialize sum
      double summation( 0.0), alpha( 0.15), beta( 0.3);
      size_t qbins( 1);

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< restraint::SasScatteringPoint>::const_iterator
          exp_data_itr( SAXS_DATA.GetExperimentalData().Begin()),
          cal_data_itr( SAXS_DATA.GetCalculatedData().Begin()),
          exp_data_itr_end( SAXS_DATA.GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // if both I values are defined
        double sigma = exp_data_itr->GetIntensity() * ( exp_data_itr->GetQvalue() + alpha) * beta;

        // add to the sum
        summation += math::Sqr( ( exp_data_itr->GetIntensity() - cal_data_itr->GetIntensity()) / sigma);
        ++qbins;
      }

      return qbins == 1 ? 0.0 : math::Sqrt( summation / ( qbins - 1));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief overloaded () operator to calculate Derivative Score or Cumulative Integral Score from two SAXS curves
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return return Derivative Score for two SAXS curves
    double SasType::operator()( const restraint::SasExperimentalAndCalculatedData &SAXS_DATA) const
    {

      double score( 0.0);

      switch( m_ScoreType)
      {
        case e_chi:
          score = CalculateChiScore( SAXS_DATA);
          break;
        case e_cumulative:
          score = CalculateCumulativeIntegralScore( SAXS_DATA);
          break;
        case e_stovgaard:
          score = CalculateStovgaardScore( SAXS_DATA);
          break;
        default:
          BCL_Assert( false, "Unknown scoring type - the accepted types are chi, cumulative, and stovgaard");
          break;

      }
      return score;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &SasType::GetAlias() const
    {
      static const std::string s_name( "SasType");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SasType::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription( "parameters for small angle x-ray/neutron scattering");
      serial.AddInitializer
      (
        "consider errors",
        "true to consider the amount of error in each measurement when computing rmsd",
        io::Serialization::GetAgent( &m_UseErrors),
        "False"
      );
      serial.AddInitializer
      (
        "score",
        "how to score differences between the profiles",
        io::Serialization::GetAgent( &m_ScoreType),
        GetFunctionDescriptor( SasType().GetScoreFunction())
      );
      return serial;
    }

  } // namespace score

} // namespace bcl
