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
#include "score/bcl_score_pofr.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_statistics.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_distance_density_point.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PofR::s_Instance
    (
      GetObjectInstances().AddInstance( new PofR())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Constructor
    PofR::PofR()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PofR
    PofR *PofR::Clone() const
    {
      return new PofR( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PofR::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief overloaded () operator to calculate Derivative Score or Cumulative Integral Score from two SAXS curves
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return return Derivative Score for two SAS curves
    double PofR::operator()( const restraint::SasExperimentalAndCalculatedDensity &SAS_DATA) const
    {
      // The score is initialized to zero
      double score( 0.0);

      // Compare dmax if model exceeds experimental boundary, reject
      double cal_dmax( SAS_DATA.GetCalculatedDensity().GetDmax());
      double exp_dmax( SAS_DATA.GetExperimentalDensity().GetDmax());

      double dmax_difference( exp_dmax - cal_dmax);

      // if the Max Dimension is violated, return a score of 10000
      if( dmax_difference < 0)
      {
        score = 10000.0;
        return score;
      }

      // Compare auc of both profiles, if model exceeds experimental auc, reject
      double auc_difference( restraint::SasAnalysis::CalculatePofRIntegralScore( SAS_DATA));

      if( auc_difference < 0)
      {
        score = 10000.0;
        return score;
      }

      // Compute Excess Integral score
       double excess_integral( restraint::SasAnalysis::CalculatePofRExcessIntegralScore( SAS_DATA));

      // Compute Oscillation score
      double cal_oscillation( restraint::SasAnalysis::CalculatePofROscillationScore( SAS_DATA.GetCalculatedDensity()));
      double exp_oscillation( restraint::SasAnalysis::CalculatePofROscillationScore( SAS_DATA.GetExperimentalDensity()));

      double oscillation_score( math::Absolute( exp_oscillation - cal_oscillation));

      if( excess_integral > 0.1)
      {
        excess_integral = 10 * excess_integral;
      }

      score = ( dmax_difference) + ( auc_difference) + ( 10 * excess_integral) + ( 10 * oscillation_score);
      //score = ( dmax_difference + auc_difference);
      BCL_MessageStd( "score: " + util::Format()( score));

      // if no violations occur, return a score
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PofR::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PofR::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {

      // return the stream
      return OSTREAM;
    }

  } // namespace score

} // namespace bcl
