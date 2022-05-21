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
#include "score/bcl_score_epr_accessibility.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_assignment.h"
#include "score/bcl_score_energy_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return the default filename for the default accessibility statistics
    const std::string &EPRAccessibility::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "SLnc-CBnc_pdbs07.histogram");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &EPRAccessibility::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "epr_accessibility");

      // end
      return s_default_scheme;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &EPRAccessibility::GetAlias() const
    {
      static const std::string s_name( "EPR Accessibility");
      return s_name;
    }

    //! @brief return parameters for member data that are not set up from the labels
    //! @return parameters for member data that are not set up from the labels
    io::Serializer EPRAccessibility::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores EPR accessibility");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to file where the statistics and in consequence the energy potentials are read from",
        io::Serialization::GetAgent( &m_HistogramFileName)
       );
      serializer.AddInitializer
      (
        "histogram lower bound",
        "the lower bounds of the bins of the histogram",
        io::Serialization::GetAgent( &m_HistogramLowerBound)
       );
      serializer.AddInitializer
      (
        "histogram upper bound",
        "the upper bounds of the bins of the histogram",
        io::Serialization::GetAgent( &m_HistogramUpperBound)
       );

      return serializer;
    }

    //////////////////////////////////<
    // construction and destruction //
    //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param SCHEME scheme to be used
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    EPRAccessibility::EPRAccessibility
    (
      const std::string &SCHEME,
      const std::string &HISTOGRAM_FILENAME
    ) :
      m_Scheme( SCHEME),
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_EnergyFunction(),
      m_HistogramLowerBound(),
      m_HistogramUpperBound()
    {
      // read the histogram file and store the energy functions
      ReadEnergyVector();
    }

    //! @brief virtual copy constructor
      EPRAccessibility *EPRAccessibility::Clone() const
    {
      return new EPRAccessibility( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EPRAccessibility::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that scores assignment of interest using pre-generated statistics
    //! @param ASSIGNMENT particular assignment of interest
    //! @return score
    double EPRAccessibility::operator()
    (
      const restraint::Assignment
      <
        storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>, double, biol::AABase
      > &ASSIGNMENT
    ) const
    {
      if( !ASSIGNMENT.GetGroupCollection().Begin()->second.Begin()->IsDefined())
      {
        return 0.0;
      }
      // create ShPtr "experimental_accessibilities" and initialize as pointing to the experimental restraints
      // in "ASSIGNMENT"
      util::ShPtr
      <
        storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>
      > experimental_accessibilities
      (
        ASSIGNMENT.GetRestraint()
      );

      // create double "score_sum" which will hold the sum of the agreement of the AABase in "ASSIGNMENT" with each of
      // the experimental measurements
      double score_sum( 0.0);

      // iterate through the experimental data to see how well each measurement agrees with the neighbor count of the
      // AABase in "ASSIGNMENT"
      for
      (
        storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>::const_iterator
          exp_data_itr( experimental_accessibilities->Begin()), exp_data_itr_end( experimental_accessibilities->End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr
      )
      {
        // create const double "current_accessibility_value" and initialize with the experimental accessibility
        // currently denoted by "exp_data_itr"
        const double current_accessibility_value( exp_data_itr->second);

        // create const double "neighbor_measure" and initialize with the neighbor measure that was calculated for the
        // AABase in "ASSIGNMENT"
        const double neighbor_measure( ASSIGNMENT.GetGroupCollection().Begin()->first);

        const double slnc_cbnc( current_accessibility_value - neighbor_measure);

        double current_score;

        // true if "sl_cb" is outside the boundaries of the SL-CB histogram
        if
        (
          slnc_cbnc < m_HistogramLowerBound || //< smaller than any value
          slnc_cbnc > m_HistogramUpperBound   //< larger than any value
        )
        {
          current_score = 0.0;
        }
        else
        {
          // get the score associated with having "neighbor_measure" calculated for AABase given an experimentally
          // measured accessibility of "current_accessibility_value"
          current_score = m_EnergyFunction( slnc_cbnc);
        }
        BCL_MessageDbg
        (
          "current score for  " +
          ( *ASSIGNMENT.GetGroupCollection().Begin()->second.Begin())->GetIdentification() + " with environment type " +
          exp_data_itr->first.GetString() + " with value " + util::Format()( current_accessibility_value) +
          " and neighbor measure value of " + util::Format()( neighbor_measure) + " is " +
          util::Format()( current_score)
        );
        score_sum += current_score;
      }

      return score_sum;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param ASSIGNMENT
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &EPRAccessibility::WriteDetailedSchemeAndValues
    (
      const restraint::Assignment
      <
        storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>, double, biol::AABase
      > &ASSIGNMENT,
      std::ostream &OSTREAM
    ) const
    {
      if( !ASSIGNMENT.GetGroupCollection().Begin()->second.Begin()->IsDefined())
      {
        return OSTREAM;
      }

      // create ShPtr "experimental_accessibilities" and initialize as pointing to the experimental restraints
      // in "ASSIGNMENT"
      util::ShPtr
      <
        storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>
      > experimental_accessibilities
      (
        ASSIGNMENT.GetRestraint()
      );

      // iterate through the experimental data to see how well each measurement agrees with the neighbor count of the
      // AABase in "ASSIGNMENT"
      for
      (
        storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>::const_iterator
          exp_data_itr( experimental_accessibilities->Begin()), exp_data_itr_end( experimental_accessibilities->End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr
      )
      {
        // create const double "current_accessibility_value" and initialize with the experimental accessibility
        // currently denoted by "exp_data_itr"
        const double current_accessibility_value( exp_data_itr->second);

        // create const double "neighbor_measure" and initialize with the neighbor measure that was calculated for the
        // AABase in "ASSIGNMENT"
        const double neighbor_measure( ASSIGNMENT.GetGroupCollection().Begin()->first);

        const double slnc_cbnc( current_accessibility_value - neighbor_measure);

        double current_score;

        // true if "sl_cb" is outside the boundaries of the SL-CB histogram
        if
        (
          slnc_cbnc < m_HistogramLowerBound || //< smaller than any value
          slnc_cbnc > m_HistogramUpperBound   //< larger than any value
        )
        {
          current_score = 0.0;
        }
        else
        {
          // get the score associated with having "neighbor_measure" calculated for AABase given an experimentally
          // measured accessibility of "current_accessibility_value"
          current_score = m_EnergyFunction( slnc_cbnc);
        }

        OSTREAM << "current score for  " <<
          ( *ASSIGNMENT.GetGroupCollection().Begin()->second.Begin())->GetIdentification() << " with environment type "
          << exp_data_itr->first.GetString() << " with value " << current_accessibility_value <<
          " and neighbor measure value of " << neighbor_measure << " is " << current_score << '\n';
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
    bool EPRAccessibility::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      ReadEnergyVector();
      return true;
    }

    //! @brief read the energy distribution for scoring EPR accessibilities
    void EPRAccessibility::ReadEnergyVector()
    {
      // initialize read
      io::IFStream read;

      // bind read to the file "m_HistogramFileName"
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      // create math::Histogram named "accessibility_histogram" which will be used to hold the accessibility histogram
      math::Histogram accessibility_histogram;

      // read in the data from "m_HistogramFileName" into "accessibility_histogram"
      read >> accessibility_histogram;

      // close and clear read stream
      io::File::CloseClearFStream( read);

      const storage::VectorND< 2, double> boundaries( accessibility_histogram.GetBoundaries());

      m_HistogramLowerBound = boundaries.First() + 0.5 * accessibility_histogram.GetBinSize();
      m_HistogramUpperBound = boundaries.Second() - 0.5 * accessibility_histogram.GetBinSize();

      // set "m_EnergyFunction" to the cubic spline that has been trained on "accessibility_histogram"
      m_EnergyFunction = EnergyDistribution::GeneratePotentialFromHistogram
      (
        accessibility_histogram, 1.0, math::e_Natural, EnergyDistribution::GetDefaultFirstDerivative(),
        true, true
      );
    }

  } // namespace score
} // namespace bcl
