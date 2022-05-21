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
#include "score/bcl_score_restraint_distance_spin_label.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "score/bcl_score_energy_distribution.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RestraintDistanceSpinLabel::s_Instance
    (
      util::Enumerated< RestraintAtomDistanceAssignment>::AddInstance( new RestraintDistanceSpinLabel())
    );

    // energy function with lower and upper bound
    storage::Triplet< math::CubicSplineDamped, double, double> RestraintDistanceSpinLabel::s_EnergyFunction;

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &RestraintDistanceSpinLabel::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "sl-cb_distances.histograms");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &RestraintDistanceSpinLabel::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "epr_distance");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual copy constructor
    RestraintDistanceSpinLabel *RestraintDistanceSpinLabel::Clone() const
    {
      return new RestraintDistanceSpinLabel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RestraintDistanceSpinLabel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() which takes an Assignment for calculating the agreement of the Assignment with the distance
    //! @param ASSIGNMENT the assignment which contains the MovableInterfaces whose distance will be scored
    //! @return return a double which is the score of the agreement of the MovableInterfaces with the distance
    double RestraintDistanceSpinLabel::operator()
    (
      const restraint::AtomDistanceAssignment &ASSIGNMENT
    ) const
    {
      if( !ASSIGNMENT.GetAtomA().AllCoordinatesDefined() || !ASSIGNMENT.GetAtomB().AllCoordinatesDefined())
      {
        return double( 0);
      }

      return ScoreDistance( ASSIGNMENT.CalculateAtomDistance(), ASSIGNMENT.GetDistance());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief ScoreDistance is the function which calculates the score a given distance should get
    //! @brief DISTANCE is the distance which will be scored - this is the distance from the protein model
    //! @brief RESTRAINT_DISTANCE is the distance that DISTANCE will be scored against
    //! @return returns the score of a distance of DISTANCE
    double RestraintDistanceSpinLabel::ScoreDistance
    (
      const double &DISTANCE, const double &RESTRAINT_DISTANCE
    ) const
    {
      // create cont double "sl_cb" and initialize with the difference between "RESTRAINT_DISTANCE" and "DISTANCE"
      double sl_cb( RESTRAINT_DISTANCE - DISTANCE);

      // create double "score" and initialize with the potential associated with "bin_number"
      const double score( GetEnergyFunction().First()( sl_cb));

      BCL_Message
      (
        util::Message::e_Verbose,
        "restraint: " + util::Format()( RESTRAINT_DISTANCE) + " actual: " + util::Format()( DISTANCE) +
        " difference: " + util::Format()( RESTRAINT_DISTANCE - DISTANCE) + " score: " + util::Format()( score)
      );

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RestraintDistanceSpinLabel::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Knowledge-based energy function score for spin labels, histogram file is " + GetDefaultHistogramFilename()
      );
      return serializer;
    }

    //! @brief gets the energy function and lower and upper bounds
    //! @return the energy function and lower and upper bounds
    storage::Triplet< math::CubicSplineDamped, double, double> &RestraintDistanceSpinLabel::GetEnergyFunction()
    {
      // if the energy function has not been set
      if( s_EnergyFunction.First().GetXValues().IsEmpty())
      {
        // read in it
        s_EnergyFunction = ReadEnergyFunction();
      }

      // end
      return s_EnergyFunction;
    }

    //! @brief reads the energy function and lower and upper bounds
    //! @return the energy function and lower and upper bounds
    storage::Triplet< math::CubicSplineDamped, double, double> RestraintDistanceSpinLabel::ReadEnergyFunction()
    {
      // create Histogram "histogram" which will be used to hold the histogram of SL-CB distances
      math::Histogram histogram;

      // create IFStream "read"
      io::IFStream read;

      // open "read" and bind it to the histogram file containing SL-CB distances
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( GetDefaultHistogramFilename()));

      // read in from "read" into "histogram"
      read >> histogram;

      // close and clear read stream
      io::File::CloseClearFStream( read);

      const storage::VectorND< 2, double> boundaries( histogram.GetBoundaries());

      const size_t start_bin( histogram.GetIndexOfFirstInformationContainingBin() - 1);

      histogram.RemoveBinsBeforeIndex( start_bin);
      const double lower_bound( histogram.GetBoundaries().First());
      const size_t last_bin( histogram.GetIndexOfLastInformationContainingBin() + 1);
      histogram.RemoveBinsAfterIndex( last_bin);
      const double upper_bound( histogram.GetBoundaries().Second());

      return storage::Triplet< math::CubicSplineDamped, double, double>
      (
        EnergyDistribution::GeneratePotentialFromHistogram
        (
          histogram, 1000000, math::e_FirstDer, EnergyDistribution::GetDefaultFirstDerivative(), true, true
        ),
        lower_bound,
        upper_bound
      );
    }

  } // namespace score
} // namespace bcl
