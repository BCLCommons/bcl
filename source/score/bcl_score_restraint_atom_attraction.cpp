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
#include "score/bcl_score_restraint_atom_attraction.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_trigonometric_transition.h"
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
    const util::SiPtr< const util::ObjectInterface> RestraintAtomAttraction::s_Instance
    (
      util::Enumerated< RestraintAtomDistanceAssignment>::AddInstance( new RestraintAtomAttraction())
    );

    //! @brief gives the default width of the cosine transition region of the peicewise function
    //! @return double which is the default width of the cosine transition region of the peicewise function
    double RestraintAtomAttraction::GetDefaultTransitionWidth()
    {
      return double( 25.0);
    }

    //! @brief gives the default range of the depth of the cosine transition region of the peicewise function
    //! @return double which gives the default range of the depth of the cosine transition region of the function
    const math::Range< double> &RestraintAtomAttraction::GetDefaultDepthRange()
    {
      static const math::Range< double> s_depth_range( -1.0, 0.0);

      return s_depth_range;
    }

    //! @brief amount to shift in the x-direction so that attraction is not level at same place as KB potential
    //! @return double amount to shift in x-direction so that attraction is not level at same place as KB potential
    double RestraintAtomAttraction::GetDefaultScoreOffset()
    {
      return 2.0;
    }

    //! @brief gives the x-coordinate of the ending point for the attraction on the left of the KB potential
    //! @return double which is the x-coordinate of the ending point for the attraction on the left of the KB potential
    double RestraintAtomAttraction::GetDefaultLeftEndWell( const math::Histogram &HISTOGRAM)
    {
      return HISTOGRAM.GetBoundaries().First() +
        double( HISTOGRAM.GetIndexOfFirstInformationContainingBin()) * HISTOGRAM.GetBinSize() +
        RestraintAtomAttraction::GetDefaultScoreOffset();
    }

    //! @brief gives the x-coordinate of the ending point for the attraction on the right of the KB potential
    //! @return double which is the x-coordinate of the ending point for the attraction on the right of the KB potential
    double RestraintAtomAttraction::GetDefaultRightEndWell( const math::Histogram &HISTOGRAM)
    {
      BCL_Assert
      (
        HISTOGRAM.GetBinning().GetSize() != HISTOGRAM.GetIndexOfLastInformationContainingBin(),
        "Your histogram has counts outside its boundaries. Increase the width of your histogram."
      );
      return HISTOGRAM.GetBoundaries().Second() -
        double( HISTOGRAM.GetBinning().GetSize() - HISTOGRAM.GetIndexOfLastInformationContainingBin() - 1) *
        HISTOGRAM.GetBinSize() - GetDefaultScoreOffset() + GetDefaultTransitionWidth();
    }

  //////////
  // data //
  //////////

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &RestraintAtomAttraction::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "atom_attraction");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintAtomAttraction::RestraintAtomAttraction() :
      m_WellPotential( 0.0, GetDefaultTransitionWidth(), GetDefaultDepthRange().GetMin(), GetDefaultDepthRange().GetMax()),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief constructor from a specified histogram file
    //! @param WELL_DEPTH range over which the function will cover
    //! @param END_OF_WELL where the well should end
    //! @param WIDTH how wide the cos portion of the scoring function should be
    //! @param MAX_TO_MIN true if the function should go from max to min value - false otherwise
    //! @param SCHEME scheme to be used
    RestraintAtomAttraction::RestraintAtomAttraction
    (
      const math::Range< double> &WELL_DEPTH,
      const double END_OF_WELL,
      const double WIDTH,
      const bool MAX_TO_MIN,
      const std::string &SCHEME
    ) :
      m_WellPotential
      (
        END_OF_WELL - WIDTH,
        END_OF_WELL,
        MAX_TO_MIN ? WELL_DEPTH.GetMax() : WELL_DEPTH.GetMin(),
        MAX_TO_MIN ? WELL_DEPTH.GetMin() : WELL_DEPTH.GetMax()
      ),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone is the virtual copy constructor
    RestraintAtomAttraction *RestraintAtomAttraction::Clone() const
    {
      return new RestraintAtomAttraction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RestraintAtomAttraction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &RestraintAtomAttraction::GetAlias() const
    {
      return m_Scheme;
    }

    //! @brief gives the peicewise function that comprises this score
    //! @return function that is this score
    const math::TrigonometricTransition &RestraintAtomAttraction::GetFunction() const
    {
      return m_WellPotential;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator scores protein model
    //! @param RESTRAINT restraint to be scored
    //! @return score
    double RestraintAtomAttraction::operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const
    {
      // calculate the distance between the two modeled atoms
      const double cb_distance( RESTRAINT.CalculateAtomDistance());

      // if the calculated distance is undefined
      if( !util::IsDefined( cb_distance))
      {
        // return default score
        return std::max( m_WellPotential.GetTransitionBeginYAxis(), m_WellPotential.GetTransitionEndYAxis());
      }

      // return the well-potential's scoring value
      return m_WellPotential( RESTRAINT.GetDistance() - cb_distance);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RestraintAtomAttraction::GetSerializer() const
    {
      io::Serializer serializer;
      if( m_Scheme == GetDefaultScheme())
      {
        serializer.SetClassDescription
        (
          "A piecewise-continuous step-like score with a cosine transition region in place of the step."
        );
        serializer.AddInitializer
        (
          "start angstroms",
          "beginning of the transition region in angstroms representing target restraint distance between atoms - "
          "model's distance between atoms",
          io::Serialization::GetAgent( &m_WellPotential.GetTransitionBeginXAxis())
        );
        serializer.AddInitializer
        (
          "end angstroms",
          "end of the transition region in angstroms representing target restraint distance between atoms - "
          "model's distance between atoms",
          io::Serialization::GetAgent( &m_WellPotential.GetTransitionEndXAxis())
        );
        serializer.AddInitializer
        (
          "start score",
          "score for distances <= start angstroms",
          io::Serialization::GetAgent( &m_WellPotential.GetTransitionBeginYAxis())
        );
        serializer.AddInitializer
        (
          "end score",
          "score for distances >= end angstroms",
          io::Serialization::GetAgent( &m_WellPotential.GetTransitionEndYAxis())
        );
      }
      else
      {
        serializer.SetClassDescription
        (
          "An piecewise-continuous step-like score with a cosine transition region in place of the step.\n "
          "f(x) = " + m_WellPotential.AsString() + " where x is the target distance of the atoms in the restraint minus "
          "their actual distance in the current model"
        );
      }
      return serializer;
    }

  } // namespace score
} // namespace bcl
