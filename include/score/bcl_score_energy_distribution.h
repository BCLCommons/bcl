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

#ifndef BCL_SCORE_ENERGY_DISTRIBUTION_H_
#define BCL_SCORE_ENERGY_DISTRIBUTION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "contact/bcl_contact_types.h"
#include "math/bcl_math_histogram_3d.h"
#include "math/bcl_math_spline_border_type.h"
#include "storage/bcl_storage_pair.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnergyDistribution
    //! @brief a class with a collection of static functions, that convert probability distributions into bayesian potentials
    //!
    //! @see @link example_score_energy_distribution.cpp @endlink
    //! @author woetzen
    //! @date Sep 03, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EnergyDistribution
    {

    public:

    //////////
    // data //
    //////////

      //! pseudo count
      static const size_t s_PseudoCount = 1;

      //! propensity for nn-observed features
      static const double s_NonObservedPropensity;

      //! minimal count for aa environment potential to be considered
      static const double s_MinimalEnvironenmentBinCount;

    ////////////////
    // operations //
    ////////////////

      //! @brief converts a given histogram into an energy distribution
      //! @param HISTOGRAM the histogram to be used
      //! @param PSEUDOCOUNT a pseudocount = 1 by default
      //! @return cubic spline with linear continuation
      static math::CubicSplineDamped EnergyfunctionFromHistogram
      (
        const math::Histogram &HISTOGRAM,
        const double PSEUDOCOUNT = double( s_PseudoCount)
      );

      //! @brief create a bicubicspline for the energy distribution derived for angle/distance statistics by symmetrizing
      //! @param ANGLE_DISTANCE_DISTRIBUTION Histogram2D that contains the counts for angle distance distribution
      //! @param DISTANCE_RANGE distance range for given contact type to be used
      //! @return BicubicSpline for given ContactType's angle distance distribution
      static math::BicubicSpline SSEPackingSymmetrize
      (
        const math::Histogram2D &ANGLE_DISTANCE_DISTRIBUTION,
        const math::Range< double> &DISTANCE_RANGE
      );

      //! @brief creates a bicubic spline for the energy distribution derived for EPR distance potentials
      //! @param HISTOGRAM histogram containing the counts for the distance exposure distributions
      //! @param DISTANCE_RANGE distance range for given contact type to be used
      //! @return bicubic spline for the given EPR distance distributions
      static math::BicubicSpline EPRDistance
      (
        const math::Histogram2D &HISTOGRAM,
        const math::Range< double> &DISTANCE_RANGE
      );

      //! @brief creates a bicubic spline for the energy distribution derived for contact energy
      //!        multiplies the entries in the associated histogram by the energy for the contact term
      //! @param HISTOGRAM histogram containing the counts for the distance exposure distributions
      //! @param ENERGY energy of the contact
      //! @return bicubic spline for the given EPR distance distributions
      static math::TricubicSpline DistanceAngleContactEnergy
      (
        const math::Histogram3D &HISTOGRAM,
        const double &ENERGY
      );

      //! @brief create a bicubicspline for the energy distribution derived for angle/distance statistics
      //! @param ANGLE_DISTANCE_DISTRIBUTION Histogram2D that contains the counts for angle distance distribution
      //! @param DISTANCE_RANGE distance range for given contact type to be used
      //! @return BicubicSpline for given ContactType's angle distance distribution
      static math::BicubicSpline SSEPacking2D
      (
        const math::Histogram2D &ANGLE_DISTANCE_DISTRIBUTION,
        const math::Range< double> &DISTANCE_RANGE
      );

      //! create a cubic spline for the potential between a pair of amino acids derived form statistics
      //! @param AA_PAIR_DISTANCE_DISTRIBUTION the distribution the energy is generated for
      //! @param BACKGROUND_DISTRIBUTION the background probability
      //! @param NORMALIZE_BY_BACKGROUND normalize by given background - otherwise each bin has a propensity purely derived from this distribution
      static math::CubicSplineDamped AAPairPotential
      (
        const math::Histogram &AA_PAIR_DISTANCE_DISTRIBUTION,
        const math::Histogram &BACKGROUND_DISTRIBUTION,
        const bool NORMALIZE_BY_BACKGROUND
      );

      //! @brief create a cubic spline for an environment potential
      //! @param AA_ENV_DISTRIBUTION the distribution of exposure measure
      //! @param BACKGROUND_DISTRIBUTION the background distribution for that exposure measure
      //! @param NORMALIZE_BY_BACKGROUND normalize by given background - otherwise each bin has a propensity purely derived from this distribution
      //! @return a CubicSpline that represents the energy distribution
      static math::CubicSplineDamped AAEnvironmentPotential
      (
        const math::Histogram &AA_ENV_DISTRIBUTION,
        const math::Histogram &BACKGROUND_DISTRIBUTION,
        const bool NORMALIZE_BY_BACKGROUND = true
      );

      //! create a Map for each aa type with its environment potential
      static storage::Map< biol::AAType, math::CubicSplineDamped>
      AAEnvironmentPotential
      (
        const storage::Vector< math::Histogram> &HISTOGRMAS_AA
      );

      //! create a map for each membrane region with a vector for each amino acid for the membrane environment potential
      static storage::Map< biol::EnvironmentType, storage::Map< biol::AAType, math::CubicSplineDamped> >
      AAMembraneEnvironmentPotential
      (
        const storage::Vector< storage::Vector< math::Histogram> > &HISTOGRMAS_MEMBRANE_AA
      );

      //! create a cubic spline for loop length potentials derived form statistics
      static math::CubicSplineDamped LoopLengthDistancePotential( const math::Histogram &LOOPLENGTHDISTRIBUTION);

      //! create a bicubic spline for loop length potentials derived form statistics
      static math::BicubicSpline LoopLengthDistancePotential
      (
        const math::Histogram2D &LOOPLENGTHDISTRIBUTION,
        const size_t MAX_RESIDUES
      );

      //! @brief create Cubic splines for loop length potentials derived form statistics for different loop lengths
      static storage::Vector< math::CubicSplineDamped> LoopLengthDistancePotential
      (
        const storage::Vector< math::Histogram> &LOOPLENGTH_DISTRIBUTIONS
      );

      //! remove the bins from the end, that are empty in all given histograms
      //! @param HISTOGRAMS list of histograms
      //! @return number of removed bins
      static void RemoveAdditionalEmptyBinsExceptOne( util::SiPtrVector< math::Histogram> &HISTOGRAMS);

      //! @brief This function reads the aa_distances histogram and parses it into a matrix and returns it
      //! matrix gives for each AAType pair, probability of observing them within the given distance cutoff
      //! @param DISTANCE_CUTOFF distance in angstrom that is the cutoff
      //! @return matrix that contains the aa pair probabilities
      static linal::Matrix< double> ReadAAPairDistanceHistogram( const double DISTANCE_CUTOFF);

      //! This function returns a matrix a probabilities for a pair of aas to be within a DISTANCE
      static linal::Matrix< double> ProbabilityOfAAPairTOBeWithinDistance( const double &DISTANCE);

      //! @brief This function converts phi psi angle histogram for an amino acid and generates a potential
      //! @param PHI_PSI_ANGLE_DISTRIBUTION  Histogram 2D that contains amino acid phi psi angle statistics
      //! @param BACKGROUND the background distribution to use
      //! @param NORMALIZE_BY_BACKGROUND normalize by given background - otherwise each bin has a propensity purely derived from this distribution
      //! @return returns a bicubic spline
      static math::BicubicSpline PhiPsiAnglePotential
      (
        const math::Histogram2D &PHI_PSI_ANGLE_DISTRIBUTION,
        const math::Histogram2D &BACKGROUND,
        const bool NORMALIZE_BY_BACKGROUND = true
      );

      //! @brief convert a histogram into energu function with a cosine background distribution
      //! @param HISTOGRAM from 0 to pi/2
      //! @return Cubic spline representing the energy
      static math::CubicSplineDamped AngleAlignmentPotential( const math::Histogram &HISTOGRAM);

      //! @brief default first derivative for CubicSplines
      //! @return pair of left and right hand side first derivative for a cubic spline
      static const storage::Pair< double, double> &GetDefaultFirstDerivative()
      {
        static const storage::Pair< double, double> s_default_first_derivative( 0.0, 0.0);

        return s_default_first_derivative;
      }

      //! @brief This function reads a histogram for an amino acid and generates a potential
      //! @param HISTOGRAM Histogram that contains statistics =
      //! @param PSEUDOCOUNT pseudocount to be added
      //! @param BORDER_FLAG Training method to be used for handling the borders
      //! @param FIRST_DERIVATIVE only to be used when e_FirstDer is being used
      //! @param NORMALIZE_ENERGIES if true normalizes the energies calculated from the histogram to be between 0 and 1
      //! @param BONUS_ONLY if true shifts the energy potential to be 0 or less (i.e. no penalties)
      //! @return returns a cubic spline
      static math::CubicSplineDamped GeneratePotentialFromHistogram
      (
        const math::Histogram &HISTOGRAM,
        const double PSEUDOCOUNT,
        const math::SplineBorderType BORDER_FLAG,
        const storage::Pair< double, double> &FIRST_DERIVATIVE = GetDefaultFirstDerivative(),
        const bool NORMALIZE_ENERGIES = false,
        const bool BONUS_ONLY = false
      );

    }; // class Energydistribution

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_ENERGY_DISTRIBUTION_H_
