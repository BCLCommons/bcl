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
#include "score/bcl_score_energy_distribution.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_bicubic_spline.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"
#include "math/bcl_math_smooth_data.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_tricubic_spline.h"
#include "score/bcl_score_aa_pair_distance.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! propensity for nn-observed features
    const double EnergyDistribution::s_NonObservedPropensity = 1.0;

    //! minimal count for aa environment potential to be considered
    const double EnergyDistribution::s_MinimalEnvironenmentBinCount = 2.5;

  ////////////////
  // operations //
  ////////////////

    //! @brief converts a given histogram into an energy distribution
    //! @param HISTOGRAM
    //! @return cubic spline with linear continuation
    math::CubicSplineDamped EnergyDistribution::EnergyfunctionFromHistogram
    (
      const math::Histogram &HISTOGRAM, const double PSEUDOCOUNT
    )
    {
      // copy the histogram
      math::Histogram current( HISTOGRAM);
      current.RemoveBinsAfterIndex( current.GetIndexOfLastInformationContainingBin());

      linal::Vector< double> counts( current.GetHistogram());
      const linal::Vector< double> bins( current.GetBinning());

      // add pseudocount
      counts += PSEUDOCOUNT;

      // calculate propensities
      counts.SetToSum( counts.Sum());

      // logarithmize
      for( double *ptr = counts.Begin(), *ptr_end( counts.End()); ptr != ptr_end; ptr++)
      {
        ( *ptr) = -log( *ptr);
      }

      math::CubicSplineDamped energydistribution;
      energydistribution.Train
                         (
                           bins( 0),
                           current.GetBinSize(),
                           counts,
                           0.0,
                           0.0
                         );

      // end
      return energydistribution;
    }

    //! @brief create a bicubicspline for the energy distribution derived for angle/distance statistics by symmetrization
    //! @param ANGLE_DISTANCE_DISTRIBUTION Histogram2D that contains the counts for angle distance distribution
    //! @param DISTANCE_RANGE distance range for given contact type to be used
    //! @return BicubicSpline for given ContactType's angle distance distribution
    math::BicubicSpline EnergyDistribution::SSEPackingSymmetrize
    (
      const math::Histogram2D &ANGLE_DISTANCE_DISTRIBUTION,
      const math::Range< double> &DISTANCE_RANGE
    )
    {
      //cutoff values for different types of sse_interaction represented by the row in the histogram - angle_distance_count
      const std::pair< size_t, size_t> lower_upper_distance_cutoff
      (
        size_t( std::floor( DISTANCE_RANGE.GetMin() / ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().Second())),
        size_t( std::ceil( DISTANCE_RANGE.GetMax() / ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().Second()))
      );

      const linal::Matrix< double> &angle_distance_counts( ANGLE_DISTANCE_DISTRIBUTION.GetHistogram());
      const storage::VectorND< 2, linal::Vector< double> > angle_distance_bins
      (
        ANGLE_DISTANCE_DISTRIBUTION.GetBinningXY()
      );

      // create vectors to sum angle counts and distance counts
      linal::Vector< double> angle_counts( angle_distance_counts.GetNumberCols());
      linal::Vector< double> distance_counts( angle_distance_counts.GetNumberRows(), double( s_PseudoCount));

      //add up all counts for each distance bin and for each angle bin in the given distance range
      for( size_t i( 0); i < angle_distance_counts.GetNumberRows(); ++i)
      {
        for( size_t j( 0); j < angle_distance_counts.GetNumberCols(); j++)
        {
          angle_counts( j) += angle_distance_counts( i, j);
          distance_counts( i) += angle_distance_counts( i, j);
        }
      }

      // normalize distance distribution by distance^2
      distance_counts /= angle_distance_bins.Second();

      // create a matrix with two rows margin to store the actual energy distribution
      linal::Matrix< double> energydistributionmatrix
      (
        lower_upper_distance_cutoff.second - lower_upper_distance_cutoff.first + 4, angle_distance_counts.GetNumberCols()
      );

      // multiply counts of distance and angle( plus pseudocount) into the energydistributionmatrix
      double energysum( 0);
      for( size_t row( 2); row < energydistributionmatrix.GetNumberRows() - 2; ++row)
      {
        for( size_t col( 0); col < energydistributionmatrix.GetNumberCols(); ++col)
        {
          energydistributionmatrix( row, col) =
            distance_counts( lower_upper_distance_cutoff.first + row - 2) * ( angle_counts( col) + 1.0);
          energysum += energydistributionmatrix( row, col);
        }
      }

      double completesum( energysum);
      //multiply counts of distance and angle( plus pseudo count) into the complete sum
      for( size_t col( 0); col < energydistributionmatrix.GetNumberCols(); ++col)
      {
        for( size_t row( lower_upper_distance_cutoff.second); row < angle_distance_counts.GetNumberRows(); ++row)
        {
          completesum += distance_counts( row) * ( angle_counts( col) + 1.0);
        }
      }

      completesum /= energydistributionmatrix.GetNumberCols() * angle_distance_counts.GetNumberRows();
      energydistributionmatrix /= completesum;

      //-log of every bin to have the energy
      for
      (
        double
          *ptr( energydistributionmatrix.Begin() + energydistributionmatrix.GetNumberCols()),
          *ptr_end( energydistributionmatrix.End() - energydistributionmatrix.GetNumberCols());
        ptr != ptr_end; ++ptr)
      {
        *ptr = -log( *ptr);
      }

      //set the second and the row bfore last to the average of the neigbouring lines
      energydistributionmatrix.ReplaceRow
      (
        0,
        linal::Vector< double>( energydistributionmatrix.GetNumberCols(), double( 3.0))
      );
      energydistributionmatrix.ReplaceRow
      (
        1,
        ( energydistributionmatrix.GetRow( 2) + energydistributionmatrix.GetRow( 0)) / 2.0
      );
      energydistributionmatrix.ReplaceRow
      (
        energydistributionmatrix.GetNumberRows() - 2,
        energydistributionmatrix.GetRow( energydistributionmatrix.GetNumberRows() - 3) / 2.0
      );

      //instantiate Spline
      math::BicubicSpline energydistribution;
      const math::SplineBorderType borderflags[ 2] = { math::e_FirstDer, math::e_Periodic};
      const double start[ 2] =
      {
          angle_distance_bins.Second()( lower_upper_distance_cutoff.first)
        - 2 * ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().Second(),
          angle_distance_bins.First()( 0)
      };
      const double binsize[ 2] =
      {
        ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().Second(), ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().First()
      };
      const bool lincont[ 2] = { true, false};
      const storage::Pair< double, double> firstorderder[ 2] =
      {
        storage::Pair< double, double>( -1, 0), storage::Pair< double, double>( 0, 0)
      };
      energydistribution.Train( borderflags, start, binsize, energydistributionmatrix, lincont, firstorderder);

      // end
      return energydistribution;
    }

    //! @brief creates a bicubic spline for the energy distribution derived for EPR distance potentials
    //! @param HISTOGRAM histogram containing the counts for the distance exposure distributions
    //! @param DISTANCE_RANGE distance range for given contact type to be used
    //! @return bicubic spline for the given EPR distance distributions
    math::BicubicSpline EnergyDistribution::EPRDistance
    (
      const math::Histogram2D &HISTOGRAM,
      const math::Range< double> &DISTANCE_RANGE
    )
    {
      const linal::Matrix< double> counts( HISTOGRAM.GetHistogram());
      const storage::VectorND< 2, linal::Vector< double> > bins( HISTOGRAM.GetBinningXY());

      // create the energy function
      const math::SplineBorderType borderflags[ 2] = { math::e_FirstDer, math::e_Periodic};
      const double start[ 2] = { bins.Second()( 2) - 2 * HISTOGRAM.GetBinSizeXY().Second(), bins.First()( 0)};
      const double binsize[ 2] = { HISTOGRAM.GetBinSizeXY().Second(), HISTOGRAM.GetBinSizeXY().First()};
      const bool lincont[ 2] = { true, false};
      const storage::Pair< double, double> firstorderder[ 2] =
        {
          storage::Pair< double, double>( -1, 0),
          storage::Pair< double, double>( 0, 0)
        };
      math::BicubicSpline function;
      function.Train( borderflags, start, binsize, counts, lincont, firstorderder);

      return function;
    }

    //! @brief creates a bicubic spline for the energy distribution derived for contact energy
    //!        multiplies the entries in the associated histogram by the energy for the contact term
    //! @param HISTOGRAM histogram containing the counts for the distance exposure distributions
    //! @param ENERGY energy of the contact
    //! @return bicubic spline for the given EPR distance distributions
    math::TricubicSpline EnergyDistribution::DistanceAngleContactEnergy
    (
      const math::Histogram3D &HISTOGRAM,
      const double &ENERGY
    )
    {
      math::Tensor< double> contact_frequency( HISTOGRAM.GetHistogram());
      if( util::IsDefined( ENERGY))
      {
        for( auto itr( contact_frequency.Begin()), itr_end( contact_frequency.End()); itr != itr_end; ++itr)
        {
          if( *itr < 0.0)
          {
            *itr = 0.0;
          }
        }
        contact_frequency *= ENERGY;
      }
      else
      {
        for( auto itr( contact_frequency.Begin()), itr_end( contact_frequency.End()); itr != itr_end; ++itr)
        {
          if( *itr < 0.0)
          {
            *itr = -1.0;
          }
          else if( *itr > 0.0)
          {
            *itr = 0.0;
          }
        }
      }

      const storage::VectorND< 3, linal::Vector< double> > bins( HISTOGRAM.GetBinningXYZ());

      // create the energy function
      const math::SplineBorderType borderflags[ 2] = { math::e_FirstDer, math::e_FirstDer};
      const double start[ 3] = { bins.First()( 0), bins.Second()( 0), bins.Third()( 0)};
      const double binsize[ 3] = { HISTOGRAM.GetBinSizeXYZ().First(), HISTOGRAM.GetBinSizeXYZ().Second(), HISTOGRAM.GetBinSizeXYZ().Third()};
      const bool lincont[ 3] = { true, true, true};
      const storage::Pair< double, double> firstorderder[ 3] =
      {
        storage::Pair< double, double>( 0, 0),
        storage::Pair< double, double>( 0, 0),
        storage::Pair< double, double>( 0, 0)
      };
      math::TricubicSpline function;
      function.Train( borderflags, start, binsize, contact_frequency, lincont, firstorderder);

      return function;
    }

    //! @brief create a bicubicspline for the energy distribution derived for angle/distance statistics
    //! @param ANGLE_DISTANCE_DISTRIBUTION Histogram2D that contains the counts for angle distance distribution
    //! @param DISTANCE_RANGE distance range for given contact type to be used
    //! @return BicubicSpline for given ContactType's angle distance distribution
    math::BicubicSpline EnergyDistribution::SSEPacking2D
    (
      const math::Histogram2D &ANGLE_DISTANCE_DISTRIBUTION,
      const math::Range< double> &DISTANCE_RANGE
    )
    {
      //cutoff values for different types of sse_interaction represented by the row in the histogram - angle_distance_count
      const std::pair< size_t, size_t> lower_upper_distance_cutoff_index
      (
        size_t( std::floor( ( DISTANCE_RANGE.GetMin() - ANGLE_DISTANCE_DISTRIBUTION.GetBoundariesY().First()) / ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().Second())),
        size_t( std::ceil( ( DISTANCE_RANGE.GetMax() - ANGLE_DISTANCE_DISTRIBUTION.GetBoundariesY().First()) / ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().Second()))
      );

      linal::Matrix< double> angle_distance_counts( ANGLE_DISTANCE_DISTRIBUTION.GetHistogram());
      const storage::VectorND< 2, linal::Vector< double> > angle_distance_bins
      (
        ANGLE_DISTANCE_DISTRIBUTION.GetBinningXY()
      );

      // variable that adds up all counts
      double integral( 0);

      //add up all counts for each distance bin and for each angle bin in the given distance range
      for( size_t i( lower_upper_distance_cutoff_index.first); i < lower_upper_distance_cutoff_index.second; ++i)
      {
        for( size_t j( 0); j < angle_distance_counts.GetNumberCols(); j++)
        {
          double &current_value( angle_distance_counts( i, j));
          // add a pseudo count
          current_value += s_PseudoCount;
          // normalize by the distance
          current_value /= angle_distance_bins.Second()( i);
          // add that value to the integral
          integral += current_value;
        }
      }

      // normalize
      angle_distance_counts *= angle_distance_counts.GetNumberCols() * angle_distance_counts.GetNumberRows() / integral;

      // initialize matrix with two row margin
      linal::Matrix< double> energydistributionmatrix
      (
        lower_upper_distance_cutoff_index.second - lower_upper_distance_cutoff_index.first + 4,
        angle_distance_bins.First().GetSize(),
        double( 0.0)
      );

      // copy the values to the matrix with the 2 row margin on each side
      for( size_t i( lower_upper_distance_cutoff_index.first); i < lower_upper_distance_cutoff_index.second; ++i)
      {
        for( size_t j( 0); j < angle_distance_counts.GetNumberCols(); j++)
        {
          energydistributionmatrix( i + 2 - lower_upper_distance_cutoff_index.first, j) = angle_distance_counts( i, j);
        }
      }

      // repulsion
      static const double s_default_repulsion( 7.0);

      //-log of every bin to get the energy
      for
      (
        double
          *ptr( energydistributionmatrix.Begin() + 2 * energydistributionmatrix.GetNumberCols()),
          *ptr_end( energydistributionmatrix.End() - 2 * energydistributionmatrix.GetNumberCols());
        ptr != ptr_end;
        ++ptr
      )
      {
        *ptr = -log( *ptr);
      }

      //set the second and the row before last to the average of the neighboring lines
      energydistributionmatrix.ReplaceRow
                               ( 0, linal::Vector< double>( energydistributionmatrix.GetNumberCols(), s_default_repulsion));
      energydistributionmatrix.ReplaceRow
                               ( 1, ( energydistributionmatrix.GetRow( 2) + energydistributionmatrix.GetRow( 0)) / 2.0);
      energydistributionmatrix.ReplaceRow
                               (
                                 energydistributionmatrix.GetNumberRows() - 1,
                                 linal::Vector< double>( energydistributionmatrix.GetNumberCols(), double( 0.0))
                               );
      energydistributionmatrix.ReplaceRow
                               (
                                 energydistributionmatrix.GetNumberRows() - 2,
                                 energydistributionmatrix.GetRow( energydistributionmatrix.GetNumberRows() - 3) / 2.0
                               );

      //instantiate Spline
      math::BicubicSpline energydistribution;
      const math::SplineBorderType           borderflags[ 2] = { math::e_FirstDer, math::e_Periodic};
      const double                                 start[ 2] =
        {
            angle_distance_bins.Second()( lower_upper_distance_cutoff_index.first)
          - 2 * ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().Second(),
            angle_distance_bins.First()( 0)
        };
      const double                               binsize[ 2] =
        {
          ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().Second(),
          ANGLE_DISTANCE_DISTRIBUTION.GetBinSizeXY().First()
        };
      const bool                                 lincont[ 2] = { true, false};
      const storage::Pair< double, double> firstorderder[ 2] =
        {
          storage::Pair< double, double>( -1, 0),
          storage::Pair< double, double>( 0, 0)
        };
      energydistribution.Train( borderflags, start, binsize, energydistributionmatrix, lincont, firstorderder);

      // end
      return energydistribution;
    }

    //! create a cubic spline for the potential between a pair of amino acids derived form statistics
    //! @param AA_PAIR_DISTANCE_DISTRIBUTION the distribution the energy is generated for
    //! @param BACKGROUND_DISTRIBUTION the background probability
    //! @param NORMALIZE_BY_BACKGROUND normalize by given background - otherwise each bin has a propensity purely derived from this distribution
    math::CubicSplineDamped EnergyDistribution::AAPairPotential
    (
      const math::Histogram &AA_PAIR_DISTANCE_DISTRIBUTION,
      const math::Histogram &BACKGROUND_DISTRIBUTION,
      const bool NORMALIZE_BY_BACKGROUND
    )
    {
      linal::Vector< double> aa_distance_counts( AA_PAIR_DISTANCE_DISTRIBUTION.GetHistogram());
      const linal::Vector< double> distance_bins( AA_PAIR_DISTANCE_DISTRIBUTION.GetBinning());

      linal::Vector< double> background_distribution( BACKGROUND_DISTRIBUTION.GetHistogram());

      const size_t usedbins( 20);

      std::fill( aa_distance_counts.Begin() + usedbins, aa_distance_counts.End(), double( 0.0));
      std::fill( background_distribution.Begin() + usedbins, background_distribution.End(), double( 0.0));
      aa_distance_counts.SetToSum( 1.0);
      if( !NORMALIZE_BY_BACKGROUND)
      {
        // fill background with r^2
        for( size_t dist( 0); dist < background_distribution.GetSize(); ++dist)
        {
          background_distribution( dist) = math::Sqr( 0.5 + dist);
        }
      }
      background_distribution.SetToSum( 1.0);

      // divide each count in bin by the background distribution for that bin
      size_t i( 0);
      for
      (
        double
          *bg_counts( background_distribution.Begin()), *bg_counts_end( background_distribution.End()),
          *counts( aa_distance_counts.Begin()), *counts_end( aa_distance_counts.End());
        counts != counts_end && bg_counts != bg_counts_end && i < usedbins;
        ++counts, ++bg_counts, ++i
      )
      {
        if( *bg_counts > double( 0))
        {
          // divide by background probability
          *counts /= *bg_counts;
        }

        // add pseudo propensity for empty bins defining the repulsion
        *counts += 0.00000001;
      }

      linal::Vector< double> energydistributionvector( usedbins, aa_distance_counts.Begin());

      // -log of every bin to have the energy
      for( double *ptr = energydistributionvector.Begin(), *ptr_end( energydistributionvector.End()); ptr != ptr_end; ptr++)
      {
        ( *ptr) = -log( *ptr);
      }

      // last value set to zero
      *( energydistributionvector.End() - 1) = 0;

      math::CubicSplineDamped energydistribution;
      energydistribution.Train
                         (
                           AA_PAIR_DISTANCE_DISTRIBUTION.GetBinning().First(),
                           AA_PAIR_DISTANCE_DISTRIBUTION.GetBinSize(),
                           energydistributionvector,
                           0,
                           0
                         );

      // end
      return energydistribution;
    }

    //! @brief create a cubic spline for an environment potential
    //! @param AA_ENV_DISTRIBUTION the distribution of exposure measure
    //! @param BACKGROUND_DISTRIBUTION the background distribution for that exposure measure
    //! @param NORMALIZE_BY_BACKGROUND normalize by given background - otherwise each bin has a propensity purely derived from this distribution
    //! @return a CubicSpline that represents the energy distribution
    math::CubicSplineDamped EnergyDistribution::AAEnvironmentPotential
    (
      const math::Histogram &AA_ENV_DISTRIBUTION,
      const math::Histogram &BACKGROUND_DISTRIBUTION,
      const bool NORMALIZE_BY_BACKGROUND
    )
    {
      // binning for the environment potential
      const linal::Vector< double> env_bins( AA_ENV_DISTRIBUTION.GetBinning());

      // actual counts
      linal::Vector< double> aa_env_counts( AA_ENV_DISTRIBUTION.GetHistogram());
      aa_env_counts.SetToSum( 1.0);

      // background
      linal::Vector< double> background_distribution( BACKGROUND_DISTRIBUTION.GetHistogram());
      background_distribution.SetToSum( 1.0);

      // background should not be used - assume each bin is equally probable
      if( !NORMALIZE_BY_BACKGROUND)
      {
        background_distribution = 1.0 / double( background_distribution.GetSize());
      }

      // divide each count in bin by the background distribution for that bin
      for
      (
        double
          *bg_counts( background_distribution.Begin()), *bg_counts_end( background_distribution.End()),
          *counts( aa_env_counts.Begin()), *counts_end( aa_env_counts.End());
        counts != counts_end && bg_counts != bg_counts_end;
        ++counts, ++bg_counts
      )
      {
        if( *counts > 0.0)
        { // divide by background probability
          *counts /= *bg_counts;
        }
        // no background, so energy should be 0 by setting propensity to 1
        else
        {
          *counts = s_NonObservedPropensity;
        }

        //-log of every bin to have the energy
        *counts = -log( *counts);
      }

      math::CubicSplineDamped energydistribution;
      energydistribution.Train
                         (
                           env_bins( 0),
                           AA_ENV_DISTRIBUTION.GetBinSize(),
                           math::SmoothData::SmoothVector( aa_env_counts, 0.5, true),
                           0,
                           0
                         );

      return energydistribution;
    }

    //! create a Map for each aa type with its environment potential
    storage::Map< biol::AAType, math::CubicSplineDamped>
    EnergyDistribution::AAEnvironmentPotential
    (
      const storage::Vector< math::Histogram> &HISTOGRMAS_AA
    )
    {
      storage::Vector< math::Histogram> histograms( HISTOGRMAS_AA);

      // sum all histograms
      math::Histogram total_sum_histograms;
      for
      (
        storage::Vector< math::Histogram>::iterator
          itr( histograms.Begin()), itr_end( histograms.Begin() + biol::AATypes::s_NumberStandardAATypes);
        itr != itr_end;
        ++itr
      )
      {
        // add pseudo count and normalize the current histogram in order to remove bias to certain aa types
        math::Histogram &current_histogram( *itr);

        // find last information containing bin
        current_histogram.GetIndexOfLastInformationContainingBin( 2.5);
        current_histogram.ResetBinsAfterIndex( current_histogram.GetIndexOfLastInformationContainingBin( 2.5));

        current_histogram.Normalize();

        BCL_Assert
        (
          total_sum_histograms.Combine( current_histogram),
          "unable to combine histograms of different parameters"
        );
      }

      storage::Map< biol::AAType, math::CubicSplineDamped> potentials;

//      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
//      {
//        io::OFStream write;
//        io::File::MustOpenOFStream( write, "aa_neighbor_count_background.gnuplot");
//        math::GnuplotHeatmap heatmap;
//        heatmap.SetFromHistogram( total_sum_histograms, false, false);
//        heatmap.SetTitleAndLabel( "amino acid neighbor count background", "neighbor count", "", "p");
//        heatmap.SetPixelAndRatio( 1080, 800, -2.0);
//        heatmap.SetFont( "arialbd", 20);
//        heatmap.SetRotationXTics( 90.0);
//        heatmap.SetFilename( "aa_neighbor_count_background");
//        heatmap.WriteScript( write);
//        io::File::CloseClearFStream( write);
//      }

      // derive potential for each amino acid
      for
      (
        biol::AATypes::const_iterator
          aa_itr( biol::GetAATypes().Begin()),
          aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        potentials[ *aa_itr] = AAEnvironmentPotential( histograms( *aa_itr), total_sum_histograms, true);
      }

      // end
      return potentials;
    }

    //! create a vector for each membrane region with a vector for each amino acid for the membrane environment potential
    storage::Map< biol::EnvironmentType, storage::Map< biol::AAType, math::CubicSplineDamped> >
    EnergyDistribution::AAMembraneEnvironmentPotential
    (
      const storage::Vector< storage::Vector< math::Histogram> > &HISTOGRAMS_MEMBRANE_AA
    )
    {
      // number of bins
      const math::Histogram &first_histogram( HISTOGRAMS_MEMBRANE_AA( 0)( 0));

      // counts of each amino acid in each region and in their bins
      storage::Vector< storage::Vector< linal::Vector< double> > >
        region_counts
        (
          HISTOGRAMS_MEMBRANE_AA.GetSize(),
          storage::Vector< linal::Vector< double> >( HISTOGRAMS_MEMBRANE_AA( 0).GetSize())
        );

      // collect the sum of all counts in each region over all aas
      linal::Vector< double> region_sum( biol::GetEnvironmentTypes().GetEnumCount());

      // total number of amino acids
      storage::Vector< double> aminoacid_counts( biol::GetAATypes().GetEnumCount(), double( 0.0));

      for
      (
        storage::Vector< biol::EnvironmentType>::const_iterator
          env_itr( biol::GetEnvironmentTypes().GetReducedTypes().Begin()),
          env_itr_end( biol::GetEnvironmentTypes().GetReducedTypes().End());
        env_itr != env_itr_end; ++env_itr
      )
      {
        for
        (
          biol::AATypes::const_iterator
            aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          region_counts( *env_itr)( *aa_itr) =
            HISTOGRAMS_MEMBRANE_AA( *env_itr)( *aa_itr).GetHistogram();
          //add pseudo count
          region_counts( *env_itr)( *aa_itr) += double( s_PseudoCount);
          //sum up sums for each aatype
          region_sum( *env_itr) += region_counts( *env_itr)( *aa_itr).Sum();
        }

        //divide by the sum in all regions by number_bins for the environment and 20 different aatypes
        for
        (
          biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          region_counts( *env_itr)( *aa_itr) /=
            ( region_sum( *env_itr) / double( first_histogram.GetNumberOfBins()) / double( 20));
        }

        //add up the sum of the aas in all regions
        for
        (
          biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          //sum up sums for each aatype
          aminoacid_counts( *aa_itr) += region_counts( *env_itr)( *aa_itr).Sum();
        }
      }

      for
      (
        storage::Vector< biol::EnvironmentType>::const_iterator
          env_itr( biol::GetEnvironmentTypes().GetReducedTypes().Begin()),
          env_itr_end( biol::GetEnvironmentTypes().GetReducedTypes().End());
        env_itr != env_itr_end; ++env_itr
      )
      {
        // Divide by the sum of the aminoacids in all regions
        for
        (
          biol::AATypes::const_iterator aa_itr( biol::GetAATypes().Begin()),
            aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          //divide by the total number of one aatype in all regions, by the number_bins and 3 for 3 regions where the counts for one aa has been summed up
          region_counts( *env_itr)( *aa_itr) /=
            aminoacid_counts( *aa_itr) / double( first_histogram.GetNumberOfBins()) / double( biol::GetEnvironmentTypes().GetReducedTypes().GetSize());

          //-log to yield the energy
          for
          (
            double *ptr( region_counts( *env_itr)( *aa_itr).Begin()),
              *ptr_end( region_counts( *env_itr)( *aa_itr).End());
            ptr != ptr_end; ++ptr
          )
          {
            // logarithmize
            *ptr = -log( *ptr);
          }
        }
      }

      //construct the cubic splines and train them for each amino acid
      storage::Map< biol::EnvironmentType, storage::Map< biol::AAType, math::CubicSplineDamped> > energy_maps;

      for
      (
        storage::Vector< biol::EnvironmentType>::const_iterator
          env_itr( biol::GetEnvironmentTypes().GetReducedTypes().Begin()),
          env_itr_end( biol::GetEnvironmentTypes().GetReducedTypes().End());
        env_itr != env_itr_end; ++env_itr
      )
      {
        for
        (
          biol::AATypes::const_iterator
            aa_itr( biol::GetAATypes().Begin()),
             aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
           aa_itr != aa_itr_end;
           ++aa_itr
        )
        {
          energy_maps[ *env_itr][ *aa_itr].Train
          (
            first_histogram.GetBinning()( 0),           // start
            first_histogram.GetBinSize(),               // bin size
            math::SmoothData::SmoothVector( region_counts( *env_itr)( *aa_itr), 0.5, true), // data
            0,      // first derivative at left border
            0       // first derivative at right border
          );
        }
      }

      // end
      return energy_maps;
    }

    math::CubicSplineDamped EnergyDistribution::LoopLengthDistancePotential( const math::Histogram &LOOPLENGTHDISTRIBUTION)
    {
      linal::Vector< double> loop_length_counts( LOOPLENGTHDISTRIBUTION.GetHistogram());

      //search for the first position where the count is zero which means this is a forbidden loop conformation
      size_t index( 0);
      for( ; index < loop_length_counts.GetSize(); ++index)
      {
        if( index > 10 && loop_length_counts( index) < 1)
        {
          break;
        }
      }

      linal::Vector< double> loop_energy_distribution( index, loop_length_counts.Begin());

      //add pseudo count
      loop_energy_distribution += double( s_PseudoCount);

      //divide by the sum over i, so that the total sum is equal to the number of bins
      loop_energy_distribution.SetToSum( index);

      //-log of every bin to have the energy
      for( double *ptr = loop_energy_distribution.Begin(); ptr != loop_energy_distribution.End(); ++ptr)
      {
        *ptr = -log( *ptr);
      }

//      BCL_MessageStd( util::Format()( loop_energy_distribution));

      math::CubicSplineDamped energydistribution;
      energydistribution.Train
      (
        LOOPLENGTHDISTRIBUTION.GetBoundaries().First() + 0.5 * LOOPLENGTHDISTRIBUTION.GetBinSize(),
        LOOPLENGTHDISTRIBUTION.GetBinSize(),
        loop_energy_distribution,
        1,
        1
      );

      return energydistribution;
    }

    //! remove the bins from the end, that are empty in all given histograms
    //! @param HISTOGRAMS list of histograms
    //! @return number of removed bins
    void EnergyDistribution::RemoveAdditionalEmptyBinsExceptOne( util::SiPtrVector< math::Histogram> &HISTOGRAMS)
    {
      size_t number_remaining_bins( 0);

      // determine the highest last information containing bin within all histograms
      for( util::SiPtrVector< math::Histogram>::const_iterator itr( HISTOGRAMS.Begin()), itr_end( HISTOGRAMS.End());
        itr != itr_end;
        ++itr)
      {
        number_remaining_bins = std::max( number_remaining_bins, ( *itr)->GetIndexOfLastInformationContainingBin());
      }

      // leave one empty bin for repulsion term
      ++number_remaining_bins;

      for
      (
        util::SiPtrVector< math::Histogram>::iterator itr( HISTOGRAMS.Begin()), itr_end( HISTOGRAMS.End());
        itr != itr_end;
        ++itr
      )
      {
        ( *itr)->RemoveBinsAfterIndex( number_remaining_bins);
      }
    }

    //! create a cubic spline for loop length potentials derived form statistics
    math::BicubicSpline
    EnergyDistribution::LoopLengthDistancePotential
    (
      const math::Histogram2D &LOOPLENGTHDISTRIBUTION,
      const size_t MAX_RESIDUES
    )
    {
      // calculate the number of bins necessary from the LOOPLENGTHDISTRIBUTION
      // x-number of residues
      // y-euclidean distance
      const size_t number_cols
                   (
                     std::min
                     (
                       size_t( MAX_RESIDUES / LOOPLENGTHDISTRIBUTION.GetBinSizeXY().First()),
                       LOOPLENGTHDISTRIBUTION.GetNumberOfBinsX()
                     ) + 2
                   );

      // get a submatrix including only the number_cols determined by the MAX_RESIDUES - numbercols-1 is necessary since
      // GetColSubMatrix includes the passed COL
      linal::Matrix< double> distribution
      (
        LOOPLENGTHDISTRIBUTION.GetHistogram().CreateSubMatrix
        (
          LOOPLENGTHDISTRIBUTION.GetHistogram().GetNumberRows(),
          number_cols
        )
      );

      // determine the number of rows by checking each column for its first 0 after a value different from 0
      // (to skip the first 0's in a column
      size_t number_rows( 0);
      //iterate over all cols
      for( size_t col( 0); col < distribution.GetNumberCols(); ++col)
      {
        //initialize current count
        double current_count( distribution( 0, col));
        //iterate over all values in this col
        for( size_t row( 1); row < distribution.GetNumberRows(); ++row)
        {
          //store previous value
          const double previous_count( current_count);
          //store current count
          current_count = distribution( row, col);
          //if current is 0 and previous was not 0 (like at the top at a column)
          if( current_count == 0.0 && previous_count != 0.0)
          {
            //check that row is larger than already found number_rows
            if( row > number_rows)
            {
              //store highest number_row so far
              number_rows = row;
            }
            //fill the rest of the column with 0.0
            while( row < distribution.GetNumberRows())
            {
              distribution( row, col) = 0.0;
              ++row;
            }
            break;
          }
        }
      }

      //add 1 to number rows, since the best row found is 1 smaller than the number for rows to be considered
      ++number_rows;

      // smoothing is neccessary since there might be certain seq distances that do not occur
      // linal::Matrix< double> smoothed_distribution( math::SmoothData::SmoothMatrix( distribution.GetRowSubMatrix( 0, number_rows - 1), 0.5, false));
      linal::Matrix< double> smoothed_distribution( distribution.CreateSubMatrix( number_rows - 1, distribution.GetNumberCols()));
      linal::Vector< double> col_sum( number_cols, double( 0));

      for( size_t col( 0); col < number_cols; ++col)
      {
        double &current_col_sum( col_sum( col));
        for( size_t row( 0); row < number_rows; ++row)
        {
          smoothed_distribution( row, col) += double( 1.0);
          current_col_sum += smoothed_distribution( row, col);
        }
      }

      for( size_t col( 0); col < number_cols; ++col)
      {
        const double current_col_sum( col_sum( col) / number_rows);
        for( size_t row( 0); row < number_rows; ++row)
        {
          //normalize and logarithmize
          smoothed_distribution( row, col) = -log( smoothed_distribution( row, col) / current_col_sum);
        }
      }

//      // replace one before last with half of the previous col and the last col with 0 - so that every number of residues
//      // Extending current energy distribution gets 0 as score
//      smoothed_distribution.ReplaceCol( number_cols - 1, linal::Vector< double>( number_rows, double( 0.0)));
//      smoothed_distribution.ReplaceCol( number_cols - 2, ( smoothed_distribution.GetCol( number_cols - 3) / 2.0));

      //instantiate Spline
      math::BicubicSpline energydistribution;
      const math::SplineBorderType             borderflags[ 2] = { math::e_FirstDer, math::e_FirstDer};
      const double                             start[ 2] = { LOOPLENGTHDISTRIBUTION.GetBoundariesY().First() + 0.5 * LOOPLENGTHDISTRIBUTION.GetBinSizeXY().Second(),
                                                             LOOPLENGTHDISTRIBUTION.GetBoundariesX().First() + 0.5 * LOOPLENGTHDISTRIBUTION.GetBinSizeXY().First()};
      const double                           binsize[ 2] = { LOOPLENGTHDISTRIBUTION.GetBinSizeXY().Second(), LOOPLENGTHDISTRIBUTION.GetBinSizeXY().First()};
      const bool                             lincont[ 2] = { true, true};
      const storage::Pair< double, double> firstorderder[ 2] = { storage::Pair< double, double>( -1, 1), storage::Pair< double, double>( 0, 0)};
      energydistribution.Train( borderflags, start, binsize, smoothed_distribution, lincont, firstorderder);

//      util::Message::SetMessageLevel( util::Message::e_Silent);
//      util::Format format;
//      format.FFP( 7, 3).W( 7);
//
//      std::stringstream stream;
//      stream << util::Format().W( 7)( "smooth") + "\t");
//      for( double residues = 0;  residues < 2 * MAX_RESIDUES; residues += 4)
//      {
//        stream << format( residues) + "\t";
//      }
//      stream << '\n';
//
//      for( double distance = 0; distance < number_rows * 4; distance += 4)
//      {
//        stream << format( distance) + "\t";
//        for( double residues = 0; residues < 2 * MAX_RESIDUES && residues < number_cols * 4; residues += 4)
//        {
//          stream << format( smoothed_distribution( size_t( distance / 4.0), size_t( residues/4.0))) << '\t';
//        }
//        stream << '\n';
//      }
//      stream << util::Format()( smoothed_distribution);
//
//      stream << util::Format().W( 7)( "test") << '\t';
//      for( double residues = 0;  residues < 3 * MAX_RESIDUES; residues += 6)
//      {
//        stream << format( residues) << '\t';
//      }
//      stream << '\n';
//
//      for( double distance = 0; distance < 120; distance +=4)
//      {
//        stream << format( distance) << '\t';
//        for( double residues = 0;  residues < 3 * MAX_RESIDUES; residues += 30)
//        {
//          stream << format( energydistribution.FFP( MakeVector( distance, residues))) << '\t';
//        }
//        stream << '\n';
//      }
//      util::Message::ResetToPreviousMessageLevel();

      return energydistribution;
    }

    //! @brief create Cubic splines for loop length potentials derived form statistics for different loop lengths
    storage::Vector< math::CubicSplineDamped> EnergyDistribution::LoopLengthDistancePotential
    (
      const storage::Vector< math::Histogram> &LOOPLENGTH_DISTRIBUTIONS
    )
    {
      // vector that will hold all final energy distributions
      storage::Vector< math::CubicSplineDamped> energy_distributions;

      // iterate over all possible loop lengths
      for
      (
        storage::Vector< math::Histogram>::const_iterator
          itr( LOOPLENGTH_DISTRIBUTIONS.Begin()), itr_end( LOOPLENGTH_DISTRIBUTIONS.End());
        itr != itr_end;
        ++itr
      )
      {
        math::Histogram current_histogram( *itr);

        // remove all 0 bins
        current_histogram.RemoveBinsAfterIndex( current_histogram.GetIndexOfLastInformationContainingBin());

        // create vector with histogram
        linal::Vector< double> distribution( current_histogram.GetHistogram());
        linal::Vector< double> binning( current_histogram.GetBinning());

        // add pseudo count
        // normalize by square distance
        for
        (
          double *bin( binning.Begin()), *dist( distribution.Begin()), *dist_end( distribution.End());
          dist != dist_end;
          ++bin, ++dist
        )
        {
          *dist += s_PseudoCount;
          *dist /= math::Sqr( *bin);
        }

        // normalize by number of bins and total sum to get propensities
        distribution.SetToSum( distribution.GetSize());

        // logarithmize
        for( double *dist( distribution.Begin()), *dist_end( distribution.End()); dist != dist_end; ++dist)
        {
          *dist = -log( *dist);
        }

        // train a cubic spline
        math::CubicSplineDamped energy;
        energy.Train
        (
          current_histogram.GetBoundaries().First() + 0.5 * current_histogram.GetBinSize(),
          current_histogram.GetBinSize(),
          distribution,
          0,
          1
        );

        // insert into energies
        energy_distributions.PushBack( energy);
      }

      // end
      return energy_distributions;
    }

    //! @brief This function reads the aa_distances histogram and parses it into a matrix and returns it
    //! matrix gives for each AAType pair, probability of observing them within the given distance cutoff
    //! @param DISTANCE_CUTOFF distance in angstrom that is the cutoff
    //! @return matrix that contains the aa pair probabilities
    linal::Matrix< double> EnergyDistribution::ReadAAPairDistanceHistogram( const double DISTANCE_CUTOFF)
    {
      linal::Matrix< double> probmatrix( biol::GetAATypes().GetEnumCount(), biol::GetAATypes().GetEnumCount());

      // read file with all histograms for each pair of sstypes
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( "aa_distances.histograms"));

      while( !read.eof())
      {
        math::Histogram dist_histogram;
        std::pair< biol::AAType, biol::AAType> aatypepair;
        std::string tmp;
        read >> tmp; aatypepair.first = biol::GetAATypes().AATypeFromOneLetterCode( tmp[ 0]);
        read >> tmp; aatypepair.second = biol::GetAATypes().AATypeFromOneLetterCode( tmp[ 0]);

        //abort if one aatype is unknown
        if( aatypepair.first == biol::GetAATypes().e_Undefined || aatypepair.second == biol::GetAATypes().e_Undefined) break;
        read >> dist_histogram;

        double probability( dist_histogram.GetCountsInBetween( 0, DISTANCE_CUTOFF) / dist_histogram.GetSumOfAllCounts());
        probmatrix( aatypepair.first, aatypepair.second) = probability;
        probmatrix( aatypepair.second, aatypepair.first) = probability;
      }

      io::File::CloseClearFStream( read);

      return probmatrix;
    }

    linal::Matrix< double> EnergyDistribution::ProbabilityOfAAPairTOBeWithinDistance( const double &DISTANCE)
    {
      // read file with all histograms for each pair of sstypes
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( "aa_distances.histograms"));

      linal::Matrix< double> probabilitymatrix
      (
        biol::GetAATypes().GetEnumCount(), biol::GetAATypes().GetEnumCount()
      );

      while( !read.eof())
      {
        math::Histogram current_aa_distance_histogram;
        std::pair< biol::AAType, biol::AAType> aatypepair;
        std::string tmp;
        read >> tmp; aatypepair.first = biol::GetAATypes().AATypeFromOneLetterCode( tmp[ 0]);
        read >> tmp; aatypepair.second = biol::GetAATypes().AATypeFromOneLetterCode( tmp[ 0]);

        //abort if one aatype is unknown
        if( aatypepair.first == biol::GetAATypes().e_Undefined || aatypepair.second == biol::GetAATypes().e_Undefined) break;
        read >> current_aa_distance_histogram;

        //calculate probability
        linal::Vector< double> aa_distance_counts( current_aa_distance_histogram.GetHistogram());
        const linal::Vector< double> distance_bins( current_aa_distance_histogram.GetBinning());

        double totalcounts( 0.0);

        //devide each count in bin by distance distribution square
        const double *bins = distance_bins.Begin();
        for
        (
          double *counts = aa_distance_counts.Begin();
          counts != aa_distance_counts.End() && bins != distance_bins.End();
          counts++, bins++
        )
        {
          //divide by distance^2
          *counts /= math::Sqr( *bins);

          totalcounts += ( *counts);
        }

        double counts_in_distance_range( 0.0);
        //add up all normalized counts till the DISTANCE
        for( size_t i( 0); i < aa_distance_counts.GetSize(); ++i)
        {
          if( distance_bins( i) > DISTANCE)
          {
            break;
          }

          counts_in_distance_range += aa_distance_counts( i);
        }

        const double probability( counts_in_distance_range / totalcounts);

        //store probability in matrix
        probabilitymatrix( aatypepair.first, aatypepair.second) = probability;
        probabilitymatrix( aatypepair.second, aatypepair.first) = probability;
      }

      io::File::CloseClearFStream( read);

      return probabilitymatrix;
    }

    //! @brief This function reads a phi psi angle histogram for the 20 amino acid types and generates a potential
    //! @param PHI_PSI_ANGLE_DISTRIBUTION  Histogram 2D that contains amino acid phi psi angle statistics
    //! @param BACKGROUND the background distribution to use
    //! @param NORMALIZE_BY_BACKGROUND normalize by given background - otherwise each bin has a propensity purely derived from this distribution
    //! @return returns a bicubic spline, energydistribution
    math::BicubicSpline EnergyDistribution::PhiPsiAnglePotential
    (
      const math::Histogram2D &PHI_PSI_ANGLE_DISTRIBUTION,
      const math::Histogram2D &BACKGROUND,
      const bool NORMALIZE_BY_BACKGROUND
    )
    {
      // get the core counts
      linal::Matrix< double> phi_psi_counts( PHI_PSI_ANGLE_DISTRIBUTION.GetHistogram());

      // background distribution
      linal::Matrix< double> background_counts( BACKGROUND.GetHistogram());

      // if the background should not be used = equal probability for each phi-psi combination
      if( !NORMALIZE_BY_BACKGROUND)
      {
        background_counts = 1.0 / double( background_counts.GetNumberOfElements());
      }

      // add pseudo count
      phi_psi_counts += 0.00001;
      background_counts += 0.00001;

      // -log of every bin to get the energy
      for( double *ptr( phi_psi_counts.Begin()), *ptr_end( phi_psi_counts.End()), *back( background_counts.Begin()); ptr != ptr_end; ++ptr, ++back)
      {
        *ptr = -log( *ptr / *back);
      }

      //invoke Spline; both borders are periodic due to spherical symmetry
      math::BicubicSpline energy_distribution;
      const math::SplineBorderType           borderflags[ 2] = { math::e_Periodic, math::e_Periodic};
      const double                                 start[ 2] =
        {
          PHI_PSI_ANGLE_DISTRIBUTION.GetBoundariesY().First() + 0.5 * PHI_PSI_ANGLE_DISTRIBUTION.GetBinSizeXY().Second(),
          PHI_PSI_ANGLE_DISTRIBUTION.GetBoundariesX().First() + 0.5 * PHI_PSI_ANGLE_DISTRIBUTION.GetBinSizeXY().First()
        };
      const double                               binsize[ 2] =
        {
          PHI_PSI_ANGLE_DISTRIBUTION.GetBinSizeXY().Second(),
          PHI_PSI_ANGLE_DISTRIBUTION.GetBinSizeXY().First()
        };
      const bool                                 lincont[ 2] = { false, false};
      const storage::Pair< double, double> firstorder[ 2] =
        { storage::Pair< double, double>( 0, 0), storage::Pair< double, double>( 0, 0)};
      energy_distribution.Train( borderflags, start, binsize, phi_psi_counts, lincont, firstorder);

      // end
      return energy_distribution;
    }

    //! @brief convert a histogram into energy function with a cosine background distribution
    //! @param HISTOGRAM from 0 to pi/2
    //! @return Cubic spline representing the energy
    math::CubicSplineDamped EnergyDistribution::AngleAlignmentPotential( const math::Histogram &HISTOGRAM)
    {
      // current histogram, binsize and counts
      const double bin_size( HISTOGRAM.GetBinSize());
      linal::Vector< double> counts( HISTOGRAM.GetHistogram());

      // initialize angle
      double angle( bin_size / double( 2.0));

      // iterate over current histogram
      for
      (
        double *ptr( counts.Begin()), *ptr_end( counts.End());
        ptr != ptr_end; ++ptr, angle += bin_size
      )
      {
        // add pseudo count
        *ptr += double( 1.0);
        // divide by background probability
        *ptr /= cos( angle);
      }

      // calculate propensities
      counts.SetToSum( counts.GetSize());

      // -log to get the energy
      for( double *ptr( counts.Begin()), *ptr_end( counts.End()); ptr != ptr_end; ++ptr)
      {
        *ptr = -log( *ptr);
      }

      // calculate spline
      return math::CubicSplineDamped().Train
      (
        HISTOGRAM.GetBoundaries().First() + 0.5 * bin_size,
        bin_size,
        counts,
        0,
        0
      );
    }

    //! @brief This function reads a histogram for an amino acid and generates a potential
    //! @param HISTOGRAM Histogram that contains statistics =
    //! @param PSEUDOCOUNT pseudocount to be added
    //! @param BORDER_FLAG Training method to be used for handling the borders
    //! @param FIRST_DERIVATIVE only to be used when e_FirstDer is being used
    //! @param NORMALIZE_ENERGIES if true normalizes the energies calculated from the histogram to be between 0 and 1
    //! @param BONUS_ONLY if true shifts the energy potential to be 0 or less (i.e. no penalties)
    //! @return returns a cubic spline
    math::CubicSplineDamped EnergyDistribution::GeneratePotentialFromHistogram
    (
      const math::Histogram &HISTOGRAM,
      const double PSEUDOCOUNT,
      const math::SplineBorderType BORDER_FLAG,
      const storage::Pair< double, double> &FIRST_DERIVATIVE,
      const bool NORMALIZE_ENERGIES,
      const bool BONUS_ONLY
    )
    {
      // create math::Vector values and initialize with the values in the histogram provided
      linal::Vector< double> values( HISTOGRAM.GetHistogram());

      // add pseudocount
      values += PSEUDOCOUNT;

      // variable that adds up all counts
      const double total_sum( values.Sum());

      // convert to propensities
      values *= values.GetSize() / total_sum;

      //-log of every bin to get the energy
      for( double *ptr( values.Begin()), *ptr_end( values.End()); ptr != ptr_end; ++ptr)
      {
        *ptr = -log( *ptr);
      }

      // true if "NORMALIZE_ENERGIES" is set to true - need to normalize the energies of "values" to be between 0 and 1
      if( NORMALIZE_ENERGIES)
      {
        // divide the energies of "values" by its range
        values /=
          (
            math::Statistics::MaximumValue( values.Begin(), values.End()) -
            math::Statistics::MinimumValue( values.Begin(), values.End())
          );
      }

      // true if "BONUS_ONLY" is set to true - need to shift the energy potential to be 0 or less (i.e. no penalties)
      if( BONUS_ONLY)
      {
        // subtract the largest value of "values" from each element of "values"
        values -= math::Statistics::MaximumValue( values.Begin(), values.End());
      }

      // create CubicSpline energy_distribution
      math::CubicSplineDamped energy_distribution;

      // create const double "" and initialize with the size of the bins in "HISTOGRAM"
      const double bin_size( HISTOGRAM.GetBinSize());

      const double start( HISTOGRAM.GetBoundaries().First() + 0.5 * bin_size);
      // train "energy_distribution"
      if( BORDER_FLAG == math::e_Natural)
      {
        energy_distribution.Train( start, bin_size, values);
      }
      else if( BORDER_FLAG == math::e_FirstDer)
      {
        energy_distribution.Train
        (
          start, //< value of the first full bin
          bin_size,
          values,
          FIRST_DERIVATIVE.First(),
          FIRST_DERIVATIVE.Second()
        );
      }
      else
      {
        BCL_Exit( "cubic spline damped does not currently support periodic or not-a-knot border types", -1);
      }

      // end
      return energy_distribution;
    }

  } // namespace score
} // namespace bcl
