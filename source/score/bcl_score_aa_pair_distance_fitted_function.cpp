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
#include "score/bcl_score_aa_pair_distance_fitted_function.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  ///////////
  // data //
  //////////

    //! default start values for the start of the repulsion term for the fitted aapairenergy distribution
    const storage::VectorND< 2, double> AAPairDistanceFittedFunction::s_DefaultRepulsionStartXY( 2.0, 10.0);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAPairDistanceFittedFunction::AAPairDistanceFittedFunction() :
      m_RepulsionStartXY( s_DefaultRepulsionStartXY)
    {
    }

    //! @brief constructor from a given histogram
    //! @param AA_PAIR_DISTANCE_DISTRIBUTION Histogram that containes aa pair distance distribution
    AAPairDistanceFittedFunction::AAPairDistanceFittedFunction
    (
      const math::Histogram &AA_PAIR_DISTANCE_DISTRIBUTION
    ) :
      m_RepulsionStartXY( s_DefaultRepulsionStartXY)
    {
      linal::Vector< double> aa_distance_counts( AA_PAIR_DISTANCE_DISTRIBUTION.GetHistogram());
      const linal::Vector< double> distance_bins( AA_PAIR_DISTANCE_DISTRIBUTION.GetBinning());

      const size_t usedbins( 20);

      //divide by the total number of all pairs
      double totalcounts( 0);

      //divide each count in bin by distance distribution square
      const double *bins( distance_bins.Begin()), *bins_end( distance_bins.End());
      size_t i( 0);
      for
      (
        double *counts( aa_distance_counts.Begin()), *counts_end( aa_distance_counts.End());
        counts != counts_end && bins != bins_end && i < usedbins; ++counts, ++bins, ++i
      )
      {
        //divide by distance^2
        *counts /= math::Sqr( *bins);

        //add pseudo count
        *counts += 1;

        totalcounts += ( *counts);
      }
//      totalcounts += AA_PAIR_DISTANCE_DISTRIBUTION.GetBoundariesCounts().Second() / Sqr( AA_PAIR_DISTANCE_DISTRIBUTION.GetBoundaries().Second());

      linal::Vector< double> energydistributionvector( 20, aa_distance_counts.Begin());

      //-log of every bin to have the energy
      for
      (
        double *ptr( energydistributionvector.Begin()), *ptr_end( energydistributionvector.End());
        ptr != ptr_end; ++ptr
      )
      {
        ( *ptr) /= totalcounts / 20;
        ( *ptr) = -log( *ptr);
      }

//      const double shift = -1.5;

      // shift all values but the last by shift, last value set to zero
//      aa_distance_counts -= shift;
      *( energydistributionvector.End() - 1) = 0;

      //determine minimal value and index;
      const size_t index_energy_minima
      (
        math::Statistics::MinimumIndex( energydistributionvector.Begin(), energydistributionvector.End())
      );
      m_AttractionMinimumXY.First() = double( distance_bins( index_energy_minima));
      m_AttractionMinimumXY.Second() = energydistributionvector( index_energy_minima);

      //determine index of last positive value which is the Repulsion End
      for( size_t i( 0); i < energydistributionvector.GetSize(); ++i)
      {
        if( energydistributionvector( i) < 0)
        {
          m_RepulsionEndXY.First() = distance_bins( i);
          break;
        }
      }
      m_RepulsionEndXY.Second() = 0.0;

      m_AttractionStartXY = storage::VectorND< 2, double>( m_RepulsionStartXY.First(), 0);
      m_AttractionEndXY = storage::VectorND< 2, double>( std::min( 20.0, 2 * m_AttractionMinimumXY.First() - m_RepulsionStartXY.First()), 0);

//        BCL_MessageStd( "m_RepulsionStartXY:    " + util::Format()( m_RepulsionStartXY   ));
//        BCL_MessageStd( "m_RepulsionEndXY:      " + util::Format()( m_RepulsionEndXY     ));
//        BCL_MessageStd( "m_AttractionStartXY:   " + util::Format()( m_AttractionStartXY  ));
//        BCL_MessageStd( "m_AttractionMinimumXY: " + util::Format()( m_AttractionMinimumXY));
//        BCL_MessageStd( "m_AttractionEndXY:     " + util::Format()( m_AttractionEndXY    ));
//
//        for( double dist = 0.0; dist < 25; dist += 1.0)
//        {
//          BCL_MessageStd( util::Format()( dist) + "\t|\t" + util::Format()( Repulsion( dist)) + "\t|\t" + util::Format()( Attraction( dist)) + "\t|\t" + util::Format()( F( dist)));
//        }
    }

    //! @brief construct from a vector of bins and a given energy distribution
    //! @param BINNING binning to be used
    //! @param ENERGY_DISTRIBUTION energy distribution to be used
    AAPairDistanceFittedFunction::AAPairDistanceFittedFunction
    (
      const linal::Vector< double> &BINNING,
      const linal::Vector< double> &ENERGY_DISTRIBUTION
    ) :
      m_RepulsionStartXY( s_DefaultRepulsionStartXY)
    {
      const size_t index_energy_minima
      (
        math::Statistics::MinimumIndex( ENERGY_DISTRIBUTION.Begin(), ENERGY_DISTRIBUTION.End())
      );
      //determine minimal value and index;
      m_AttractionMinimumXY.First() = double( BINNING( index_energy_minima));
      m_AttractionMinimumXY.Second() = ENERGY_DISTRIBUTION( index_energy_minima);

      //determine index of last positive value which is the Repulsion End
      for( size_t i( 0); i < ENERGY_DISTRIBUTION.GetSize(); ++i)
      {
        if( ENERGY_DISTRIBUTION( i) < 0)
        {
          m_RepulsionEndXY.First() = BINNING( i);
          break;
        }
      }

      m_AttractionStartXY = storage::VectorND< 2, double>( m_RepulsionStartXY.First(), 0);
      m_AttractionEndXY = storage::VectorND< 2, double>( std::min( 20.0, 2 * m_AttractionMinimumXY.First() - m_RepulsionStartXY.First()), 0);
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAPairDistanceFittedFunction copied from this one
    AAPairDistanceFittedFunction *AAPairDistanceFittedFunction::Clone() const
    {
      return new AAPairDistanceFittedFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairDistanceFittedFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the energy according to the distance
    //! @param DISTANCE distance to be used
    //! @return the energy calculated for the given distance
    double AAPairDistanceFittedFunction::operator()( const double &DISTANCE) const
    {
//      if( DISTANCE < 15.0)
//      {
//        BCL_Message
//        (
//          util::Message::e_Standard,
//          util::Format()( DISTANCE) + "|" + util::Format()( Repulsion( DISTANCE)) +
//            "|" + util::Format()( Attraction( DISTANCE))
//        );
//      }
      return Repulsion( DISTANCE) + Attraction( DISTANCE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAPairDistanceFittedFunction::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_RepulsionStartXY, ISTREAM);
      io::Serialize::Read( m_RepulsionEndXY, ISTREAM);
      io::Serialize::Read( m_AttractionStartXY, ISTREAM);
      io::Serialize::Read( m_AttractionMinimumXY, ISTREAM);
      io::Serialize::Read( m_AttractionEndXY, ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAPairDistanceFittedFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write member
      io::Serialize::Write( m_RepulsionStartXY, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RepulsionEndXY, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AttractionStartXY, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AttractionMinimumXY, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AttractionEndXY, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates the repulsion for the given distance
    //! @param DISTANCE distance to be used
    //! @return the repulsion for the given distance
    double AAPairDistanceFittedFunction::Repulsion( const double &DISTANCE) const
    {
      if( DISTANCE <= m_RepulsionStartXY.First())
      {
        return m_RepulsionStartXY.Second();
      }

      else if( DISTANCE >= m_RepulsionEndXY.First())
      {
        return m_RepulsionEndXY.Second();
      }

      else
      {
        return
        0.5 * ( m_RepulsionStartXY.Second() - m_RepulsionEndXY.Second()) *
        std::cos
        (
          math::g_Pi *
          ( DISTANCE - m_RepulsionStartXY.First())
          / ( m_RepulsionEndXY.First() - m_RepulsionStartXY.First())
        ) + 0.5 * ( m_RepulsionStartXY.Second() - m_RepulsionEndXY.Second());
      }
    }

    //! @brief calculates the attraction for the given distance
    //! @param DISTANCE distance to be used
    //! @return the attraction for the given distance
    double AAPairDistanceFittedFunction::Attraction( const double &DISTANCE) const
    {
      //if argument is smaller than begin of attraction return start
      if( DISTANCE <= m_AttractionStartXY.First())
      {
        return m_AttractionStartXY.Second();
      }

      //if argument is larger than begin of attraction return end
      else if( DISTANCE >= m_AttractionEndXY.First())
      {
        return m_AttractionEndXY.Second();
      }

      else
      {
        return
        m_AttractionStartXY.Second() +
          ( m_AttractionStartXY.Second() - m_AttractionMinimumXY.Second()) *
          0.5 *
          (
            std::cos
            (
              2 * math::g_Pi *
              (
                ( DISTANCE - m_AttractionStartXY.First()) /
                ( m_AttractionEndXY.First() - m_AttractionStartXY.First())
              )
             ) - 1
          );
      }
    }

  } // namespace score
} // namespace bcl
