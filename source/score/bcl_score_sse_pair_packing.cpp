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
#include "score/bcl_score_sse_pair_packing.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_angle.h"
#include "math/bcl_math_histogram_2d.h"
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
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &SSEPairPacking::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "sse_angle_distance.histograms2D");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &SSEPairPacking::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "ssepack");

      // end
      return s_default_scheme;

    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param SCHEME scheme to be used
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param DISTANCE_RANGE_MAP map of distance ranges for valid contact types
    SSEPairPacking::SSEPairPacking
    (
      const std::string &SCHEME,
      const std::string &HISTOGRAM_FILENAME,
      const storage::Map< contact::Type, math::Range< double> > &DISTANCE_RANGE_MAP
    ) :
      m_Scheme( SCHEME),
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_DistanceRangeMap( DISTANCE_RANGE_MAP),
      m_MinimalInterfaceLength(),
      m_EnergyFunctions()
    {
      // read the histogram file and store the energy functions
      ReadEnergyVector();
    }

    //! @brief virtual copy constructor
    SSEPairPacking *SSEPairPacking::Clone() const
    {
      return new SSEPairPacking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPairPacking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that score packing of a pair of SSEs of interest
    //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest
    //! @return potential of interaction
    double SSEPairPacking::operator()
    (
      const assemble::SSEGeometryPacking &SSE_PACK
    ) const
    {
      //if packing weight is small just return no score
      if( SSE_PACK.GetInteractionWeight() < std::numeric_limits< double>::epsilon())
      {
        return 0.0;
      }

      // store contact type
      contact::Type contact_type( SSE_PACK.GetContactType());

      // switch over contact type
      // if UNKNOWN
      if( contact_type == contact::GetTypes().e_Undefined)
      {
        return 0.0;
      }
      // if HELIX_HELIX
      else if
      (
        contact_type == contact::GetTypes().HELIX_HELIX
      )
      {
        return
        (
          SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
          *
          (
            m_EnergyFunctions.Find( SSE_PACK.GetContactType())->second->F
            (
              linal::MakeVector( SSE_PACK.GetDistance(), SSE_PACK.GetTwistAngle())
            )
          )
        );
      }
      // for any helix and strand (sheet) combination
      else if
      (
        contact_type == contact::GetTypes().HELIX_SHEET            ||
        contact_type == contact::GetTypes().SHEET_HELIX            ||
        contact_type == contact::GetTypes().UNDEFINED_HELIX_STRAND ||
        contact_type == contact::GetTypes().UNDEFINED_STRAND_HELIX ||
        contact_type == contact::GetTypes().HELIX_STRAND           ||
        contact_type == contact::GetTypes().STRAND_HELIX
      )
      {
        return
        (
          SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
          *
          (
            m_EnergyFunctions.Find( contact::GetTypes().HELIX_SHEET)->second->F
            (
              linal::MakeVector( SSE_PACK.GetDistance(), SSE_PACK.GetTwistAngle())
            )
          )
          +
          SSE_PACK.GetInteractionWeight() * ( 1.0 - SSE_PACK.GetRelativePositionWeight())
          *
          (
            m_EnergyFunctions.Find( contact::GetTypes().HELIX_STRAND)->second->F
            (
              linal::MakeVector( SSE_PACK.GetDistance(), SSE_PACK.GetTwistAngle())
            )
          )
        );
      }
      // if STRAND_STRAND, SHEET_SHEET or UNDEFINED_STRAND_STRAND
      else if
      (
        contact_type == contact::GetTypes().STRAND_STRAND ||
        contact_type == contact::GetTypes().SHEET_SHEET   ||
        contact_type == contact::GetTypes().UNDEFINED_STRAND_STRAND
      )
      {
        //apply transition function
        return
        (
          SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
          *
          (
            m_EnergyFunctions.Find( contact::GetTypes().SHEET_SHEET)->second->F
            (
              linal::MakeVector( SSE_PACK.GetDistance(), SSE_PACK.GetTwistAngle())
            )
          )
        );
      }
      // if none above
      else
      {
        return util::GetUndefined< double>();
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &SSEPairPacking::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_HistogramFileName, ISTREAM);
      io::Serialize::Read( m_DistanceRangeMap, ISTREAM);

      // read the histogram file and store the energy functions
      m_EnergyFunctions.Reset();
      ReadEnergyVector();

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &SSEPairPacking::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HistogramFileName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceRangeMap, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &SSEPairPacking::WriteDetailedSchemeAndValues
    (
      const assemble::SSEGeometryPacking &SSE_PACK,
      std::ostream &OSTREAM
    ) const
    {
      // write sstype, Angle to membrane plane to the STREAM
      OSTREAM << SSE_PACK.GetContactType().GetName() << '\t'
              << SSE_PACK.GetInteractionWeight() << '\t'
              << math::Angle::Degree( SSE_PACK.GetRelativePosition()) << '\t'
              << SSE_PACK.GetRelativePositionWeight() << '\t'
              << SSE_PACK.GetDistance() << '\t'
              << math::Angle::Degree( SSE_PACK.GetTwistAngle()) << '\t'
              << operator()( SSE_PACK) << '\n';

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read energy distribution for scoring pairs of sses
    void SSEPairPacking::ReadEnergyVector()
    {
      // read file with all histograms for each pair of sstypes
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      // read minimal interface length and set the member variable to that one
      read >> m_MinimalInterfaceLength;

      // while reading contact type and a histogram not reaching the end of file
      while( read.good())
      {
        // initialize contact type
        contact::Type contact_type;

        // initialize histogram to read
        math::Histogram2D current_angle_distance_histogram;

        // read the contact type and the histogram
        read >> contact_type >> current_angle_distance_histogram;

        // calculate energydistribution, write it in spline and store it in map - also with swapped types

        util::ShPtr< math::BicubicSpline> current_spline;

        // for fragment packing, use the 2d energy distributions
        if( m_Scheme.find( "_fr") != std::string::npos)
        {
          current_spline = util::ShPtr< math::BicubicSpline>
          (
            new math::BicubicSpline
            (
              EnergyDistribution::SSEPacking2D( current_angle_distance_histogram, m_DistanceRangeMap.Find( contact_type)->second)
            )
          );
        }
        // use the symmetrize
        else
        {
          current_spline = util::ShPtr< math::BicubicSpline>
          (
            new math::BicubicSpline
            (
              EnergyDistribution::SSEPackingSymmetrize( current_angle_distance_histogram, m_DistanceRangeMap.Find( contact_type)->second)
            )
          );
        }

        m_EnergyFunctions[ contact_type] = current_spline;
        if( contact_type == contact::GetTypes().HELIX_SHEET)
        {
          m_EnergyFunctions[ contact::GetTypes().SHEET_HELIX] = current_spline;
        }
        else if( contact_type == contact::GetTypes().HELIX_STRAND)
        {
          m_EnergyFunctions[ contact::GetTypes().STRAND_HELIX] = current_spline;
        }

        // end of all histograms data
        if( contact_type == contact::GetTypes().SHEET_SHEET)
        {
          break;
        }
      }

      //close input stream
      io::File::CloseClearFStream( read);
    }

  } // namespace score
} // namespace bcl
