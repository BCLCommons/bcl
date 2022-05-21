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
#include "score/bcl_score_strand_pairing.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "score/bcl_score_energy_distribution.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> StrandPairing::s_Instance
    (
      util::Enumerated< SSEPackInterface>::AddInstance( new StrandPairing())
    );

  //////////
  // data //
  //////////

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &StrandPairing::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "strand_angle_distance.histograms2D");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &StrandPairing::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "strand");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param SCHEME scheme to be used
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param DISTANCE_RANGE distance range
    StrandPairing::StrandPairing
    (
      const std::string &SCHEME,
      const std::string &HISTOGRAM_FILENAME,
      const math::Range< double> &DISTANCE_RANGE
    ) :
      m_Scheme( SCHEME),
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_DistanceRange( DISTANCE_RANGE),
      m_MinimalInterfaceLength(),
      m_EnergyFunction()
    {
      if( !command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        // read the histogram file and store the energy functions
        ReadEnergyVector();
      }
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new StrandPairing object that is copied from this one
    StrandPairing *StrandPairing::Clone() const
    {
      return new StrandPairing( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StrandPairing::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }
    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &StrandPairing::GetAlias() const
    {
      static const std::string s_name( "StrandPairing");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer StrandPairing::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "A Function derived class for scoring pairing of strands within one sheet");
      serializer.AddInitializer
      (
        "histogram file name",
        "path to file where the statistics and in consequence the energy potentials are read from",
        io::Serialization::GetAgent( &m_HistogramFileName),
        util::Format()( GetDefaultHistogramFilename())
      );
      serializer.AddInitializer
      (
        "distance range",
        "distance range",
        io::Serialization::GetAgent( &m_DistanceRange),
        contact::GetTypes().STRAND_STRAND->GetDistanceRange().GetString()
      );
      serializer.AddInitializer
      (
        "Minimal interface length",
        "minimal interface length for packing to get a full weight",
        io::Serialization::GetAgent( &m_MinimalInterfaceLength),
        "0"
      );
      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns whether the given pair of SSEs are valid
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return whether the given pair of SSEs are valid
    bool StrandPairing::AreValidSSEs( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      // return true only if both SSEs are strands
      return SSE_A.GetType() == biol::GetSSTypes().STRAND && SSE_B.GetType() == biol::GetSSTypes().STRAND;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates strand pair packing potential
    //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest - if not strands just returns 0
    //! @return strand pairing potential
    double
    StrandPairing::operator()
    (
      const assemble::SSEGeometryPacking &SSE_PACK
    ) const
    {
      // if packign weight is small just return no score
      if( SSE_PACK.GetInteractionWeight() < std::numeric_limits< double>::epsilon())
      {
        return 0.0;
      }

      // store the contact type
      contact::Type contact_type( SSE_PACK.GetContactType());

      // switch over contacttype
      if
      (
          contact_type == contact::GetTypes().e_Undefined            ||
          contact_type == contact::GetTypes().HELIX_HELIX            ||
          contact_type == contact::GetTypes().HELIX_SHEET            ||
          contact_type == contact::GetTypes().SHEET_HELIX            ||
          contact_type == contact::GetTypes().HELIX_STRAND           ||
          contact_type == contact::GetTypes().STRAND_HELIX           ||
          contact_type == contact::GetTypes().UNDEFINED_HELIX_STRAND ||
          contact_type == contact::GetTypes().UNDEFINED_STRAND_HELIX
      )
      {
        return double( 0);
      }
      // if strand strand related contact type
      else if
      (
          contact_type == contact::GetTypes().SHEET_SHEET   ||
          contact_type == contact::GetTypes().STRAND_STRAND ||
          contact_type == contact::GetTypes().UNDEFINED_STRAND_STRAND
      )
      {
        // calculate weight for transition region and multiply it with score
        return
            (
                SSE_PACK.GetStrandStrandPairingWeight() * SSE_PACK.GetInteractionWeight()
                *
                (
                    m_EnergyFunction->F( linal::MakeVector( SSE_PACK.GetDistance(), SSE_PACK.GetTwistAngle()))
                )
            );
      }

      // return undefined
      return util::GetUndefined< double>();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &StrandPairing::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_HistogramFileName, ISTREAM);
      io::Serialize::Read( m_DistanceRange, ISTREAM);

      // read the histogram file and store the energy functions
      ReadEnergyVector();

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &StrandPairing::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HistogramFileName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceRange, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    StrandPairing::WriteDetailedSchemeAndValues
    (
      const assemble::SSEGeometryPacking &SSE_PACK,
      std::ostream &OSTREAM
    ) const
    {
      //write sstype, Angel to membrane plane first and last aminoacid and value to the STREAM
      OSTREAM << SSE_PACK.GetContactType().GetName() << '\t'

          // output interaction angleas and weight. it shows how parallel main axes of SSEs are. optimal is 0, 0, 1
          //              << math::Degree( SSE_PACK.AngleSSEConnection( *SSE_PAIR.First(), SSE_PACK.GetShortestConnection()))
          //              << '\t'
          //              << math::Degree( SSE_PACK.AngleSSEConnection( *SSE_PAIR.Second(), SSE_PACK.GetShortestConnection()))
          //              << '\t'
          << SSE_PACK.GetInteractionWeight() << '\t'

          // output planarity within sheet and associated weight, perfect strand strand pair would be 90, 90, 1
          //              << math::Degree( math::Absolute( math::g_Pi / 2 - math::Absolute( math::g_Pi / 2 - linal::ProjAngle( SSE_PACK.GetShortestConnection().GetEndPoint(), SSE_PACK.GetShortestConnection().GetStartPoint(), SSE_PACK.GetShortestConnection().GetEndPoint() + SSE_PAIR.Second()->GetAxis( coord::GetAxes().e_X))))) << '\t'
          //              << math::Degree( math::Absolute( math::g_Pi / 2 - math::Absolute( math::g_Pi / 2 - linal::ProjAngle( SSE_PACK.GetShortestConnection().GetStartPoint(), SSE_PACK.GetShortestConnection().GetEndPoint(), SSE_PACK.GetShortestConnection().GetStartPoint() + SSE_PAIR.First()->GetAxis( coord::GetAxes().e_X))))) << '\t'
          << SSE_PACK.GetStrandStrandPairingWeight() << '\t'

          // output shortest distance and twis angle
          //               << ( 1.0 - math::WeightBetweenZeroAndPi_ThreeSections( SSE_PACK.GetRelativePosition())) << '\t'
          << SSE_PACK.GetDistance() << '\t'
          << math::Angle::Degree( SSE_PACK.GetTwistAngle()) << '\t'
          << operator()( SSE_PACK) << '\n';

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read energy distribution for scoring pairs of  of strands within one sheet
    void StrandPairing::ReadEnergyVector()
    {
      //read file with all histograms for each pair of sstypes
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      // read minimal interface length and set the member variable to that one
      read >> m_MinimalInterfaceLength;

      //read the sstype pair and according histogram of counts
      contact::Type contact_type;
      read >> contact_type;

      // assert the histogram contains the correct contact type (STRAND-STRAND) histogram
      BCL_Assert
      (
        contact_type == contact::GetTypes().STRAND_STRAND,
        "file is not a histogram for strand-strand packing"
      );

      // read the histogram
      math::Histogram2D current_angle_distance_histogram;
      read >> current_angle_distance_histogram;

      //calculate energydistribution, write it in spline and store it in map - also with swapped types
      m_EnergyFunction =
          util::ShPtr< math::BicubicSpline>
      (
        new math::BicubicSpline
        (
          EnergyDistribution::SSEPacking2D( current_angle_distance_histogram, m_DistanceRange)
        )
      );

      //close input stream
      io::File::CloseClearFStream( read);
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool StrandPairing::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      ReadEnergyVector();
      return true;
    }

  } // namespace score
} // namespace bcl

