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
#include "score/bcl_score_sse_pair_clash.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &SSEPairClash::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "sseclash");

      // end
      return s_default_scheme;

    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEPairClash::s_Instance
    (
      util::Enumerated< SSEPackInterface>::AddInstance( new SSEPairClash)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a SCHEME
    //! @param SIGMOID_WIDTH width of the sigmoid repulsive term, before it reaches 1.0
    //! @param SCHEME scheme to be used
    SSEPairClash::SSEPairClash
    (
      const double SIGMOID_WIDTH,
      const std::string &SCHEME
    ) :
      m_Scheme( SCHEME),
      m_MinimalInterfaceLength( 0.0),
      m_SigmoidWidth( SIGMOID_WIDTH)
    {
    }

    //! @brief virtual copy constructor
    SSEPairClash *SSEPairClash::Clone() const
    {
      return new SSEPairClash( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPairClash::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SSEPairClash::GetAlias() const
    {
      static const std::string s_name( "SSEPairClash");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEPairClash::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Evaluates clashes between SSEs.");
      serializer.AddInitializer
      (
        "min interface length",
        "minimum length of the interface between SSEs to receive full weight",
        io::Serialization::GetAgent( &m_MinimalInterfaceLength),
        "0.0"
      );
      serializer.AddInitializer
      (
        "sigmoid width",
        "width of the sigmoidal transition function",
        io::Serialization::GetAgent( &m_SigmoidWidth),
        "1.0"
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that score packing of a pair of SSEs of interest
    //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest
    //! @return potential of interaction
    double SSEPairClash::operator()
    (
      const assemble::SSEGeometryPacking &SSE_PACK
    ) const
    {
      // if packing weight is small just return no score
      if( SSE_PACK.GetInteractionWeight() < std::numeric_limits< double>::epsilon())
      {
        return 0.0;
      }

      // store contact type
      const contact::Type contact_type( SSE_PACK.GetContactType());

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
            CalculateRepulsiveTerm( SSE_PACK.GetDistance(), contact_type->GetMinimalSSEDistance())
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
            CalculateRepulsiveTerm( SSE_PACK.GetDistance(), contact::GetTypes().HELIX_SHEET->GetMinimalSSEDistance())
          )
          +
          SSE_PACK.GetInteractionWeight() * ( 1.0 - SSE_PACK.GetRelativePositionWeight())
          *
          (
            CalculateRepulsiveTerm( SSE_PACK.GetDistance(), contact::GetTypes().HELIX_STRAND->GetMinimalSSEDistance())
          )
        );
      }
      // if strand strand related contact type
      else if
      (
        contact_type == contact::GetTypes().STRAND_STRAND ||
        contact_type == contact::GetTypes().SHEET_SHEET   ||
        contact_type == contact::GetTypes().UNDEFINED_STRAND_STRAND
      )
      {
        // apply transition function
        return
        (
          SSE_PACK.GetInteractionWeight() * SSE_PACK.GetRelativePositionWeight()
          *
          (
            CalculateRepulsiveTerm( SSE_PACK.GetDistance(), contact::GetTypes().SHEET_SHEET->GetMinimalSSEDistance())
          )
          +
          SSE_PACK.GetInteractionWeight() * SSE_PACK.GetStrandStrandPairingWeight()
          *
          (
            CalculateRepulsiveTerm( SSE_PACK.GetDistance(), contact::GetTypes().STRAND_STRAND->GetMinimalSSEDistance())
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
    std::istream &SSEPairClash::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_MinimalInterfaceLength, ISTREAM);
      io::Serialize::Read( m_SigmoidWidth, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &SSEPairClash::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinimalInterfaceLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SigmoidWidth, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &SSEPairClash::WriteDetailedSchemeAndValues
    (
      const assemble::SSEGeometryPacking &SSE_PACK,
      std::ostream &OSTREAM
    ) const
    {
      //write sstype, Angle to membrane plane to the STREAM
      OSTREAM << SSE_PACK.GetContactType().GetName() << '\t'
              << SSE_PACK.GetInteractionWeight() << '\t'
              << math::Angle::Degree( SSE_PACK.GetRelativePosition()) << '\t'
              << SSE_PACK.GetRelativePositionWeight() << '\t'
              << SSE_PACK.GetDistance() << '\t'
              << math::Angle::Degree( SSE_PACK.GetTwistAngle()) << '\t'
              << operator()( SSE_PACK) << '\n';

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief repulsive term from distance and shortest observed distance
    //! @param DISTANCE the actual distance between sses of interest
    //! @param SHORTEST_OBSERVED_DISTANCE shortest distance observed
    double SSEPairClash::CalculateRepulsiveTerm
    (
      const double DISTANCE,
      const double SHORTEST_OBSERVED_DISTANCE
    ) const
    {
      const double difference( SHORTEST_OBSERVED_DISTANCE - DISTANCE);

      // no repulsion
      if( difference <= 0.0)
      {
        return 0.0;
      }
      // full repulsion
      else if( difference >= m_SigmoidWidth)
      {
        return 1.0;
      }

      return math::WeightBetweenZeroAndPi( ( ( m_SigmoidWidth - difference) / m_SigmoidWidth) * math::g_Pi);
    }

  } // namespace score
} // namespace bcl
