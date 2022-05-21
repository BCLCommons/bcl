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
#include "score/bcl_score_body_extent_agreement.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> BodyExtentAgreement::s_Instance
    (
      GetObjectInstances().AddInstance( new BodyExtentAgreement())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BodyExtentAgreement::BodyExtentAgreement() :
      m_LowerTolerance(),
      m_LowerTransitionWidth(),
      m_UpperTolerance(),
      m_UpperTransitionWidth(),
      m_EnergyWellDepth(),
      m_Axis()
    {
    }

    //! @brief constructor
    //! @param LOWER_TOLERANCE the amount an extent can be shorter than the restraint and still be considered perfect
    //! @param LOWER_TRANSITION_WIDTH width of range over which a too short extent goes from perfect to non-agreeing
    //! @param UPPER_TOLERANCE the amount an extent can be longer than the restraint and still be considered perfect
    //! @param UPPER_TRANSITION_WIDTH width of range over which a too long extent goes from perfect to non-agreeing
    //! @param ENERGY_WELL_DEPTH the amount that an agreeing restraint will give a bonus to the energy
    //! @param AXIS defines which extent (x, y, z) is being checked for restraint agreement
    BodyExtentAgreement::BodyExtentAgreement
    (
      const double LOWER_TOLERANCE, const double LOWER_TRANSITION_WIDTH, const double UPPER_TOLERANCE,
      const double UPPER_TRANSITION_WIDTH, const double ENERGY_WELL_DEPTH, const coord::Axis &AXIS
    ) :
      m_LowerTolerance( LOWER_TOLERANCE),
      m_LowerTransitionWidth( LOWER_TRANSITION_WIDTH),
      m_UpperTolerance( UPPER_TOLERANCE),
      m_UpperTransitionWidth( UPPER_TRANSITION_WIDTH),
      m_EnergyWellDepth( ENERGY_WELL_DEPTH),
      m_Axis( AXIS)
    {
    }

    //! @brief Clone is the virtual copy constructor
    BodyExtentAgreement *BodyExtentAgreement::Clone() const
    {
      return new BodyExtentAgreement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BodyExtentAgreement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the scheme
    //! @return string that describes that score
    const std::string &BodyExtentAgreement::GetScheme() const
    {
      static const std::string s_scheme( "body_agreement");
      return s_scheme;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &BodyExtentAgreement::GetAlias() const
    {
      static const std::string s_name( "BodyExtentAgreement");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer BodyExtentAgreement::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores agreement between two coord::Bodies");
      serializer.AddInitializer
      (
        "lower tolerance",
        "the amount an extent can be shorter than the restraint and still be considered perfect",
        io::Serialization::GetAgent( &m_LowerTolerance)
      );
      serializer.AddInitializer
      (
        "lower transition width",
        "range over which a too short extent goes from perfect to non-agreeing",
        io::Serialization::GetAgent( &m_LowerTransitionWidth)
      );
      serializer.AddInitializer
      (
        "upper tolerance",
        "the amount an extent can be longer than the restraint and still be considered perfect",
        io::Serialization::GetAgent( &m_UpperTolerance)
      );
      serializer.AddInitializer
      (
        "upper transition width",
        "range over which a too long extent goes from perfect to non-agreeing",
        io::Serialization::GetAgent( &m_UpperTransitionWidth)
      );
      serializer.AddInitializer
      (
        "energy well depth",
        "amount that an agreeing restraint will give a bonus to the energy",
        io::Serialization::GetAgent( &m_EnergyWellDepth)
      );
      serializer.AddInitializer
      (
        "axis",
        "which extent is being checked for restraint agreement",
        io::Serialization::GetAgent( &m_Axis)
      );

      return serializer;
    }
 
  ////////////////
  // operations //
  ////////////////

    //! @brief operator() which takes a density rod and an assigned SSE
    //! @param BODY the density rod
    //! @param SSE the secondary structure element assigned to that body
    //! @return return a double which is the agreement the body with the sse
    double BodyExtentAgreement::operator()( const assemble::SSEGeometryInterface &BODY, const assemble::SSE &SSE) const
    {
      // create double "deviation" initialize with the difference between extents of the restraint and the test body
      const double deviation( GetExtentDeviation( BODY, SSE));

      // true if extent completely outside of restraint
      if
      (
        deviation < -m_LowerTolerance - m_LowerTransitionWidth || deviation > m_UpperTolerance + m_UpperTransitionWidth
      )
      {
        BCL_MessageDbg
        (
          "completely outside of BodyExtentAgreement restraint : deviation : " + util::Format()( deviation)
        );
        BCL_MessageDbg
        (
          "-m_LowerTolerance - m_LowerTransitionWidth: " + util::Format()( -m_LowerTolerance - m_LowerTransitionWidth)
        );
        BCL_MessageDbg
        (
          "m_UpperTolerance + m_UpperTransitionWidth : " + util::Format()( m_UpperTolerance + m_UpperTransitionWidth)
        );

//        // no energy bonus given
//        return 0.0;

        // give penalty for large extent deviation (penalty is linearly proportional to extent deviation)
        // in case the length of the restraint body (density) is shorter than length of the test body (helix)
        if( deviation < -m_LowerTolerance - m_LowerTransitionWidth)
        {
          return ( math::Absolute( deviation) - m_LowerTolerance - m_LowerTransitionWidth);
        }
        // in case the length of the restraint body (density) is longer than length of the test body (helix)
        else if( deviation > m_UpperTolerance + m_UpperTransitionWidth)
        {
          return ( deviation - m_UpperTolerance - m_UpperTransitionWidth);
        }
      }
      // true if extent completely agrees with restraint
      else if( deviation >= -m_LowerTolerance && deviation <= m_UpperTolerance)
      {
        BCL_MessageDbg
        (
          "completely agreeing with BodyExtentAgreement restraint : deviation : " + util::Format()( deviation)
        );
        BCL_MessageDbg( "-m_LowerTolerance : " + util::Format()( -m_LowerTolerance));
        BCL_MessageDbg( "m_UpperTolerance  : " + util::Format()( m_UpperTolerance));

        // entire possible energy bonus given
        return m_EnergyWellDepth;
      }
      // true if extent is slightly less than restraint
      else if( deviation >= -m_LowerTolerance - m_LowerTransitionWidth && deviation < -m_LowerTolerance)
      {
        BCL_MessageDbg
        (
          "slightly less than restraint : deviation : " + util::Format()( deviation)
        );
        BCL_MessageDbg
        (
          "-m_LowerTolerance - m_LowerTransitionWidth: " + util::Format()( -m_LowerTolerance - m_LowerTransitionWidth)
        );
        BCL_MessageDbg
        (
          "m_LowerTolerance : " + util::Format()( m_LowerTolerance)
        );

        // create double "violation" and initialize with the difference between "m_LowerTolerance" and "deviation"
        double violation( -m_LowerTolerance - deviation);
        BCL_MessageDbg( "violation : " + util::Format()( violation));
        return 0.5 * m_EnergyWellDepth * cos( math::g_Pi * violation / m_LowerTransitionWidth) + 0.5 * m_EnergyWellDepth;
      }
      // true if slightly larger than restraint
      else if( deviation <= m_UpperTolerance + m_UpperTransitionWidth && deviation > m_UpperTolerance)
      {
        BCL_MessageDbg
        (
          "slightly larger than restraint : deviation : " + util::Format()( deviation)
        );
        BCL_MessageDbg
        (
          "m_UpperTolerance + m_UpperTransitionWidth : " + util::Format()( m_UpperTolerance + m_UpperTransitionWidth)
        );
        BCL_MessageDbg( "m_UpperTolerance : " + util::Format()( m_UpperTolerance));
        double violation( deviation - m_UpperTolerance);
        BCL_MessageDbg( "violation : " + util::Format()( violation));
        return 0.5 * m_EnergyWellDepth * cos( math::g_Pi * violation / m_UpperTransitionWidth) + 0.5 * m_EnergyWellDepth;
      }
      else
      {
        BCL_Exit( "could not determine the agreement of the to bodies", -1);
      }

      return 0.0;
    }

    //! @brief GetExtentDeviation gives the value of the extent position difference for an axis
    //!        between two bodies i.e. restraint body and test body
    //! @param BODY_A is the first coord::GeometryInterface used in the calculation
    //! @param BODY_B is the second coord::GeometryInterface used in the calculation
    //! @return returns the value of the difference in position of the extent of an axis
    double BodyExtentAgreement::GetExtentDeviation
    (
      const coord::GeometryInterface &BODY_A, const coord::GeometryInterface &BODY_B
    ) const
    {
      double body_a_extent( 2 * BODY_A.GetExtent( m_Axis));
      double body_b_extent( 2 * BODY_B.GetExtent( m_Axis));
      double difference( body_a_extent - body_b_extent);

      BCL_MessageDbg( "body_a_extent : " + util::Format()( body_a_extent));
      BCL_MessageDbg( "body_b_extent : " + util::Format()( body_b_extent));
      BCL_MessageDbg
      (
        "difference of extents: " + util::Format()( difference)
      );
      return difference;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &BodyExtentAgreement::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_LowerTolerance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LowerTransitionWidth, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UpperTolerance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UpperTransitionWidth, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EnergyWellDepth, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Axis, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &BodyExtentAgreement::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_LowerTolerance, ISTREAM);
      io::Serialize::Read( m_LowerTransitionWidth, ISTREAM);
      io::Serialize::Read( m_UpperTolerance, ISTREAM);
      io::Serialize::Read( m_UpperTransitionWidth, ISTREAM);
      io::Serialize::Read( m_EnergyWellDepth, ISTREAM);
      io::Serialize::Read( m_Axis, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace score
} // namespace bcl
