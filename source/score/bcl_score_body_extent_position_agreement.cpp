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
#include "score/bcl_score_body_extent_position_agreement.h"

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
    const util::SiPtr< const util::ObjectInterface> BodyExtentPositionAgreement::s_Instance
    (
      GetObjectInstances().AddInstance( new BodyExtentPositionAgreement())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BodyExtentPositionAgreement::BodyExtentPositionAgreement() :
      m_LowerTolerance( 1.5),
      m_LowerTransitionWidth( 3.0),
      m_UpperTolerance( 1.5),
      m_UpperTransitionWidth( 3.0),
      m_EnergyWellDepth( -1.0)
    {
    }

    //! @brief Clone is the virtual copy constructor
    BodyExtentPositionAgreement *BodyExtentPositionAgreement::Clone() const
    {
      return new BodyExtentPositionAgreement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BodyExtentPositionAgreement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object
    //! @return the name of the object
    const std::string &BodyExtentPositionAgreement::GetAlias() const
    {
      static const std::string s_name( "BodyExtentPositionAgreement");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer BodyExtentPositionAgreement::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores the agreement of two coord::Bodies based on the deviation between their x, y , and z extents");
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

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator() which takes a density rod and an assigned SSE
    //! @param BODY the density rod
    //! @param SSE the secondary structure element assigned to that body
    //! @return return a double which is the agreement the body with the sse
    double BodyExtentPositionAgreement::operator()( const assemble::SSEGeometryInterface &BODY, const assemble::SSE &SSE) const
    {
      BCL_Assert( BODY.GetCenter() == SSE.GetCenter(), "origin of bodies is not the same");
//      cost double x_deviation( GetAbsoluteExtentDeviation( BODIES.First(), BODIES.Second(), coord::GetAxes().e_X));
//      const double y_deviation( GetAbsoluteExtentDeviation( BODIES.First(), BODIES.Second(), coord::GetAxes().e_Y));
      const double z_deviation( GetAbsoluteExtentDeviation( BODY, SSE, coord::GetAxes().e_Z));

      // true if completely outside of restraint
      if( z_deviation > m_LowerTolerance + m_LowerTransitionWidth || z_deviation > m_UpperTolerance + m_UpperTransitionWidth)
      {
        return 0;
      }
      // true if completely agreeing with restraint
      else if( z_deviation < m_LowerTolerance && z_deviation < m_UpperTolerance)
      {
        return m_EnergyWellDepth;
      }
      // true if slightly less than restraint
      else if( z_deviation < m_LowerTolerance + m_LowerTransitionWidth)
      {
        double violation( m_LowerTolerance - z_deviation);
        return 0.5 * m_EnergyWellDepth * cos( math::g_Pi * violation / m_LowerTransitionWidth) + 0.5 * m_EnergyWellDepth;
      }
      // true if slightly larger than restraint
      else if( z_deviation < m_UpperTolerance + m_UpperTransitionWidth)
      {
        double violation( z_deviation - m_UpperTolerance);
        return 0.5 * m_EnergyWellDepth * cos( math::g_Pi * violation / m_UpperTransitionWidth) + 0.5 * m_EnergyWellDepth;
      }
      else
      {
        BCL_Exit( "could not determine the agreement of the to bodies", -1);
      }

      return 0.0;
    }

    //! @brief GetAbsoluteExtentDeviation gives the absolute value of the extent position difference for an axis
    //! @param BODY_A is the first coord::GeometryInterface used in the calculation
    //! @param BODY_B is the second coord::GeometryInterface used in the calculation
    //! @param AXIS is the RotationAxis for which the calculation will be performed
    //! @return returns the absolute value of the difference in position of the extent of an axis
    double BodyExtentPositionAgreement::GetAbsoluteExtentDeviation
    (
      const coord::GeometryInterface &BODY_A, const coord::GeometryInterface &BODY_B, const coord::Axis &AXIS
    ) const
    {
      double body_a_extent_position( GetExtentPosition( BODY_A, AXIS));
      double body_b_extent_position( GetExtentPosition( BODY_B, AXIS));
      double difference( math::Absolute( body_a_extent_position - body_b_extent_position));
      BCL_MessageDbg( "body_a_extent_position : " + util::Format()( body_a_extent_position));
      BCL_MessageDbg( "body_b_extent_position : " + util::Format()( body_b_extent_position));
      BCL_MessageDbg
      (
        "BodyExtentPositionAgreement::GetAbsoluteExtentDeviation difference: "
        + util::Format()( difference)
      );
      return math::Absolute( body_a_extent_position - body_b_extent_position);
    }

    //! @brief GetExtentPosition calculates the position of an extent for a defined axis for a coord::GeometryInterface
    //! @param BODY the coord::GeometryInterface for which the position of its extent in a given axis will be calculated
    //! @param AXIS is the RotationAxis for which the extent position will be calculated
    //! @return double which is the value of the origin in the axis "XYZ" plus one half the extent of BODY in that
    //!         axis direction
    double BodyExtentPositionAgreement::GetExtentPosition( const coord::GeometryInterface &BODY, const coord::Axis &AXIS) const
    {
//      return BODY.GetOrigin()( XYZ) + BODY.GetExtent( XYZ) / 2;
        return BODY.GetExtent( AXIS);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &BodyExtentPositionAgreement::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &BodyExtentPositionAgreement::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

  } // namespace score
} // namespace bcl
