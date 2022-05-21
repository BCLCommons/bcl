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
#include "biol/bcl_biol_chi_angle.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    //! @brief conversion to a string from a Chi
    //! @param CHI the chi angle to get a string for
    //! @return a string representing that chi
    const std::string &ChiAngle::GetChiName( const ChiAngle::Chi &CHI)
    {
      static const std::string s_descriptors[] =
      {
        "e_One",
        "e_Two",
        "e_Three",
        "e_Four",
        "e_Five",
        "e_Undefined",
        GetStaticClassName< ChiAngle>()
      };
      return s_descriptors[ size_t( CHI)];
    }

      //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ChiAngle::s_Instance
    (
      GetObjectInstances().AddInstance( new ChiAngle())
    );

    //! @brief gives the maximum magnitude that a chi angle could have
    //! @param UNIT the unit the max angle should be given in
    //! @return double which is the maximum magnitude that a chi angle could have
    double ChiAngle::GetMaxAngle( const math::Angle::Unit &UNIT)
    {
      if( UNIT == math::Angle::e_Degree)
      {
        return 180;
      }
      if( UNIT == math::Angle::e_Radian)
      {
        return math::g_Pi;
      }

      return util::GetUndefinedDouble();
    }

    //! @brief gives the total number of angular units in a circle
    //! @param UNIT the unit the angular units are in
    //! @return double the total number of angular units in a circle
    double ChiAngle::GetCircle( const math::Angle::Unit &UNIT)
    {
      return GetMaxAngle( UNIT) * 2.0;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ChiAngle::ChiAngle() :
      m_Chi( ChiAngle::e_Undefined),
      m_Angle( util::GetUndefinedDouble()),
      m_Unit( math::Angle::e_Undefined)
    {
    }

    //! @brief constructor from member variable parameters
    //! @param CHI the chi this will correspond to
    ChiAngle::ChiAngle( const Chi &CHI) :
      m_Chi( CHI),
      m_Angle( util::GetUndefinedDouble()),
      m_Unit( math::Angle::e_Undefined)
    {
    }

    //! @brief constructor from member variable parameters
    //! @param CHI the chi this will correspond to
    //! @param ANGLE the angle of this chi
    //! @param UNIT the unit type the angle is measured in
    ChiAngle::ChiAngle( const Chi &CHI, const double ANGLE, const math::Angle::Unit &UNIT) :
      m_Chi( CHI),
      m_Angle( ANGLE),
      m_Unit( UNIT)
    {
      BCL_Assert
      (
        math::Absolute( GetAngle( m_Unit) <= GetMaxAngle( m_Unit)) || !util::IsDefined( ANGLE),
        "chi angle in units of " + m_Unit.GetString() +
        " should not be larger in magnitude than " + util::Format()( GetMaxAngle( m_Unit))
      );
    }

    //! @brief Clone function
    //! @return pointer to new ChiAngle
    ChiAngle *ChiAngle::Clone() const
    {
      return new ChiAngle( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ChiAngle::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the chi this corresponds to
    //! @return Chi which is the chi this corresponds to
    const ChiAngle::Chi &ChiAngle::GetChi() const
    {
      return m_Chi;
    }

    //! @brief gives the angular value of this angle
    //! @param ANGLE_UNIT the unit the angle should be given in
    //! @return double which is the angular value of this angle in the units desired according ANGLE_UNIT
    double ChiAngle::GetAngle( const math::Angle::Unit &ANGLE_UNIT) const
    {
      // true if angle stored as degrees and degrees are desired
      if( m_Unit == math::Angle::e_Degree && ANGLE_UNIT == math::Angle::e_Degree)
      {
        // return the angle
        return m_Angle;
      }
      // true if angle stored as radians and degrees are desired
      if( m_Unit == math::Angle::e_Radian && ANGLE_UNIT == math::Angle::e_Degree)
      {
        // convert the radians to degrees and return value
        return math::Angle::Degree( m_Angle);
      }
      // true if angle stored as degrees and radians are desired
      if( m_Unit == math::Angle::e_Degree && ANGLE_UNIT == math::Angle::e_Radian)
      {
        // convert the degrees to radian and return value
        return math::Angle::Radian( m_Angle);
      }
      // true if angle stored as radians and radians are desired
      if( m_Unit == math::Angle::e_Radian && ANGLE_UNIT == math::Angle::e_Radian)
      {
        // return the angle
        return m_Angle;
      }

      return util::GetUndefinedDouble();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the absolute angular difference between this chi angle and a provided chi angle
    //! @param CHI_ANGLE the chi angle whose difference will be calculated from this
    //! @param ANGLE_UNIT the unit the angular difference should be given in
    //! @return double which is the absolute angular difference between this and the given ChiAngle in desired units
    double ChiAngle::CalculateAngleDifference( const ChiAngle &CHI_ANGLE, const math::Angle::Unit &ANGLE_UNIT) const
    {
      // true if the chi angles are not for the same chi or either of the chi angles are not defined
      if
      (
        CHI_ANGLE.GetChi() != GetChi() || !util::IsDefined( CHI_ANGLE.GetAngle( ANGLE_UNIT))
        || !util::IsDefined( GetAngle( ANGLE_UNIT))
      )
      {
        return util::GetUndefinedDouble();
      }

      double absolute_difference( math::Absolute( CHI_ANGLE.GetAngle( ANGLE_UNIT) - GetAngle( ANGLE_UNIT)));

      // need to take into account that if one angle is -175 and the other is
      // 175 the difference is 10 degrees, not 350
      // so if the difference is >= than 180 then subtract it from 360
      // in order to get the smaller portion of the circle between the angles
      if( absolute_difference > GetMaxAngle( ANGLE_UNIT))
      {
        absolute_difference = GetCircle( ANGLE_UNIT) - absolute_difference;
      }

      return absolute_difference;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ChiAngle::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Chi, ISTREAM);
      io::Serialize::Read( m_Angle, ISTREAM);
      io::Serialize::Read( m_Unit, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ChiAngle::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Chi, OSTREAM, INDENT);
      io::Serialize::Write( m_Angle, OSTREAM, INDENT);
      io::Serialize::Write( m_Unit, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief reads stream formatted in way easier for user to create
    //!        format is space separated
    //!        <chi enum> <angle value> <angle unit>
    //!        An example is
    //!        e_Two 90 degree
    //!        Another example is
    //!        e_Three 1.1 radian
    //! @return istream the ChiAngle was read from
    std::istream &ChiAngle::ReadSimple( std::istream &ISTREAM)
    {
      std::string chi;
      std::getline( ISTREAM, chi);
      BCL_MessageDbg( "chi line is |" + chi + "|");
      if( !chi.empty())
      {
        std::stringstream read( chi);
        io::Serialize::Read( m_Chi, read);
        io::Serialize::Read( m_Angle, read);
        io::Serialize::Read( m_Unit, read);
      }

      return ISTREAM;
    }

    //! @brief gives description of this in format as read by ReadSimple
    //! @return std::string gives description of this in ReadSimple format
    std::string ChiAngle::WriteSimple( const math::Angle::Unit &ANGLE_UNIT) const
    {
      return m_Chi.GetString() + " " +
      util::Format().W( 8).FFP( 3)( GetAngle( ANGLE_UNIT)) + " " + math::Angle::UnitEnum( ANGLE_UNIT).GetString();
    }

  } // namespace biol
} // namespace bcl
