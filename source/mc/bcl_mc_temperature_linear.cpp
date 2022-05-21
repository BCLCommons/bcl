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
#include "mc/bcl_mc_temperature_linear.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> TemperatureLinear::s_Instance
    (
      util::Enumerated< TemperatureInterface>::AddInstance( new TemperatureLinear())
    );

  ////////////////////////////////////
  //  construction and destruction  //
  ////////////////////////////////////

    //! @brief default constructor
    TemperatureLinear::TemperatureLinear() :
      m_LastIterationNumber(),
      m_Delta(),
      m_StartTemperature()
    {
    }

    //! @brief constructor from a counter, starting and ending temperature and total number of iterations
    //! @param ITERATION_COUNTER iteration counter
    //! @param START_TEMPERATURE starting temperature
    //! @param END_TEMPERATURE ending temperature
    //! @param NUMBER_OF_ITERATIONS number of iterations
    TemperatureLinear::TemperatureLinear
    (
      const double START_TEMPERATURE,
      const double END_TEMPERATURE,
      const size_t NUMBER_OF_ITERATIONS
    ) :
      m_LastIterationNumber( 0),
      m_Delta( ( END_TEMPERATURE - START_TEMPERATURE) / NUMBER_OF_ITERATIONS),
      m_StartTemperature( START_TEMPERATURE)
    {
    }

    //! @brief virtual copy constructor
    TemperatureLinear *TemperatureLinear::Clone() const
    {
      return new TemperatureLinear( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemperatureLinear::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns last calculated temperature without updating
    //! @return last calculated temperature without updating
    double TemperatureLinear::GetLastCalculatedTemperature() const
    {
      return m_StartTemperature + m_LastIterationNumber * m_Delta;
    }

    //! @brief return current temperature
    //! @param TRACKER the current tracker
    //! @return current temperature
    double TemperatureLinear::GetTemperature( const opti::TrackerBase &TRACKER) const
    {
      // update the last iteration counter
      m_LastIterationNumber = TRACKER.GetIteration();

      // return temperature
      return GetLastCalculatedTemperature();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TemperatureLinear::GetAlias() const
    {
      static const std::string s_name( "LinearTemparatureControl");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TemperatureLinear::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Temperature control that allows linear adjustment of the temperature.");
      serializer.AddInitializer
      (
        "step size of temperature change",
        "delta by which the temperature is changed",
        io::Serialization::GetAgent( &m_Delta)
      );
      serializer.AddInitializer
      (
        "starting temperature",
        "temperature at the beginning of the simulation",
        io::Serialization::GetAgent( &m_StartTemperature)
      );
      serializer.AddDataMember
      (
        "last iteration number",
        io::Serialization::GetAgent( &m_LastIterationNumber)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset this temperature
    void TemperatureLinear::Reset()
    {
      // reset members
      m_LastIterationNumber = 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TemperatureLinear::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LastIterationNumber, ISTREAM);
      io::Serialize::Read( m_Delta, ISTREAM);
      io::Serialize::Read( m_StartTemperature, ISTREAM);

      // end
      return ISTREAM;
    };

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &TemperatureLinear::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LastIterationNumber, OSTREAM) << '\n';
      io::Serialize::Write( m_Delta, OSTREAM) << '\n';
      io::Serialize::Write( m_StartTemperature, OSTREAM);

      // end
      return OSTREAM;
    };

  } // namespace mc
} // namespace bcl
