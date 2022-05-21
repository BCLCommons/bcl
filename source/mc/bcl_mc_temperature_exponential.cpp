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
#include "mc/bcl_mc_temperature_exponential.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math.h"
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
    const util::SiPtr< const util::ObjectInterface> TemperatureExponential::s_Instance
    (
      util::Enumerated< TemperatureInterface>::AddInstance( new TemperatureExponential())
    );

  ////////////////////////////////////
  //  construction and destruction  //
  ////////////////////////////////////

    //! @brief default constructor
    TemperatureExponential::TemperatureExponential() :
      m_LastIterationNumber( 0),
      m_Scale( 0),
      m_StartTemperature( 0)
    {
    }

    //! @brief constructor from a starting and ending temperature and total number of iterations
    //! @param START_TEMPERATURE starting temperature
    //! @param END_TEMPERATURE ending temperature
    //! @param NUMBER_OF_ITERATIONS number of iterations
    TemperatureExponential::TemperatureExponential
    (
      const double START_TEMPERATURE,
      const double END_TEMPERATURE,
      const size_t NUMBER_OF_ITERATIONS
    ) :
      m_LastIterationNumber( 0),
      m_Scale( 1.0),
      m_StartTemperature( START_TEMPERATURE)
    {
      // if the starting temperature is not 0 and the start is larger than the end
      if( START_TEMPERATURE != 0.0 && START_TEMPERATURE > END_TEMPERATURE)
      {
        // calculate the scale
        m_Scale = pow( END_TEMPERATURE / START_TEMPERATURE, 1.0 / double( NUMBER_OF_ITERATIONS));
      }
    }

    //! @brief virtual copy constructor
    TemperatureExponential *TemperatureExponential::Clone() const
    {
      return new TemperatureExponential( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemperatureExponential::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TemperatureExponential::GetAlias() const
    {
      static const std::string s_name( "ExponentialTemparatureControl");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TemperatureExponential::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Temperature control that allows exponential change of the simulated temperature."
      );
      serializer.AddInitializer
      (
        "scaling factor for temperature adjustment",
        "scaling factor for the temperature adjustment",
        io::Serialization::GetAgent( &m_Scale)
      );
      serializer.AddInitializer
      (
        "starting temperature",
        "temperature at the start of the simulation",
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
    void TemperatureExponential::Reset()
    {
      // reset members
      m_LastIterationNumber = 0;
    }

    //! @brief returns last calculated temperature without updating
    //! @return last calculated temperature without updating
    double TemperatureExponential::GetLastCalculatedTemperature() const
    {
      return m_StartTemperature * math::Pow( m_Scale, double( m_LastIterationNumber));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TemperatureExponential::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LastIterationNumber, ISTREAM);
      io::Serialize::Read( m_Scale, ISTREAM);
      io::Serialize::Read( m_StartTemperature, ISTREAM);

      // end
      return ISTREAM;
    };

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &TemperatureExponential::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LastIterationNumber, OSTREAM) << '\n';
      io::Serialize::Write( m_Scale, OSTREAM) << '\n';
      io::Serialize::Write( m_StartTemperature, OSTREAM);

      // end
      return OSTREAM;
    };

  } // namespace mc
} // namespace bcl
