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
#include "mc/bcl_mc_temperature_default.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> TemperatureDefault::s_Instance
    (
      util::Enumerated< TemperatureInterface>::AddInstance( new TemperatureDefault())
    );

  ////////////////////////////////////
  //  construction and destruction  //
  ////////////////////////////////////

    //! @brief default constructor
    TemperatureDefault::TemperatureDefault() :
      m_Temperature( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from a starting and ending temperature
    //! @param TEMPERATURE temperature
    TemperatureDefault::TemperatureDefault( const double TEMPERATURE) :
      m_Temperature( TEMPERATURE)
    {
    }

    //! @brief virtual copy constructor
    TemperatureDefault *TemperatureDefault::Clone() const
    {
      return new TemperatureDefault( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemperatureDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TemperatureDefault::GetAlias() const
    {
      static const std::string s_alias( "DefaultTemperatureControl");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TemperatureDefault::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Keeps the temperature constant at the defined value.");
      serializer.AddInitializer
      (
        "temperature",
        "temperature is kept constant at this value",
        io::Serialization::GetAgent( &m_Temperature)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset temperature
    void TemperatureDefault::Reset()
    {
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TemperatureDefault::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Temperature, ISTREAM);

      // end
      return ISTREAM;
    };

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &TemperatureDefault::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Temperature, OSTREAM);

      // end
      return OSTREAM;
    };

  } // namespace mc
} // namespace bcl
