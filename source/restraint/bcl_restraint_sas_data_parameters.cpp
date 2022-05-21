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

// include header of this class
#include "restraint/bcl_restraint_sas_data_parameters.h"
// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasDataParameters::s_Instance
    (
      GetObjectInstances().AddInstance( new SasDataParameters())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasDataParameters::SasDataParameters() :
      m_UseSansImplementation( false),
      m_Qvalue( util::GetUndefinedDouble()),
      m_Sasa( util::GetUndefinedDouble()),
      m_ExcludedVolume( util::GetUndefinedDouble()),
      m_HydrationShell( util::GetUndefinedDouble()),
      m_DeuteriumExchangeRate( util::GetUndefinedDouble())
    {
    }

    SasDataParameters::SasDataParameters( const double &Q_VALUE) :
        m_UseSansImplementation( false),
        m_Qvalue( Q_VALUE),
        m_Sasa( util::GetUndefinedDouble()),
        m_ExcludedVolume( util::GetUndefinedDouble()),
        m_HydrationShell( util::GetUndefinedDouble()),
        m_DeuteriumExchangeRate( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from given input data
    //! @param QVALUE the scattering angle from SAXS curve ( x-axis)
    //! @param INTENSITY the wavelength intensity for a given scattering angle ( y-axis)
    //! @param ERROR - the variation in measurement for a given intensity
    SasDataParameters::SasDataParameters
    (
       const bool &SANS_IMPLEMENTATION,
       const double &Q_VALUE,
       const double &SASA_VALUE,
       const double &EXCLUDED_VOLUME,
       const double &HYDRATION_SHELL,
       const double &DEUTERIUM_EXCHANGE_RATE
    ) :
       m_UseSansImplementation( SANS_IMPLEMENTATION),
       m_Qvalue( Q_VALUE),
       m_Sasa( SASA_VALUE),
       m_ExcludedVolume( EXCLUDED_VOLUME),
       m_HydrationShell( HYDRATION_SHELL),
       m_DeuteriumExchangeRate( DEUTERIUM_EXCHANGE_RATE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SaxsScatteringPoint
    SasDataParameters *SasDataParameters::Clone() const
    {
      return new SasDataParameters( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasDataParameters::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasDataParameters::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_UseSansImplementation, ISTREAM);
      io::Serialize::Read( m_Qvalue, ISTREAM);
      io::Serialize::Read( m_Sasa, ISTREAM);
      io::Serialize::Read( m_ExcludedVolume, ISTREAM);
      io::Serialize::Read( m_HydrationShell, ISTREAM);
      io::Serialize::Read( m_DeuteriumExchangeRate, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasDataParameters::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_UseSansImplementation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Qvalue, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Sasa, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExcludedVolume, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HydrationShell, OSTREAM, INDENT)<< '\n';
      io::Serialize::Write( m_DeuteriumExchangeRate, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
