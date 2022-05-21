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
#include "restraint/bcl_restraint_sas_distance_density_point.h"

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
    const util::SiPtr< const util::ObjectInterface> SasDistanceDensityPoint::s_Instance
    (
      GetObjectInstances().AddInstance( new SasDistanceDensityPoint())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasDistanceDensityPoint::SasDistanceDensityPoint() :
      m_Rvalue( util::GetUndefinedDouble()),
      m_Density( util::GetUndefinedDouble()),
      m_Error( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from given input data
    //! @param QVALUE the scattering angle from SAXS curve ( x-axis)
    //! @param INTENSITY the wavelength intensity for a given scattering angle ( y-axis)
    //! @param ERROR - the variation in measurement for a given intensity
    SasDistanceDensityPoint::SasDistanceDensityPoint
    (
      const double RVALUE,
      const double DENSITY,
      const double MEASUREMENT_ERROR
    ) :
       m_Rvalue( RVALUE),
       m_Density( DENSITY),
       m_Error( MEASUREMENT_ERROR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasDistanceDensityPoint
    SasDistanceDensityPoint *SasDistanceDensityPoint::Clone() const
    {
      return new SasDistanceDensityPoint( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasDistanceDensityPoint::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasDistanceDensityPoint::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Rvalue, ISTREAM);
      io::Serialize::Read( m_Density, ISTREAM);
      io::Serialize::Read( m_Error, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasDistanceDensityPoint::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Rvalue, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Density, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Error, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
