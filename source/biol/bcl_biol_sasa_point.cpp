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
#include "biol/bcl_biol_sasa_point.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasaPoint::s_Instance
    (
      GetObjectInstances().AddInstance( new SasaPoint())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasaPoint::SasaPoint() :
      m_AtomNumber( util::GetUndefinedSize_t()),
      m_SolventExcludedSurface( util::GetUndefinedDouble()),
      m_SolventAccessibleSurface( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from given input data
    //! @param ATOMNUM - The atom number in the PDB
    //! @param SES     - The Solvent Excluded Surface
    //! @param SAS     - The Solvent Accessible Surface
    SasaPoint::SasaPoint
    (
      const size_t ATOMNUM,
      const double SES,
      const double SAS
    ) :
       m_AtomNumber( ATOMNUM),
       m_SolventExcludedSurface( SES),
       m_SolventAccessibleSurface( SAS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasaPoint
    SasaPoint *SasaPoint::Clone() const
    {
      return new SasaPoint( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasaPoint::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasaPoint::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AtomNumber, ISTREAM);
      io::Serialize::Read( m_SolventExcludedSurface, ISTREAM);
      io::Serialize::Read( m_SolventAccessibleSurface, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasaPoint::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AtomNumber, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SolventExcludedSurface, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SolventAccessibleSurface, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
