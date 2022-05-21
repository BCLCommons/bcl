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
#include "sdf/bcl_sdf_mdl_header.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_entry_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! @brief default constructor
    MdlHeader::MdlHeader() :
      m_NumberAtoms( 0),
      m_NumberBonds( 0),
      m_IsValid( false)
    {
    }

    //! @brief constructor from # atoms, # bonds
    MdlHeader::MdlHeader
    (
      const size_t &NUMBER_ATOMS,
      const size_t &NUMBER_BONDS
    ) :
      m_NumberAtoms( NUMBER_ATOMS),
      m_NumberBonds( NUMBER_BONDS),
      m_IsValid( true)
    {
    }

    //! @brief virtual copy constructor
    MdlHeader *MdlHeader::Clone() const
    {
      return new MdlHeader( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MdlHeader::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the string of the header
    //! @return the string of the header
    std::string MdlHeader::ToMdlLine() const
    {
      std::string mdl_header( GetDefaultLine( e_HeaderLine));
      GetMdlEntryTypes().Header_NumberAtomLines->Set( mdl_header, m_NumberAtoms);
      GetMdlEntryTypes().Header_NumberBondLines->Set( mdl_header, m_NumberBonds);
      return mdl_header;
    }

    //! @brief return the string of the header
    //! @return the string of the header
    void MdlHeader::SetFromMdlLine( const std::string &MDL_HEADER, const size_t NUMBER_DESC_LINES)
    {
      // reinitialize
      m_IsValid = false;
      m_NumberAtoms = m_NumberBonds = 0;

      // first check that the size is correct
      static const size_t s_min_size( GetDefaultLine( e_HeaderLine).size());

      if( MDL_HEADER.size() >= s_min_size)
      {
        // mdl V3000 format is not supported, so check that the header is == the default (V2000)
        m_IsValid =
          GetMdlEntryTypes().Header_CtabVersion->GetTrimmedString( MDL_HEADER)
          == GetMdlEntryTypes().Header_CtabVersion->GetDefault();
        if( m_IsValid)
        {
          m_NumberAtoms = GetMdlEntryTypes().Header_NumberAtomLines->GetUnsignedInt( MDL_HEADER);
          m_NumberBonds = GetMdlEntryTypes().Header_NumberBondLines->GetUnsignedInt( MDL_HEADER);
        }
      }
      else if
      (
        NUMBER_DESC_LINES == size_t( 3) && MDL_HEADER.size() >= 6
        && isdigit( MDL_HEADER[ 2]) && isdigit( MDL_HEADER[ 5])
        && ( isdigit( MDL_HEADER[ 1]) || ( MDL_HEADER[ 1] == ' ' && MDL_HEADER[ 0] == ' '))
        && ( isdigit( MDL_HEADER[ 0]) || MDL_HEADER[ 0] == ' ')
        && ( isdigit( MDL_HEADER[ 4]) || ( MDL_HEADER[ 4] == ' ' && MDL_HEADER[ 3] == ' '))
        && ( isdigit( MDL_HEADER[ 3]) || MDL_HEADER[ 3] == ' ')
      )
      {
        // we've already read 3 description lines, test whether it looks enough like a header line to proceed
        m_IsValid = true;
        m_NumberAtoms = GetMdlEntryTypes().Header_NumberAtomLines->GetUnsignedInt( MDL_HEADER);
        m_NumberBonds = GetMdlEntryTypes().Header_NumberBondLines->GetUnsignedInt( MDL_HEADER);
      }
    }

    //! @brief read MdlHeader object from std::istream
    //! @param ISTREAM istream that contains MdlHeader object
    //! @return istream after MdlHeader object was extracted
    std::istream &MdlHeader::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_NumberAtoms, ISTREAM);
      io::Serialize::Read( m_NumberBonds, ISTREAM);
      io::Serialize::Read( m_IsValid,     ISTREAM);
      // return
      return ISTREAM;
    }

    //! @brief write MdlHeader into std::ostream
    //! @param OSTREAM ostream that gets MdlHeader object
    //! @param INDENT indentation
    //! @return ostream after MdlHeader object was inserted
    std::ostream &MdlHeader::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_NumberAtoms, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_NumberBonds, OSTREAM, 0) << '\t';
      io::Serialize::Write( m_IsValid,     OSTREAM, 0);
      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MdlHeader::s_Instance
    (
      GetObjectInstances().AddInstance( new MdlHeader())
    );

  } // namespace sdf
} // namespace bcl

