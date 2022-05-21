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
#include "nmr/bcl_nmr_star_tag_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from tag category and data type
    //! @param CATEGORY category for this tag
    //! @param DESCRIPTION description string for this tag
    //! @param DATA_TYPE data type for this tag
    StarTagData::StarTagData
    (
      const StarTagCategory &CATEGORY,
      const std::string &DESCRIPTION,
      const util::CPPDataTypes::Types DATA_TYPE
    ) :
      m_TagCategory( CATEGORY),
      m_Description( DESCRIPTION),
      m_DataType( DATA_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StarTagData
    StarTagData *StarTagData::Clone() const
    {
      return new StarTagData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarTagData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StarTagData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_TagCategory, ISTREAM);
      io::Serialize::Read( m_Description, ISTREAM);
      io::Serialize::Read( m_DataType, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StarTagData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_TagCategory, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DataType, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace nmr
} // namespace bcl
