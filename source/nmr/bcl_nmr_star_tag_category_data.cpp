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
#include "nmr/bcl_nmr_star_tag_category_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

    //! @brief SaveFrame as string
    //! @param FRAME the save frame
    //! @return the string for FRAME
    const std::string &StarTagCategoryData::GetSaveFrameName( const SaveFrame &FRAME)
    {
      static const std::string s_names[] =
      {
        "DistanceConstraint",
        "RDCConstraint",
        GetStaticClassName< SaveFrame>()
      };
      return s_names[ FRAME];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from StarSaveFrame
    //! @param SAVE_FRAME save frame to assign to this tag category
    //! @param DESCRIPTION string used to recognize this category in an NMR-STAR file
    StarTagCategoryData::StarTagCategoryData( const SaveFrame &SAVE_FRAME, const std::string &DESCRIPTION) :
      m_SaveFrame( SAVE_FRAME),
      m_Description( DESCRIPTION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StarTagCategoryData
    StarTagCategoryData *StarTagCategoryData::Clone() const
    {
      return new StarTagCategoryData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarTagCategoryData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StarTagCategoryData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SaveFrame, ISTREAM);
      io::Serialize::Read( m_Description, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StarTagCategoryData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SaveFrame, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace nmr
} // namespace bcl
