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
#include "pdb/bcl_pdb_line_criterium.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_line.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LineCriterium::s_Instance
    (
      GetObjectInstances().AddInstance( new LineCriterium())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LineCriterium::LineCriterium() :
      m_MeetAllCriteria( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LineCriterium
    LineCriterium *LineCriterium::Clone() const
    {
      return new LineCriterium( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LineCriterium::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief just add a line type
    //! @param LINE_TYPE the line needs to have a specific line type
    void LineCriterium::AddCriterium( const LineType &LINE_TYPE)
    {
      m_Criteria[ LINE_TYPE];
    }

    //! @brief set meet all critiera
    //! is set to true, all critiera need be met to be true, if false, only one has to be true
    //! @param MEET_ALL_CRITIERA
    void LineCriterium::SetMeetAllCriteria( const bool MEET_ALL_CRITERIA)
    {
      m_MeetAllCriteria = MEET_ALL_CRITERIA;
    }

    //! @brief reset the criteria
    void LineCriterium::Reset()
    {
      m_Criteria.Reset();
    }

    //! @brief filter all lines that match the criterium from a given set of lines
    //! @param CRITERIUM the criterium that each copied line has to fullfill
    //! @return the list of lines that match the given criterium
    util::ShPtrList< Line>
    LineCriterium::Filter( const util::ShPtrList< Line> &LINES, const util::FunctionInterface< Line, bool> &CRITERIUM)
    {
      util::ShPtrList< Line> result;

      // iterate through all lines
      for( util::ShPtrList< Line>::const_iterator itr( LINES.Begin()), itr_end( LINES.End()); itr != itr_end; ++itr)
      {
        if( CRITERIUM( **itr))
        {
          result.PushBack( *itr);
        }
      }

      return result;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that checks of line meets all criteria
    //! @param LINE the line to check
    //! @return true, if the criteria are met (for intersect, if all are met, for union if one is met)
    bool LineCriterium::operator()( const Line &LINE) const
    {
      // check if there are criteria for that line
      const storage::Map< LineType, storage::Vector< storage::Pair< EntryType, std::string> > >::const_iterator
        line_type_itr( m_Criteria.Find( LINE.GetType()));

      if( line_type_itr == m_Criteria.End())
      {
        return false;
      }

      if( line_type_itr->second.IsEmpty())
      {
        return true;
      }

      if( m_MeetAllCriteria)
      {
        // iterate through criteria
        for
        (
          storage::Vector< storage::Pair< EntryType, std::string> >::const_iterator
            crit_itr( line_type_itr->second.Begin()), crit_itr_end( line_type_itr->second.End());
          crit_itr != crit_itr_end;
          ++crit_itr
        )
        {
          if( !LINE.MatchesCriteria( *crit_itr))
          {
            return false;
          }
        }

        return true;
      }

      // meet one criteria
      // iterate through criteria
      for
      (
        storage::Vector< storage::Pair< EntryType, std::string> >::const_iterator
          crit_itr( line_type_itr->second.Begin()), crit_itr_end( line_type_itr->second.End());
        crit_itr != crit_itr_end;
        ++crit_itr
      )
      {
        if( LINE.MatchesCriteria( *crit_itr))
        {
          return true;
        }
      }

      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LineCriterium::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Criteria       , ISTREAM);
      io::Serialize::Read( m_MeetAllCriteria, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LineCriterium::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Criteria       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MeetAllCriteria, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
