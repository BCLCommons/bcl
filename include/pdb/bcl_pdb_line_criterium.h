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

#ifndef BCL_PDB_LINE_CRITERIUM_H_
#define BCL_PDB_LINE_CRITERIUM_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_pdb_entry_type_data.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LineCriterium
    //! @brief is a unary predicate, that defines criteria for contents of a line
    //! @details multiple criteria can be initialized. A line is passed to the operator. If it matches all or one
    //!          (as desired) the operator returns true, flase otherwise.
    //!
    //! @see @link example_pdb_line_criterium.cpp @endlink
    //! @author woetzen
    //! @date Feb 20, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LineCriterium :
      public util::FunctionInterface< Line, bool>
    {

    private:

    //////////
    // data //
    //////////

      //! meet al criteria (false - only one needs to be met)
      bool m_MeetAllCriteria;

      //! map of line types and entrytypes with desired values
      storage::Map< LineType, storage::Vector< storage::Pair< EntryType, std::string> > > m_Criteria;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LineCriterium();

      //! @brief Clone function
      //! @return pointer to new LineCriterium
      LineCriterium *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief just add a line type
      //! @param LINE_TYPE the line needs to have a specific line type
      void AddCriterium( const LineType &LINE_TYPE);

      //! @brief add criteria
      //! @param ENTRY_TYPE the entry type to check
      //! @param DATA the data for that entry that needs to match
      template< typename t_DataType>
      void AddCriterium( const EntryType &ENTRY_TYPE, const t_DataType &DATA)
      {
        m_Criteria[ ENTRY_TYPE->GetLineType()].PushBack
        (
          storage::Pair< EntryType, std::string>( ENTRY_TYPE, ENTRY_TYPE->GetFormat()( DATA))
        );
      }

      //! @brief set meet all critiera
      //! is set to true, all critiera need be met to be true, if false, only one has to be true
      //! @param MEET_ALL_CRITIERA
      void SetMeetAllCriteria( const bool MEET_ALL_CRITERIA);

      //! @brief reset the criteria
      void Reset();

      //! @brief filter all lines that match the criterium from a given set of lines
      //! @param CRITERIUM the criterium that each copied line has to fullfill
      //! @return the list of lines that match the given criterium
      static util::ShPtrList< Line>
      Filter( const util::ShPtrList< Line> &LINES, const util::FunctionInterface< Line, bool> &CRITERIUM);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that checks of line meets all criteria
      //! @param LINE the line to check
      //! @return true, if the criteria are met (for intersect, if all are met, for union if one is met)
      bool operator()( const Line &LINE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LineCriterium

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_LINE_CRITERIUM_H_ 
