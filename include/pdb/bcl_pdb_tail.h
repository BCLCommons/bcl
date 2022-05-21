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

#ifndef BCL_PDB_TAIL_H_
#define BCL_PDB_TAIL_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_pdb_line.h"
#include "bcl_pdb_line_group_interface.h"
#include "bcl_pdb_line_type_data.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Tail
    //! @brief holds all lines in a pdb file after the coordinate section.
    //! @details access to the CONECT entries and the MASTER record
    //!
    //! @see @link example_pdb_tail.cpp @endlink
    //! @author woetzen
    //! @date Feb 21, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Tail :
      public LineGroupInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Data, containing LineType and the content for each pdb line
      util::ShPtrList< Line> m_ConectLines;

      //! master record
      mutable
      util::ShPtr< Line>     m_MasterRecord;

      //! end
      util::ShPtr< Line>     m_End;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Tail();

      //! @brief Clone function
      //! @return pointer to new Tail
      Tail *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief linetypes within group
      //! @return set of line types
      const storage::Set< LineType> &GetTypesOfLines() const;

      //! @brief access to lines of given type
      //! @param LINE_TYPE the desire line type
      //! @return lines of given type
      util::ShPtrList< Line> GetLines( const LineType &LINE_TYPE) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief locate lines of given criterium
      //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
      //! @return lines that are considered by criterium
      util::ShPtrList< Line> CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const;

      //! @brief pushback a new line into that group
      //! @param ShPtr to the line
      //! @return true, if it fits into that group (line type is eligible)
      bool PushBack( const util::ShPtr< Line> &LINE);

      //! @brief reset the line group
      void Reset();

      //! @brief connections defined through atom serials in the CONNECT lines
      //! @return Map of atom serial center atoms and a set of serials for all connected atoms
      storage::Map< size_t, storage::Set< size_t> > GetConnections() const;

      //! @brief update the MASTER record for header information
      //! @param HEAD a pdb head section will all lines before the coordinate section
      void UpdateMasterRecord( const Head &HEAD) const;

      //! @brief update the MASTER record for coordinate information
      //! @param FIRST_MODEL only use the first model
      void UpdateMasterRecord( const Model &FIRST_MODEL) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @return outputstream which was written to
      std::ostream &WriteLines( std::ostream &OSTREAM) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief initialize the master record with zeros
      void InitializeMasterRecord();

    }; // class Tail

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_TAIL_H_ 
