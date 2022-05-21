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
#include "pdb/bcl_pdb_tail.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_head.h"
#include "pdb/bcl_pdb_line_criterium.h"
#include "pdb/bcl_pdb_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Tail::s_Instance
    (
      GetObjectInstances().AddInstance( new Tail())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Tail::Tail() :
      m_MasterRecord( new Line( GetLineTypes().MASTER)),
      m_End( new Line( GetLineTypes().END))
    {
      InitializeMasterRecord();
    }

    //! @brief Clone function
    //! @return pointer to new Tail
    Tail *Tail::Clone() const
    {
      return new Tail( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Tail::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief linetypes within group
    //! @return set of line types
    const storage::Set< LineType> &Tail::GetTypesOfLines() const
    {
      static const storage::Set< LineType> s_line_types
      (
        GetLineTypes().CONECT.GetIterator(),
        GetLineTypes().End()
      );

      // end
      return s_line_types;
    }

    //! @brief access to lines of given type
    //! @param LINE_TYPE the desire line type
    //! @return lines of given type
    util::ShPtrList< Line> Tail::GetLines( const LineType &LINE_TYPE) const
    {
      util::ShPtrList< Line> lines;

      if( LINE_TYPE == m_MasterRecord->GetType())
      {
        lines.PushBack( m_MasterRecord);
      }
      else if( LINE_TYPE == m_End->GetType())
      {
        lines.PushBack( m_End);
      }
      else if( LINE_TYPE == GetLineTypes().CONECT)
      {
        lines.Append( m_ConectLines);
      }

      // end
      return lines;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @return lines that are considered by criterium
    util::ShPtrList< Line> Tail::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const
    {
      util::ShPtrList< Line> lines;

      // check conect lines
      lines.Append( LineCriterium::Filter( m_ConectLines, CRITERIUM));

      // master
      if( CRITERIUM( *m_MasterRecord))
      {
        lines.PushBack( m_MasterRecord);
      }

      // end
      if( CRITERIUM( *m_End))
      {
        lines.PushBack( m_End);
      }

      // end
      return lines;
    }

    //! @brief pushback a new line into that group
    //! @param ShPtr to the line
    //! @return true, if it fits into that group (line type is eligible)
    bool Tail::PushBack( const util::ShPtr< Line> &LINE)
    {
      if( LINE->GetType() == GetLineTypes().CONECT)
      {
        m_ConectLines.PushBack( LINE);
        return true;
      }
      else if( LINE->GetType() == m_MasterRecord->GetType())
      {
        m_MasterRecord = LINE;
        return true;
      }
      else if( LINE->GetType() == m_End->GetType())
      {
        m_End = LINE;
        return true;
      }

      return false;
    }

    //! @brief reset the line group
    void Tail::Reset()
    {
      m_ConectLines.Reset();
      InitializeMasterRecord();
    }

    //! @brief connections defined through atom serials in the CONNECT lines
    //! @return Map of atom serial center atoms and a set of serials for all connected atoms
    storage::Map< size_t, storage::Set< size_t> > Tail::GetConnections() const
    {
      const util::ShPtrList< Line> &lines( GetLines( GetLineTypes().CONECT));

      storage::Map< size_t, storage::Set< size_t> > connections;

      // iterate over the connection lines
      for( util::ShPtrList< Line>::const_iterator itr( lines.Begin()), itr_end( lines.End()); itr != itr_end; ++itr)
      {
        const size_t center( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().CONECTAtomSerial));
        storage::Set< size_t> &connected( connections[ center]);

        // iterate through all connected serials
        for
        (
          EntryTypes::const_iterator
            entry_itr( GetEntryTypes().CONECTBondedAtomSerial1.GetIterator()),
            entry_itr_end( GetEntryTypes().CONECTBondedAtomSerial4.GetIterator() + 1);
          entry_itr != entry_itr_end;
          ++entry_itr
        )
        {
          const size_t serial( ( *itr)->GetNumericalValue< size_t>( *entry_itr));
          if( !util::IsDefined( serial))
          {
            continue;
          }
          connected.Insert( connected.End(), serial);
        }
      }

      // end
      return connections;
    }

    //! @brief update the MASTER record for header information
    //! @param HEAD a pdb head section will all lines before the coordinate section
    void Tail::UpdateMasterRecord( const Head &HEAD) const
    {
      // Number of REMARK records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumRemark, HEAD.Count( GetLineTypes().REMARK));

      // Number of HET records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumHet, HEAD.Count( GetLineTypes().HET));

      // Number of HELIX records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumHelix, HEAD.Count( GetLineTypes().HELIX));

      // Number of SHEET records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumSheet, HEAD.Count( GetLineTypes().SHEET));

      // Number of SITE records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumSite, HEAD.Count( GetLineTypes().SITE));

      // Number of coordinate transformation records  (ORIGX+SCALE+MTRIX)
      {
        size_t number_trans( 0);
        for
        (
          LineTypes::const_iterator
            itr( GetLineTypes().ORIGX1.GetIterator()), itr_end( GetLineTypes().MTRIX3.GetIterator() + 1);
          itr != itr_end;
          ++itr
        )
        {
          number_trans += HEAD.Count( *itr);
        }
        m_MasterRecord->Put( GetEntryTypes().MASTERNumXform, number_trans);
      }

      // Number of CONECT records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumConect, m_ConectLines.GetSize());

      // Number of SEQRES records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumSeqres, HEAD.Count( GetLineTypes().SEQRES));
    }

    //! @brief update the MASTER record for coordinate information
    //! @param FIRST_MODEL only use the first model
    void Tail::UpdateMasterRecord( const Model &FIRST_MODEL) const
    {
      // Number of atomic coordinate records (ATOM+HETATM)
      m_MasterRecord->Put
      (
        GetEntryTypes().MASTERNumCoord,
        FIRST_MODEL.Count( GetLineTypes().ATOM) + FIRST_MODEL.Count( GetLineTypes().HETATM)
      );

      // Number of TER records
      m_MasterRecord->Put( GetEntryTypes().MASTERNumTer, FIRST_MODEL.Count( GetLineTypes().TER));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &Tail::WriteLines( std::ostream &OSTREAM) const
    {
      // iterate through all lines
      for
      (
        util::ShPtrList< Line>::const_iterator itr( m_ConectLines.Begin()), itr_end( m_ConectLines.End());
        itr != itr_end;
        ++itr
      )
      {
        OSTREAM << ( *itr)->GetString() << '\n';
      }

      // master
      OSTREAM << m_MasterRecord->GetString() << '\n';

      // end
      OSTREAM << m_End->GetString() << '\n';

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Tail::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ConectLines , ISTREAM);
      io::Serialize::Read( m_MasterRecord, ISTREAM);
      io::Serialize::Read( m_End         , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Tail::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ConectLines , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MasterRecord, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_End         , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize the master record with zeros
    void Tail::InitializeMasterRecord()
    {
      // iterate through all entries
      for
      (
        EntryTypes::const_iterator
          itr( GetEntryTypes().MASTERNumRemark.GetIterator()),
          itr_end( GetEntryTypes().MASTERNumSeqres.GetIterator() + 1);
        itr != itr_end;
        ++itr
      )
      {
        m_MasterRecord->Put( *itr, 0);
      }
    }

  } // namespace pdb
} // namespace bcl
