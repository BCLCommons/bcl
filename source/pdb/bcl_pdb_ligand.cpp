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
#include "pdb/bcl_pdb_ligand.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Ligand::s_Instance
    (
      GetObjectInstances().AddInstance( new Ligand())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Ligand::Ligand()
    {
    }

    //! @brief constructor from lines
    Ligand::Ligand( const util::ShPtrList< Line> &LINES) :
      m_Lines( LINES)
    {
      if( LINES.IsEmpty())
      {
        return;
      }

      const Line &first_line( *m_Lines.FirstElement());
      BCL_Assert( first_line.GetType() == GetLineTypes().HETATM, "ligands can only be constructed from HETATM lines");

      m_ResidueName = first_line.GetString(               GetEntryTypes().HETATMResidueName);
      m_ChainID     = first_line.GetChar(                 GetEntryTypes().HETATMChainID);
      m_PDBID       = first_line.GetNumericalValue< int>( GetEntryTypes().HETATMResidueSequenceID);
      m_ICode       = first_line.GetChar(                 GetEntryTypes().HETATMInsertionCode);
    }

    //! @brief Clone function
    //! @return pointer to new Ligand
    Ligand *Ligand::Clone() const
    {
      return new Ligand( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Ligand::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return Lines
    util::ShPtrList< Line> const &Ligand::GetLines() const
    {
      return m_Lines;
    }

    //! change Lines
    util::ShPtrList< Line> &Ligand::ChangeLines()
    {
      return m_Lines;
    }

    //! @brief access to the fullname
    //! @return the full ligand name
    const std::string &Ligand::GetFullname() const
    {
      return m_Fullname;
    }

    //! @brief access to the formula
    //! @return the ligand formula
    const std::string &Ligand::GetFormula() const
    {
      return m_Formula;
    }

    //! @brief set full name
    //! @param FULL_NAME as reported for that residue name in the HETNAME line
    void Ligand::SetFullname( const std::string &FULL_NAME)
    {
      m_Fullname = FULL_NAME;
    }

    //! @brief set the formula
    //! @brief FORMULA as reported in the FORMUL lines for the residue name
    void Ligand::SetFormula( const std::string &FORMULA)
    {
      m_Formula = FORMULA;
    }

    //! @brief add connections for an atom serial
    //! @param CONNECTIONS map of atom serials with sets of all atoms that are connected; all connections are present
    //!        twice and inverted, as they come from the pdb
    void Ligand::AddConnections( const storage::Map< size_t, storage::Set< size_t> > &CONNECTIONS)
    {
      // get the atom serial associated with that residues
      const storage::Set< size_t> serials( AtomSerials());

      // iterate over the connections
      for
      (
        storage::Map< size_t, storage::Set< size_t> >::const_iterator
          center_itr( CONNECTIONS.Begin()), center_itr_end( CONNECTIONS.End());
        center_itr != center_itr_end;
        ++center_itr
      )
      {
        if( !serials.Contains( center_itr->first))
        {
          continue;
        }

        // collect all internal connections
        {
          storage::Set< size_t> connected;
          std::set_intersection
          (
            center_itr->second.Begin(), center_itr->second.End(),
            serials.Begin(), serials.End(),
            std::inserter( connected.InternalData(), connected.End())
          );

          if( !connected.IsEmpty())
          {
            m_InternalConnections[ center_itr->first] = connected;
          }
        }

        // collect external connections
        {
          storage::Set< size_t> connected;
          std::set_difference
          (
            center_itr->second.Begin(), center_itr->second.End(),
            serials.Begin(), serials.End(),
            std::inserter( connected.InternalData(), connected.End())
          );

          if( !connected.IsEmpty())
          {
            m_ExternalConnections[ center_itr->first] = connected;
          }
        }
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom serials for that residue
    //! @return set of atom serials
    storage::Set< size_t> Ligand::AtomSerials() const
    {
      storage::Set< size_t> serials;

      // iterate over all atom lines
      for( util::ShPtrList< Line>::const_iterator itr( m_Lines.Begin()), itr_end( m_Lines.End()); itr != itr_end; ++itr)
      {
        serials.Insert( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().HETATMSerial));
      }

      // end
      return serials;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Ligand::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ResidueName        , ISTREAM);
      io::Serialize::Read( m_ChainID            , ISTREAM);
      io::Serialize::Read( m_PDBID              , ISTREAM);
      io::Serialize::Read( m_ICode              , ISTREAM);
      io::Serialize::Read( m_Lines              , ISTREAM);
      io::Serialize::Read( m_Fullname           , ISTREAM);
      io::Serialize::Read( m_Formula            , ISTREAM);
      io::Serialize::Read( m_InternalConnections, ISTREAM);
      io::Serialize::Read( m_ExternalConnections, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Ligand::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ResidueName        , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ChainID            , OSTREAM) << '\t';
      io::Serialize::Write( m_PDBID              , OSTREAM) << '\t';
      io::Serialize::Write( m_ICode              , OSTREAM) << '\n';
      io::Serialize::Write( m_Lines              , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Fullname           , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Formula            , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_InternalConnections, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExternalConnections, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace pdb
  
} // namespace bcl
