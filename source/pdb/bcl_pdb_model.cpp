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
#include "pdb/bcl_pdb_model.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "pdb/bcl_pdb_line_criterium.h"
#include "pdb/bcl_pdb_residue_simple.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Model::s_Instance
    (
      GetObjectInstances().AddInstance( new Model())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Model::Model() :
      m_ChainLines()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Model
    Model *Model::Clone() const
    {
      return new Model( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Model::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief linetypes within group
    //! @return set of line types
    const storage::Set< LineType> &Model::GetTypesOfLines() const
    {
      static const storage::Set< LineType> s_line_types
      (
        GetLineTypes().ATOM.GetIterator(), GetLineTypes().HETATM.GetIterator() + 1
      );

      return s_line_types;
    }

    //! @brief get hetatm lines
    //! @return reference to HETATM lines
    const storage::Map< char, util::ShPtrList< Line> > &Model::GetHETATMLines() const
    {
      return m_HetatmLines;
    }

    //! @brief access to lines of given type
    //! @param LINE_TYPE the desire line type
    //! @return lines of given type
    util::ShPtrList< Line> Model::GetLines( const LineType &LINE_TYPE) const
    {
      util::ShPtrList< Line> lines;

      // iterate through all chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Begin()), chain_itr_end( m_ChainLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->GetType() == LINE_TYPE)
          {
            lines.PushBack( *itr);
          }
        }
      }

      // iterate through all hetatms
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_HetatmLines.Begin()), chain_itr_end( m_HetatmLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->GetType() == LINE_TYPE)
          {
            lines.PushBack( *itr);
          }
        }
      }

      // end
      return lines;
    }

    //! @brief count the number of lines for a given line type
    //! @param LINE_TYPE the desire line type
    //! @return the number of lines for the given line type in the model
    size_t Model::Count( const LineType &LINE_TYPE) const
    {
      size_t count( 0);

      // iterate through all chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Begin()), chain_itr_end( m_ChainLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->GetType() == LINE_TYPE)
          {
            ++count;
          }
        }
      }

      // iterate through all hetatm
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_HetatmLines.Begin()), chain_itr_end( m_HetatmLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          if( ( *itr)->GetType() == LINE_TYPE)
          {
            ++count;
          }
        }
      }

      // end
      return count;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @return lines that are considered by criterium
    util::ShPtrList< Line> Model::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const
    {
      // collect lines for each chain
      util::ShPtrList< Line> lines;

      // iterate through chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator itr( m_ChainLines.Begin()), itr_end( m_ChainLines.End()); itr != itr_end; ++itr)
      {
        lines.Append( LineCriterium::Filter( itr->second, CRITERIUM));
      }

      // iterate through hetatms
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator itr( m_HetatmLines.Begin()), itr_end( m_HetatmLines.End()); itr != itr_end; ++itr)
      {
        lines.Append( LineCriterium::Filter( itr->second, CRITERIUM));
      }

      // end
      return lines;
    }

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @param CHAIN_ID only lines for that chain if know apriori
    //! @return lines that are considered by criterium
    util::ShPtrList< Line>
    Model::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM, const char CHAIN_ID) const
    {
      util::ShPtrList< Line> collected_lines;

      // from chain
      {
        const storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Find( CHAIN_ID));

        // no such chain
        if( chain_itr != m_ChainLines.End())
        {
          collected_lines = LineCriterium::Filter( chain_itr->second, CRITERIUM);
        }
      }

      // from hetatm
      {
        const storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_HetatmLines.Find( CHAIN_ID));

        // no such chain
        if( chain_itr != m_HetatmLines.End())
        {
          collected_lines.Append( LineCriterium::Filter( chain_itr->second, CRITERIUM));
        }
      }

      return collected_lines;
    }

    //! @brief pushback a new line into that group
    //! @param ShPtr to the line
    //! @return true, if it fits into that group (line type is eligible)
    bool Model::PushBack( const util::ShPtr< Line> &LINE)
    {
      char chain_id( '\0');

      // atom line
      if( LINE->GetType() == GetLineTypes().ATOM)
      {
        chain_id = LINE->GetChar( GetEntryTypes().ATOMChainID);
      }
      // hetatm line
      else if( LINE->GetType() == GetLineTypes().HETATM)
      {
        chain_id = LINE->GetChar( GetEntryTypes().HETATMChainID);
      }
      // anisou line
      else if( LINE->GetType() == GetLineTypes().ANISOU)
      {
        chain_id = LINE->GetChar( GetEntryTypes().ANISOUChainID);
      }
      // ter line
      else if( LINE->GetType() == GetLineTypes().TER)
      {
        chain_id = LINE->GetChar( GetEntryTypes().TERChainID);
        if( chain_id == ' ')
        {
          if( !m_ChainLines.IsEmpty())
          {
            storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.End());
            --chain_itr;
            if( !chain_itr->second.IsEmpty())
            {
              const Line &last_chain_line( *chain_itr->second.LastElement());
              if( last_chain_line.GetType() == GetLineTypes().ATOM)
              {
                chain_id = last_chain_line.GetChar( GetEntryTypes().ATOMChainID);
              }
              else if( last_chain_line.GetType() == GetLineTypes().ATOM)
              {
                chain_id = last_chain_line.GetChar( GetEntryTypes().HETATMChainID);
              }
              else
              {
                BCL_MessageCrt( "TER line does not make sense here!");
                return false;
              }
            }
            else
            {
              BCL_MessageCrt( "TER line does not make sense here!");
              return false;
            }
          }
          else
          {
            BCL_MessageCrt( "TER line does not make sense here!");
            return false;
          }
          BCL_MessageVrb
          (
            "TER line might not be complete: " + LINE->GetString() + "\nassuming it is terminating the highest chainid "
                "read so far: " + chain_id
          );
        }
      }
      // did not find any of the considered line types
      else
      {
        return false;
      }

      // insert the line for the correct chain id
      // check if the TER for the chain in chain lines was given already
      storage::Map< char, util::ShPtrList< Line> >::iterator chain_line_itr( m_ChainLines.Find( chain_id));
      if( chain_line_itr == m_ChainLines.End())
      {
        m_ChainLines[ chain_id].PushBack( LINE);
      }
      else if( chain_line_itr->second.IsEmpty() || chain_line_itr->second.LastElement()->GetType() != GetLineTypes().TER)
      {
        chain_line_itr->second.PushBack( LINE);
      } // TER was already given for that chain id; insert into the hetatm lines
      else
      {
        // only insert hetatm lines
        if( LINE->GetType() != GetLineTypes().HETATM)
        {
          return false;
        }
        m_HetatmLines[ chain_id].PushBack( LINE);
      }

      return true;
    }

    //! @brief reset the line group
    void Model::Reset()
    {
      m_ChainLines.Reset();
      m_HetatmLines.Reset();
      m_StructuredChains.Reset();
    }

    //! @brief create the chains as they were given by the atoms
    //! @return map with chainid as key and list of residues as data
    storage::Map< char, storage::List< ResidueSimple> > Model::GetChains() const
    {
      storage::Map< char, storage::List< ResidueSimple> > chains;

      // iterate over chains
      // iterate through all chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Begin()), chain_itr_end( m_ChainLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
        )
        {
          // current line
          const Line &current_line( **itr);

          // normally, that should be an atom line
          if( current_line.GetType() == GetLineTypes().ATOM)
          {
            // current residue from first atom line
            const ResidueSimple current_residue
            (
              current_line.GetString( GetEntryTypes().ATOMResidueName),
              current_line.GetChar( GetEntryTypes().ATOMChainID),
              current_line.GetNumericalValue< int>( GetEntryTypes().ATOMResidueSequenceID),
              current_line.GetChar( GetEntryTypes().ATOMInsertionCode)
            );
            chains[ current_residue.GetChainID()].PushBack( current_residue);

            // skip all other atom lines belonging to that residue
            while( ++itr != itr_end)
            {
              const Line &next_line( **itr);
              if
              (
                   next_line.GetType() != GetLineTypes().ATOM
                || current_residue.GetResidueName() != next_line.GetString( GetEntryTypes().ATOMResidueName)
                || current_residue.GetChainID()     != next_line.GetChar( GetEntryTypes().ATOMChainID)
                || current_residue.GetPDBID()       != next_line.GetNumericalValue< int>( GetEntryTypes().ATOMResidueSequenceID)
                || current_residue.GetICode()       != next_line.GetChar( GetEntryTypes().ATOMInsertionCode)
              )
              {
                break;
              }
            } // advance iterator to next residue
          } // if ATOM line
          // sometime residues are not standard
          else if( current_line.GetType() == GetLineTypes().HETATM)
          {
            // current residue from first HETATM line
            const ResidueSimple current_residue
            (
              current_line.GetString( GetEntryTypes().HETATMResidueName),
              current_line.GetChar( GetEntryTypes().HETATMChainID),
              current_line.GetNumericalValue< int>( GetEntryTypes().HETATMResidueSequenceID),
              current_line.GetChar( GetEntryTypes().HETATMInsertionCode)
            );
            chains[ current_residue.GetChainID()].PushBack( current_residue);

            // skip all other hetatm lines
            while( ++itr != itr_end)
            {
              const Line &next_line( **itr);
              if
              (
                   next_line.GetType() != GetLineTypes().HETATM
                || current_residue.GetResidueName() != next_line.GetString( GetEntryTypes().HETATMResidueName)
                || current_residue.GetChainID()     != next_line.GetChar( GetEntryTypes().HETATMChainID)
                || current_residue.GetPDBID()       != next_line.GetNumericalValue< int>( GetEntryTypes().HETATMResidueSequenceID)
                || current_residue.GetICode()       != next_line.GetChar( GetEntryTypes().HETATMInsertionCode)
              )
              {
                break;
              }
            } // advance iterator to next residue
          } // if HETATM line
          else // other line type
          {
            ++itr;
          }
        } // iterate over pdb lines
      } // iterate through chains

      // end
      return chains;
    }

    namespace
    {
      //! @brief get the longest contiguous stretch of AAs from a given position
      //! @param ITR, ITR_END iterator positions
      size_t GetNextLongestContiguousStretch
      (
        storage::List< storage::Map< char, Residue> >::const_iterator ITR,
        const storage::List< storage::Map< char, Residue> >::const_iterator &ITR_END
      )
      {
        if( ITR == ITR_END)
        {
          return 0;
        }
        int sz( 1);
        const int initial_pdb_id( ITR->Begin()->second.GetPDBID());
        for( ++ITR; ITR != ITR_END && ITR->Begin()->second.GetPDBID() == initial_pdb_id + sz; ++ITR, ++sz)
        {
        }
        return sz;
      }
    }

    //! @brief initialize the structure
    //! @param CHAIN_ID
    //! @param SEQRES sequences for the chain
    //! @param MISSING_RESIDUES residues that are missing in the ATOM section (REMARK 465)
    void Model::InitializeStructuredChain
    (
      const char CHAIN_ID,
      const storage::List< ResidueSimple> &SEQRES,
      const storage::List< ResidueSimple> &MISSING_RESIDUES
    )
    {
      // initialize the chain from the seqres
      storage::List< Residue> &structure( m_StructuredChains[ CHAIN_ID]);
      for( storage::List< ResidueSimple>::const_iterator res_itr( SEQRES.Begin()), res_itr_end( SEQRES.End()); res_itr != res_itr_end; ++res_itr)
      {
        structure.PushBack( Residue( *res_itr));
      }

      // locate all atom lines for the current chain
      // all atom lines are collected in residues, which have a Map of alternate locations for one and the same residue
      storage::List< storage::Map< char, Residue> >
        atom_line_residues( AtomLineResiduesFromChainID( CHAIN_ID));

      // merge the atom line residues with the missing residues
      const storage::List< storage::Map< char, Residue> >
        merged_atom_line_residues( MergeAtomLinesAndMissingResidueLine( atom_line_residues, MISSING_RESIDUES));

      // iterator on list of residues from atom lines
      storage::List< storage::Map< char, Residue> >::const_iterator
        atom_itr( merged_atom_line_residues.Begin()), atom_itr_end( merged_atom_line_residues.End());

      // iterator for residues in that chain from seqres
      storage::List< Residue>::iterator res_itr( structure.Begin()), res_itr_end( structure.End());

      // align atom line and residue gaps
      while( atom_itr != atom_itr_end && res_itr != res_itr_end)
      {
        // make copies of "res_itr" and "atom_itr"
        // they point to residues and atom lines, respectively, with matching residue types
        storage::List< storage::Map< char, Residue> >::const_iterator atom_itr2( atom_itr);
        storage::List< Residue>::iterator res_itr2( res_itr);

        // create size_t "number_aligned_residues" this will hold the number of aligned residues
        size_t number_aligned_residues( 0);

        // determine the maximum # of residues that could possibly be aligned by AlignResiduesAndAtomLines
        size_t dist_res_end( std::distance( res_itr2, res_itr_end));
        size_t dist_atom_end( GetNextLongestContiguousStretch( atom_itr2, atom_itr_end));
        size_t max_residues_aligned
        (
          std::min
          (
            size_t( s_NumberRequiredAlignedResidues),
            std::min( dist_res_end, dist_atom_end)
          )
        );
        storage::List< Residue>::iterator last_aligned_residue( res_itr);

        // match up first residue from Structure with first by res name corresponding atom line
        // and check if at least max_residues_aligned can be aligned if name matches
        while
        (
             !CompareAtomLineResWithResByName( *atom_itr, *res_itr, false)
          || ( number_aligned_residues = AlignResiduesAndAtomLines( atom_itr2, res_itr2, atom_itr_end, res_itr_end, false, max_residues_aligned)) < max_residues_aligned
        )
        {
          ++res_itr;
          if( res_itr == res_itr_end)
          {
            // outer loop will also break, since res_itr points to end
            break;
          }
          else
          {
            BCL_MessageDbg
            (
              "could align " + util::Format()( number_aligned_residues) + " starting from residue\n" +
              util::Format()( *res_itr)
            );
          }
          atom_itr2 = atom_itr;
          res_itr2 = res_itr;
          dist_res_end = std::distance( res_itr2, res_itr_end);
          dist_atom_end = GetNextLongestContiguousStretch( atom_itr2, atom_itr_end);
          max_residues_aligned =
            std::min
            (
              size_t( s_NumberRequiredAlignedResidues),
              std::min( dist_res_end, dist_atom_end)
            );
        }

        // align the residues that could be aligned so far by assigning atom line to residues
        if( !AlignResiduesAndAtomLines( atom_itr, res_itr, atom_itr_end, res_itr_end, true, util::GetUndefined< size_t>()))
        {
          BCL_MessageStd
          (
            "could not find any residue in m_Structure that matches with given atom line residue:\n" +
            util::Format()( *atom_itr)
          );
          ++atom_itr;
          res_itr = last_aligned_residue;
        }
      }

      // add the pdb id to residues in the beginning and end, if they do not have some yet
      res_itr = structure.Begin();
      size_t first_res_with_defined_pdb_id( 0);
      // find the first res_itr with a defined pdb id
      while( res_itr != res_itr_end && !util::IsDefined( res_itr->GetPDBID()))
      {
        ++res_itr;
        ++first_res_with_defined_pdb_id;
      }

      // all residues had undefined pdbid
      if( res_itr == res_itr_end)
      {
        BCL_MessageCrt
        (
          "there are no alignable atom lines for that chain in the pdb: " + util::Format()( CHAIN_ID)
        );
        return;
      }

      // residues in front had undefined pdb ids
      if( res_itr != structure.Begin())
      {
        // get the current pdbid
        int current_res_id( res_itr->GetPDBID());
        --res_itr;
        --current_res_id;

        //while( res_itr != res_itr_begin)
        for
        (
          size_t current( 0); current != first_res_with_defined_pdb_id; ++current
        )
        {
          // for pdbs there is no pdb residue id of 0
          if( current_res_id == 0)
          {
            --current_res_id;
          }
          res_itr->SetPDBID( current_res_id);

          // this has to be checked since VS does not allow --Begin on forward iterators
          // and this can't be asserted in the for loop sentinel, since we still have to fix the first residue
          if( res_itr == structure.Begin())
          {
            break;
          }

          BCL_MessageDbg( "did set the resid of that residue\n" + util::Format()( *res_itr));
          --res_itr;
          --current_res_id;
        }
      }

      // iterator to rev begin of chain
      storage::List< Residue>::reverse_iterator
        res_itr_rev( structure.ReverseBegin()),
        res_itr_rev_end( structure.ReverseEnd());

      size_t number_skipped_to_first_defined_pdb_id( 0);
      // find the first res_itr_rev with a defined pdb id
      while( res_itr_rev != res_itr_rev_end && !util::IsDefined( res_itr_rev->GetPDBID()))
      {
        ++res_itr_rev;
        ++number_skipped_to_first_defined_pdb_id;
      }

      // all residues had undefined pdbid
      if( res_itr_rev == res_itr_rev_end)
      {
        BCL_MessageCrt
        (
          "there are no alignable atom lines for that chain in the pdb: " + util::Format()( CHAIN_ID)
        );
        return;
      }

      // residues in tail had undefined pdb ids
      if( res_itr_rev != structure.ReverseBegin())
      {
        // get the current pdbid
        int current_res_id( res_itr_rev->GetPDBID());
        --res_itr_rev;
        ++current_res_id;

        //while( res_itr_rev != res_itr_rev_begin)
        for( size_t current( 0); current < number_skipped_to_first_defined_pdb_id; ++current)
        {
          // for pdbs there is no pdb residue id of 0
          if( current_res_id == 0)
          {
            ++current_res_id;
          }

          res_itr_rev->SetPDBID( current_res_id);

          BCL_MessageDbg( "did set the resid of that residue\n" + util::Format()( *res_itr_rev));

          // this has to be checked since VS does not allow --Begin on forward iterators
          // and this can't be asserted in the for loop sentinel, since we still have to fix the first residue
          if( res_itr_rev == structure.ReverseBegin())
          {
            break;
          }

          --res_itr_rev;
          ++current_res_id;
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &Model::WriteLines( std::ostream &OSTREAM) const
    {
      // iterate through all chains
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Begin()), chain_itr_end( m_ChainLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          OSTREAM << ( *itr)->GetString() << '\n';
        }
      }

      // iterate through all hetatm lines
      for( storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_HetatmLines.Begin()), chain_itr_end( m_HetatmLines.End()); chain_itr != chain_itr_end; ++chain_itr)
      {
        // iterate through all lines
        for
        (
          util::ShPtrList< Line>::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          OSTREAM << ( *itr)->GetString() << '\n';
        }
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Model::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ChainLines , ISTREAM);
      io::Serialize::Read( m_HetatmLines, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Model::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChainLines , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HetatmLines, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

    //! @brief list that stores for each pdbid and insertion code a map of alternate location ids their residue
    //! HETATM lines are also considered, but are only stored as hard copies within the residues as ATOM lines
    //! @param CHAINID chain for which this list is generated
    //! @return list that stores for each pdbid and insertion code a map of alternate location ids their residue
    storage::List< storage::Map< char, Residue> >
    Model::AtomLineResiduesFromChainID( const char CHAINID) const
    {
      util::ShPtrList< Line> chain_atom_lines;
      // from chain
      {
        const storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( m_ChainLines.Find( CHAINID));

        // no such chain
        if( chain_itr != m_ChainLines.End())
        {
          // search for all atom and het atom lines for that chain
          LineCriterium atomline_chain_criteria;
          atomline_chain_criteria.SetMeetAllCriteria( false);
          atomline_chain_criteria.AddCriterium( GetLineTypes().ATOM);
          atomline_chain_criteria.AddCriterium( GetLineTypes().HETATM);

          chain_atom_lines = LineCriterium::Filter( chain_itr->second, atomline_chain_criteria);
        }
      }

      // list that stores for each pdbid and insertion code a map of alternate location ids their residue
      storage::List< storage::Map< char, Residue> > atom_line_residues;

      // iterator for atom lines
      util::ShPtrList< Line>::const_iterator
        atom_itr( chain_atom_lines.Begin()), atom_itr_end( chain_atom_lines.End());

      // list that stores for each pdbid and insertion code a map of alternate location ids their residue
      storage::List< storage::Map< char, Residue> >::iterator alr_res_itr( atom_line_residues.End());

      // create string "current_res_id" initialize with sequence id of atom line currently denoted by "ATOM_LINE_ITR"
      int previous_res_id( util::GetUndefined< int>());

      // create string "current_insert_id" initialize with insertion code of atom line denoted by "ATOM_LINE_ITR"
      char previous_insert_icode( util::GetUndefined< char>());

      while( atom_itr != atom_itr_end)
      {
        // shptr to current atom line
        util::ShPtr< Line> current_line( *atom_itr);
        if( current_line->GetType() == GetLineTypes().HETATM)
        {
          current_line = util::ShPtr< Line>( Line::CopyHetatmToAtomLine( *current_line).Clone());
        }

        // create string "current_res_id" initialize with sequence id of atom line currently denoted by "ATOM_LINE_ITR"
        const int current_res_id( current_line->GetNumericalValue< int>( GetEntryTypes().ATOMResidueSequenceID));

        // create string "current_insert_id" initialize with insertion code of atom line denoted by "ATOM_LINE_ITR"
        const char current_insert_icode( current_line->GetChar( GetEntryTypes().ATOMInsertionCode));

        // hit new residue in atom lines
        if( current_res_id != previous_res_id || previous_insert_icode != current_insert_icode)
        {
          // check that first insertion id has all side chain atom lines
          if( alr_res_itr != atom_line_residues.End())
          {
            CompleteFirstAlternateResidue( *alr_res_itr);
          }

          // insert new empty residue map
          atom_line_residues.PushBack( storage::Map< char, Residue>());

          // set iterator to the new residue map
          alr_res_itr = atom_line_residues.Last();
        }

        // alternate location id
        const char current_alternate_location_id( current_line->GetChar( GetEntryTypes().ATOMAlternateLocationID));

        // point to current residue for that alternate_location_id
        storage::Map< char, Residue>::iterator current_residue_itr( alr_res_itr->Find( current_alternate_location_id));

        // if there is no residue with this current_alternate_location_id
        if( current_residue_itr == alr_res_itr->End())
        {
          // make a new one
          alr_res_itr->operator[]( current_alternate_location_id) =
              Residue
              (
                current_line->GetString( GetEntryTypes().ATOMResidueName),
                CHAINID,
                current_res_id,
                current_insert_icode
              );
          current_residue_itr = alr_res_itr->Find( current_alternate_location_id);
        }

        // assert that residue name is identical for that atom line and that residue
        BCL_Assert
        (
          current_residue_itr->second.GetResidueName() == current_line->GetString( GetEntryTypes().ATOMResidueName),
          "two atom line with same alternate location id, name pdbid and insertion code:\n" +
          ( *atom_itr)->GetString() + "\n" + util::Format()( current_residue_itr->second)
        );

        // insert current line to the residue
        current_residue_itr->second.ChangeLines().PushBack( current_line);

        // remember pdb id and i code
        previous_res_id = current_res_id;
        previous_insert_icode = current_insert_icode;

        // go to next line
        ++atom_itr;
      }

      // end
      return atom_line_residues;
    }

    //! @brief complete first residue in a map of alternate residues
    //! @param RESIDUES map of residues where first one is to be completed
    void
    Model::CompleteFirstAlternateResidue( storage::Map< char, Residue> &RESIDUES)
    {
      // nothing to do if there is only one alternate residue
      if( RESIDUES.GetSize() <= 1)
      {
        return;
      }

      // append atoms lines from the second alternative location
      if( RESIDUES.Begin()->first == ' ')
      {
        Residue &correct_residues( RESIDUES.Begin()->second);
        correct_residues.ChangeLines().Append( ( ++RESIDUES.Begin())->second.GetLines());
      }
    }

    //! @brief add missing residues to the atom line residues
    //! @param ATOM_LINE_RESIDUES residues acquired form residues
    //! @param MISSING_RESIDUES residues that are not located in experiment REMARK 465
    //! @return list of maps (alternate location residues) merged with the missing residues
    storage::List< storage::Map< char, Residue> >
    Model::MergeAtomLinesAndMissingResidueLine
    (
      const storage::List< storage::Map< char, Residue> > &ATOM_LINE_RESIDUES,
      const storage::List< ResidueSimple> &MISSING_RESIDUES
    )
    {
      // assuming that the missing residues and the atom line residues have consistent numbering, the missing residues
      // are just added in between
      storage::List< storage::Map< char, Residue> > merged_residues;

      // iterators for the atom line residues
      storage::List< storage::Map< char, Residue> >::const_iterator
        atom_itr( ATOM_LINE_RESIDUES.Begin()), atom_itr_end( ATOM_LINE_RESIDUES.End());

      // iterator on the missing residues
      storage::List< ResidueSimple>::const_iterator
        mis_itr( MISSING_RESIDUES.Begin()), mis_itr_end( MISSING_RESIDUES.End());

      // try to inserted sorted
      while( atom_itr != atom_itr_end && mis_itr != mis_itr_end)
      {
        bool nothing_happened( true);
        // insert atom line residues till the first missing residue
        for( ; atom_itr != atom_itr_end && atom_itr->Begin()->second < *mis_itr; ++atom_itr, nothing_happened &= false)
        {
          merged_residues.PushBack( *atom_itr);
        }

        // end of atom line residues
        if( atom_itr == atom_itr_end)
        {
          break;
        }

        // insert missing residues till the first atom line
        for( ; mis_itr != mis_itr_end && *mis_itr < atom_itr->Begin()->second; ++mis_itr, nothing_happened &= false)
        {
          merged_residues.InsertElement
          (
            storage::Map< char, Residue>::Create
            (
              std::pair< char, Residue>( ' ', Residue( *mis_itr))
            )
          );
        }

        if( nothing_happened)
        {
          BCL_MessageCrt( "unable to align missing residues to ATOM section residues! Ignoring Missing residues!");
          return ATOM_LINE_RESIDUES;
        }
      }

      // fill up to the end with atom line residues
      for( ; atom_itr != atom_itr_end; ++atom_itr)
      {
        merged_residues.InsertElement( *atom_itr);
      }

      // fill up to the end with missing residues
      for( ; mis_itr != mis_itr_end; ++mis_itr)
      {
        merged_residues.InsertElement
        (
          storage::Map< char, Residue>::Create
          (
            std::pair< char, Residue>( ' ', Residue( *mis_itr))
          )
        );
      }

      // check that merged list is as long as sum of both lists
      BCL_Assert
      (
        merged_residues.GetSize() == ATOM_LINE_RESIDUES.GetSize() + MISSING_RESIDUES.GetSize(),
        "was not able to merge missing residues with residues from atom lines"
      )

      // end
      return merged_residues;
    }

    //! @brief compare a map of alternate residues to a given residue if one matches by name
    //! @param ATOM_RES map of alternate residues
    //! @param RES residue to be matched
    //! @param SET_RESIDUE if match was found, set the given RES to the residue
    //! @return true, if match was found
    bool Model::CompareAtomLineResWithResByName
    (
      const storage::Map< char, Residue> &ATOM_RES,
      Residue &RESIDUE,
      const bool SET_RESIDUE
    )
    {
      // check for match to every alternate residue
      for
      (
        storage::Map< char, Residue>::const_iterator
          atom_res_itr( ATOM_RES.Begin()), atom_res_itr_end( ATOM_RES.End());
        atom_res_itr != atom_res_itr_end;
        ++atom_res_itr
      )
      {
        // there was one residue of the ones with the alternate location id, that did match
        if( biol::GetAATypes().HaveSameParent( atom_res_itr->second.GetResidueName(), RESIDUE.GetResidueName()))
        {
          if( SET_RESIDUE)
          {
            RESIDUE = atom_res_itr->second;
          }
          return true;
        }
      }

      // no match found
      return false;
    }

    //! @brief align range of residues and atom lines
    //! @param ATOM_LINE_ITR itr to atomline that matches up by name with the residue
    //! @param RES_ITR itr to residue that matches up by name with the atom line
    //! @param ATOM_LINE_ITR_END total end of atom lines
    //! @param RES_ITR_END total end of residues
    //! @param SET_RESIDUE_ATOM_LINES if true, atom lines are inserted into the residues
    //! @param MAX_NUMBER_RES_TO_ALIGN stop aligning after this number is exceeded
    //! @return number fo residues that have been aligned
    size_t Model::AlignResiduesAndAtomLines
    (
      storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR,
      storage::List< Residue>::iterator    &RES_ITR,
      const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR_END,
      const storage::List< Residue>::iterator    &RES_ITR_END,
      const bool SET_RESIDUE_ATOM_LINES,
      const size_t MAX_NUMBER_RES_TO_ALIGN
    )
    {
      // create size_t "number_residue_matches" to hold the number of sequential residues that can be aligned
      size_t number_residue_matches( 0);

      // require that all following residues with the same type as the current type be aligned
      size_t max_number_res_to_align( MAX_NUMBER_RES_TO_ALIGN);
      if( util::IsDefined( MAX_NUMBER_RES_TO_ALIGN))
      {
        max_number_res_to_align
          = std::min
            (
              std::max
              (
                MAX_NUMBER_RES_TO_ALIGN,
                GetLengthHomopolymerChain( ATOM_LINE_ITR, ATOM_LINE_ITR_END, RES_ITR, RES_ITR_END)
              ),
              size_t( std::min( std::distance( RES_ITR, RES_ITR_END), std::distance( ATOM_LINE_ITR, ATOM_LINE_ITR_END)))
            );
      }

      // see how many sequential residues match up
      while
      (
        // try to align only MAX_NUMBER_RES_TO_ALIGN
          number_residue_matches < max_number_res_to_align
        // check that "ATOM_LINE_ITR" has not reached "ATOM_LINE_ITR_END"
        && ATOM_LINE_ITR != ATOM_LINE_ITR_END
        // check that "RES_ITR" has not reached "RES_ITR_END"
        && RES_ITR != RES_ITR_END
        // check that residues match up and evtl set the res to the matched one (SET_RESIDUE_ATOM_LINES flag)
        && CompareAtomLineResWithResByName( *ATOM_LINE_ITR, *RES_ITR, SET_RESIDUE_ATOM_LINES)
      )
      {
        const Residue &current_res( ATOM_LINE_ITR->Begin()->second);
        // Retrieve current res id for sequence distance check
        const int current_res_id( current_res.GetPDBID());

        // got to next atom line residue
        ++ATOM_LINE_ITR;

        // no residues in atom lines left, nothing to match up anymore
        if( ATOM_LINE_ITR == ATOM_LINE_ITR_END)
        {
          // increase the number of residue matches
          ++number_residue_matches;

          // move "RES_ITR" to next residue
          ++RES_ITR;

          break;
        }

        const Residue &next_res( ATOM_LINE_ITR->Begin()->second);

        // different between current and next atom line res id
        const int res_id_diff( next_res.GetPDBID() - current_res_id);

        // if the res id difference is unexpected, there could be a deletion, not mentioned in the REMARK 465
        // since the atom lines also contain the residues with deletions, those residues would not be empty, if they actually
        // come from the ATOM lines
        if( res_id_diff > 1 && !current_res.GetLines().IsEmpty() && !next_res.GetLines().IsEmpty())
        {
          storage::List< Residue>::iterator peek_res_itr( RES_ITR);
          BCL_MessageDbg
          (
            "chain " + util::Format()( RES_ITR->GetChainID()) + ": "
            "two following atom line residues do have difference in res id larger than 1: " +
            util::Format()( res_id_diff) + ". Will try to skip that many residues! till atom line residue line:\n" +
            util::Format()( *ATOM_LINE_ITR)
          );

          // try to advance iterator by the number of residues that are defined by the given difference
          // peek to the next residue and confirm that it would agree at least by its residue name
          storage::AdvanceIterator( peek_res_itr, RES_ITR_END, res_id_diff);

          // if peek hits end of residues, assume that difference does not mean anything
          if( peek_res_itr == RES_ITR_END)
          {
            BCL_MessageDbg
            (
              "chain " + util::Format()( RES_ITR->GetChainID()) + ": "
              "could not skip the residues, since the resulting residue does not exist"
            );

            // increase the number of residue matches
            ++number_residue_matches;

            // move "RES_ITR" to next residue
            ++RES_ITR;

            // next matching pair
            continue;
          }

          // if the peek residue is not at the end
          storage::List< storage::Map< char, Residue> >::const_iterator atom_itr2( ATOM_LINE_ITR);
          storage::List< Residue>::iterator res_itr2( peek_res_itr);

          const size_t number_remaining_residues
          (
            std::min
            (
              std::distance( atom_itr2, ATOM_LINE_ITR_END),
              std::distance( res_itr2, RES_ITR_END)
            )
          );
          const size_t min_aligned_peek_residues
          (
            std::min
            (
              number_remaining_residues,
              size_t( s_NumberRequiredAlignedResidues)
            )
          );

          const size_t number_aligned_peek_residues
          (
            AlignResiduesAndAtomLines
            (
              atom_itr2, res_itr2, ATOM_LINE_ITR_END, RES_ITR_END, false, min_aligned_peek_residues
            )
          );
          if( number_aligned_peek_residues < min_aligned_peek_residues)
          {
            BCL_MessageDbg
            (
              "chain " + util::Format()( RES_ITR->GetChainID()) + ": "
              "could not skip the residues, since the resulting pair of residue did not match: " +
              util::Format()( *peek_res_itr) + "\n!=\n" + util::Format()( *ATOM_LINE_ITR)
            );

            // increase the number of residue matches
            ++number_residue_matches;

            // move "RES_ITR" to next residue
            ++RES_ITR;

            // next iteration
            continue;
          }

          // if aligned more than 3 peek residues
          // set the pdbid of the residues that do not have an atom line
          if( SET_RESIDUE_ATOM_LINES)
          {
            // go to the start of the "gapped" residues
            ++RES_ITR;
            for
            (
              int number_skipped_residues( 1);
              RES_ITR != peek_res_itr && number_skipped_residues < res_id_diff;
              ++RES_ITR, ++number_skipped_residues
            )
            {
              RES_ITR->SetPDBID( current_res_id + number_skipped_residues);
            }
          }

          BCL_MessageDbg
          (
            "chain " + util::Format()( RES_ITR->GetChainID()) + ": "
            "successful skipping the lines, since the residues did match: " +
            util::Format()( *peek_res_itr) + "\n!=\n" + util::Format()( *ATOM_LINE_ITR) +
            "\nand " + util::Format()( number_aligned_peek_residues) + " peek residues could be aligned"
          );

          // set RES_ITR to peek res_itr position
          RES_ITR = peek_res_itr;

          // increase the number of matched residues by the difference in residue ids
          number_residue_matches += res_id_diff - 1;

          // next iteration
          continue;
        }
        // increase the number of residue matches
        ++number_residue_matches;

        // move "RES_ITR" to next residue
        ++RES_ITR;
      }
      return number_residue_matches;
    }

    //! @brief Get the number of residues before getting to the next type
    //! @param ATOM_LINE_ITR itr to atomline that matches up by name with the residue
    //! @param ATOM_LINE_ITR_END total end of atom lines
    //! @param RES_ITR itr to residue that matches up by name with the atom line
    //! @param RES_ITR_END total end of residues
    //! @return min(number next residues in atom lines that have the same res type,
    //!             number next residues in seqres lines that have the same res type)
    size_t Model::GetLengthHomopolymerChain
    (
      const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR,
      const storage::List< storage::Map< char, Residue> >::const_iterator &ATOM_LINE_ITR_END,
      const storage::List< Residue>::iterator &RES_ITR,
      const storage::List< Residue>::iterator &RES_ITR_END
    )
    {
      if( RES_ITR == RES_ITR_END)
      {
        return 0;
      }
      // determine the # of AAs till the next AA type in the atom lines
      const Residue &first_residue( *RES_ITR);
      const char chainid( first_residue.GetChainID());
      size_t number_seqresidues_till_next_type( 0);
      for( storage::List< Residue>::const_iterator res_itr( RES_ITR); res_itr != RES_ITR_END; ++res_itr)
      {
        if( res_itr->GetResidueName() == first_residue.GetResidueName() && res_itr->GetChainID() == chainid)
        {
          ++number_seqresidues_till_next_type;
        }
        else
        {
          break;
        }
      }

      // determine the # of AAs till the next AA type in the seqres lines
      storage::List< storage::Map< char, Residue> >::const_iterator atom_line_itr( ATOM_LINE_ITR);
      const storage::Map< char, Residue> &first_seqres( *atom_line_itr);
      size_t number_residues_till_next_type( 1);
      for
      (
        ++atom_line_itr;
        atom_line_itr != ATOM_LINE_ITR_END;
        ++atom_line_itr, ++number_residues_till_next_type
      )
      {
        bool does_match( false);
        for
        (
          storage::Map< char, Residue>::const_iterator
            first_res_itr( first_seqres.Begin()), first_res_itr_end( first_seqres.End());
          first_res_itr != first_res_itr_end && !does_match;
          ++first_res_itr
        )
        {
          // check for match to every alternate residue
          for
          (
            storage::Map< char, Residue>::const_iterator
              atom_res_itr( atom_line_itr->Begin()), atom_res_itr_end( atom_line_itr->End());
            atom_res_itr != atom_res_itr_end;
            ++atom_res_itr
          )
          {
            // there was one residue of the ones with the alternate location id, that did match
            if
            (
              biol::GetAATypes().HaveSameParent
              (
                atom_res_itr->second.GetResidueName(),
                first_res_itr->second.GetResidueName()
              )
            )
            {
              does_match = true;
              break;
            }
          }
        }
        if( !does_match)
        {
          break;
        }
      }
      return std::max( number_residues_till_next_type, number_seqresidues_till_next_type);
    }

    //! @brief Inserts missing resid informations if they are in HETATM lines instead of ATOM lines
    //! @param CHAIN_ID the chain if of that chain
    //! @param CHAIN the residues for that chain
    void Model::InsertLinesFromHETATMLines
    (
      const char CHAIN_ID,
      storage::List< Residue> &CHAIN
    ) const
    {
      // instantiate an iterator pointing to the first residue of current chain (called res_itr,
      // because the chain structure has been build with the atom line information
      for
      (
        storage::List< Residue>::iterator res_itr( CHAIN.Begin()), res_itr_end( CHAIN.End());
        res_itr != res_itr_end;
        ++res_itr
      )
      {
        //checks for emptiness of current residue - if it is not empty -> go to next residue
        if( !res_itr->GetLines().IsEmpty())
        {
          continue;
        }

        // construct criteria vector to search for current residue to find HETATM lines
        LineCriterium residue_criteria;
        residue_criteria.SetMeetAllCriteria( true);
        residue_criteria.AddCriterium( GetEntryTypes().HETATMResidueSequenceID, res_itr->GetPDBID());
        residue_criteria.AddCriterium( GetEntryTypes().HETATMInsertionCode, res_itr->GetICode());

        //set new type of each line to be ATOM line
        util::ShPtrList< Line> hetatm_reslines( LineCriterium::Filter( m_ChainLines.Find( CHAIN_ID)->second, residue_criteria));
        for
        (
          util::ShPtrList< Line>::iterator
            hetatom_itr( hetatm_reslines.Begin()), hetatom_itr_end( hetatm_reslines.End());
          hetatom_itr != hetatom_itr_end; ++hetatom_itr
        )
        {
          if
          (
            biol::GetAATypes().HaveSameParent
            (
              res_itr->GetResidueName(),
              ( *hetatom_itr)->GetString( GetEntryTypes().HETATMResidueName)
            )
          )
          {
            util::ShPtr< Line> atom_line( Line::CopyHetatmToAtomLine( **hetatom_itr).Clone());
            res_itr->ChangeLines().PushBack( atom_line);
          }
        }
      }
    }

  } // namespace pdb
} // namespace bcl
