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
#include "pdb/bcl_pdb_handler.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "pdb/bcl_pdb_line_criterium.h"
#include "pdb/bcl_pdb_site.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    //! Flag to choose which helix classes to process
    //! default is helix class 1 (right handed alpha helix) if there are no classes given in the list
    util::ShPtr< command::FlagInterface> &Handler::GetFlagHelixClasses()
    {
      static util::ShPtr< command::FlagInterface> s_flag_helix_class
      (
        new command::FlagDynamic
        (
          "helix_classes",
          "list of helix classes to be read in, by default, only 1-right handed helix will be considered",
          command::Parameter
          (
            "helix_class_code",
            "helix class form the pdb format",
            command::ParameterCheckRanged< size_t>( 1, 10),
            "1"
          ),
          0, 10
        )
      );

      return s_flag_helix_class;
    }

    //! Flag to trigger the merging of overlapping SSE definitions - although they should not occur, there are still
    //! many pdb files that do not adhere to the format standard
    util::ShPtr< command::FlagInterface> &Handler::GetFlagMergeOverlappingSSEs()
    {
      static util::ShPtr< command::FlagInterface> s_flag_merge
      (
        new command::FlagStatic
        (
          "merge_overlapping_sses",
          "merge overlapping sses given in a pdb file. Although overlapping SSEs are not pdb format conform, they occur"
          " quite often. If two SSE definitions are of the same type and they overlap, they will be merged into one."
        )
      );

      return s_flag_merge;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Handler::Handler( const bool IGNORE_CLASH) :
      m_IgnoreClash( IGNORE_CLASH),
      m_HelixClasses( GetHelixClassesFromCommandLine())
    {
    }

    //! construct PDBRreader from the std::string, containing the pdb file. It reads all information and performs a
    //! SequenceCheck. NOT appropiate when there is no SEQRES information in the pdb file!!!
    //! In that case, use default constructor and call function Read().
    Handler::Handler( std::istream &ISTREAM, const bool IGNORE_CLASH) :
      m_IgnoreClash( IGNORE_CLASH),
      m_HelixClasses( GetHelixClassesFromCommandLine())
    {
      Read( ISTREAM);
    }

    //! copy constructor
    Handler::Handler( const Handler &PDBREADER, const bool IGNORE_CLASH) :
      m_Head( PDBREADER.m_Head),
      m_Models( PDBREADER.m_Models),
      m_Tail( PDBREADER.m_Tail),
      m_SEQRES( PDBREADER.m_SEQRES),
      m_SSEStructure( PDBREADER.m_SSEStructure),
      m_IgnoreClash( IGNORE_CLASH),
      m_HelixClasses( GetHelixClassesFromCommandLine()),
      m_NonProteinChainsIDs( PDBREADER.m_NonProteinChainsIDs)
    {
    }

    //! copy constructor
    Handler *Handler::Clone() const
    {
      return new Handler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief linetypes within group
    //! @return set of line types
    const storage::Set< LineType> &Handler::GetTypesOfLines() const
    {
      static const storage::Set< LineType> s_line_types( GetLineTypes().Begin(), GetLineTypes().End());
      return s_line_types;
    }

    //! @brief return the Lines for a given line type
    //! @param LINE_TYPE the line type of interest
    //! @return reference to the lines - can be empty if no lines are available
    util::ShPtrList< Line> Handler::GetLines( const LineType &LINE_TYPE) const
    {
      util::ShPtrList< Line> lines;

      // from head
      lines.Append( m_Head.GetLines( LINE_TYPE));

      // from all models
      for( storage::Vector< Model>::const_iterator itr( m_Models.Begin()), itr_end( m_Models.End()); itr != itr_end; ++itr)
      {
        lines.Append( itr->GetLines( LINE_TYPE));
      }

      // from tail
      lines.Append( m_Tail.GetLines( LINE_TYPE));

      // end
      return lines;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Handler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @return lines that are considered by criterium
    util::ShPtrList< Line> Handler::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const
    {
      util::ShPtrList< Line> lines;

      // from head
      lines.Append( m_Head.CollectLines( CRITERIUM));

      // from all models
      for( storage::Vector< Model>::const_iterator itr( m_Models.Begin()), itr_end( m_Models.End()); itr != itr_end; ++itr)
      {
        lines.Append( itr->CollectLines( CRITERIUM));
      }

      // from tail
      lines.Append( m_Tail.CollectLines( CRITERIUM));

      // end
      return lines;
    }

    // returns the entire sequence as one letter code of a chain CHAINID
    std::string Handler::GetSequence( const char CHAINID) const
    {
      // initialize empty sequence
      std::string sequence( "");

      // get chain for given CHAINID
      const storage::Map< char, storage::List< ResidueSimple> >::const_iterator chain_itr( m_SEQRES.Find( CHAINID));
      if( chain_itr == m_SEQRES.End())
      {
        return sequence;
      }

      // add one-letter-code of each residue in chain to sequence
      for
      (
        storage::List< ResidueSimple>::const_iterator res_itr( chain_itr->second.Begin()), res_itr_end( chain_itr->second.End());
        res_itr != res_itr_end;
        ++res_itr
      )
      {
        sequence += biol::GetAATypes().AATypeFromThreeLetterCode( res_itr->GetResidueName())->GetOneLetterCode();
      }

      // end
      return sequence;
    }

    //! @brief get the helix classes to be considered as given in the command line
    //! @return set of helix classes
    storage::Set< biol::SSType> Handler::GetHelixClassesFromCommandLine()
    {
      // get all classes specified on the commandline
      const storage::Vector< size_t> class_list( GetFlagHelixClasses()->GetNumericalList< size_t>());

      // the set of classes
      storage::Set< biol::SSType> classes;

      for
      (
        storage::Vector< size_t>::const_iterator itr( class_list.Begin()), itr_end( class_list.End());
        itr != itr_end;
        ++itr
      )
      {
        classes.Insert( biol::GetSSTypes().SSTypeFromPDBHelixClass( *itr));
      }

      // if empty, default is right handed helix
      if( classes.IsEmpty())
      {
        classes.Insert( biol::GetSSTypes().HELIX);
      }

      // end
      return classes;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all chain ids
    std::string Handler::ChainIDs() const
    {
      storage::Set< char> keys( m_SEQRES.GetKeys());
      keys.InsertElements( m_NonProteinChainsIDs);
      return std::string( keys.Begin(), keys.End());
    }

    //! @brief get chain ids for protein chains
    std::string Handler::GetProteinChains() const
    {
      storage::Set< char> all_chainids( m_SEQRES.GetKeys());

      std::string protein_chainids;
      std::set_difference
      (
        all_chainids.Begin(), all_chainids.End(),
        m_NonProteinChainsIDs.Begin(), m_NonProteinChainsIDs.End(),
        std::inserter( protein_chainids, protein_chainids.begin())
      );

      // end
      return protein_chainids;
    }

    //! @brief pushback a new line into that group
    //! @param LINE ShPtr to the line
    //! @return true, if it fits into that group (line type is eligible)
    bool Handler::PushBack( const util::ShPtr< Line> &LINE)
    {
      if( m_Head.PushBack( LINE))
      {
        return true;
      }
      if( m_Tail.PushBack( LINE))
      {
        return true;
      }

      if( m_Models.IsEmpty())
      {
        m_Models.PushBack( Model());
      }

      if( m_Models.GetSize() == 1 && m_Models.FirstElement().PushBack( LINE))
      {
        return true;
      }

      return false;
    }

    //! @brief reset the line group
    void Handler::Reset()
    {
      m_Head.Reset();
      m_Models.Reset();
      m_Tail.Reset();
      m_SEQRES.Reset();
      m_SSEStructure.Reset();
      m_NonProteinChainsIDs.Reset();
    }

    //! @brief Appends a list of lines
    //! @param PDBLINES list of lines
    //! @return true, if all lines could be inserted
    bool Handler::AppendLines( const util::ShPtrList< Line> &PDBLINES)
    {
      // iterate through argument lines
      for( util::ShPtrList< Line>::const_iterator itr( PDBLINES.Begin()), itr_end( PDBLINES.End()); itr != itr_end; ++itr)
      {
        if( !PushBack( *itr))
        {
          return false;
        }
      }

      // end
      return true;
    }

    //! @brief extracts the path and tag from the provided full path
    storage::VectorND< 2, std::string> Handler::ExtractPathAndPDBTag( const std::string &FULL_PATH)
    {
      const size_t pos_slash( FULL_PATH.rfind( PATH_SEPARATOR));
      const size_t pos_point( FULL_PATH.rfind( "."));

      std::string path, pdb_tag;

      // if a dot exists
      if( pos_point != std::string::npos)
      {
        // if no slash found then just filename without any path is given
        if( pos_slash == std::string::npos)
        {
          path = ".";
          pdb_tag = FULL_PATH.substr( 0, pos_point);
        }
        // if a slash is found that is before the point then take whatever in between as the pdb tag
        else if( pos_point > pos_slash)
        {
          path = FULL_PATH.substr( 0, pos_slash);
          pdb_tag = FULL_PATH.substr( pos_slash + 1, pos_point - pos_slash - 1);
        }
        else
        {
          BCL_Exit( "Unable to deduct the path and pdb tag from the provided full path \"" + FULL_PATH + "\"", -1);
        }
      }

      return storage::VectorND< 2, std::string>( path, pdb_tag);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &Handler::WriteLines( std::ostream &OSTREAM) const
    {
      return Write( OSTREAM, 0);
    }

    //! reads Handler from std::istream
    std::istream &Handler::Read( std::istream &ISTREAM)
    {
      bool unknown_linetype( false);

      // previously stored pdb lines are cleaned
      Reset();

      util::SiPtr< LineGroupInterface> current_group( &m_Head);

      while( !ISTREAM.eof())
      {
        std::string buffer;
        std::getline( ISTREAM, buffer);
        if( buffer.empty())
        {
          continue;
        }

        // create line and get line type
        const util::ShPtr< Line> current_line( new Line( buffer));
        const LineType &line_type( current_line->GetType());

        if( line_type == GetLineTypes().e_Undefined && !unknown_linetype)
        {
          BCL_MessageVrb( "File contains at least one unknown linetype! " + buffer);
          unknown_linetype = true;
          continue;
        }

        // end of pdb file
        if( line_type == GetLineTypes().END)
        {
          break;
        }

        // store the line
        if( !current_group->PushBack( current_line))
        {
          // if the current group is tail already, than the current line is found at the wrong place
          if( current_group == util::SiPtr< LineGroupInterface>( &m_Tail))
          {
            BCL_MessageCrt( "pdb file does not appear to have proper order of records!");
            BCL_MessageVrb( "ignoring line: " + current_line->GetString());
            continue;
          }

          // is it the start of a new model
          if( line_type == GetLineTypes().MODEL)
          {
            BCL_MessageStd
            (
              "reading model nr: " + current_line->GetString( GetEntryTypes().MODELSerial)
            );
            // create a new model
            m_Models.PushBack( Model());
            current_group = &m_Models.LastElement();
            continue;
          }

          // end of model - set to head
          if( line_type == GetLineTypes().ENDMDL)
          {
            current_group = &m_Head;
            continue;
          }

          // start of a model without explicit indication
          if( Model().GetTypesOfLines().Contains( line_type))
          {
            BCL_MessageDbg( "start reading model");
            // create a new model
            m_Models.PushBack( Model());
            current_group = &m_Models.LastElement();
            current_group->PushBack( current_line);
            continue;
          }

          // is it a line that should go to the tail
          if( m_Tail.GetTypesOfLines().Contains( line_type))
          {
            current_group = &m_Tail;
            m_Tail.PushBack( current_line);
            continue;
          }

          // try to insert the line back into the head group
          if( m_Head.PushBack( current_line))
          {
            BCL_MessageStd
            (
              "Header line following atom or other model-specific lines merged with header lines: "
              + current_line->GetString()
            );
            continue;
          }

          BCL_MessageVrb
          (
            "order of records in pdb is not standard conform! ignoring current line:\n" + current_line->GetString()
          );
          current_group = &m_Head;
        } // change the group
      }

      // if no model was read, insert an empty model
      if( m_Models.IsEmpty())
      {
        m_Models.PushBack( Model());
      }

      m_SEQRES = m_Head.GetSEQRESProteinChains();
      m_SSEStructure = m_Head.SSEDefinitions( m_HelixClasses, GetFlagMergeOverlappingSSEs()->GetFlag());

      // no seqres given in pdb header
      if( m_SEQRES.IsEmpty())
      {
        BCL_MessageStd
        (
          "no SEQRES information given, try to retrieve sequence from ATOM section, which requires that every residue "
          "is given in there!"
        )
        m_SEQRES = m_Models.FirstElement().GetChains();
      }

      // initialize the models
      {
        size_t model_nr( 1);
        for( storage::Vector< Model>::iterator itr( m_Models.Begin()), itr_end( m_Models.End()); itr != itr_end; ++itr, ++model_nr)
        {
          const storage::Map< char, storage::List< ResidueSimple> > missing_residues( m_Head.GetMissingResidues( model_nr));
          // iterate through chains
          for
          (
            storage::Map< char, storage::List< ResidueSimple> >::const_iterator
              seqres_itr( m_SEQRES.Begin()), seqres_itr_end( m_SEQRES.End());
            seqres_itr != seqres_itr_end;
            ++seqres_itr
          )
          {

            const storage::Map< char, storage::List< ResidueSimple> >::const_iterator mis_itr( missing_residues.Find( seqres_itr->first));
            itr->InitializeStructuredChain
            (
              seqres_itr->first,
              seqres_itr->second,
              mis_itr == missing_residues.End() ? storage::List< ResidueSimple>() : mis_itr->second
            );
          }
        }
      }

      //return
      return ISTREAM;
    }

    //! Writes all existing lines into std::ostrem STREAM
    std::ostream &Handler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write header
      m_Head.WriteLines( OSTREAM);

      // write all models
      const bool write_model_separators( m_Models.GetSize() != 1);
      size_t model_count( 1);

      // iterate through models
      for( storage::Vector< Model>::const_iterator itr( m_Models.Begin()), itr_end( m_Models.End()); itr != itr_end; ++itr)
      {
        if( write_model_separators)
        {
          Line model_line( GetLineTypes().MODEL);
          model_line.Put( GetEntryTypes().MODELSerial, model_count);
          OSTREAM << model_line.GetString() << '\n';
        }
        itr->WriteLines( OSTREAM);
        if( write_model_separators)
        {
          Line end_line( GetLineTypes().ENDMDL);
          OSTREAM << end_line.GetString() << '\n';
        }
      }

      // tail
      m_Tail.UpdateMasterRecord( m_Head);
      if( !m_Models.IsEmpty())
      {
        m_Tail.UpdateMasterRecord( m_Models.FirstElement());
      }
      m_Tail.WriteLines( OSTREAM);

      // return
      return OSTREAM;
    }

    //! @brief return a map of hetero residues
    //! @brief return all ligands
    //! @return List of ligands
    util::ShPtrList< Ligand> Handler::GetLigands() const
    {
      util::ShPtrList< Ligand> ligands;

      // connections
      const storage::Map< size_t, storage::Set< size_t> > connections( m_Tail.GetConnections());

      // fullnames
      const storage::Map< std::string, std::string> fullnames( m_Head.GetHetFullname());

      // formula
      const storage::Map< std::string, std::string> formulas( m_Head.GetHetFormula());

      // search all HET entries
      const util::ShPtrList< Line> het_lines( GetLines( GetLineTypes().HET));

      // iterate through all HET lines
      for( util::ShPtrList< Line>::const_iterator itr( het_lines.Begin()), itr_end( het_lines.End()); itr != itr_end; ++itr)
      {
        // new residues for that HET entry
        // select HETATM lines corresponding to the HET entries
        LineCriterium criteria;
        criteria.SetMeetAllCriteria( true);
        criteria.AddCriterium( GetEntryTypes().HETATMResidueSequenceID, ( *itr)->GetString( GetEntryTypes().HETSequenceID));
        criteria.AddCriterium( GetEntryTypes().HETATMInsertionCode    , ( *itr)->GetString( GetEntryTypes().HETInsertionCode));

        // locate the lines
        const size_t expected_nr_lines( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().HETNumberAtoms));
        const util::ShPtrList< Line> hetatm_lines( m_Models.FirstElement().CollectLines( criteria, ( *itr)->GetChar( GetEntryTypes().HETChainID)));

        if( expected_nr_lines != hetatm_lines.GetSize())
        {
          BCL_MessageCrt
          (
            "different number of HETATM lines found than expected: " + ( *itr)->GetString() +
            util::Format()( hetatm_lines.GetSize())
          );
        }

        // new ligand
        util::ShPtr< Ligand> sp_ligand( new Ligand( hetatm_lines));

        // connections
        sp_ligand->AddConnections( connections);

        // name
        {
          const storage::Map< std::string, std::string>::const_iterator name_itr( fullnames.Find( sp_ligand->GetResidueName()));
          if( name_itr == fullnames.End())
          {
            BCL_MessageVrb( "no fullname for ligand available: " + sp_ligand->GetResidueName());
            sp_ligand->SetFullname( "NO_HETNAM");
          }
          else
          {
            sp_ligand->SetFullname( name_itr->second);
          }
        }
        // formula
        {
          const storage::Map< std::string, std::string>::const_iterator formul_itr( formulas.Find( sp_ligand->GetResidueName()));
          if( formul_itr == formulas.End())
          {
            BCL_MessageVrb( "no formula for ligand available: " + sp_ligand->GetResidueName());
          }
          else
          {
            sp_ligand->SetFormula( formul_itr->second);
          }
        }

        ligands.PushBack( sp_ligand);
      }

      // end
      return ligands;
    }

    //! @brief return all sites
    //! @return list of all sites
    util::ShPtrList< Site> Handler::GetSites() const
    {
      static const std::string s_site_identifier_string(    "SITE_IDENTIFIER:");
      static const std::string s_site_evidence_code_string( "EVIDENCE_CODE:");
      static const std::string s_site_description_string(   "SITE_DESCRIPTION:");

      util::ShPtrList< Site> sites;

      // ligands
      const util::ShPtrList< Ligand> ligands( GetLigands());

      // locate all remark 800 lines
      LineCriterium criterium_remark;
      criterium_remark.AddCriterium( GetEntryTypes().REMARK_Number, 800);
      const util::ShPtrList< Line> remark_lines( m_Head.CollectLines( criterium_remark, GetLineTypes().REMARK));

      // iterate through all lines
      for
      (
        util::ShPtrList< Line>::const_iterator itr( remark_lines.Begin()), itr_end( remark_lines.End());
        itr != itr_end;
        ++itr
      )
      {
        // first line needs to be a site identifier string
        if( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteIdentifierString) != s_site_identifier_string)
        {
          continue;
        }

        // site identifier / name
        const std::string site_identifier( util::TrimString( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteIdentifier)));

        // next line
        ++itr;
        if( itr == itr_end) break;
        if( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteEvidenceCodeString) != s_site_evidence_code_string)
        {
          BCL_MessageCrt( "missing REMARK 800 entries for site identifier: " + site_identifier);
          continue;
        }

        // evidence code
        const Site::EvidenceCodeEnum evidence_code
        (
          util::TrimString( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteEvidenceCode))
        );

        // next line
        ++itr;
        if( itr == itr_end) break;
        if( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteDescriptionString) != s_site_description_string)
        {
          BCL_MessageCrt( "missing REMARK 800 entries for site identifier: " + site_identifier);
          continue;
        }

        std::string site_description( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteDescription));

        // site description over multiple lines
        ++itr;
        while( itr != itr_end)
        {
          const std::string current_description( ( *itr)->GetString( GetEntryTypes().REMARK_String));
          if( util::TrimString( current_description).empty())
          {
            break;
          }

          // sometimes there is no empty REMARK 800 line between two SITE definitions, so they have to be put back
          if( ( *itr)->GetString( GetEntryTypes().REMARK_800_SiteIdentifierString) == s_site_identifier_string)
          {
            --itr;
            break;
          }
          site_description += ' ';
          site_description += current_description;
          ++itr;
        }

        // insert new site
        util::ShPtr< Site> sp_site( new Site( site_identifier, evidence_code, site_description));
        sp_site->FindLigand( ligands);

        // find all SITE residues
        LineCriterium criterium_site;
        criterium_site.AddCriterium( GetEntryTypes().SITEName, site_identifier);

        // locate site
        const util::ShPtrList< Line> site_lines( m_Head.CollectLines( criterium_site, GetLineTypes().SITE));
        if( site_lines.IsEmpty())
        {
          BCL_MessageCrt( "required SITE lines are not found: " + site_identifier);

          // last line of REMARK 800
          if( itr == itr_end) break;

          continue;
        }
        // number of residues
        const size_t number_residues( site_lines.FirstElement()->GetNumericalValue< size_t>( GetEntryTypes().SITENumberResidues));
        const size_t max_number_residues_per_site( 4);

        // only soft assert, since the SiteNumberResidues can only have two digits
        if( ( number_residues - 1) / max_number_residues_per_site + 1 != site_lines.GetSize())
        {
          BCL_MessageCrt
          (
            "number residues in SITE " + site_identifier + " does not match the number of lines: " +
            util::Format()( number_residues) + " do not fit in nr lines: " + util::Format()( site_lines.GetSize())
          );
        }

        size_t res_nr( 0);
        size_t previous_line_nr( 0);
        bool record_error( false);
        // iterate through all site lines
        for
        (
          util::ShPtrList< Line>::const_iterator sitel_itr( site_lines.Begin()), sitel_itr_end( site_lines.End());
          sitel_itr != sitel_itr_end;
          ++sitel_itr
        )
        {
          const size_t current_line_nr( ( *sitel_itr)->GetNumericalValue< size_t>( GetEntryTypes().SITESequenceNumber));
          if( previous_line_nr + 1 != current_line_nr)
          {
            record_error = true;
            break;
          }

          // update previours line nr
          previous_line_nr = current_line_nr;

          for
          (
            EntryTypes::const_iterator ent_res_itr( GetEntryTypes().SITEResidueName1.GetIterator()), ent_res_itr_last( GetEntryTypes().SITEResidueName4.GetIterator()),
              ent_chain_itr( GetEntryTypes().SITEChainID1.GetIterator()),
              ent_seq_itr( GetEntryTypes().SITEResidueSequenceID1.GetIterator()),
              ent_icode_itr( GetEntryTypes().SITEResidueInsertionCode1.GetIterator());
            ent_res_itr <= ent_res_itr_last && res_nr < number_residues;
            ent_res_itr += 4, ent_chain_itr += 4, ent_seq_itr += 4, ent_icode_itr += 4, ++res_nr
          )
          {
            // add residue if it is part of the chain to the chain residues
            const ResidueSimple site_residue
            (
              ( *sitel_itr)->GetString( *ent_res_itr), ( *sitel_itr)->GetChar( *ent_chain_itr),
              ( *sitel_itr)->GetNumericalValue< int>( *ent_seq_itr), ( *sitel_itr)->GetChar( *ent_icode_itr)
            );

            if( IsInChain( site_residue))
            {
              sp_site->AddChainResidue( site_residue);
            }
            else
            {
              sp_site->AddHetatmResidue( site_residue);
            }
          }
        }

        // record number not continuous, i.e. if the site name is not unique
        if( record_error)
        {
          BCL_MessageCrt
          (
            "problems parsing SITE lines for " + site_identifier + " since the sequence numbering is off"
          )
          // last line of REMARK 800
          if( itr == itr_end) break;
          continue;
        }

        // insert
        sites.PushBack( sp_site);

        // last line of REMARK 800
        if( itr == itr_end) break;
      }

      return sites;
    }

    //! @brief check if given residue is in one of the chains (if not might be in the HETATM section)
    //! @param RESIDUE the residue to check for
    //! @return true, if that residue in found in one of the chains, false oterhwise
    bool Handler::IsInChain( const ResidueInterface &RESIDUE) const
    {
      // first model
      const Model &first_model( m_Models.FirstElement());

      const storage::Map< char, util::ShPtrList< Line> >::const_iterator chain_itr( first_model.GetHETATMLines().Find( RESIDUE.GetChainID()));

      // if it is not within HETATM lines
      if( chain_itr == first_model.GetHETATMLines().End())
      {
        return true;
      }

      // collect all lines in hetatm section for that residue - if it is empty, residue must be part of chain
      return LineCriterium::Filter( chain_itr->second, RESIDUE.GetCriterium( GetLineTypes().HETATM)).IsEmpty();
    }

  } // namespace pdb
} // namespace bcl
