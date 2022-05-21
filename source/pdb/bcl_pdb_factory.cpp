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
#include "pdb/bcl_pdb_factory.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_dssp.h"
#include "biol/bcl_biol_membrane.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_handler.h"
#include "pdb/bcl_pdb_printer_biomatrix.h"
#include "sspred/bcl_sspred_pdb.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  //////////
  // data //
  //////////

    //! static bool for whether all command line defaults have been saved
    bool Factory::s_HaveDefaults( Factory::ResetFlagDefaults());

    //! Flag to switch between numerated AtomIDs and the original PDBAtomIDs
    //! - default will numerate the atomIDs as they are written
    util::ShPtr< command::FlagInterface> &
    Factory::GetFlagWritePDBAtomID()
    {
      static util::ShPtr< command::FlagInterface> s_flag_write_pdb_atom_id
      (
        new command::FlagStatic
        (
          "write_pdb_atom_ids",
          "switch between numerated AtomIDs and the original PDBAtomIDs - default will numerate the atom IDs as they are written"
        )
      );

      return s_flag_write_pdb_atom_id;
    }

    //! Flag to switch between numerated ResIDs( SeqIDs) and the original PDBResIDs and PDBICodes
    //! default will write SeqIDs
    util::ShPtr< command::FlagInterface> &
    Factory::GetFlagWritePDBResID()
    {
      static util::ShPtr< command::FlagInterface> s_flag_write_pdb_res_id
      (
        new command::FlagStatic
        (
          "write_pdb_res_ids",
          "switch between numerated ResidueIDs( SeqIDs) and the original PDBResIDs and PDBICodes - if flag is not set, sequence IDs (starting with 1 for first residue) will be written"
        )
      );

      return s_flag_write_pdb_res_id;
    }

    //! command line flag to be used to change the Defaut AAClass that is used to generate AASequences/Chains/Proteins
    util::ShPtr< command::FlagInterface> &
    Factory::GetFlagAAClass()
    {
      static util::ShPtr< command::FlagInterface> s_flag_aa_class
      (
        new command::FlagStatic
        (
          "aaclass",
          "choose the amino acid class that will be written",
          command::Parameter
          (
            "amino_acid_class",
            "choice of amino acid class that is different by the amino acids used",
            command::ParameterCheckEnumerate< biol::AAClasses>(),
            biol::GetAAClasses().e_AABackBone.GetName()
          )
        )
      );

      return s_flag_aa_class;
    }

    //! command line flag to be used to change the default minimal size of the sse types tp be included in protein models
    util::ShPtr< command::FlagInterface> &
    Factory::GetFlagMinSSESize()
    {
      static util::ShPtr< command::ParameterInterface> s_min_helix_size_parameter
      (
        new command::Parameter
        (
          "min_helix_size",
          "minimal size of helices to be added to protein model",
          command::ParameterCheckRanged< size_t>( 0, 999),
          "0"
        )
      );

      static util::ShPtr< command::ParameterInterface> s_min_strand_size_parameter
      (
        new command::Parameter
        (
          "min_strand_size",
          "minimal size of strand to be added to protein model",
          command::ParameterCheckRanged< size_t>( 0, 999),
          "0"
        )
      );

      static util::ShPtr< command::ParameterInterface> s_min_loop_size_parameter
      (
        new command::Parameter
        (
          "min_loop_size",
          "minimal size of loops to be added to protein model",
          command::ParameterCheckRanged< size_t>( 0, 999),
          "0"
        )
      );

      static util::ShPtr< command::FlagInterface> s_min_sse_size_flag
      (
        new command::FlagStatic
        (
          "min_sse_size",
          "change the default minimal size of each secondary structure type that is to be included when protein models "
          "are constructed from a pdb file"
        )
      );

      // insert parameters if it was constructed for the first time
      if( s_min_sse_size_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> min_sse_size_flag( s_min_sse_size_flag);
        min_sse_size_flag->PushBack( s_min_helix_size_parameter);
        min_sse_size_flag->PushBack( s_min_strand_size_parameter);
        min_sse_size_flag->PushBack( s_min_loop_size_parameter);
      }

      // end
      return s_min_sse_size_flag;
    }

    //! command line flag to be used to allow conversion of unnatural aminoacids to their parent amino acid
    util::ShPtr< command::FlagInterface> &Factory::GetFlagConvertToNaturalAAType()
    {
      static util::ShPtr< command::FlagInterface> s_flag_convert_to_natural_aa_type
      (
        new command::FlagStatic
        (
          "convert_to_natural_aa_type",
          "when creating amino acids from a pdb, unnatural amino acid types are converted to their parent amino acid: "
          "e.g. MSE -> MET"
        )
      );

      return s_flag_convert_to_natural_aa_type;
    }

    //! command line flag to be used, if pdb does not have any sse definition to use the backbone conformation to determine that
    util::ShPtr< command::FlagInterface> &Factory::GetFlagSSEsFromBackBone()
    {
      static util::ShPtr< command::FlagInterface> s_sse_from_back_bone
      (
        new command::FlagStatic
        (
          "sse_from_backbone",
          "only if no pdb SSE definitions are given for certain chain within pdb, use the backbone conformation to calculate them"
        )
      );

      //end
      return s_sse_from_back_bone;
    }

    //! command line flag to be used, if pdb SSE definitions should be reassigned by DSSP
    util::ShPtr< command::FlagInterface> &Factory::GetFlagDSSP()
    {
      static util::ShPtr< command::FlagInterface> s_sse_from_dssp
      (
        new command::FlagStatic
        (
          "sse_from_dssp",
          "reassign pdb SSE definitions using dssp",
          command::Parameter
          (
            "max_hbond_energy",
            "maximal H-Bond energy to be considered an H-Bond",
            command::ParameterCheckRanged< double>( -10, 0),
            util::Format()( biol::DSSP::GetDefaultHBondMaxEnergy())
          )
        )
      );

      //end
      return s_sse_from_dssp;
    }

    //! command line flag to be used to write zero coordinates for undefined amino acids
    util::ShPtr< command::FlagInterface> &Factory::GetFlagWriteZeroCoordinatesForUndefinedAminoAcids()
    {
      static util::ShPtr< command::FlagInterface> s_write_zero_coordinates
      (
        new command::FlagStatic
        (
          "write_zero_coordinates",
          "If amino acids are undefined, write them to the pdb with zero coordinates"
        )
      );

      //end
      return s_write_zero_coordinates;
    }

    // command line flag to be used to write hydrogen atoms if available
    util::ShPtr< command::FlagInterface> &Factory::GetFlagWriteHydrogens()
    {
      static util::ShPtr< command::FlagInterface> s_write_hdydrogens
      (
        new command::FlagStatic
        (
          "write_hydrogens",
          "output hydrogen atoms if available"
        )
      );

      //end
      return s_write_hdydrogens;
    }

    //! command line flag to switch between IUPAC (standard) and pdb atom name output (when flag is used)
    util::ShPtr< command::FlagInterface> &Factory::GetFlagPDBAtomName()
    {
      static util::ShPtr< command::FlagInterface> s_write_pdb_atom_names
      (
        new command::FlagStatic
        (
          "write_pdb_atom_names",
          "instead of IUPAC atom names, pdb atom nameing will be used (e.g. HH21 -> 1HH2 for ARGININE)"
        )
      );

      //end
      return s_write_pdb_atom_names;
    }

    //! command line flag for specifying which biomolecule number to read in for generating multimers
    util::ShPtr< command::FlagInterface> &Factory::GetFlagBiomolecule()
    {
      static util::ShPtr< command::FlagInterface> s_flag_biomolecule
      (
        new command::FlagStatic
        (
          "biomolecule",
          "use the REMARK 350  BIOMT entries to generate the biomolecule",
          command::Parameter
          (
            "BIOMOLECULE",
            "biomolecule number",
            command::ParameterCheckRanged< size_t>( 1, 99),
            "1"
          )
        )
      );

      // end
      return s_flag_biomolecule;
    }

    //! function to return all flags for the pdb factory
    util::ShPtrVector< command::FlagInterface> &Factory::GetAllFlags()
    {
      static util::ShPtrVector< command::FlagInterface> s_all_flags
      (
        util::ShPtrVector< command::FlagInterface>::Create
        (
          GetFlagWritePDBAtomID(),
          GetFlagWritePDBResID(),
          GetFlagAAClass(),
          GetFlagMinSSESize(),
          GetFlagSSEsFromBackBone(),
          GetFlagDSSP(),
          GetFlagConvertToNaturalAAType(),
          GetFlagWriteZeroCoordinatesForUndefinedAminoAcids(),
          GetFlagWriteHydrogens(),
          GetFlagPDBAtomName(),
          GetFlagBiomolecule()
        )
      );
      return s_all_flags;
    }

    //! @brief reset defaults on all flags; to ensure a consistent set of defaults
    //! This function should be called whenever a given example or app may depend on the normal default parameters
    //! @return true on success
    //! @note should be called during static initialization to acquire the actual defaults
    bool Factory::ResetFlagDefaults()
    {
      static util::ShPtrVector< command::FlagInterface> &s_flags( GetAllFlags());

      // record of the original defaults for all flags
      static storage::Vector< storage::Vector< std::string> > s_original_defaults;

      if( s_original_defaults.IsEmpty())
      {
        // fill up the vector with the original defaults
        for
        (
          util::ShPtrVector< command::FlagInterface>::const_iterator
            itr( s_flags.Begin()), itr_end( s_flags.End());
          itr != itr_end;
          ++itr
        )
        {
          s_original_defaults.PushBack( storage::Vector< std::string>());
          for
          (
            util::ShPtrVector< command::ParameterInterface>::const_iterator
              itr_param( ( *itr)->GetParameterList().Begin()), itr_param_end( ( *itr)->GetParameterList().End());
            itr_param != itr_param_end;
            ++itr_param
          )
          {
            // store the default value
            s_original_defaults.LastElement().PushBack( ( *itr_param)->GetDefaultValue());
          }
        }
      }
      else
      {
        // perform the reset on all the flags and reset defaults on parameters
        storage::Vector< storage::Vector< std::string> >::const_iterator itr_defaults( s_original_defaults.Begin());
        for
        (
          util::ShPtrVector< command::FlagInterface>::iterator
            itr( s_flags.Begin()), itr_end( s_flags.End());
          itr != itr_end;
          ++itr, ++itr_defaults
        )
        {
          // reset the flag itself
          ( *itr)->ResetFlag();

          // reset the defaults
          storage::Vector< std::string>::const_iterator itr_default( itr_defaults->Begin());
          for
          (
            util::ShPtrVector< command::ParameterInterface>::iterator
              itr_param( ( *itr)->GetParameterList().Begin()), itr_param_end( ( *itr)->GetParameterList().End());
            itr_param != itr_param_end;
            ++itr_param, ++itr_default
          )
          {
            // test whether there is a default value
            if( ( *itr_param)->GetWasDefaultGiven())
            {
              // yes, so set the default back to its original value
              ( *itr_param)->SetDefaultParameter( *itr_default);
            }
          }
        }
      }
      return true;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Factory::Factory() :
      m_AAClass( GetFlagAAClass()->GetFirstParameter()->GetValue()),
      m_Printers
      (
        1,
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > >( new PrinterBiomatrix())
      )
    {
    }

    //! @brief constructor from AAClass
    //! @param AA_CLASS AAClass of interest
    Factory::Factory( const biol::AAClass &AA_CLASS) :
      m_AAClass( AA_CLASS),
      m_Printers
      (
        1,
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > >( new PrinterBiomatrix())
      )
    {
    }

    //! @brief virtual copy constructor
    Factory *Factory::Clone() const
    {
      return new Factory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Factory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief resets the printers
    void Factory::ResetPrinters()
    {
      m_Printers.Reset();
      m_Printers.PushBack
      (
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > >( new PrinterBiomatrix())
      );
    }

    //! @brief appends printer to list
    //! @param SP_PRINTER printer to be appended
    void Factory::AppendPrinter
    (
      const util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > > &SP_PRINTER
    )
    {
      m_Printers.Append( SP_PRINTER);
    }

  ///////////////////////
  // operations - read //
  ///////////////////////

    //! @brief builds biol::Atom from Line
    //! @param ATOM_TYPE type of atom in line
    //! @param LINE Line which contains Atom information
    //! @return ShPtr to Atom read from Line
    util::ShPtr< biol::Atom>
    Factory::AtomFromLine
    (
      const biol::AtomType &ATOM_TYPE,
      const Line &LINE
    )
    {
      // return util::ShPtr for biol::Atom, constructed from Atom Name, PDB Atom Serial and its Position
      return util::ShPtr< biol::Atom>
      (
        new biol::Atom
        (
          LINE.RetrieveCoordinates(),
          ATOM_TYPE,
          LINE.GetNumericalValue< size_t>( GetEntryTypes().ATOMSerial),
          LINE.GetNumericalValue< double>( GetEntryTypes().ATOMTempFactor)
        )
      );
    }

    //! @brief builds biol::Atom from Line
    //! @param LINE Line which contains Atom information
    //! @return ShPtr to Atom read from Line
    util::ShPtr< biol::Atom>
    Factory::AtomFromLine( const Line &LINE)
    {
      const biol::AtomType type( biol::GetAtomTypes().TypeFromPDBAtomName( util::TrimString( LINE.GetString( GetEntryTypes().ATOMName))));

      // if the AtomName in the Line is undefined
      if( !type.IsDefined())
      {
        BCL_MessageStd
        (
          " This atom name leads to undefined biol::AtomType |" +
            util::TrimString( LINE.GetString( GetEntryTypes().ATOMName)) + "|"
        );
      }

      // return util::ShPtr for biol::Atom, constructed from Atom Name, PDB Atom Serial and its Position
      return AtomFromLine( type, LINE);
    }

    //! @brief builds all Atoms for provided AtomTypes from Lines
    //! @param LINES List of Lines which contains Atom information
    //! @param AA_TYPE aa type for which atoms are desired
    //! @return ShPtrVector of Atoms read from given Lines
    util::ShPtrVector< biol::Atom>
    Factory::AtomsFromLines
    (
      const util::ShPtrList< Line> &LINES,
      const biol::AAType &AA_TYPE,
      const biol::AAType &AA_TYPE_IN_FILE
    ) const
    {
      util::ShPtrVector< biol::Atom> atoms;

      if( LINES.IsEmpty())
      {
        return atoms;
      }

      storage::Set< biol::AtomType> missing_types( AA_TYPE->GetAllowedAtomTypes());

      // find extra atom types that the aa type in the file may have (for unnatural aa types) that will not be
      // represented in the bcl
      storage::Set< biol::AtomType> extra_types;
      if( AA_TYPE != AA_TYPE_IN_FILE)
      {
        std::set_difference
        (
          AA_TYPE_IN_FILE->GetAllowedAtomTypes().Begin(),
          AA_TYPE_IN_FILE->GetAllowedAtomTypes().End(),
          missing_types.Begin(),
          missing_types.End(),
          std::inserter( extra_types.InternalData(), extra_types.Begin())
        );
      }

      // iterate over atoms in residue lines to match this type
      for
      (
        util::ShPtrList< Line>::const_iterator line_itr( LINES.Begin()), line_itr_end( LINES.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        // get the atom type from the pdb atom name
        const biol::AtomType &current_atom_type
        (
          AA_TYPE->GetAtomTypeFromAtomName( ( *line_itr)->GetString( GetEntryTypes().ATOMName))
        );

        // if the types match, store the atom and update the found status
        if( !current_atom_type.IsDefined())
        {
          const biol::AtomType real_atom_type
          (
            biol::GetAtomTypes().TypeFromPDBAtomName( ( *line_itr)->GetString( GetEntryTypes().ATOMName))
          );

          // if the atom type is in the atom in the file but is not in the aa type used in the bcl (e.g. when using
          // convert to natural aa type), quietly skip the aa
          if( extra_types.Erase( real_atom_type) == 0)
          {
            BCL_MessageCrt
            (
              "found line with atom type, that is not compatible with the residue type: " + ( *line_itr)->GetString()
            );
          }

          continue;
        }

        // try to delete the type
        if( !biol::GetAtomTypes().GetTerminalExtraAtomTypes().Has( current_atom_type) && missing_types.Erase( current_atom_type) == 0)
        {
          BCL_MessageCrt
          (
            "found additional atom types for residue with line: " + ( *line_itr)->GetString()
          );
          continue;
        }

        // if type could have been removed, atom was not seen yet, insert
        atoms.PushBack( AtomFromLine( current_atom_type, **line_itr));
      }

      // iterate over missing types and check that all non-hydrogens atoms could be constructed
      for
      (
        storage::Set< biol::AtomType>::const_iterator itr( missing_types.Begin()), itr_end( missing_types.End());
        itr != itr_end;
        ++itr
      )
      {
        // hydrogen can be missed
        if( ( *itr)->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          continue;
        }

        //
        BCL_MessageVrb
        (
          "at least one non hydrogen atom line is missing for residue starting with line: " +
          LINES.FirstElement()->GetString()
        );

        break;
      }

      // return all atoms
      return atoms;
    }

    //! @brief build amino acid from pdb residue
    //! @param RESIDUE Residue that contains amino acid information
    //! @param SEQ_ID sequence id of amino acid to be read
    //! @return ShPtr to amino acid that was read from given Residue
    util::ShPtr< biol::AABase>
    Factory::AminoAcidFromResidue( const Residue &RESIDUE, const int SEQ_ID) const
    {
      // aa type for current aa
      const biol::AAType aa_type_in_file( biol::GetAATypes().AATypeFromThreeLetterCode( RESIDUE.GetResidueName()));
      const biol::AAType current_aa_type
      (
        aa_type_in_file.IsDefined() && GetFlagConvertToNaturalAAType()->GetFlag()
        ? aa_type_in_file->GetParentType()
        : aa_type_in_file
      );

      // aa data for new amino acid
      util::ShPtr< biol::AAData> sp_aa_data
      (
        new biol::AAData
        (
          current_aa_type,
          SEQ_ID,
          RESIDUE.GetPDBID(),
          RESIDUE.GetICode(),
          RESIDUE.GetChainID()
        )
      );

      // construct an amino acid without any coordinates of specified AAType
      util::ShPtr< biol::AABase> amino_acid( ( *m_AAClass)->Empty( sp_aa_data));

      // set atoms for that amino acid
      //if no atom information is available return amino acid without atoms
      if( RESIDUE.GetLines().IsEmpty())
      {
        // construct an amino acid without any coordinates
        return amino_acid;
      };

      // get the atom types for this AA
      const storage::Set< biol::AtomType> atom_types( current_aa_type->GetAllowedAtomTypes());

      // get all atoms used for this
      util::ShPtrVector< biol::Atom> atoms( AtomsFromLines( RESIDUE.GetLines(), current_aa_type, aa_type_in_file));

      // for glycine, construct a HA2, which corresponds to CB in carbon side chain, for any other AA construct CB if missing
      if( !biol::Atom::FindAtom( atoms, current_aa_type->GetFirstSidechainAtomType())->GetType().IsDefined())
      {
        // first side chain atom coordinate
        const linal::Vector3D first_sc_atom_coord
        (
          FirstSidechainAtomCoordinateFromBackboneAtoms
          (
            atoms,
            biol::GetAtomTypes().CA->GetBondLength( current_aa_type->GetFirstSidechainAtomType())
          )
        );

        // if information to construct pseudo side chain atom is not sufficient
        if( !first_sc_atom_coord.IsDefined())
        {
          // give warning
          BCL_MessageCrt
          (
            "Insufficient data to construct first side chain atom for residue with leading Line:\n" +
            RESIDUE.GetLines().FirstElement()->GetString()
           );
        }
        else
        {
          // insert first side chain atom
          atoms.PushBack( util::ShPtr< biol::Atom>( new biol::Atom( first_sc_atom_coord, current_aa_type->GetFirstSidechainAtomType())));
        }
      }

      // first side chain atom is is in place - set the atoms and return the amino acid
      amino_acid->SetAtoms( atoms);

      // end
      return amino_acid;
    }

    //! @brief builds AASequence from CHAINID and util::ShPtrVector of Residues
    //! @param PDB_ID PdbID of the sequence to be read
    //! @param CHAIN_ID ChainID of the sequence to be read
    //! @param RESIDUES List of Residues that contain amino acid information for the sequence
    //! @return ShPtr to AASequence that was read from Residues
    util::ShPtr< biol::AASequence>
    Factory::AASequenceFromResidues
    (
      const std::string &PDB_ID,
      const char CHAIN_ID,
      const storage::List< Residue> &RESIDUES
    ) const
    {
      //instantiate new aasequences
      const std::string pdb_id( ">" + PDB_ID + ":" + CHAIN_ID + "|PDBID|CHAIN|SEQUENCE");
      util::ShPtr< biol::AASequence> aasequence( new biol::AASequence( m_AAClass, 0, CHAIN_ID, pdb_id));

      //give residues increasing number for later sequence-distance calculations
      int seqid( 1);

      //loop over all pdb residues and construt util::ShPtr< AMINOACID_TYPE> of amino acids
      for
      (
        storage::List< Residue>::const_iterator residue_itr( RESIDUES.Begin()),
          residue_itr_end( RESIDUES.End());
        residue_itr != residue_itr_end;
        ++residue_itr, ++seqid
      )
      {
        aasequence->PushBack( AminoAcidFromResidue( *residue_itr, seqid));
        // check that chain id matches
        BCL_Assert
        (
          aasequence->GetLastAA()->GetChainID() == CHAIN_ID,
          std::string( "amino acid has wrong chain id, should be: ") + CHAIN_ID +
          " but is " + aasequence->GetLastAA()->GetChainID()
        );
      }

      //end
      return aasequence;
    }

    //! @brief cut secondary structure elements from a given sequence
    //! @param SEQEUENCE the amino acid sequence
    //! @param SSE_RESIDUE secondary structure element defined by starting and end residue
    //! @return ShPtr to SSE, undefined if error occurs (no such amino acids, wrong chain, wrong amino acid order etc.)
    util::ShPtr< assemble::SSE> Factory::SecondaryStructureElementFromSequence
    (
      const biol::AASequence &SEQUENCE,
      const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_RESIDUE
    )
    {
      biol::AASequence::const_iterator aa_itr_first( SEQUENCE.Begin()), aa_itr_last( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());

      // find first amino acid of sse
      while( aa_itr_first != aa_itr_end && SSE_RESIDUE.Second() != **aa_itr_first)
      {
        ++aa_itr_first;
      }

      // find last amino acid of sse
      while( aa_itr_last != aa_itr_end && SSE_RESIDUE.Third() != **aa_itr_last)
      {
        ++aa_itr_last;
      }

      // if either could not be found, error
      if( aa_itr_first == aa_itr_end || aa_itr_last == aa_itr_end)
      {
        BCL_MessageVrb
        (
          "unable to find SSE in given SEQUENCE: " + SSE_RESIDUE.First().GetName() + " " +
          SSE_RESIDUE.Second().GetIdentification() + " " + SSE_RESIDUE.Third().GetIdentification()
        );

        return util::ShPtr< assemble::SSE>();
      }

      // make last the end
      ++aa_itr_last;

      if( std::distance( aa_itr_first, aa_itr_last) < 0)
      {
        BCL_MessageStd( "given sse residues are of no valid sse definition!");

        BCL_MessageStd( "residues: " + util::Format()( SSE_RESIDUE));
        BCL_MessageStd( "first: " + ( *aa_itr_first)->GetIdentification());
        BCL_MessageStd( "last: " + ( *aa_itr_last)->GetIdentification());
        return util::ShPtr< assemble::SSE>();
      }

      // return new sse from subsequence
      util::ShPtr< assemble::SSE> sp_sse
        (
          new assemble::SSE
          (
            biol::AASequence( util::ShPtrVector< biol::AABase>( aa_itr_first, aa_itr_last), SEQUENCE.GetChainID(), SEQUENCE.GetFastaHeader()),
            SSE_RESIDUE.First()
          )
        );

      // iterate over the residues in the sequence
      for
      (
        biol::AASequence::iterator
          aa_itr( sp_sse->Begin()), aa_itr_end( sp_sse->End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // set the data
        ( *aa_itr)->SetSSPrediction
        (
          sspred::GetMethods().e_PDB,
          sspred::PDB( SSE_RESIDUE.First(), biol::GetEnvironmentTypes().e_Solution)
        );
      }

      // end
      return sp_sse;
    }

    //! @brief cut secondary structure elements from a given sequence
    //! @param SEQEUENCE the amino acid sequence
    //! @param SSE_RESIDUES secondary structure elements defined by starting and end residue
    //! @return ShPtrList of SSEs
    util::ShPtrVector< assemble::SSE> Factory::SecondaryStructureElementsFromSequence
    (
      const biol::AASequence &SEQUENCE,
      const storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > &SSE_RESIDUES
    )
    {
      util::ShPtrVector< assemble::SSE> sses;

      // iterate through secondary structure element definitions
      for
      (
        storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> >::const_iterator
           sse_itr( SSE_RESIDUES.Begin()), sse_itr_end( SSE_RESIDUES.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        const util::ShPtr< assemble::SSE> sp_sse( SecondaryStructureElementFromSequence( SEQUENCE, *sse_itr));
        if( sp_sse.IsDefined())
        {
          sses.PushBack( sp_sse);
        }
      }

      // end
      return sses;
    }

    //! @brief builds a AASequence> from a pdb file read by the given handler
    //! @param HANDLER Handler that contains pdb file
    //! @param PDB_ID Optional PDB ID to be used in the fasta header, otherwise, it is retrieved from the handler
    //! @return ShPtrVector of AASequences that were read from the given Handler
    util::ShPtrVector< biol::AASequence>
    Factory::AASequencesFromPDB( const Handler &HANDLER, const std::string &PDB_ID) const
    {
      // instantiate ShPtrVector of chains
      util::ShPtrVector< biol::AASequence> sequences;

      // get the first model
      const Model &first_model( HANDLER.GetModels().FirstElement());

      //get iterator on pdb structure
      //loop over all chains and instantiate AASequences
      for
      (
        storage::Map< char, storage::List< Residue> >::const_iterator
          chain_itr( first_model.GetStructuredChains().Begin()), chain_itr_end( first_model.GetStructuredChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        sequences.PushBack
        (
          AASequenceFromResidues( PDB_ID == "" ? HANDLER.GetHead().GetPDBID() : PDB_ID, chain_itr->first, chain_itr->second)
        );
      }

      //return sequneces
      return sequences;
    }

    //! @brief merge two overlapping sequences into one
    //! @param SEQUENCE_A sequence a
    //! @param SEQUENCE_B sequence b
    //! @return the new sequence - empty if sses are not overlapping/consecutive
    biol::AASequence Factory::MergeOverlappingSequences
    (
      const biol::AASequence &SEQUENCE_A,
      const biol::AASequence &SEQUENCE_B
    )
    {
      if( !biol::DoOverlap( SEQUENCE_A, SEQUENCE_B))
      {
        return biol::AASequence();
      }

      util::SiPtr< const biol::AASequence> first_seq( SEQUENCE_A);
      util::SiPtr< const biol::AASequence> second_seq( SEQUENCE_B);

      if( first_seq->GetFirstAA()->GetSeqID() <= second_seq->GetLastAA()->GetSeqID())
      {
        std::swap( first_seq, second_seq);
      }

      biol::AASequence new_seq( *first_seq);
      // search for last amino acid of first sequence in the second sequence
      const util::ShPtr< biol::AAData> &last_aa_data( new_seq.GetLastAA()->GetData());
      biol::AASequence::const_iterator aa_itr( second_seq->Begin()), aa_itr_end( second_seq->End());
      for( ; aa_itr != aa_itr_end && ( *aa_itr)->GetData() != last_aa_data; ++aa_itr);

      if( aa_itr == aa_itr_end)
      {
        return new_seq;
      }

      // append the rest of second sequence
      ++aa_itr;
      for( ; aa_itr != aa_itr_end; ++aa_itr)
      {
        new_seq.PushBack( aa_itr->HardCopy());
      }

      // end
      return new_seq;
    }

    //! @brief remove the amino acids that are common to both sses
    //! @param SSE_A overlapping sse a
    //! @param SSE_B overlapping sse b
    //! @return two sses that do not have the amino acids that were common to the two argument sses
    storage::VectorND< 2, util::ShPtr< assemble::SSE> > Factory::RemoveOverlappingAAs
    (
      const util::ShPtr< assemble::SSE> &SSE_A,
      const util::ShPtr< assemble::SSE> &SSE_B
    )
    {
      if( !biol::DoOverlap( *SSE_A, *SSE_B))
      {
        return storage::VectorND< 2, util::ShPtr< assemble::SSE> >( SSE_A, SSE_B);
      }

      util::ShPtr< assemble::SSE> sse_first( SSE_A);
      util::ShPtr< assemble::SSE> sse_second( SSE_B);
      if( sse_second->GetFirstAA()->GetSeqID() <= sse_first->GetFirstAA()->GetSeqID())
      {
        std::swap( sse_first, sse_second);
      }

      biol::AASequence::const_iterator aa_itr1( sse_first->Begin()), aa_itr1_end( sse_first->End());
      biol::AASequence::const_iterator aa_itr2( sse_second->Begin()), aa_itr2_end( sse_second->End());

      // fill first
      util::ShPtrVector< biol::AABase> new_first;
      while( aa_itr1 != aa_itr1_end && ( *aa_itr1)->GetData() != ( *aa_itr2)->GetData())
      {
        new_first.PushBack( aa_itr1->HardCopy());
        ++aa_itr1;
      }

      // skip common amino acids
      for
      (
        ;
        aa_itr1 != aa_itr1_end && aa_itr2 != aa_itr2_end && ( *aa_itr1)->GetData() == ( *aa_itr2)->GetData();
        ++aa_itr1, ++aa_itr2
      );

      // fill second
      util::ShPtrVector< biol::AABase> new_second;
      while( aa_itr2 != aa_itr2_end)
      {
        new_second.PushBack( aa_itr2->HardCopy());
        ++aa_itr2;
      }

      // end
      return storage::VectorND< 2, util::ShPtr< assemble::SSE> >
      (
        util::ShPtr< assemble::SSE>( new_first.IsEmpty() ? NULL : new assemble::SSE( biol::AASequence( new_first, sse_first->GetChainID()), sse_first->GetType())),
        util::ShPtr< assemble::SSE>( new_second.IsEmpty() ? NULL : new assemble::SSE( biol::AASequence( new_second, sse_second->GetChainID()), sse_second->GetType()))
      );
    }

    //! @brief process overlapping secondary structure elements
    //! @param SECONDARY_STRUCTURE_ELEMENTS possibly overlapping secondary structure elements
    //! @param MERGE_OVERLAPPING merge overlapping sses, if they are of the same type
    //! @return storage set of non overlapping SSEs
    storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
    Factory::ProcessOverlappingSSEs
    (
      const util::ShPtrVector< assemble::SSE> &SECONDARY_STRUCTURE_ELEMENTS,
      const bool MERGE_OVERLAPPING
    )
    {
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> non_overlapping_sses;

      // iterate over sses
      for
      (
        util::ShPtrVector< assemble::SSE>::const_iterator
          sse_itr( SECONDARY_STRUCTURE_ELEMENTS.Begin()), sse_itr_end( SECONDARY_STRUCTURE_ELEMENTS.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> sp_sse_insert( *sse_itr);

        // search overlapping sse
        std::set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::iterator overlap_itr;

        // modify the sse to insert until there is no overlapping sse left
        while( sp_sse_insert.IsDefined() && ( overlap_itr = non_overlapping_sses.Find( sp_sse_insert)) != non_overlapping_sses.End())
        {
          // copy overlapping sse and remove from set
          util::ShPtr< assemble::SSE> sp_sse_overlap( *overlap_itr);

          // identical sses
          if
          (
               sp_sse_insert->GetFirstAA()->GetData() == sp_sse_overlap->GetFirstAA()->GetData()
            && sp_sse_insert->GetLastAA()->GetData() == sp_sse_overlap->GetLastAA()->GetData()
          )
          {
            sp_sse_insert = util::ShPtr< assemble::SSE>();
            break;
          }
          non_overlapping_sses.RemoveElement( overlap_itr);

          // merge overlapping
          if( MERGE_OVERLAPPING)
          {
            // check if overlapping sse has different type
            if( sp_sse_overlap->GetType() != sp_sse_insert->GetType())
            {
              // remove the sse
              BCL_MessageCrt
              (
                "found overlapping sses of different type; cannot merge them -> will not consider either\n:" +
                sp_sse_insert->GetIdentification() + "\n" + sp_sse_overlap->GetIdentification()
              );
              sp_sse_insert = util::ShPtr< assemble::SSE>();
              break;
            }
            else
            {
              sp_sse_insert = util::ShPtr< assemble::SSE>
              (
                new assemble::SSE( MergeOverlappingSequences( *sp_sse_insert, *sp_sse_overlap), sp_sse_insert->GetType())
              );
            }
          }
          // split overlapping by removing all overlapping amino acids
          else
          {
            // remove overlapping amino acids
            const storage::VectorND< 2, util::ShPtr< assemble::SSE> > new_sses( RemoveOverlappingAAs( sp_sse_overlap, sp_sse_insert));
            if( new_sses.First().IsDefined() && non_overlapping_sses.Find( new_sses.First()) == non_overlapping_sses.End())
            {
              non_overlapping_sses.Insert( new_sses.First());
              sp_sse_insert = new_sses.Second();
            }
            else if( new_sses.Second().IsDefined() && non_overlapping_sses.Find( new_sses.Second()) == non_overlapping_sses.End())
            {
              non_overlapping_sses.Insert( new_sses.Second());
              sp_sse_insert = new_sses.First();
            }
            else
            {
              BCL_MessageCrt
              (
                "unable to merge SSEs:\n" +
                sp_sse_insert->GetIdentification() + '\n' + sp_sse_overlap->GetIdentification()
              );
              continue;
            }
          }
        }

        // insert sse
        if( sp_sse_insert.IsDefined())
        {
          non_overlapping_sses.Insert( sp_sse_insert);
        }
      }

      // end
      return non_overlapping_sses;
    }

    //! @brief add loops to a set of secondary structure elements
    //! @param SECONDARY_STRUCTURE_ELEMENTS secondary structure elements
    //! @param SEQUENCE the full sequence; coils will be subsequences
    //! @return the set of secondary structure elements with coils
    storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
    Factory::AddLoops
    (
      const storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> &SECONDARY_STRUCTURE_ELEMENTS,
      const biol::AASequence &SEQUENCE
    )
    {
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> all_sses;

      // check that there is a least one sse
      if( SECONDARY_STRUCTURE_ELEMENTS.IsEmpty())
      {
        // create ShPtr to SSE "new_sse"
        util::ShPtr< assemble::SSE> new_sse( new assemble::SSE( SEQUENCE, biol::GetSSTypes().COIL));

        // iterate over amino acids and set the prediction
        for
        (
          biol::AASequence::iterator aa_itr( new_sse->Begin()), aa_itr_end( new_sse->End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // set the data
          ( *aa_itr)->SetSSPrediction
          (
            sspred::GetMethods().e_PDB,
            sspred::PDB( biol::GetSSTypes().COIL, biol::GetEnvironmentTypes().e_Solution)
          );
        }

        all_sses.Insert( new_sse);

        return all_sses;
      }

      // create ShPtrVector iterator "seq_itr" and "seq_itr_end" for iterating over the aa sequence of this chain
      biol::AASequence::const_iterator seq_itr( SEQUENCE.Begin()), seq_itr_end( SEQUENCE.End());

      // loop over all SSEs in this chain
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( SECONDARY_STRUCTURE_ELEMENTS.Begin()), sse_itr_end( SECONDARY_STRUCTURE_ELEMENTS.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // create size_t "sse_begin_seq_id" initialize with the SeqID of the first aa in the SSE denoted by "sse_itr"
        const int sse_begin_seq_id( ( *sse_itr)->GetFirstAA()->GetSeqID());

        // create ShPtr to SSE "new_sse" and initialize with a new SSE of type COIL
        biol::AASequence new_seq;

        // set the chain ID
        new_seq.SetChainID( SEQUENCE.GetChainID());

        // set the fasta header
        new_seq.SetFastaHeader( SEQUENCE.GetFastaHeader());

        // while seq ID of aa denoted by "seq_itr" has not reached "sse_begin_seq_id" and the end of the aa sequence
        // has not been reached
        while( ( *seq_itr)->GetSeqID() != sse_begin_seq_id && seq_itr != seq_itr_end)
        {
          // create a new amino acid of the corresponding aa class (should set coordinates to undefined)
          util::ShPtr< biol::AABase> new_aa( ( *seq_itr)->Clone());

          // set the data
          new_aa->SetSSPrediction
          (
            sspred::GetMethods().e_PDB,
            sspred::PDB( biol::GetSSTypes().COIL, biol::GetEnvironmentTypes().e_Solution)
          );

          // add the new aa to "new_sse"
          new_seq.PushBack( new_aa);

          // move "seq_itr" to the next aa
          ++seq_itr;
        }

        // true if "new_sse" has some aas in it
        if( new_seq.GetSize() > 0)
        {
          // insert "new_sse" into "new_sses"
          all_sses.Insert( util::ShPtr< assemble::SSE>( new assemble::SSE( new_seq, biol::GetSSTypes().COIL)));
        }

        // move "seq_itr" to the place of the next non-sse region
        const size_t sse_size( ( *sse_itr)->GetSize());
        for( size_t sse_index( 0); sse_index != sse_size && seq_itr != seq_itr_end; ++sse_index)
        {
          ++seq_itr;
        }

        // add the current sse also to the all sses
        all_sses.Insert( all_sses.End(), *sse_itr);
      }

      // we have to make sure to add the loop after the last SSE
      const size_t last_sse_end_id( ( *SECONDARY_STRUCTURE_ELEMENTS.ReverseBegin())->GetLastAA()->GetSeqID());
      const size_t last_aa_id( SEQUENCE.GetLastAA()->GetSeqID());

      // if there are any residues at the end of the sequence that are not in the last SSE
      if( last_sse_end_id < last_aa_id)
      {
        // create ShPtr to SSE "new_sse" and initialize with a new SSE of type COIL
        biol::AASequence new_seq;

        // set the chain ID
        new_seq.SetChainID( SEQUENCE.GetChainID());

        // set the fasta header
        new_seq.SetFastaHeader( SEQUENCE.GetFastaHeader());

        // iterate until you hit the end of amino acids in the sequecne
        while( seq_itr != seq_itr_end)
        {
          // create a new amino acid of the corresponding aa class (should set coordinates to 0.000)
          util::ShPtr< biol::AABase> new_aa( ( *seq_itr)->Clone());

          // set the data
          new_aa->SetSSPrediction
          (
            sspred::GetMethods().e_PDB,
            sspred::PDB( biol::GetSSTypes().COIL, biol::GetEnvironmentTypes().e_Solution)
          );

          // add the new aa to "new_sse"
          new_seq.PushBack( new_aa);

          // move "seq_itr" to the next aa
          ++seq_itr;
        }

        // push the last loop
        all_sses.Insert( util::ShPtr< assemble::SSE>( new assemble::SSE( new_seq, biol::GetSSTypes().COIL)));
      }

      // end
      return all_sses;
    }

    //! @brief builds a Chain with sequence and SSEs from given Residues and SSE Residues
    //! builds a chain of template amino acids from a chainID and a util::ShPtrVector of SSEs represented as
    //! util::ShPtrVector of PDB Residues - each secondary structure element will have a minimal length of SSE_MINSIZE
    //! @param CHAIN_ID ChainID of the chain of interest
    //! @param SSE_RESIDUES Residues belonging to SSEs
    //! @param SEQUENCE_RESIDUES Residues from the sequence
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @param PDB_ID PdbID of the protein
    //! @return Chain with the full sequence and SSEs read from given pdb information
    util::ShPtr< assemble::Chain>
    Factory::ChainFromPDBSSElementsAndSequence
    (
      const char CHAIN_ID,
      const storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > &SSE_RESIDUES,
      const storage::List< Residue> &SEQUENCE_RESIDUES,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE,
      const std::string &PDB_ID
    ) const
    {
      // create the current chain with the sequence
      util::ShPtr< biol::AASequence> sp_sequence( AASequenceFromResidues( PDB_ID, CHAIN_ID, SEQUENCE_RESIDUES));
      util::ShPtr< assemble::Chain> current_chain;

      // assign secondary structure usign dssp
      if( GetFlagDSSP()->GetFlag())
      {
        assemble::ProteinModel model;
        util::ShPtr< assemble::Chain> tmp_chain( new assemble::Chain( sp_sequence));
        tmp_chain->Insert( util::ShPtr< assemble::SSE>( new assemble::SSE( *sp_sequence, biol::GetSSTypes().COIL)));
        model.Insert( tmp_chain);

        const biol::DSSP dssp( GetFlagDSSP()->GetFirstParameter()->GetNumericalValue< double>());
        math::MutateResult< assemble::ProteinModel> result( dssp( model));
        if( result.GetArgument().IsDefined())
        {
          current_chain = result.GetArgument()->GetChain( sp_sequence->GetChainID());
          BCL_MessageVrb( "Used DSSP to assign SSEs; found " + util::Format()( current_chain->GetData().GetSize()) + " SSEs");
        }
      }
      // if the sses are supposed to be created from the backbone conformation ignoring and no sses are given in pdb for that chain
      else if( GetFlagSSEsFromBackBone()->GetFlag())
      {
        // create chain with sses from the conformation
        current_chain = util::ShPtr< assemble::Chain>
        (
          assemble::ConstructChainWithSSEsFromConformation( sp_sequence).Clone()
        );
      }

      else
      {
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> non_overlapping_sses
        (
          ProcessOverlappingSSEs( SecondaryStructureElementsFromSequence( *sp_sequence, SSE_RESIDUES), false)
        );

        non_overlapping_sses = AddLoops( non_overlapping_sses, *sp_sequence);

        current_chain = util::ShPtr< assemble::Chain>
        (
          new assemble::Chain( sp_sequence, util::ShPtrVector< assemble::SSE>( non_overlapping_sses.Begin(), non_overlapping_sses.End()))
        );
      }

      // for each sse type remove sses that are too short are not even considered
      current_chain->FilterByMinSSESizes( SSE_TYPE_MINSIZE);

      // return Chain
      return current_chain;
    }

    //! builds chain without sses
    //! @param CHAIN_ID chain id of given fasta
    //! @param ISTREAM fasta file stream
    //! @return chain without SSEs read from the given fasta stream
    util::ShPtr< assemble::Chain>
    Factory::ChainFromFastaStream( const char CHAIN_ID, std::istream &ISTREAM) const
    {
      // create sequence from fasta
      util::ShPtr< biol::AASequence> sp_aa_seq
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( ISTREAM, m_AAClass, CHAIN_ID))
      );

      //! create chain
      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( sp_aa_seq));

      // return
      return sp_chain;
    }

    //! @brief builds model of with one or more chains, each formed of full sequence and SSEs from given Handler
    //! The types and sizes of SSEs to be considered are determined by the additional parameters
    //! @param HANDLER Handler that contains pdb information
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @return Protein Model read from given Handler
    assemble::ProteinModel
    Factory::ProteinModelFromPDB
    (
      const Handler &HANDLER,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE
    ) const
    {
      return ProcessModel( HANDLER, SSE_TYPE_MINSIZE, HANDLER.GetModels().FirstElement());
    }

    //! @brief builds an ensemble from the models stored in the handler
    //! @param HANDLER Handler that contains pdb information
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @return Protein ensemble read from given Handler
    assemble::ProteinEnsemble Factory::ProteinEnsembleFromPDB
    (
      const Handler &HANDLER,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE
    ) const
    {
      // initialize ensemble
      assemble::ProteinEnsemble ensemble;

      // iterate over the models
      for
      (
        storage::Vector< Model>::const_iterator model_itr( HANDLER.GetModels().Begin()),
          model_itr_end( HANDLER.GetModels().End());
        model_itr != model_itr_end; ++model_itr
      )
      {
        // build the model and add it to the ensemble
        ensemble.InsertElement( util::CloneToShPtr( ProcessModel( HANDLER, SSE_TYPE_MINSIZE, *model_itr)));
      }

      // end
      return ensemble;
    }

    //! @brief creates and returns a protein model based on a PDB filename
    //! @param PDB_FILENAME is the pdb filename from which the protein model will be created
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @param IGNORE_CLASH whether to ignore clashes or not
    //! @return ProteinModel which was created from "PDB_FILENAME"
    assemble::ProteinModel Factory::ProteinModelFromPDBFilename
    (
      const std::string &PDB_FILENAME,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE,
      const bool IGNORE_CLASH
    ) const
    {
      // initialize write and read stream objects
      io::IFStream read;

      // open pdb file
      io::File::MustOpenIFStream( read, PDB_FILENAME);

      // true is used to advice handler to ignore clashes in the structure and insert residues as they are suggested by
      // the numbering without regard for the bond information
      Handler pdb( read, IGNORE_CLASH);
      io::File::CloseClearFStream( read);

      // create a protein model from the handler
      assemble::ProteinModel model( ProteinModelFromPDB( pdb, SSE_TYPE_MINSIZE));

      // add filename to the protein model data
      util::ShPtr< assemble::ProteinModelData> sp_data( model.GetProteinModelData().HardCopy());
      sp_data->Insert
      (
        assemble::ProteinModelData::e_PDBFile,
        util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( PDB_FILENAME))
      );
      model.SetProteinModelData( sp_data);
      return model;
    }

    //! @brief returns the map of SSTypes and their corresponding minimal size
    //! @param FLAG the flag to get ss type min sizes from
    //! @return map of ss type and associated value
    storage::Map< biol::SSType, size_t>
    Factory::GetCommandlineSSETypeMinSizes( const command::FlagInterface &FLAG)
    {
      return
        GetSSETypeMinSizes
        (
          FLAG.GetParameterList()( biol::GetSSTypes().HELIX)->GetNumericalValue< size_t>(),
          FLAG.GetParameterList()( biol::GetSSTypes().STRAND)->GetNumericalValue< size_t>(),
          FLAG.GetParameterList()( biol::GetSSTypes().COIL)->GetNumericalValue< size_t>()
        );
    }

    //! @brief returns the map of SSTypes and their corresponding minimal size
    //! @param MIN_HELIX_SIZE the minimum helix size
    //! @param MIN_STRAND_SIZE the minimum strand size
    //! @param MIN_COIL_SIZE the minimum coil size
    //! @return map of ss type and associated value
    storage::Map< biol::SSType, size_t>
    Factory::GetSSETypeMinSizes
    (
      const size_t &MIN_HELIX_SIZE,
      const size_t &MIN_STRAND_SIZE,
      const size_t &MIN_COIL_SIZE
    )
    {
      storage::Map< biol::SSType, size_t> sse_type_min_sizes;
      sse_type_min_sizes[ biol::GetSSTypes().HELIX] = MIN_HELIX_SIZE;
      sse_type_min_sizes[ biol::GetSSTypes().STRAND] = MIN_STRAND_SIZE;
      sse_type_min_sizes[ biol::GetSSTypes().COIL] = MIN_COIL_SIZE;

      const storage::Set< biol::SSType> helix_classes( Handler::GetHelixClassesFromCommandLine());

      // iterate over the helix types
      for
      (
        storage::Set< biol::SSType>::const_iterator
          helix_class_itr( helix_classes.Begin()), helix_class_itr_end( helix_classes.End());
        helix_class_itr != helix_class_itr_end;
        ++helix_class_itr
      )
      {
        sse_type_min_sizes[ *helix_class_itr] = MIN_HELIX_SIZE;
      }

      // end
      return sse_type_min_sizes;
    }

  ////////////////////////
  // operations - write //
  ////////////////////////

    //! @brief write header information to Line
    //! @param CLASSIFICATION string for header classification
    //! @param TIME time (date) information to be placed in the header
    //! @return ShPtr to Line that contains Header information
    util::ShPtr< Line> Factory::WriteHeaderToLine( const std::string &CLASSIFICATION, const util::Time &TIME)
    {
      // instantiate new header line
      util::ShPtr< Line> header_line( new Line( GetLineTypes().HEADER));

      // add the entries
      header_line->Put( GetEntryTypes().HEADERClassification, CLASSIFICATION);
      header_line->Put( GetEntryTypes().HEADERDate, ConvertTimeToPDBDate( TIME));
      header_line->Put( GetEntryTypes().HEADERIDCode, Head::GetBCLPdbID());

      // return the Line
      return header_line;
    }

    //! @brief write Atom to Line of a certain Residue and the CHAINID and the atom SERIAL
    //! @param ATOM Atom of interest
    //! @param AMINO_ACID Amino acid which atom belongs to
    //! @param CHAIN_ID chain id of the sequence atom belongs to
    //! @param SERIAL serial id of the atom
    //! @return ShPtr to Line that contains Atom information
    util::ShPtr< Line> Factory::WriteAtomToLine
    (
      const biol::Atom &ATOM,
      const biol::AABase &AMINO_ACID,
      const char CHAIN_ID,
      const size_t SERIAL
    )
    {
      util::ShPtr< Line> atom_line;

      if( AMINO_ACID.GetType()->IsNaturalAminoAcid())
      {
        //instantiate new atom line
        atom_line = util::ShPtr< Line>( new Line( GetLineTypes().ATOM));

        // write defined coordinates
        if( ATOM.GetCoordinates().IsDefined())
        {
          atom_line->PutCoordinates( ATOM.GetCoordinates());

          // write default occupancy of 1.00 with the precision 2
          atom_line->Put( GetEntryTypes().ATOMOccupancy, double( 1));
        }
        // true if need to write zero coordinates for undefined atoms
        else if( GetFlagWriteZeroCoordinatesForUndefinedAminoAcids()->GetFlag())
        {
          atom_line->PutCoordinates( linal::Vector3D( 0.0, 0.0, 0.0));

          // write default occupancy of -1.00 with the precision 2 for missing residues, compatible to rosetta
          atom_line->Put( GetEntryTypes().ATOMOccupancy, double( -1));
        }
        else
        {
          // return empty line
          return util::ShPtr< Line>();
        }

        // write atom ID
        if( GetFlagWritePDBAtomID()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().ATOMSerial, ATOM.GetPdbID());
        }
        else
        {
          atom_line->Put( GetEntryTypes().ATOMSerial, SERIAL);
        }

        // use the pdb atom name instead of the IUPAC name
        if( GetFlagPDBAtomName()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().ATOMName, AMINO_ACID.GetType()->GetPDBAtomName( ATOM.GetType()));
        }
        else
        {
          atom_line->Put( GetEntryTypes().ATOMName, ATOM.GetType()->AtomTypeData::GetName());
        }

        // residue name and chain id
        atom_line->Put( GetEntryTypes().ATOMResidueName, AMINO_ACID.GetType()->GetThreeLetterCode());
        atom_line->Put( GetEntryTypes().ATOMChainID, CHAIN_ID);

        // write res ID
        if( GetFlagWritePDBResID()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().ATOMResidueSequenceID, AMINO_ACID.GetPdbID());
          atom_line->Put( GetEntryTypes().ATOMInsertionCode, AMINO_ACID.GetPdbICode());
        }
        else
        {
          atom_line->Put( GetEntryTypes().ATOMResidueSequenceID, AMINO_ACID.GetSeqID());
        }

        // write temp factor with the precision 2
        atom_line->Put
        (
          GetEntryTypes().ATOMTempFactor,
          util::IsDefined( ATOM.GetBFactor()) ? ATOM.GetBFactor() : double( 0.0)
        );

        // element type
        atom_line->Put( GetEntryTypes().ATOMElement, std::string( ATOM.GetType()->GetElementType()->GetChemicalSymbol()));
      }
      // non natural amino acid as HETATM line
      else
      {
        //instantiate new atom line
        atom_line = util::ShPtr< Line>( new Line( GetLineTypes().HETATM));

        // write defined coordinates
        if( ATOM.GetCoordinates().IsDefined())
        {
          atom_line->PutCoordinates( ATOM.GetCoordinates());

          // write default occupancy of 1.00 with the precision 2
          atom_line->Put( GetEntryTypes().HETATMOccupancy, double( 1));
        }
        // true if need to write zero coordinates for undefined atoms
        else if( GetFlagWriteZeroCoordinatesForUndefinedAminoAcids()->GetFlag())
        {
          atom_line->PutCoordinates( linal::Vector3D( 0.0, 0.0, 0.0));

          // write default occupancy of -1.00 with the precision 2 for missing residues, compatible to rosetta
          atom_line->Put( GetEntryTypes().HETATMOccupancy, double( -1));
        }
        else
        {
          // return empty line
          return util::ShPtr< Line>();
        }

        // write atom ID
        if( GetFlagWritePDBAtomID()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().HETATMSerial, ATOM.GetPdbID());
        }
        else
        {
          atom_line->Put( GetEntryTypes().HETATMSerial, SERIAL);
        }

        // use the pdb atom name instead of the IUPAC name
        if( GetFlagPDBAtomName()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().HETATMName, AMINO_ACID.GetType()->GetPDBAtomName( ATOM.GetType()));
        }
        else
        {
          atom_line->Put( GetEntryTypes().HETATMName, ATOM.GetType()->AtomTypeData::GetName());
        }

        // residue name and chain id
        atom_line->Put( GetEntryTypes().HETATMResidueName, AMINO_ACID.GetType()->GetThreeLetterCode());
        atom_line->Put( GetEntryTypes().HETATMChainID, CHAIN_ID);

        // write res ID
        if( GetFlagWritePDBResID()->GetFlag())
        {
          atom_line->Put( GetEntryTypes().HETATMResidueSequenceID, AMINO_ACID.GetPdbID());
          atom_line->Put( GetEntryTypes().HETATMInsertionCode, AMINO_ACID.GetPdbICode());
        }
        else
        {
          atom_line->Put( GetEntryTypes().HETATMResidueSequenceID, AMINO_ACID.GetSeqID());
        }

        // write temp factor with the precision 2
        atom_line->Put
        (
          GetEntryTypes().HETATMTempFactor,
          util::IsDefined( ATOM.GetBFactor()) ? ATOM.GetBFactor() : double( 0.0)
        );

        // element type
        atom_line->Put( GetEntryTypes().HETATMElement, std::string( ATOM.GetType()->GetElementType()->GetChemicalSymbol()));
      }

      //return
      return atom_line;
    }

    //! @brief write Residue to Lines of a certain CHAINID and the atom SERIAL
    //! @param AMINO_ACID Amino acid of interest
    //! @param CHAIN_ID chain id of the sequence amino acid belongs to
    //! @param SERIAL serial id of the amino acid
    //! @return ShPtr to Line that contains amino acid information
    util::ShPtrList< Line> Factory::WriteResiduesToLines
    (
      const biol::AABase &AMINO_ACID,
      const char CHAIN_ID,
      size_t &SERIAL
    )
    {
      //instantiate residue lines of amino acids
      util::ShPtrList< Line> residue_lines;

      //loop over all atoms in residue
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator atom_itr( AMINO_ACID.GetAtoms().Begin()),
          atom_itr_end( AMINO_ACID.GetAtoms().End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        // in non debug mode
        if( !util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
        {
          // by default, write all atom atoms that are not hydrogens or have defined coordinates
          // write hydrogen atoms only if the corresponding flag is set
          if
          (
               !GetFlagWriteHydrogens()->GetFlag()
            && ( *atom_itr)->GetType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen
          )
          {
            continue;
          }
        }

        // ShPtr to new line
        const util::ShPtr< Line> sp_line( WriteAtomToLine( **atom_itr, AMINO_ACID, CHAIN_ID, SERIAL));
        if( sp_line.IsDefined())
        {
          // write the actual line
          residue_lines.PushBack( sp_line);
          ++SERIAL;
        }
      }

      //return
      return residue_lines;
    }

    //! @brief write all atoms of a sequence to a list of lines to be used eventually for writing to a pdb
    //! @param AA_SEQUENCE AASequence of interest
    //! @param SERIAL serial id of the sequence
    //! @return ShPtrList of Lines that contain AASequence information
    util::ShPtrList< Line> Factory::WriteAASequenceToLines( const biol::AASequence &AA_SEQUENCE, size_t &SERIAL)
    {
      //instantiate lines for the entire sequence
      util::ShPtrList< Line> sequence_lines;

      //loop over all residues
      for
      (
        biol::AASequence::const_iterator residue_itr( AA_SEQUENCE.GetData().Begin()),
          residue_itr_end( AA_SEQUENCE.GetData().End());
        residue_itr != residue_itr_end;
        ++residue_itr
      )
      {
        sequence_lines.Append( WriteResiduesToLines( **residue_itr, AA_SEQUENCE.GetChainID(), SERIAL));
      }

      //return
      return sequence_lines;
    }

    //! @brief write a Chain to a SharedpointerVector of Lines
    //! @param CHAIN Chain of interest
    //! @param SERIAL serial id
    //! @return ShPtrList of Lines that contain Chain information
    util::ShPtrList< Line> Factory::WriteChainToLines( const assemble::Chain &CHAIN, size_t &SERIAL)
    {
      //instantiate lines for the entire sequence
      util::ShPtrList< Line> chain_lines;

      //loop over all sses in chain
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( CHAIN.GetData().Begin()), sse_itr_end( CHAIN.GetData().End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        chain_lines.Append( WriteAASequenceToLines( ( **sse_itr), SERIAL));
      }

      // TER line
      if( chain_lines.IsEmpty())
      {
        return chain_lines;
      }

      const Line &last_line( *chain_lines.LastElement());

      util::ShPtr< Line> ter_line( new Line( GetLineTypes().TER));
      ter_line->Put( GetEntryTypes().TERSerial, SERIAL);
      ter_line->Put( GetEntryTypes().TERChainID, CHAIN.GetChainID());
      if( last_line.GetType() == GetLineTypes().ATOM)
      {
        ter_line->Put( GetEntryTypes().TERInsertionCode, last_line.GetChar( GetEntryTypes().ATOMInsertionCode));
        ter_line->Put( GetEntryTypes().TERResidueName, last_line.GetString( GetEntryTypes().ATOMResidueName));
        ter_line->Put( GetEntryTypes().TERResidueSequenceID, last_line.GetString( GetEntryTypes().ATOMResidueSequenceID));
      }
      else
      {
        ter_line->Put( GetEntryTypes().TERInsertionCode, last_line.GetChar( GetEntryTypes().HETATMInsertionCode));
        ter_line->Put( GetEntryTypes().TERResidueName, last_line.GetString( GetEntryTypes().HETATMResidueName));
        ter_line->Put( GetEntryTypes().TERResidueSequenceID, last_line.GetString( GetEntryTypes().HETATMResidueSequenceID));
      }
      chain_lines.PushBack( ter_line);
      ++SERIAL;

      //return
      return chain_lines;
    }

    //! @brief write a single Protein Model to pdb lines
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param WRITE_BODY_INFORMATION flag to determine whether to write body information
    //! @param CLASSIFICATION string for header classification
    //! @param TIME time (date) information to be placed in the header
    //! @return lines to which ProteinModel was written to
    util::ShPtrList< Line> Factory::WriteCompleteModelToPDBLines
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const bool WRITE_BODY_INFORMATION,
      const std::string &CLASSIFICATION,
      const util::Time &TIME
    ) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // the HEADER section
      if( CLASSIFICATION.empty())
      {
        lines.PushBack( util::ShPtr< Line>( new Line( std::string( "HEADER  BCL MODEL"))));
      }
      else
      {
        lines.PushBack( WriteHeaderToLine( CLASSIFICATION, TIME));
      }

      // the SEQRES section
      lines.Append( WriteSeqResToLines( PROTEIN_MODEL));

      size_t serial( 1);
      // the HELIX section
      lines.Append( WriteHelixDefinitionsToLines( PROTEIN_MODEL.GetSSEs( biol::GetSSTypes().GetHelixTypes()), serial));

      serial = 1;
      // the STRAND section
      lines.Append( WriteStrandDefinitionsToLines( PROTEIN_MODEL.GetSSEs( biol::GetSSTypes().STRAND), serial));

      serial = 1;
      // write the residues and atoms to pdblines
      lines.Append( WriteProteinModelToLines( PROTEIN_MODEL, serial));

      //write bodyinforamtions for each sse
      if( WRITE_BODY_INFORMATION)
      {
        lines.Append( WriteBodyInformationToLines( PROTEIN_MODEL));
      }

      // iterate over printers
      for
      (
        util::ShPtrList< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> > >::const_iterator
          itr( m_Printers.Begin()), itr_end( m_Printers.End());
        itr != itr_end; ++itr
      )
      {
        lines.Append( ( *itr)->operator ()( PROTEIN_MODEL));
      }

      // end
      return lines;
    }

    //! @brief write a ProteinModel to a SharedpointerVector of Lines
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param SERIAL serial id
    //! @return ShPtrList of Lines that contain ProteinModel information
    util::ShPtrList< Line> Factory::WriteProteinModelToLines
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      size_t &SERIAL
    )
    {
      //instantiate lines for the entire sequence
      util::ShPtrList< Line> proteinmodel_lines;

      //loop over all chains of proteinmodel
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        proteinmodel_lines.Append( WriteChainToLines( **chain_itr, SERIAL));
      }

      //return
      return proteinmodel_lines;
    }

    //! @brief write SSE definition for one SSE to Line
    //! @param THIS_SSE SSE of interest
    //! @param SERIAL SSE serial id
    //! @return ShPtr to Line containing SSE definition for given SSE
    util::ShPtr< Line> Factory::WriteSSEDefinitionToLine( const assemble::SSE &THIS_SSE, const size_t &SERIAL)
    {
      // if Helix type over SSType
      if( biol::GetSSTypes().GetHelixTypes().Contains( THIS_SSE.GetType()))
      {
        return WriteHelixDefinitionToLine( THIS_SSE, SERIAL);
      }
      // else if STRAND
      else if( THIS_SSE.GetType() == biol::GetSSTypes().STRAND)
      {
        return WriteStrandDefinitionToLine( THIS_SSE, SERIAL);
      }
      else
      {
        return util::ShPtr< Line>();
      }
    }

    //! @brief write helix definition for one SSE to Line
    //! @param THIS_SSE SSE of interest
    //! @param HELIX_SERIAL helix serial id
    //! @return ShPtr to Line containing helix SSE definition for given SSE
    util::ShPtr< Line> Factory::WriteHelixDefinitionToLine
    (
      const assemble::SSE &THIS_SSE,
      const size_t HELIX_SERIAL
    )
    {
      // empty sse
      if( THIS_SSE.GetSize() == 0)
      {
        return util::ShPtr< Line>();
      }

      util::ShPtr< Line> helix_line( new Line( GetLineTypes().HELIX));

      //write all information for helix definiton line
      //helix serial
      helix_line->Put( GetEntryTypes().HELIXSerial, HELIX_SERIAL);
      //helix id
      helix_line->Put( GetEntryTypes().HELIXID, HELIX_SERIAL);

      // first residue name
      helix_line->Put
      (
        GetEntryTypes().HELIXResidueName_Initial, THIS_SSE.GetFirstAA()->GetType()->GetThreeLetterCode()
      );
      // last residue name
      helix_line->Put
      (
        GetEntryTypes().HELIXResidueName_Terminal, THIS_SSE.GetLastAA()->GetType()->GetThreeLetterCode()
      );

      // for res id use either the original or the bcl numerated seqid
      if( GetFlagWritePDBResID()->GetFlag())
      {
        helix_line->Put( GetEntryTypes().HELIXSequenceID_Initial, THIS_SSE.GetFirstAA()->GetPdbID());
        helix_line->Put( GetEntryTypes().HELIXInsertionCode_Initial, THIS_SSE.GetFirstAA()->GetPdbICode());
        helix_line->Put( GetEntryTypes().HELIXSequenceID_Terminal, THIS_SSE.GetLastAA()->GetPdbID());
        helix_line->Put( GetEntryTypes().HELIXInsertionCode_Terminal, THIS_SSE.GetLastAA()->GetPdbICode());
      }
      else
      {
        helix_line->Put( GetEntryTypes().HELIXSequenceID_Initial, THIS_SSE.GetFirstAA()->GetSeqID());
        helix_line->Put( GetEntryTypes().HELIXSequenceID_Terminal, THIS_SSE.GetLastAA()->GetSeqID());
      }

      // chain ID
      helix_line->Put( GetEntryTypes().HELIXChainID_Initial, THIS_SSE.GetChainID());
      helix_line->Put( GetEntryTypes().HELIXChainID_Terminal, THIS_SSE.GetChainID());

      // helix class - currently right handed alpha helix
      helix_line->Put< size_t>( GetEntryTypes().HELIXClass, biol::GetSSTypes().PDBHelixClassFromSSType( THIS_SSE.GetType()));

      // length of helix
      helix_line->Put( GetEntryTypes().HELIXLength, THIS_SSE.GetSize());

      // return
      return helix_line;
    }

    //! @brief write Helix definitions to Lines
    //! @param SSE_VECTOR SiPtrVector of SSEs of interest
    //! @param HELIX_SERIAL helix serial id
    //! @return ShPtrList of Lines containing helix SSE defintions for given vector of SSEs
    util::ShPtrList< Line> Factory::WriteHelixDefinitionsToLines
    (
      const util::SiPtrVector< const assemble::SSE> &SSE_VECTOR,
      size_t &HELIX_SERIAL
    )
    {
      util::ShPtrList< Line> helixlines;

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          helix_itr( SSE_VECTOR.Begin()), helix_itr_end( SSE_VECTOR.End());
        helix_itr != helix_itr_end; ++helix_itr)
      {
        helixlines.PushBack( WriteHelixDefinitionToLine( **helix_itr, HELIX_SERIAL++));
      }

      return helixlines;
    }

    //! @brief write Strand definition for one SSE to Line
    //! @param THIS_SSE SSE of interest
    //! @param STRAND_SERIAL strand serial id
    //! @return ShPtr to Line containing strand SSE definition for given SSE
    util::ShPtr< Line> Factory::WriteStrandDefinitionToLine
    (
      const assemble::SSE &THIS_SSE,
      const size_t STRAND_SERIAL
    )
    {
      // empty sse
      if( THIS_SSE.GetSize() == 0)
      {
        return util::ShPtr< Line>();
      }

      util::ShPtr< Line> strand_line( new Line( GetLineTypes().SHEET));

      //write all information for strand definiton line
      //strand ID
      strand_line->Put( GetEntryTypes().SHEETStrandID, STRAND_SERIAL);
      //first residue name
      strand_line->Put
      (
        GetEntryTypes().SHEETResidueName_Initial, THIS_SSE.GetFirstAA()->GetType()->GetThreeLetterCode()
      );
      // last residue name
      strand_line->Put
      (
        GetEntryTypes().SHEETResidueName_Terminal, THIS_SSE.GetLastAA()->GetType()->GetThreeLetterCode()
      );

      // for res id use either the original or the bcl numerated seqid
      if( GetFlagWritePDBResID()->GetFlag())
      {
        strand_line->Put( GetEntryTypes().SHEETSequenceID_Initial, THIS_SSE.GetFirstAA()->GetPdbID());
        strand_line->Put( GetEntryTypes().SHEETInsertionCode_Initial, THIS_SSE.GetFirstAA()->GetPdbICode());
        strand_line->Put( GetEntryTypes().SHEETSequenceID_Terminal, THIS_SSE.GetLastAA()->GetPdbID());
        strand_line->Put( GetEntryTypes().SHEETInsertionCode_Terminal, THIS_SSE.GetLastAA()->GetPdbICode());
      }
      else
      {
        strand_line->Put( GetEntryTypes().SHEETSequenceID_Initial, THIS_SSE.GetFirstAA()->GetSeqID());
        strand_line->Put( GetEntryTypes().SHEETSequenceID_Terminal, THIS_SSE.GetLastAA()->GetSeqID());
      }

      //chain ID
      strand_line->Put( GetEntryTypes().SHEETChainID_Initial, THIS_SSE.GetChainID());
      strand_line->Put( GetEntryTypes().SHEETChainID_Terminal, THIS_SSE.GetChainID());

      // end
      return strand_line;
    }

    //! @brief write strand definitions to Lines
    //! @param SSE_VECTOR SiPtrVector of SSEs of interest
    //! @param STRAND_SERIAL strand serial id
    //! @return ShPtrList of Lines containing strand SSE definitions for given vector of SSEs
    util::ShPtrList< Line> Factory::WriteStrandDefinitionsToLines
    (
      const util::SiPtrVector< const assemble::SSE> &SSE_VECTOR,
      size_t &STRAND_SERIAL
    )
    {
      util::ShPtrList< Line> strandlines;

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          strand_itr( SSE_VECTOR.Begin()), strand_itr_end( SSE_VECTOR.End());
        strand_itr != strand_itr_end; ++strand_itr
      )
      {
        strandlines.PushBack( WriteStrandDefinitionToLine( **strand_itr, STRAND_SERIAL++));
      }

      return strandlines;
    }

    //! @brief write assemble::SSEGeometry Information to Lines of one SSE
    //! @param THIS_SSE SSE of interest
    //! @return ShPtrList of Lines containing body information for the given SSE
    util::ShPtrList< Line> Factory::WriteBodyInformationToLines( const assemble::SSE &THIS_SSE)
    {
      //instantiate lines for current body
      util::ShPtrList< Line> body_lines;

      //instantiate std::string storing the temporary line
      std::string tmp;

      //write extension
      tmp =
        std::string( "REMARK  ") + std::string( "Extension in Z direction: ") +
        util::Format()( THIS_SSE.GetExtent( coord::GetAxes().e_Z));
      body_lines.PushBack( util::ShPtr< Line>( new Line( tmp)));

      // insert atoms at the beginning and end of bodies
      util::ShPtr< Line> begin_atom( new Line( GetLineTypes().ATOM));
      util::ShPtr< Line> end_atom( new Line( GetLineTypes().ATOM));
      begin_atom->PutCoordinates( THIS_SSE.BeginOfZ());
      end_atom->PutCoordinates( THIS_SSE.EndOfZ());
      body_lines.PushBack( begin_atom);
      body_lines.PushBack( end_atom);

      //return
      return body_lines;
    }

    //! @brief write assemble::SSEGeometry Information to Lines for given ProteinModel
    //! writes TransformationMatrix and assemble::SSEGeometry extension of all SSES to Lines
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return ShPtrList of Lines containing body information for the given ProteinModel
    util::ShPtrList< Line> Factory::WriteBodyInformationToLines( const assemble::ProteinModel &PROTEIN_MODEL)
    {
      //instantiate lines for the entire sequence
      util::ShPtrList< Line> body_lines;

      //instantiate Lines for sheet and helix definitons
      util::ShPtrList< Line> sse_helix_definition_lines;
      util::ShPtrList< Line> sse_strand_definition_lines;

      //loop over all chains of proteinmodel
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        //loop over all sses
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          //according to type write sse line
          // the HELIX lines
          if( biol::GetSSTypes().GetHelixTypes().Contains( ( *sse_itr)->GetType()))
          {
            sse_helix_definition_lines.PushBack
            (
              WriteHelixDefinitionToLine( **sse_itr, sse_helix_definition_lines.GetSize())
            );
            sse_helix_definition_lines.Append( WriteBodyInformationToLines( **sse_itr));
          }
          // the STRAND lines
          else if( ( *sse_itr)->GetType() == biol::GetSSTypes().STRAND)
          {
            sse_strand_definition_lines.PushBack
            (
              WriteStrandDefinitionToLine( **sse_itr, sse_strand_definition_lines.GetSize())
            );
            sse_strand_definition_lines.Append( WriteBodyInformationToLines( **sse_itr));
          }
        }
      }

      util::ShPtrList< Line> sse_definition_lines( sse_helix_definition_lines);
      sse_definition_lines.Append( sse_strand_definition_lines);

      //end
      return sse_definition_lines;
    }

    //! @brief write SSE definitions to Lines for a util::ShPtrVector of models
    //! @param PROTEIN_MODELS ShPtrVector of ProteinModels of interest
    //! @return ShPtrList of Lines containing SSE definition information for given ShPtrVector of ProteinModels
    util::ShPtrList< Line> Factory::WriteSSEDefinitionsToLines
    (
      const util::ShPtrVector< assemble::ProteinModel> &PROTEIN_MODELS
    )
    {
      //instantiate Lines for sheet and helix definitons
      util::ShPtrList< Line> sse_helix_definition_lines;
      util::ShPtrList< Line> sse_strand_definition_lines;
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          model_itr( PROTEIN_MODELS.Begin()), model_itr_end( PROTEIN_MODELS.End());
        model_itr != model_itr_end; ++model_itr)
      {
        size_t serial( sse_helix_definition_lines.GetSize() + 1);
        // the HELIX lines
        sse_helix_definition_lines.Append
        (
          WriteHelixDefinitionsToLines( ( *model_itr)->GetSSEs( biol::GetSSTypes().GetHelixTypes()), serial)
        );

        serial = 1;
        // the STRAND lines
        sse_helix_definition_lines.Append
        (
          WriteStrandDefinitionsToLines( ( *model_itr)->GetSSEs( biol::GetSSTypes().STRAND), serial)
        );
      }

      //put Helix and strand definition lines in one util::ShPtrVector of Lines
      util::ShPtrList< Line> sse_definition_lines( sse_helix_definition_lines);
      sse_definition_lines.Append( sse_strand_definition_lines);

      //return helix and strand definiton lines
      return sse_definition_lines;
    }

    //! @brief write SEQRES lines for a sequence
    //! @param SEQUENCE AASequence of interest
    //! @return ShPtrList of Lines containing SEQRES information for given AASequence
    util::ShPtrList< Line> Factory::WriteSeqResToLines( const biol::AASequence &SEQUENCE)
    {
      const char chainid( SEQUENCE.GetChainID());

      size_t number_aas( 0);
      std::stringstream string_sequence;

      //iterate over all aas in chain an write the threelettercodes in a string and count the aas
      for
      (
        biol::AASequence::const_iterator
          aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        string_sequence << ( *aa_itr)->GetType()->GetThreeLetterCode() << ' ';
        number_aas++;
      }

      //string to read in the threelettercodes for the string_sequence
      std::string threelettercode;
      //the seqreslines that will contain the complete seqrespart for the pdb in Lines
      util::ShPtrList< Line> seqreslines;
      //iterate over all seqserials each containing 13 AAs in the end
      for( size_t seqserial( 1); seqserial <= ( number_aas / 13) + 1; ++seqserial)
      {
        //new SEQRES line
        util::ShPtr< Line> seqresline( new Line( GetLineTypes().SEQRES));
        //write seqserial to Line
        seqresline->Put( GetEntryTypes().SEQRESSerial, seqserial);
        //write chainid to seqresline
        seqresline->Put( GetEntryTypes().SEQRESChainID, chainid);
        //write number of residues in complete chain
        seqresline->Put( GetEntryTypes().SEQRESNrOfResiduesInChain, number_aas);
        //write 13 aas in each line
        for( size_t seqres_pos( 0); seqres_pos < 13 && seqres_pos < ( number_aas - ( seqserial - 1) * 13); ++seqres_pos)
        {
          //read the threelettercode
          string_sequence >> threelettercode;
          //write the current three letter code to the Line
          seqresline->Put( EntryType( GetEntryTypes().SEQRESName_1 + seqres_pos), threelettercode);
        }
        //push the line to the lines
        seqreslines.PushBack( seqresline);
      }

      return seqreslines;
    }

    //! @brief write SEQRES information of a Chain into Lines
    //! @param CHAIN Chain of interest
    //! @return ShPtrList of Lines containing SEQRES information for given Chain
    util::ShPtrList< Line> Factory::WriteSeqResToLines( const assemble::Chain &CHAIN)
    {
      //if there is a sequence given with the model write the seqresinformation from the complete sequence
      if( CHAIN.GetSequence().IsDefined())
      {
        return WriteSeqResToLines( *CHAIN.GetSequence());
      }

      //if there are no sses in the Proteinmodel return an empty Seqreslines vector
      if( CHAIN.GetData().IsEmpty())
      {
        return util::ShPtrList< Line>();
      }

      //store the chainid
      const char chainid( CHAIN.GetChainID());

      size_t number_aas( 0);
      std::stringstream string_sequence;

      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( CHAIN.GetData().Begin()), sse_itr_end( CHAIN.GetData().End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        for
        (
          biol::AASequence::const_iterator
            aa_itr( ( *sse_itr)->Begin()), aa_itr_end( ( *sse_itr)->End());
          aa_itr != aa_itr_end; ++aa_itr)
        {
          number_aas++;
          while( int( number_aas) < ( *aa_itr)->GetSeqID())
          {
            string_sequence << biol::GetAATypes().e_Undefined->GetThreeLetterCode() << ' ';
            number_aas++;
          }
          string_sequence << ( *aa_itr)->GetType()->GetThreeLetterCode() << ' ';
        }
      }

      std::string threelettercode;
      util::ShPtrList< Line> seqreslines;
      for( size_t seqserial( 1); seqserial <= ( number_aas / 13) + 1; ++seqserial)
      {
        util::ShPtr< Line> seqresline( new Line( GetLineTypes().SEQRES));
        seqresline->Put( GetEntryTypes().SEQRESSerial, seqserial);
        seqresline->Put( GetEntryTypes().SEQRESChainID, chainid);
        seqresline->Put( GetEntryTypes().SEQRESNrOfResiduesInChain, number_aas);
        for( size_t seqres_pos( 0); seqres_pos < 13 && seqres_pos < ( number_aas - ( seqserial - 1) * 13); ++seqres_pos)
        {
          string_sequence >> threelettercode;
          seqresline->Put( EntryType( GetEntryTypes().SEQRESName_1 + seqres_pos), threelettercode);
        }
        seqreslines.PushBack( seqresline);
      }

      //return
      return seqreslines;
    }

    //! @brief write SEQRES information of a ProteinModel into Lines
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return ShPtrList of Lines containing SEQRES information for given ProteinModel
    util::ShPtrList< Line> Factory::WriteSeqResToLines( const assemble::ProteinModel &PROTEIN_MODEL)
    {
      //accumulate seqreslines for each chain
      util::ShPtrList< Line> seqreslines;

      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr)
      {
        seqreslines.Append( WriteSeqResToLines( **chain_itr));
      }

      //return
      return seqreslines;
    }

    //! @brief write a single Protein Model to a pdb file stream
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param OSTREAM output stream to which ProteinModel will be written to
    //! @param WRITE_BODY_INFORMATION flag to determine whether to write body information
    //! @param CLASSIFICATION string for header classification
    //! @param TIME time (date) information to be placed in the header
    //! @return std::ostream to which ProteinModel was written to
    std::ostream &Factory::WriteModelToPDB
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM,
      const bool WRITE_BODY_INFORMATION,
      const std::string &CLASSIFICATION,
      const util::Time &TIME
    ) const
    {
      // create handler
      Handler newpdb;

      // add lines
      newpdb.AppendLines
      (
        WriteCompleteModelToPDBLines( PROTEIN_MODEL, WRITE_BODY_INFORMATION, CLASSIFICATION, TIME)
      );

      // write lines
      return newpdb.WriteLines( OSTREAM);
    }

    //! @brief write a Chain to a pdb file stream
    //! @param CHAIN Chain of interest
    //! @param OSTREAM output stream to which Chain will be written to
    //! @return std::ostream to which Chain was written to
    std::ostream &Factory::WriteChainToPDB( const assemble::Chain &CHAIN, std::ostream &OSTREAM)
    {
      Handler newpdb;
      newpdb.PushBack
      (
        util::ShPtr< Line>
        (
          new Line( std::string( "HEADER  BCL MODEL"))// + util::Time::GetCurrent().GetTimeAsDate()))
        )
      );

      // the SEQRES section
      newpdb.AppendLines( WriteSeqResToLines( CHAIN));

      size_t serial( 1);
      // the HELIX section
      newpdb.AppendLines( WriteHelixDefinitionsToLines( CHAIN.GetSSEs( biol::GetSSTypes().GetHelixTypes()), serial));

      serial = 1;
      // the STRAND section
      newpdb.AppendLines( WriteStrandDefinitionsToLines( CHAIN.GetSSEs( biol::GetSSTypes().STRAND), serial));

      serial = 1;
      // write the residues and atoms to pdblines
      newpdb.AppendLines( WriteChainToLines( CHAIN, serial));

      return newpdb.WriteLines( OSTREAM);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Factory from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Factory::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AAClass, ISTREAM);
      io::Serialize::Read( m_Printers, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write Factory to given OSTREAM
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return std::ostream to which was written
    std::ostream &Factory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write members
      io::Serialize::Write( m_AAClass, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Printers, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief construct the first side chain atom from given atoms
    //! those atoms need to contain the CA, N and O
    //! @param ATOMS CA, N and C atom
    //! @param BOND_LENGTH the distance of the constructed side chain atom coordinate to the CA atom
    //! @return position of first side chain atom for L-amino acid
    linal::Vector3D Factory::FirstSidechainAtomCoordinateFromBackboneAtoms
    (
      const util::SiPtrVector< const biol::Atom> &ATOMS,
      const double BOND_LENGTH
    )
    {
      // initialize pointers to atoms required for building a pseudo HA2
      // CA
      const util::SiPtr< const biol::Atom> ca_ptr( biol::Atom::FindAtom( ATOMS, biol::GetAtomTypes().CA));
      if( !ca_ptr->GetCoordinates().IsDefined())
      {
        return linal::Vector3D( util::GetUndefined< double>());
      }

      // N
      const util::SiPtr< const biol::Atom> n_ptr( biol::Atom::FindAtom( ATOMS, biol::GetAtomTypes().N));
      if( !n_ptr->GetCoordinates().IsDefined())
      {
        return linal::Vector3D( util::GetUndefined< double>());
      }

      // C
      const util::SiPtr< const biol::Atom> c_ptr( biol::Atom::FindAtom( ATOMS, biol::GetAtomTypes().C));
      if( !c_ptr->GetCoordinates().IsDefined())
      {
        return linal::Vector3D( util::GetUndefined< double>());
      }

      // angle for tetrahedral geometry
      static const double s_tetrahedral_angle( 109.5 / 180.0 * math::g_Pi);

      // construct atom coordinate
      return linal::CoordinatesAngle
             (
               ca_ptr->GetCoordinates(),
               n_ptr->GetCoordinates(),
               c_ptr->GetCoordinates(),
               BOND_LENGTH, s_tetrahedral_angle, s_tetrahedral_angle
             );
    }

    //! @brief returns the current Date using the PDB format
    //! @param TIME time (date) information to be placed in the header
    //! @return string of the current Date using the PDB format
    std::string Factory::ConvertTimeToPDBDate( const util::Time &TIME)
    {
      // split the time as date string by " "
      storage::Vector< std::string> split_date( util::SplitString( TIME.GetTimeAsDate(), " "));

      // make the month uppercase
      for( size_t i( 0); i != 3; ++i)
      {
        split_date( 1)[ i] = toupper( split_date( 1)[ i]);
      }

      // return the formatted date
      return split_date( 2) + "-" + split_date( 1) + "-" + split_date( 4)[ 2] + split_date( 4)[ 3];
    }

    //! @brief builds model of with one or more chains, each formed of full sequence and SSEs from given Handler
    //! The types and sizes of SSEs to be considered are determined by the additional parameters
    //! @param HANDLER Handler that contains pdb information
    //! @param SSE_TYPE_MINSIZE map of considered sstypes with their minimal size
    //! @param MODEL model in the handler to build
    //! @return Protein Model read from given Handler
    assemble::ProteinModel Factory::ProcessModel
    (
      const Handler &HANDLER,
      const storage::Map< biol::SSType, size_t> &SSE_TYPE_MINSIZE,
      const Model &MODEL
    ) const
    {
      // model
      assemble::ProteinModel protein_model;

      // structured chain
      const storage::Map< char, storage::List< Residue> > structured_protein_chains( MODEL.GetStructuredChains());

      // secondary structure element definitions
      const storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
        sse_definitions( HANDLER.GetSSEStructure());

      // get iterator on pdb ssestructure
      for
      (
        storage::Map< char, storage::List< Residue> >::const_iterator
          chain_itr( structured_protein_chains.Begin()), chain_itr_end( structured_protein_chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        util::ShPtr< assemble::Chain> sp_chain;

        const storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >::const_iterator
          sse_itr( sse_definitions.Find( chain_itr->first));
        if( sse_itr == sse_definitions.End())
        {
          BCL_MessageVrb( "cannot find sse definitions for chain: " + util::Format()( chain_itr->first));
          //create each chain in pdb
          sp_chain = util::ShPtr< assemble::Chain>
          (
            ChainFromPDBSSElementsAndSequence
            (
              chain_itr->first, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> >(),
              chain_itr->second,
              SSE_TYPE_MINSIZE,
              HANDLER.GetHead().GetPDBID()
            )
          );
        }
        else
        {
          //create each chain in pdb
          sp_chain = util::ShPtr< assemble::Chain>
          (
            ChainFromPDBSSElementsAndSequence
            (
              chain_itr->first, sse_itr->second,
              chain_itr->second,
              SSE_TYPE_MINSIZE,
              HANDLER.GetHead().GetPDBID()
            )
          );
        }

        // insert chain into model
        protein_model.Insert( sp_chain);
      }

      // try to get the membrane
      const util::ShPtr< biol::Membrane> sp_membrane( HANDLER.GetHead().GetMembrane());

      // if a membrane is defined
      if( sp_membrane.IsDefined())
      {
        // add to protein model data
        util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
        sp_data->Insert( assemble::ProteinModelData::e_Membrane, sp_membrane);
        protein_model.SetProteinModelData( sp_data);

        // set pdb environments
        sspred::PDB::SetEnvironmentTypes( protein_model);
      }

      // if the biomolecule flag was set
      if( GetFlagBiomolecule()->GetFlag())
      {
        // get any bio transformation matrices
        const storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > bio_matrices
        (
          HANDLER.GetHead().GetBioTransformationMatrices( HANDLER.GetProteinChains())
          [
            GetFlagBiomolecule()->GetFirstParameter()->GetNumericalValue< size_t>()
          ]
        );

        // if matrices were found
        if( !bio_matrices.IsEmpty())
        {
          // construct the protein model multiplier from the first biomolecule
          util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
          (
            new assemble::ProteinModelMultiplier( bio_matrices, protein_model)
          );

          // add to protein model data
          util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
          sp_data->Insert( assemble::ProteinModelData::e_Multiplier, sp_multiplier);
          protein_model.SetProteinModelData( sp_data);
        }
        else
        {
          BCL_MessageCrt
          (
            "no biomatrices in REMARK 350 in pdb given, assuming that all atoms in pdb belong to biological molecule"
          );
        }
      }

      //return model
      return protein_model;
    }

  } // namespace pdb
} // namespace bcl
