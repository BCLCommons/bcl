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
#include "sspred/bcl_sspred_palsse.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Palsse::Palsse() :
      m_SSType()
    {
    }

    //! @brief constructor from an sse type
    //! @param SS_TYPE type of sse
    Palsse::Palsse( const biol::SSType &SS_TYPE) :
      m_SSType( SS_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Palsse
    Palsse *Palsse::Clone() const
    {
      return new Palsse( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Palsse::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &Palsse::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".palsse");

      // end
      return s_file_extension;

    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D Palsse::GetThreeStatePrediction() const
    {
      return m_SSType->GetThreeStatePrediction();
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> Palsse::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( GetThreeStatePrediction(), biol::GetEnvironmentTypes().e_Solution);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Palsse::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      // issue warning
      BCL_MessageCrt( "Palsse ReadPredictionsForAA should not have been called");

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Palsse::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence( GetMethods().e_PALSSE, AA_SEQUENCE, **GetMethods().e_PALSSE);

      // initialize necessary variables
      std::string line;
      const char sequence_chain_id( AA_SEQUENCE.GetChainID());

      // track all chain ids seen in dssp that are not == the sequence chain id
      std::string palsse_chain_ids;

      // get the minimum pdb id in the sequence
      const int min_pdb_id( AA_SEQUENCE.GetFirstAA()->GetPdbID());
      const int max_pdb_id( AA_SEQUENCE.GetLastAA()->GetPdbID());

      if( max_pdb_id - min_pdb_id < 4)
      {
        BCL_MessageStd
        (
          "Too short of an AA_Sequence to have any secondary structure, skipping for chain "
          + util::Format()( AA_SEQUENCE.GetChainID())
        );
        return ISTREAM;
      }

      // create a vector, index is pdb id to ss type
      storage::Vector< biol::SSType> aa_ss_types( max_pdb_id - min_pdb_id + 1, biol::GetSSTypes().COIL);

      // track the number of ss types found for the chain
      size_t number_aas_found( 0);

      // chain ID position for helices and strands
      const size_t chain_id_pos_helix( 19);
      const size_t chain_id_pos_strand( 21);
      const size_t chain_id_epos_helix( 31);
      const size_t chain_id_epos_strand( 32);

      // start and end locations for pdb start and end id
      const size_t pdbid_startpos_helix( chain_id_pos_helix + 1);
      const size_t pdbid_startpos_strand( chain_id_pos_strand + 1);
      const size_t pdbid_endpos_helix( chain_id_epos_helix + 1);
      const size_t pdbid_endpos_strand( chain_id_epos_strand + 1);
      const size_t pdbid_size_strand( 4);
      const size_t pdbid_size_helix( 5);

      BCL_Assert( ISTREAM.good(), "PALSSE file did not exist!");

      // read through the file
      while( ISTREAM.good())
      {
        std::getline( ISTREAM, line);

        // does the line contain SS summary information?
        const bool is_helix_line( util::StartsWith( line, "HELIX "));
        const bool is_strand_line( util::StartsWith( line, "SHEET "));
        if( !is_helix_line && !is_strand_line)
        {
          // no, continue
          continue;
        }

        const char chain_id_str( line[ is_helix_line ? chain_id_pos_helix : chain_id_pos_strand]);
        const char chain_id_estr( line[ is_helix_line ? chain_id_epos_helix : chain_id_epos_strand]);

        // check that the chain IDs are identical
        if( chain_id_str != chain_id_estr)
        {
          BCL_MessageCrt( "SS's should not span chains in palsse file: " + line);
          continue;
        }

        // get the chain id.  Stride automatically changes blank chain ids to - for easier parsing, so map them back
        const char chain_id( chain_id_str == '-' ? ' ' : chain_id_str);

        // skip undesired chain ids
        if( chain_id != sequence_chain_id)
        {
          if( palsse_chain_ids.find( chain_id) == std::string::npos)
          {
            // new chain id
            palsse_chain_ids += chain_id;
          }
          continue;
        }

        const std::string pdb_id_start_str
        (
          is_helix_line
          ? line.substr( pdbid_startpos_helix, pdbid_size_helix)
          : line.substr( pdbid_startpos_strand, pdbid_size_strand)
        );
        const std::string pdb_id_end_str
        (
          is_helix_line
          ? line.substr( pdbid_endpos_helix, pdbid_size_helix)
          : line.substr( pdbid_endpos_strand, pdbid_size_strand)
        );

        // check that the pdb ids are numeric
        int pdb_start( util::GetUndefined< int>()), pdb_end( util::GetUndefined< int>());
        if( !util::TryConvertFromString( pdb_start, pdb_id_start_str, util::GetLogger()))
        {
          BCL_MessageCrt( "Non-numeric PDB start id in " + line);
          continue;
        }
        if( !util::TryConvertFromString( pdb_end, pdb_id_end_str, util::GetLogger()))
        {
          BCL_MessageCrt( "Non-numeric PDB start id in " + line);
          continue;
        }

        // determine SS-type
        biol::SSType ss_type( is_helix_line ? biol::GetSSTypes().HELIX : biol::GetSSTypes().STRAND);

        const size_t ss_size( pdb_end - pdb_start + 1);
        number_aas_found += ss_size;

        // set the types
        // get the first pdb id that is in the AA_Sequence
        int pdb_id( std::max( pdb_start, min_pdb_id));

        if( pdb_id > min_pdb_id && aa_ss_types( pdb_id - min_pdb_id - 1) != biol::GetSSTypes().COIL)
        {
          // adjacent SSE types (not separated by any residues)
          // Reassign the starting residue of this sse and the ending residue of the last sse to be coil
          aa_ss_types( pdb_id - min_pdb_id - 1) = biol::GetSSTypes().COIL;
          ++pdb_id;
        }
        const int pdb_end_id( std::min( pdb_end, max_pdb_id));
        for( ; pdb_id <= pdb_end_id; ++pdb_id)
        {
          if( aa_ss_types( pdb_id - min_pdb_id) != biol::GetSSTypes().COIL)
          {
            // make AAs that connect two SSEs of coil type; necessary to have reasonable elements for folding
            aa_ss_types( pdb_id - min_pdb_id) = biol::GetSSTypes().COIL;
          }
          else
          {
            aa_ss_types( pdb_id - min_pdb_id) = ss_type;
          }
        }
        // check whether the ending residues are both parts of sses not separated by coils
        if( pdb_end_id < max_pdb_id && aa_ss_types( pdb_end_id - min_pdb_id + 1) != biol::GetSSTypes().COIL)
        {
          // adjacent SSE types (not separated by any residues)
          // Reassign the ending residue of this sse and the starting residue of the next sse to be coil
          // to avoid confusing it for one long sse, and to maintain the relative sse sizes
          aa_ss_types( pdb_end_id - min_pdb_id + 1) = biol::GetSSTypes().COIL;
          aa_ss_types( pdb_end_id - min_pdb_id) = biol::GetSSTypes().COIL;
        }
      }

      if( !number_aas_found && palsse_chain_ids.empty())
      {
        BCL_MessageStd
        (
          "PALSSE file empty; assuming that the protein with size: "
          + util::Format()( AA_SEQUENCE.GetSize()) + " is unstructured"
        );
      }
      else if( !number_aas_found)
      {
        BCL_MessageStd
        (
          "PALSSE file did not contain a chain " + util::Format()( sequence_chain_id)
          + " but did for " + palsse_chain_ids
          + "; assuming that the chain with size "
          + util::Format()( AA_SEQUENCE.GetSize()) + " is unstructured"
        );
      }

      // set the ss type prediction for each AA
      for
      (
        biol::AASequence::iterator itr_aa( AA_SEQUENCE.Begin()), itr_aa_end( AA_SEQUENCE.End());
        itr_aa != itr_aa_end;
        ++itr_aa
      )
      {
        if( !util::IsDefined( ( *itr_aa)->GetPdbID()))
        {
          BCL_MessageStd
          (
            "Undefined PDB id for " + AA_SEQUENCE.GetSequenceId()
            + " at position " + util::Format()( std::distance( AA_SEQUENCE.Begin(), itr_aa))
            + " it is very likely that pdb::Model::s_NumberRequiredAlignedResidues needs to be increased!"
          );
          // no pdb id; likely contains residues that lacked coordinates, set to coil
        }
        else if( ( *itr_aa)->GetPdbID() < min_pdb_id || ( *itr_aa)->GetPdbID() > max_pdb_id)
        {
          BCL_MessageStd
          (
            "Major issue with numbering: " + util::Format()( min_pdb_id) + "-" + util::Format()( max_pdb_id)
            + " vs " + util::Format()( ( *itr_aa)->GetPdbID())
            + " at distance " + util::Format()( std::distance( AA_SEQUENCE.Begin(), itr_aa))
          );
          ( *itr_aa)->SetSSPrediction( GetMethods().e_PALSSE, Palsse( biol::GetSSTypes().COIL));
        }
        else
        {
          // get the type out of the vector
          ( *itr_aa)->SetSSPrediction( GetMethods().e_PALSSE, Palsse( aa_ss_types( ( *itr_aa)->GetPdbID() - min_pdb_id)));
        }
      }

      // end
      return ISTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Palsse::Read( std::istream &ISTREAM)
    {
      // read members
      ISTREAM >> m_SSType;

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Palsse::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_SSType;

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
