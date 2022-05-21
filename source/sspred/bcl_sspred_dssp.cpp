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
#include "sspred/bcl_sspred_dssp.h"

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
    Dssp::Dssp() :
      m_SSType()
    {
    }

    //! @brief constructor from an sse type
    //! @param SS_TYPE type of sse
    Dssp::Dssp( const biol::SSType &SS_TYPE) :
      m_SSType( SS_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Dssp
    Dssp *Dssp::Clone() const
    {
      return new Dssp( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Dssp::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &Dssp::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".dssp");

      // end
      return s_file_extension;

    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D Dssp::GetThreeStatePrediction() const
    {
      return m_SSType->GetThreeStatePrediction();
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> Dssp::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( GetThreeStatePrediction(), biol::GetEnvironmentTypes().e_Solution);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Dssp::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      // issue warning
      BCL_MessageCrt( "Dssp ReadPredictionsForAA should not have been called");

      // end
      return ISTREAM;
    }

    //! @brief helper function that builds the map from Dssp type to biol type
    storage::Map< char, biol::SSType> GetDsspTypeToBiolSSTypeMap()
    {
      // this is a complete listing of all the types returned by Dssp.
      storage::Map< char, biol::SSType> mapping;
      mapping[ 'H'] =  biol::GetSSTypes().HELIX;
      mapping[ 'G'] =  biol::GetSSTypes().COIL;
      mapping[ 'I'] =  biol::GetSSTypes().COIL;
      mapping[ 'B'] =  biol::GetSSTypes().STRAND; // beta bridge = 1 residue strand; can be set to coil by changing the min strand length in the descriptor
      mapping[ 'E'] =  biol::GetSSTypes().STRAND; // conventional strand
      mapping[ 'T'] =  biol::GetSSTypes().COIL;   // Turn ~ Coil
      mapping[ 'S'] =  biol::GetSSTypes().COIL;   // Curved ~ Coil
      mapping[ ' '] =  biol::GetSSTypes().COIL;   // N/A ~ Coil
      return mapping;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Dssp::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // static map from dssp type to biol::SSType
      static const storage::Map< char, biol::SSType> s_dssp_to_biol_ss_type( GetDsspTypeToBiolSSTypeMap());

      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence( GetMethods().e_DSSP, AA_SEQUENCE, **GetMethods().e_DSSP);

      // initialize necessary variables
      std::string line;
      const char sequence_chain_id( AA_SEQUENCE.GetChainID());

      // track all chain ids seen in dssp that are not == the sequence chain id
      std::string dssp_chain_ids;

      // get the minimum pdb id in the sequence
      const int min_pdb_id( AA_SEQUENCE.GetFirstAA()->GetPdbID());
      const int max_pdb_id( AA_SEQUENCE.GetLastAA()->GetPdbID());

      // create a vector, index is pdb id to ss type
      storage::Vector< biol::SSType> aa_ss_types( max_pdb_id - min_pdb_id + 1, biol::GetSSTypes().COIL);

      // track the number of ss types found for the chain
      size_t number_aas_found( 0);

      // read through the file until the residue lines are reached
      while( ISTREAM.good())
      {
        std::getline( ISTREAM, line);

        // does the line contain SS summary information?
        if( util::StartsWith( line, "  #  RESIDUE"))
        {
          // Reached residue lines
          break;
        }
      }

      // read through residue lines
      while( ISTREAM.good())
      {
        std::getline( ISTREAM, line);

        // skip lines that are too short
        if( line.size() < 17)
        {
          continue;
        }

        // the split string contains items in the following order
        // characters 0-4: Residue ID
        // characters 5-9: PDB_ID
        // character   11: Chain ID
        // character   13: 1-letter AA code
        // character   16: 1-letter SS type
        // followed by a bunch of irrelevant (for our purposes) information

        const char ss_type_chr( line[ 16]);
        const std::string pdb_id_str( line.substr( 5, 5));

        if( pdb_id_str[ 4] == ' ' && pdb_id_str[ 3] == ' ')
        {
          // empty pdb id, happens with gaps, continue
          continue;
        }

        // get the chain id.  Dssp automatically changes blank chain ids to - for easier parsing, so map them back
        const char chain_id( line[ 11] == '-' ? ' ' : line[ 11]);

        // skip undesired chain ids
        if( chain_id != sequence_chain_id)
        {
          if( dssp_chain_ids.find( chain_id) == std::string::npos)
          {
            // new chain id
            dssp_chain_ids += chain_id;
          }
          continue;
        }

        // check that the pdb ids are numeric
        int pdb_id( util::GetUndefined< int>());
        if( !util::TryConvertFromString( pdb_id, pdb_id_str, util::GetLogger()))
        {
          BCL_MessageCrt( "Non-numeric PDB id in " + line);
          continue;
        }
        ++number_aas_found;

        // determine SS-type
        storage::Map< char, biol::SSType>::const_iterator itr( s_dssp_to_biol_ss_type.Find( ss_type_chr));
        biol::SSType ss_type;
        if( itr == s_dssp_to_biol_ss_type.End())
        {
          ss_type = biol::GetSSTypes().COIL;
        }
        else
        {
          ss_type = itr->second;
        }
        if( ss_type == biol::GetSSTypes().COIL)
        {
          // coil type is set by default, so continue
          continue;
        }
        if( pdb_id >= min_pdb_id && size_t( pdb_id - min_pdb_id) < aa_ss_types.GetSize())
        {
          aa_ss_types( pdb_id - min_pdb_id) = ss_type;
        }
      }

      // check whether any predictions were found for the chain id
      BCL_Assert
      (
        number_aas_found || AA_SEQUENCE.CountDefinedAACoordinates() < size_t( 6),
        "Dssp file contained no SS information for chain " + util::Format()( sequence_chain_id)
        + ". Chains present: " + dssp_chain_ids + " Rerun dssp and try again"
      );

      if( !number_aas_found && dssp_chain_ids.empty())
      {
        BCL_MessageStd
        (
          "DSSP file empty; assuming that the protein with size: "
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
          ( *itr_aa)->SetSSPrediction( GetMethods().e_DSSP, Dssp( biol::GetSSTypes().COIL));
        }
        else if( ( *itr_aa)->GetPdbID() < min_pdb_id || ( *itr_aa)->GetPdbID() > max_pdb_id)
        {
          BCL_MessageStd
          (
            "Major issue with numbering: " + util::Format()( min_pdb_id) + "-" + util::Format()( max_pdb_id)
            + " vs " + util::Format()( ( *itr_aa)->GetPdbID())
            + " at distance " + util::Format()( std::distance( AA_SEQUENCE.Begin(), itr_aa))
          );
          ( *itr_aa)->SetSSPrediction( GetMethods().e_DSSP, Dssp( biol::GetSSTypes().COIL));
        }
        else
        {
          // get the type out of the vector
          ( *itr_aa)->SetSSPrediction( GetMethods().e_DSSP, Dssp( aa_ss_types( ( *itr_aa)->GetPdbID() - min_pdb_id)));
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
    std::istream &Dssp::Read( std::istream &ISTREAM)
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
    std::ostream &Dssp::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_SSType;

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
