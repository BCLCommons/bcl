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
#include "sspred/bcl_sspred_stride.h"

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
    Stride::Stride() :
      m_SSType()
    {
    }

    //! @brief constructor from an sse type
    //! @param SS_TYPE type of sse
    Stride::Stride( const biol::SSType &SS_TYPE) :
      m_SSType( SS_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Stride
    Stride *Stride::Clone() const
    {
      return new Stride( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Stride::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &Stride::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".stride");

      // end
      return s_file_extension;

    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D Stride::GetThreeStatePrediction() const
    {
      return m_SSType->GetThreeStatePrediction();
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> Stride::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( GetThreeStatePrediction(), biol::GetEnvironmentTypes().e_Solution);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Stride::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      // issue warning
      BCL_MessageCrt( "Stride ReadPredictionsForAA should not have been called");

      // end
      return ISTREAM;
    }

    //! @brief helper function that builds the map from Stride type to biol type
    storage::Map< std::string, biol::SSType> GetStrideTypeToBiolSSTypeMap()
    {
      // this is a complete listing of all the types returned by Stride.
      storage::Map< std::string, biol::SSType> mapping;
      mapping[ "AlphaHelix"  ] =  biol::GetSSTypes().e_HelixLeftAlpha;
      mapping[ "310Helix"    ] =  biol::GetSSTypes().e_HelixRight310;
      mapping[ "PiHelix"     ] =  biol::GetSSTypes().e_HelixRightPi;
      mapping[ "Strand"      ] =  biol::GetSSTypes().STRAND;
      mapping[ "Bridge"      ] =  biol::GetSSTypes().COIL; // essentially a 1-residue strand
      mapping[ "Coil"        ] =  biol::GetSSTypes().COIL; // Region without any hydrogen bonding, highly unstructured

      // Various turn types with well-defined phi-psi angles
      // see http://en.wikipedia.org/wiki/Turn_%28biochemistry%29
      // TODO: Adding these as enums; it may be useful to train a secondary structure algorithm that could
      // recongize standard turn types, which would be useful for loop building
      mapping[ "TurnI"       ] =  biol::GetSSTypes().COIL;
      mapping[ "TurnI'"      ] =  biol::GetSSTypes().COIL;
      mapping[ "TurnII"      ] =  biol::GetSSTypes().COIL;
      mapping[ "TurnII'"     ] =  biol::GetSSTypes().COIL;
      mapping[ "TurnVIa"     ] =  biol::GetSSTypes().COIL;
      mapping[ "TurnVIb"     ] =  biol::GetSSTypes().COIL;
      mapping[ "TurnVIII"    ] =  biol::GetSSTypes().COIL;
      mapping[ "TurnIV"      ] =  biol::GetSSTypes().COIL;
      mapping[ "GammaClassic"] =  biol::GetSSTypes().COIL;
      mapping[ "GammaInv"    ] =  biol::GetSSTypes().COIL;
      mapping[ "Turn"        ] =  biol::GetSSTypes().COIL;
      mapping[ "Unknown"     ] =  biol::GetSSTypes().COIL; // Probably never comes up, but its in the source code
      return mapping;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Stride::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // static map from stride type to biol::SSType
      static const storage::Map< std::string, biol::SSType> s_stride_to_biol_ss_type( GetStrideTypeToBiolSSTypeMap());

      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence( GetMethods().e_Stride, AA_SEQUENCE, **GetMethods().e_Stride);

      // initialize necessary variables
      std::string line;
      const char sequence_chain_id( AA_SEQUENCE.GetChainID());

      // track all chain ids seen in stride that are not == the sequence chain id
      std::string stride_chain_ids;

      // get the minimum pdb id in the sequence
      const int min_pdb_id( AA_SEQUENCE.GetFirstAA()->GetPdbID());
      const int max_pdb_id( AA_SEQUENCE.GetLastAA()->GetPdbID());

      // create a vector, index is pdb id to ss type
      storage::Vector< biol::SSType> aa_ss_types( max_pdb_id - min_pdb_id + 1, biol::GetSSTypes().COIL);

      // track the number of ss types found for the chain
      size_t number_aas_found( 0);

      // read through the file
      while( ISTREAM.good())
      {
        std::getline( ISTREAM, line);

        // does the line contain SS summary information?
        if( !util::StartsWith( line, "LOC "))
        {
          if( !number_aas_found && line.size() > 12 && util::StartsWith( line, "ASG "))
          {
            if( line[ 9] == sequence_chain_id)
            {
              // valid line for chain; just happens to be a coil
              number_aas_found += 1;
            }
          }
          // no, continue
          continue;
        }

        // split the line
        storage::Vector< std::string> split_string( util::SplitString( line, " "));

        // the split string contains items in the following order
        // item 0: LOC[[
        // item 1: SS-type
        // item 2: 3-letter code start residue
        // item 3: PDB-id start residue
        // item 4: chain-id start residue
        // item 5: 3-letter code end residue
        // item 6: PDB-id end residue
        // item 7: chain-id end residue (should always be == item 4)
        // item 8: PDB ID (if PDB ID is known)
        if( split_string.GetSize() != 9 && split_string.GetSize() != 8)
        {
          BCL_MessageCrt( "Bad line in stride file: " + line);
          continue;
        }

        const std::string &ss_type_str( split_string( 1));
        const std::string &chain_id_str( split_string( 4));
        const std::string &pdb_id_start_str( split_string( 3));
        const std::string &pdb_id_end_str( split_string( 6));

        // check that the chain IDs are identical
        if( chain_id_str != split_string( 7))
        {
          BCL_MessageCrt( "SS's should not span chains in stride file: " + line);
          continue;
        }

        // check that the chain id is a single character.  Spaces are not allowed for chain id
        if( chain_id_str.size() != size_t( 1))
        {
          BCL_MessageCrt( "Multi-character chain-ids not allowed: " + line);
          continue;
        }

        // get the chain id.  Stride automatically changes blank chain ids to - for easier parsing, so map them back
        const char chain_id( chain_id_str[ 0] == '-' ? ' ' : chain_id_str[ 0]);

        // skip undesired chain ids
        if( chain_id != sequence_chain_id)
        {
          if( stride_chain_ids.find( chain_id) == std::string::npos)
          {
            // new chain id
            stride_chain_ids += chain_id;
          }
          continue;
        }

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
        storage::Map< std::string, biol::SSType>::const_iterator itr( s_stride_to_biol_ss_type.Find( ss_type_str));
        biol::SSType ss_type;
        if( itr == s_stride_to_biol_ss_type.End())
        {
          ss_type = biol::GetSSTypes().COIL;
        }
        else
        {
          ss_type = itr->second;
        }
        const size_t ss_size( pdb_end - pdb_start + 1);
        number_aas_found += ss_size;
        if( ss_type == biol::GetSSTypes().COIL)
        {
          // coil type is set by default, so continue
          continue;
        }
        // ignore tiny strands
        if( ss_type == biol::GetSSTypes().STRAND)
        {
          if( ss_size == 0)
          {
            continue;
          }
        }
        // ignore tiny helices
        else
        {
          if( ss_size < 2)
          {
            continue;
          }
        }

        // set the types
        // get the first pdb id that is in the AA_Sequence
        for
        (
          int pdb_id( std::max( pdb_start, min_pdb_id)), pdb_end_id( std::min( pdb_end, max_pdb_id));
          pdb_id <= pdb_end_id;
          ++pdb_id
        )
        {
          aa_ss_types( pdb_id - min_pdb_id) = ss_type;
        }
      }

      size_t number_defined_aatypes( 0);
      if( !number_aas_found)
      {
        // iterate over the residues
        for
        (
          biol::AASequence::const_iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // if there are undefined coords
          if( ( *aa_itr)->GetType()->GetParentType().GetIndex() < size_t( biol::AATypes::s_NumberStandardAATypes))
          {
            ++number_defined_aatypes;
          }
        }
      }

      // check whether any predictions were found for the chain id
      // For proteins of size < 7, it is possible for the stride file to be empty because it finds no hydrogen bonding
      BCL_Assert
      (
        number_aas_found || AA_SEQUENCE.CountDefinedAACoordinates() < size_t( 6) || number_defined_aatypes < size_t( 6),
        "Stride file contained no SS information for chain " + util::Format()( sequence_chain_id)
        + ". Chains present: " + stride_chain_ids + " Rerun stride and try again"
      );

      if( !number_aas_found && stride_chain_ids.empty())
      {
        BCL_MessageStd
        (
          "Stride file empty; assuming that the protein with size: "
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
          ( *itr_aa)->SetSSPrediction( GetMethods().e_Stride, Stride( biol::GetSSTypes().COIL));
        }
        else if( ( *itr_aa)->GetPdbID() < min_pdb_id || ( *itr_aa)->GetPdbID() > max_pdb_id)
        {
          BCL_MessageStd
          (
            "Major issue with numbering: " + util::Format()( min_pdb_id) + "-" + util::Format()( max_pdb_id)
            + " vs " + util::Format()( ( *itr_aa)->GetPdbID())
            + " at distance " + util::Format()( std::distance( AA_SEQUENCE.Begin(), itr_aa))
          );
          ( *itr_aa)->SetSSPrediction( GetMethods().e_Stride, Stride( biol::GetSSTypes().COIL));
        }
        else
        {
          // get the type out of the vector
          ( *itr_aa)->SetSSPrediction( GetMethods().e_Stride, Stride( aa_ss_types( ( *itr_aa)->GetPdbID() - min_pdb_id)));
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
    std::istream &Stride::Read( std::istream &ISTREAM)
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
    std::ostream &Stride::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_SSType;

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
