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
#include "contact/bcl_contact_prediction_map.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "contact/bcl_contact_ann.h"
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////
  // data //
  //////////

    //! How many residues to exclude from begin and end
    static const size_t s_NumberSkippedResiduesInEnds = size_t( 4);

    //! The minimum number of residues a sequence should have contact prediction
    static const size_t s_MinimumNumberResidues = size_t( 20);

    //! minimal sequence distance for contacts
    static const size_t s_MinSequenceDistance = size_t( 4);

    //! The number of values used for describing each residue for prediction (7 properties + 3 JUFO + 20 Blast)
    static const size_t s_NumberOfValuesPerResidue = size_t( 30);

    //! identifier string to be used for prediction maps
    static const std::string s_Identifier = "CHAINS";

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PredictionMap::PredictionMap() :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< double>, double> >(),
      m_Sequences()
    {
    }

    //! @brief construct prediction map ( from ANNs) from one chain( intra-sequence contacts)
    //! @param CHAIN chain
    PredictionMap::PredictionMap( const assemble::Chain &CHAIN) :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< double>, double> >(),
      m_Sequences( 1, CHAIN.GetSequence())
    {
      // make sure the sequence matches the size
      BCL_Assert
      (
        CHAIN.GetSequence()->GetSize() >= s_MinimumNumberResidues,
        "The given sequence needs to be at least " + util::Format()( s_MinimumNumberResidues) +
        " residues long, but it has only " + util::Format()( CHAIN.GetSequence()->GetSize()) + " residues"
      );
      // fill the map
      FillMap();
    }

    //! @brief construct predion map from ProteinModel
    //! @param PROTEIN_MODEL Protein Model for which prediction map is going to be constructed
    PredictionMap::PredictionMap( const assemble::ProteinModel &PROTEIN_MODEL) :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< double>, double> >(),
      m_Sequences()
    {

      // iterate over the chains in the model
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // get the sequence
        const util::ShPtr< biol::AASequence> &sp_sequence( ( *chain_itr)->GetSequence());

        // make sure the sequence matches the size
        BCL_Assert
        (
          sp_sequence->GetSize() >= s_MinimumNumberResidues,
          "The given sequence needs to be at least " + util::Format()( s_MinimumNumberResidues) +
          " residues long, but it has only " + util::Format()( sp_sequence->GetSize()) + " residues"
        );

        // insert into sequences.
        m_Sequences.PushBack( sp_sequence);
      }

      // fill the map
      FillMap();
    }

    //! @brief virtual copy constructor
    PredictionMap *PredictionMap::Clone() const
    {
      return new PredictionMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PredictionMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief iterates over the stored sequences and finds the one with matching CHAIN_ID
    //! @param CHAIN_ID chain id of the sequence that is beings searched
    //! @return AASequence with the searched chain id
    const util::ShPtr< biol::AASequence> &PredictionMap::GetSequence( const char CHAIN_ID) const
    {
      // iterate over every sequence stored
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator sequence_itr( m_Sequences.Begin()),
          sequence_itr_end( m_Sequences.End());
        sequence_itr != sequence_itr_end;
        ++sequence_itr
      )
      {
        // if the sequence id matches return it
        if( ( *sequence_itr)->GetChainID() == CHAIN_ID)
        {
          return *sequence_itr;
        }
      }

      // else Exit
      BCL_Exit( "No sequence with the provided chain id \'" + util::Format()( CHAIN_ID) + "\' is not stored!!", -1);
      return m_Sequences.FirstElement();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Returns predictions for the provided AA Data pointers
    //! @param AA_DATA_POINTERS pair of pointers to AAData of amino acids of interest
    //! @return the pair of vector of predictions and the merged prediction
    const storage::Pair< linal::Vector< double>, double> &PredictionMap::GetPredictions
    (
      const storage::VectorND< 2, util::SiPtr< const biol::AAData> > &AA_DATA_POINTERS
    ) const
    {
      // initialize undefined predictions vector
      static const storage::Pair< linal::Vector< double>, double> s_undefined_predictions
      (
        linal::Vector< double>( Types::s_NumberValidTypes, util::GetUndefined< double>()),
        util::GetUndefined< double>()
      );

      // search in the hash map
      storage::HashMap< size_t, storage::Pair< linal::Vector< double>, double> >::const_iterator itr
      (
        Find( AA_DATA_POINTERS)
      );

      // if not found return undefined
      if( itr == End())
      {
        return s_undefined_predictions;
      }

      // end
      return itr->second;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PredictionMap::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &PredictionMap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

    //! @brief helper function to read predictions from a file
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PredictionMap::ReadPredictionMap
    (
      std::istream &ISTREAM
    )
    {
      // header that is going to be searched
      std::string header;

      // while reading header and not end of file yet
      while( ISTREAM >> header && !ISTREAM.eof())
      {
        // if header matches the identifier
        if( header == s_Identifier)
        {

          // read the two chain ids
          char chain_id_a, chain_id_b;
          ISTREAM >> chain_id_a >> chain_id_b;

          // check the chain_id_a exists
          BCL_Assert
          (
            GetSequence( chain_id_a).IsDefined(),
            "No sequence is stored in the map with given chain id: " + std::string( 1, chain_id_a)
          );

          // check the chain_id_b exists
          BCL_Assert
          (
            GetSequence( chain_id_b).IsDefined(),
            "No sequence is stored in the map with given chain id: " + std::string( 1, chain_id_b)
          );

          // read the predictions for two sequences with the provided chain ids
          ReadPredictions
          (
            ISTREAM,
            *GetSequence( chain_id_a),
            *GetSequence( chain_id_b)
          );
        }
      }

      // end
      return ISTREAM;
    }

    //! @brief helper function to write predictions to a file
    //! @param OSTREAM output stream to write to
    //! @param MERGED_THRESHOLD if merged predictions are below this threshold do not print them out
    //! @return output stream which was written to
    std::ostream &PredictionMap::WritePredictionMap
    (
      std::ostream &OSTREAM,
      const double MERGED_THRESHOLD
    ) const
    {

      // iterate over every sequence
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator seq_itr_a( m_Sequences.Begin()),
          seq_itr_end( m_Sequences.End());
        seq_itr_a != seq_itr_end;
        ++seq_itr_a
      )
      {
        // versus every other sequence
        for
        (
          util::ShPtrVector< biol::AASequence>::const_iterator seq_itr_b( m_Sequences.Begin());
          seq_itr_b != seq_itr_end;
          ++seq_itr_b
        )
        {
          // output the header
          OSTREAM << s_Identifier << ' '
                  << util::Format()( (*seq_itr_a)->GetChainID()) << ' '
                  << util::Format()( (*seq_itr_b)->GetChainID()) << '\n';

          // write predictions
          WritePredictions( OSTREAM, **seq_itr_a, **seq_itr_b, MERGED_THRESHOLD);
        }
      }
      return OSTREAM;
    }

    //! @brief Reads predictions for the provided sequences
    //! @param ISTREAM input stream
    //! @param SEQUENCE_A first sequence
    //! @param SEQUENCE_B second sequence
    //! @return istream which was read from
    std::istream &PredictionMap::ReadPredictions
    (
      std::istream &ISTREAM,
      const biol::AASequence &SEQUENCE_A,
      const biol::AASequence &SEQUENCE_B
    )
    {

      // initialize sequence id and AAType pairs and the vector hold predictions
      int seq_id_a, seq_id_b;
      char type_a, type_b;
      linal::Vector< double> predictions( Types::s_NumberValidTypes, util::GetUndefined< double>());

      // while it is possible to read
      while( ISTREAM >> seq_id_a >> type_a >> seq_id_b >> type_b && !ISTREAM.eof())
      {
        // make the sure sequences match
        BCL_Assert
        (
          SEQUENCE_A.GetAA( seq_id_a - 1)->GetType()->GetOneLetterCode() == type_a &&
          SEQUENCE_B.GetAA( seq_id_b - 1)->GetType()->GetOneLetterCode() == type_b,
          " The provided contact pair does not match the provided sequences" +
          util::Format()( SEQUENCE_A.GetAA( seq_id_a - 1)->GetType()->GetOneLetterCode()) +
          " vs " + util::Format()( type_a) + " and " +
          util::Format()( SEQUENCE_B.GetAA( seq_id_b - 1)->GetType()->GetOneLetterCode()) +
          " vs " + util::Format()( type_b)
        );

        // read the predictions vector
        ISTREAM >> predictions( GetTypes().HELIX_HELIX)
                >> predictions( GetTypes().HELIX_SHEET)
                >> predictions( GetTypes().SHEET_HELIX)
                >> predictions( GetTypes().STRAND_STRAND)
                >> predictions( GetTypes().SHEET_SHEET);

        // read the merged prediction from the file
        double merged_prediction;
        ISTREAM >> merged_prediction;

        // set prediction for unhandled HELIX_STRAND and STRAND_HELIX to 0
        predictions( GetTypes().HELIX_STRAND) = 0.0;
        predictions( GetTypes().STRAND_HELIX) = 0.0;

        // insert the predictions into the predictionmap
        InsertPredictions( *SEQUENCE_A.GetAA( seq_id_a - 1), *SEQUENCE_B.GetAA( seq_id_b - 1), predictions);
      }

      // end
      return ISTREAM;
    }

    //! @brief Writes predictions for the provided sequences
    //! @param OSTREAM output stream
    //! @param SEQUENCE_A first sequence
    //! @param SEQUENCE_B second sequence
    //! @param MERGED_THRESHOLD if merged predictions are below this threshold do not print them out
    //! @return ostream which was written to
    std::ostream &PredictionMap::WritePredictions
    (
      std::ostream &OSTREAM,
      const biol::AASequence &SEQUENCE_A,
      const biol::AASequence &SEQUENCE_B,
      const double MERGED_THRESHOLD
    ) const
    {

      // iterate over first sequence
      for
      (
        biol::AASequence::const_iterator aa_itr_a( SEQUENCE_A.Begin()),
          aa_itr_a_end( SEQUENCE_A.End());
        aa_itr_a != aa_itr_a_end;
        ++aa_itr_a
      )
      {
        // iterate over second sequence
        for
        (
          biol::AASequence::const_iterator aa_itr_b( SEQUENCE_B.Begin()),
            aa_itr_b_end( SEQUENCE_B.End());
          aa_itr_b != aa_itr_b_end;
          ++aa_itr_b
        )
        {
          // get the predictions for the two amino acids
          const storage::Pair< linal::Vector< double>, double> predictions
          (
            GetPredictions
            (
              storage::VectorND< 2, util::SiPtr< const biol::AAData> >
              (
                ( *aa_itr_a)->GetData(),
                ( *aa_itr_b)->GetData()
              )
            )
          );

          // if predictions are not defined, meaning they were not calculated, skip this residue pair
          if( !predictions.First().IsDefined())
          {
            continue;
          }
          // if the merged prediction is below the threshold, skip this residue pair
          if( predictions.Second() < MERGED_THRESHOLD)
          {
            continue;
          }

          static const util::Format s_format( util::Format().W( 5).FFP( 3));
          // output values
          OSTREAM << util::Format().W( 4)( ( *aa_itr_a)->GetSeqID()) << ' '
                  << ( *aa_itr_a)->GetType()->GetOneLetterCode() << ' '
                  << util::Format().W( 4)( ( *aa_itr_b)->GetSeqID()) << ' '
                  << ( *aa_itr_b)->GetType()->GetOneLetterCode() << ' '
                  << s_format( predictions.First()( GetTypes().HELIX_HELIX)) << ' '
                  << s_format( predictions.First()( GetTypes().HELIX_SHEET)) << ' '
                  << s_format( predictions.First()( GetTypes().SHEET_HELIX)) << ' '
                  << s_format( predictions.First()( GetTypes().STRAND_STRAND)) << ' '
                  << s_format( predictions.First()( GetTypes().SHEET_SHEET)) << ' '
                  << s_format( predictions.Second()) << '\n';
        }
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief This function generates the description for a single amino acid to be used for Contact ANNs
    //! @param AMINO_ACID amino acid of interest for which description will be generated
    //! @return the descriptor vector composed of 7(properties)+3(jufo)+20(blast)=30 values for a single amino acid
    linal::Vector< double> PredictionMap::GenerateDescription( const biol::AABase &AMINO_ACID)
    {
      // initialize an empty description vector
      linal::Vector< double> description_vector( s_NumberOfValuesPerResidue);

      //insert 7 amino acid properties
      description_vector.ReplaceElements( 0, AMINO_ACID.GetType()->GetPropertiesForANN());

      // insert jufo
      description_vector.ReplaceElements
      (
        biol::AATypeData::s_NumberPropertyTypesForANN,
        AMINO_ACID.GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction()
      );

      // insert blast profile
      description_vector.ReplaceElements
      (
        biol::AATypeData::s_NumberPropertyTypesForANN + 3,
        AMINO_ACID.GetBlastProfile().GetProfile()
      );

      // end
      return description_vector;
    }

    //! @brief This function generates the descriptions for all amino acids in the provided SEQUENCE
    //! @brief to used for contact ANNS. It aims to prevent repetitive generation of descriptions
    //! @param SEQUENCE AASequence of interest for which descriptions will be generated
    //! @return the descriptor matrix which include descriptions for each individiual amino acid
    linal::Matrix< double> PredictionMap::GenerateDescription( const biol::AASequence &SEQUENCE)
    {
      // initialize empty matrix
      linal::Matrix< double> description_matrix( SEQUENCE.GetSize(), s_NumberOfValuesPerResidue);

      // initialize row number
      size_t row( 0), row_end( description_matrix.GetNumberRows());

      // iterate over every amino acid in the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( SEQUENCE.Begin()),
          aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end && row < row_end;
        ++aa_itr, ++row
      )
      {
        // insert description for amino acid behind aa_itr
        description_matrix.ReplaceRow( row, GenerateDescription( **aa_itr));
      }

      // end
      return description_matrix;
    }

    //! @brief Merges prediction from 5 contact types into one
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param PREDICTIONS vector of predictions for 5 contact types
    //! @return merged value of 5 predictions
    double PredictionMap::MergePredictions
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const linal::Vector< double> &PREDICTIONS
    ) const
    {

      // store sspredictions for A and B
      linal::Vector3D sspredictions_a
      (
        AMINO_ACID_A.GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction()
      );
      linal::Vector3D sspredictions_b
      (
        AMINO_ACID_B.GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction()
      );

      // calculate the merged value
      const double merged_value
      (
        (
          PREDICTIONS( GetTypes().HELIX_HELIX) *
          sspredictions_a( biol::GetSSTypes().HELIX) *
          sspredictions_b( biol::GetSSTypes().HELIX)
        ) +
        (
          PREDICTIONS( GetTypes().HELIX_SHEET) *
          sspredictions_a( biol::GetSSTypes().HELIX) *
          sspredictions_b( biol::GetSSTypes().STRAND)
        ) +
        (
          PREDICTIONS( GetTypes().SHEET_HELIX) *
          sspredictions_a( biol::GetSSTypes().STRAND) *
          sspredictions_b( biol::GetSSTypes().HELIX)
        ) +
        (
          PREDICTIONS( GetTypes().STRAND_STRAND) *
          sspredictions_a( biol::GetSSTypes().STRAND) *
          sspredictions_b( biol::GetSSTypes().STRAND)
        ) +
        (
          PREDICTIONS( GetTypes().SHEET_SHEET) *
          sspredictions_a( biol::GetSSTypes().STRAND) *
          sspredictions_b( biol::GetSSTypes().STRAND)
        )
      );

      // return minimum of the merged value or 1.0
      return ( merged_value > 1.0 ? 1.0 : merged_value);
    }

    //! @brief for the given aa data pair, insert the predictions into the ObjectNDHashmap
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param PREDICTIONS vector of predictions for 5 contact types
    void PredictionMap::InsertPredictions
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const linal::Vector< double> &PREDICTIONS
    )
    {
      linal::Vector< double> this_prediction( PREDICTIONS);

      // calculate the merged predictions from predictions
      const double merged_prediction( MergePredictions( AMINO_ACID_A, AMINO_ACID_B, PREDICTIONS));

      // insert h,i
      Insert
      (
        storage::VectorND< 2, util::SiPtr< const biol::AAData> >
        (
          util::SiPtr< const biol::AAData>( AMINO_ACID_A.GetData()),
          util::SiPtr< const biol::AAData>( AMINO_ACID_B.GetData())
        ),
        storage::Pair< linal::Vector< double>, double>( this_prediction, merged_prediction)
      );

    }

    //! @brief reverse the predictions so HELIX_SHEET and SHEET_HELIX values are switched
    //! @param PREDICTIONS vector of predictions which is going to be reversed
    //! @return vector where the predictions so HELIX_SHEET and SHEET_HELIX values are switched
    linal::Vector< double> PredictionMap::ReversePredictions( const linal::Vector< double> &PREDICTIONS)
    {
      // reverse the predictions
      linal::Vector< double> predictions_swapped( PREDICTIONS);
      std::swap( predictions_swapped( GetTypes().HELIX_SHEET), predictions_swapped( GetTypes().SHEET_HELIX));

      // end
      return predictions_swapped;
    }

    //! @brief Given a specific contact type and two sequence windows, this functions puts together the input
    //! @brief data and calls the corresponding ANN and returns the single prediction
    //! @param CONTACT_TYPE contact type two amino acids are in
    //! @param DESCRIPTION_A description for the first amino acid window
    //! @param DESCRIPTION_B description for the second amino acid window
    //! @param POSITIONS position descriptors for the pair of amino acids
    //! @return prediction from ANN specific to the provided contact type for the provided descriptions and positions
    double PredictionMap::PredictSingleContact
    (
      const Type &CONTACT_TYPE,
      const linal::Vector< double> &DESCRIPTION_A,
      const linal::Vector< double> &DESCRIPTION_B,
      const linal::Vector3D &POSITIONS
    )
    {
      //instantiate a vector for the input data
      linal::Vector< double> input_data( DESCRIPTION_A.GetSize() + DESCRIPTION_B.GetSize() + POSITIONS.GetSize());

      if( CONTACT_TYPE == GetTypes().SHEET_HELIX)
      {
        //fill the vector with the 3 inputdatas
        input_data.ReplaceElements( 0, DESCRIPTION_B);
        input_data.ReplaceElements( DESCRIPTION_B.GetSize(), DESCRIPTION_A);
        input_data.ReplaceElements( DESCRIPTION_A.GetSize() + DESCRIPTION_B.GetSize(), POSITIONS);
      }
      else
      {
        //fill the vector with the 3 inputdatas
        input_data.ReplaceElements( 0, DESCRIPTION_A);
        input_data.ReplaceElements( DESCRIPTION_A.GetSize(), DESCRIPTION_B);
        input_data.ReplaceElements( DESCRIPTION_A.GetSize() + DESCRIPTION_B.GetSize(), POSITIONS);
      }

      //switch over the contact types and call the according neural network for the prediction
      double value( util::GetUndefined< double>());

      if( CONTACT_TYPE == GetTypes().HELIX_HELIX)
      {
        value = ANN_CONTACT_HELIX_HELIX( input_data);
      }
      else if( CONTACT_TYPE == GetTypes().HELIX_SHEET)
      {
        value = ANN_CONTACT_HELIX_SHEET( input_data);
      }
      else if( CONTACT_TYPE == GetTypes().SHEET_HELIX)
      {
        value = ANN_CONTACT_HELIX_SHEET( input_data);
      }
      else if( CONTACT_TYPE == GetTypes().STRAND_STRAND)
      {
        value = ANN_CONTACT_STRAND_STRAND( input_data);
      }
      else if( CONTACT_TYPE == GetTypes().SHEET_SHEET)
      {
        value = ANN_CONTACT_SHEET_SHEET( input_data);
      }

      // no networks available for that currently
      else if( CONTACT_TYPE == GetTypes().HELIX_STRAND)
      {
        value = double( 0);
      }
      else if( CONTACT_TYPE == GetTypes().STRAND_HELIX)
      {
        value = double( 0);
      }

      else
      {
        return util::GetUndefined< double>();
      }

      if( util::IsDefined( value))
      {
        value = std::min( 1.0, std::max( 0.0, value));
      }
      return value;
    }

    //! @brief This function fills m_Map with prediction for each residue couple from each contact type ANN for
    //! @brief sequences stored in m_Sequences
    void PredictionMap::FillMap()
    {
      // check jufo and blast for profiles are stored for the provided sequences
      CheckBlastAndJufoForSequences();

      // iterate over every sequence stored
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator seq_itr_a( m_Sequences.Begin()),
          seq_itr_end( m_Sequences.End());
        seq_itr_a != seq_itr_end;
        ++seq_itr_a
      )
      {
        FillMap( **seq_itr_a);

        // iterate over every other sequence
        for
        (
          util::ShPtrVector< biol::AASequence>::const_iterator seq_itr_b( seq_itr_a + 1);
          seq_itr_b != seq_itr_end;
          ++seq_itr_b
        )
        {
          FillMap( **seq_itr_a, **seq_itr_b);
        }
      }
    }

    //! @brief This function fills m_Map with prediction for each residue couple from each contacttype ANN for
    //! @brief two different sequences
    //! @param SEQUENCE AASequence to be used for filling the map with predictions
    void PredictionMap::FillMap( const biol::AASequence &SEQUENCE)
    {
      // fill up the matrix for data for every single residue
      linal::Matrix< double> data( GenerateDescription( SEQUENCE));

      // temporary vector for storing 5 predictions between 2 AA
      linal::Vector< double> predictions( Types::s_NumberValidTypes);

      // iterate over every residue h
      for
      (
        size_t aa_ctr_a( s_NumberSkippedResiduesInEnds);
        aa_ctr_a < SEQUENCE.GetSize() - s_NumberSkippedResiduesInEnds - s_MinSequenceDistance;
        ++aa_ctr_a
      )
      {
        // and every other following residue i
        for
        (
          size_t aa_ctr_b( aa_ctr_a + s_MinSequenceDistance);
          aa_ctr_b < SEQUENCE.GetSize() - s_NumberSkippedResiduesInEnds;
          ++aa_ctr_b
        )
        {
          // for every contact type form inputs and collect network results
          for( size_t type( GetTypes().HELIX_HELIX); type < Types::s_NumberValidTypes; ++type)
          {
            Type contact_type( type);

            predictions( contact_type) =
              PredictSingleContact
              (
                contact_type,
                linal::Vector< double>
                (
                  s_NumberOfValuesPerResidue * ( contact_type->GetWindowLengths().First()),
                  data[ aa_ctr_a - contact_type->GetWindowRadii().First()]
                ),
                linal::Vector< double>
                (
                  s_NumberOfValuesPerResidue * ( contact_type->GetWindowLengths().Second()),
                  data[ aa_ctr_b - contact_type->GetWindowRadii().Second()]
                ),
                linal::Vector3D( aa_ctr_a, aa_ctr_b - aa_ctr_a, SEQUENCE.GetSize() - aa_ctr_b)
              );
          }

          // insert the predictions for the selected residue pair
          InsertPredictions
          (
            *SEQUENCE.GetAA( aa_ctr_a), *SEQUENCE.GetAA( aa_ctr_b), predictions
          );

          // insert the predictions for the reversed residue pair
          InsertPredictions
          (
            *SEQUENCE.GetAA( aa_ctr_b), *SEQUENCE.GetAA( aa_ctr_a), ReversePredictions( predictions)
          );

        } // loop over i
      } // loop over h
    }

    //! @brief This function fills m_Map with prediction for each residue couple from each contacttype ANN for
    //! @brief two different sequences
    //! @param SEQUENCE_A First AASequence to be used for filling the map with predictions
    //! @param SEQUENCE_B Second AASequence to be used for filling the map with predictions
    void PredictionMap::FillMap
    (
      const biol::AASequence &SEQUENCE_A,
      const biol::AASequence &SEQUENCE_B
    )
    {
      // fill up the matrix for data for every single residue
      const linal::Matrix< double> data_a( GenerateDescription( SEQUENCE_A));
      const linal::Matrix< double> data_b( GenerateDescription( SEQUENCE_B));

      // temporary vector for storing 5 predictions between 2 AA
      linal::Vector< double> predictions( Types::s_NumberValidTypes);

      // iterate over every residue h
      for
      (
        size_t aa_ctr_a( s_NumberSkippedResiduesInEnds);
        aa_ctr_a < SEQUENCE_A.GetSize() - s_NumberSkippedResiduesInEnds;
        ++aa_ctr_a
      )
      {
        // and every other following residue i
        for
        (
          size_t aa_ctr_b( s_NumberSkippedResiduesInEnds);
          aa_ctr_b < SEQUENCE_B.GetSize() - s_NumberSkippedResiduesInEnds;
          ++aa_ctr_b
        )
        {
          // accumulate the predictions for each contacttype
          for( size_t type = GetTypes().HELIX_HELIX; type < Types::s_NumberValidTypes; ++type)
          {
            Type contact_type( type);

            predictions( contact_type) =
              PredictSingleContact
              (
                contact_type,
                linal::Vector< double>
                (
                  s_NumberOfValuesPerResidue * ( contact_type->GetWindowLengths().First()),
                  data_a[ aa_ctr_a - contact_type->GetWindowRadii().First()]
                ),
                linal::Vector< double>
                (
                  s_NumberOfValuesPerResidue * ( contact_type->GetWindowLengths().Second()),
                  data_b[ aa_ctr_b - contact_type->GetWindowRadii().Second()]
                ),
                linal::Vector3D( 40, 80, 40)
              );
          }

          // insert the predictions for the selected residue pair
          InsertPredictions
          (
            *SEQUENCE_A.GetAA( aa_ctr_a), *SEQUENCE_B.GetAA( aa_ctr_b), predictions
          );

          // insert the predictions for the reversed residue pair
          InsertPredictions
          (
            *SEQUENCE_B.GetAA( aa_ctr_b), *SEQUENCE_A.GetAA( aa_ctr_a), ReversePredictions( predictions)
          );

        } // loop over i
      } // loop over h
    }

    //! @brief this functions checks that jufo and blast profiles are stored for every sequence provided
    void PredictionMap::CheckBlastAndJufoForSequences() const
    {
      // iterate over every sequence stored
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator seq_itr( m_Sequences.Begin()),
          seq_itr_end( m_Sequences.End());
        seq_itr != seq_itr_end;
        ++seq_itr
      )
      {
         //check that blast profile is stored
         BCL_Assert
         (
           ( *seq_itr)->GetFirstAA()->GetBlastProfilePtr().IsDefined(),
           "Blast Profile is not available for chain " + util::Format()( ( *seq_itr)->GetChainID())
         );

         //check that jufo is stored
         BCL_Assert
         (
           ( *seq_itr)->GetFirstAA()->GetSSPrediction( sspred::GetMethods().e_JUFO).IsDefined(),
           "Jufo Profile is not available for chain " + util::Format()( ( *seq_itr)->GetChainID())
         );
      }
    }

  } // namespace contact
} // namespace bcl
