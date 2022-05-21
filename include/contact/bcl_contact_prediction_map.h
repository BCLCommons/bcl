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

#ifndef BCL_CONTACT_PREDICTION_MAP_H_
#define BCL_CONTACT_PREDICTION_MAP_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_contact_types.h"
#include "biol/bcl_biol_aa_data.h"
#include "storage/bcl_storage_object_nd_hash_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PredictionMap
    //! @brief This class, derived from ObjectNDHashMap, is designed to store predictions  in the map from ANNs
    //! @details For each amino acid pair in the given amino acid sequences, ANNs are used to predict the contact probability
    //! for this amino acid for each 5 contact types. These predictions are then stored in a hash-map structure
    //! where retrieval by passing the pair of amino acids of interest is possible. In addition to 5 predictions for
    //! each amino acid pair, the merged predictions are also stored.
    //!
    //! @see @link example_contact_prediction_map.cpp @endlink
    //! @author karakam, woetzen
    //! @date 15.08.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PredictionMap :
      public storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< double>, double> >
    {

    private:

    //////////
    // data //
    //////////

      //! this is ShPtrVector of sequences
      util::ShPtrVector< biol::AASequence> m_Sequences;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PredictionMap();

      //! @brief construct prediction map ( from ANNs) from one chain( intra-sequence contacts)
      //! @param CHAIN chain
      PredictionMap( const assemble::Chain &CHAIN);

      //! @brief construct predion map from ProteinModel
      //! @param PROTEIN_MODEL Protein Model for which prediction map is going to be constructed
      PredictionMap( const assemble::ProteinModel &PROTEIN_MODEL);

      //! @brief virtual copy constructor
      PredictionMap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns ShPtrVector of AASequences
      //! @return ShPtrVector of AASequences
      const util::ShPtrVector< biol::AASequence> &GetSequences() const
      {
        return m_Sequences;
      }

      //! @brief iterates over the stored sequences and finds the one with matching CHAIN_ID
      //! @param CHAIN_ID chain id of the sequence that is beings searched
      //! @return AASequence with the searched chain id
      const util::ShPtr< biol::AASequence> &GetSequence( const char CHAIN_ID) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Returns predictions for the provided AA Data pointers
      //! @param AA_DATA_POINTERS pair of pointers to AAData of amino acids of interest
      //! @return the pair of vector of predictions and the merged prediction
      const storage::Pair< linal::Vector< double>, double> &GetPredictions
      (
        const storage::VectorND< 2, util::SiPtr< const biol::AAData> > &AA_DATA_POINTERS
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief helper function to read predictions from a file
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadPredictionMap
      (
        std::istream &ISTREAM
      );

      //! @brief helper function to write predictions to a file
      //! @param OSTREAM output stream to write to
      //! @param MERGED_THRESHOLD if merged predictions are below this threshold do not print them out
      //! @return output stream which was written to
      std::ostream &WritePredictionMap
      (
        std::ostream &OSTREAM,
        const double MERGED_THRESHOLD = 0.0
      ) const;

    private:

      //! @brief Reads predictions for the provided sequences
      //! @param ISTREAM input stream
      //! @param SEQUENCE_A first sequence
      //! @param SEQUENCE_B second sequence
      //! @return istream which was read from
      std::istream &ReadPredictions
      (
        std::istream &ISTREAM,
        const biol::AASequence &SEQUENCE_A,
        const biol::AASequence &SEQUENCE_B
      );

      //! @brief Writes predictions for the provided sequences
      //! @param OSTREAM output stream
      //! @param SEQUENCE_A first sequence
      //! @param SEQUENCE_B second sequence
      //! @param MERGED_THRESHOLD if merged predictions are below this threshold do not print them out
      //! @return ostream which was written to
      std::ostream &WritePredictions
      (
        std::ostream &OSTREAM,
        const biol::AASequence &SEQUENCE_A,
        const biol::AASequence &SEQUENCE_B,
        const double MERGED_THRESHOLD
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief This function generates the description for a single amino acid to be used for Contact ANNs
      //! @param AMINO_ACID amino acid of interest for which description will be generated
      //! @return the descriptor vector composed of 7(properties)+3(jufo)+20(blast)=30 values for a single amino acid
      linal::Vector< double> GenerateDescription( const biol::AABase &AMINO_ACID);

      //! @brief This function generates the descriptions for all amino acids in the provided SEQUENCE
      //! @brief to used for contact ANNS. It aims to prevent repetitive generation of descriptions
      //! @param SEQUENCE AASequence of interest for which descriptions will be generated
      //! @return the descriptor matrix which include descriptions for each individiual amino acid
      linal::Matrix< double> GenerateDescription( const biol::AASequence &SEQUENCE);

      //! @brief Merges prediction from 5 contact types into one
      //! @param AMINO_ACID_A first amino acid
      //! @param AMINO_ACID_B second amino acid
      //! @param PREDICTIONS vector of predictions for 5 contact types
      //! @return merged value of 5 predictions
      double MergePredictions
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const linal::Vector< double> &PREDICTIONS
      ) const;

      //! @brief for the given aa data pair, insert the predictions into the ObjectNDHashmap
      //! @param AMINO_ACID_A first amino acid
      //! @param AMINO_ACID_B second amino acid
      //! @param PREDICTIONS vector of predictions for 5 contact types
      void InsertPredictions
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const linal::Vector< double> &PREDICTIONS
      );

      //! @brief reverse the predictions so HELIX_SHEET and SHEET_HELIX values are switched
      //! @param PREDICTIONS vector of predictions which is going to be reversed
      //! @return vector where the predictions so HELIX_SHEET and SHEET_HELIX values are switched
      static linal::Vector< double> ReversePredictions( const linal::Vector< double> &PREDICTIONS);

      //! @brief Given a specific contact type and two sequence windows, this functions puts together the input
      //! @brief data and calls the corresponding ANN and returns the single prediction
      //! @param CONTACT_TYPE contact type two amino acids are in
      //! @param DESCRIPTION_A description for the first amino acid window
      //! @param DESCRIPTION_B description for the second amino acid window
      //! @param POSITIONS position descriptors for the pair of amino acids
      //! @return prediction from ANN specific to the provided contact type for the provided descriptions and positions
      double PredictSingleContact
      (
        const Type &CONTACT_TYPE,
        const linal::Vector< double> &DESCRIPTION_A,
        const linal::Vector< double> &DESCRIPTION_B,
        const linal::Vector3D &POSITIONS
      );

      //! @brief This function fills m_Map with prediction for each residue couple from each contacttype ANN for
      //! @brief sequences stored in m_Sequences
      void FillMap();

      //! @brief This function fills m_Map with prediction for each residue couple from each contacttype ANN for
      //! @brief two different sequences
      //! @param SEQUENCE AASequence to be used for filling the map with predictions
      void FillMap( const biol::AASequence &SEQUENCE);

      //! @brief This function fills m_Map with prediction for each residue couple from each contacttype ANN for
      //! @brief two different sequences
      //! @param SEQUENCE_A First AASequence to be used for filling the map with predictions
      //! @param SEQUENCE_B Second AASequence to be used for filling the map with predictions
      void FillMap
      (
        const biol::AASequence &SEQUENCE_A,
        const biol::AASequence &SEQUENCE_B
      );

      //! @brief this functions checks that jufo and blast profiles are stored for every sequence provided
      //! @brief so that the predictions can be generated correctly
      void CheckBlastAndJufoForSequences() const;

    }; // class Prediction Map

  } // namespace contact
} // namespace bcl

#endif //BCL_CONTACT_PREDICTION_MAP_H_

