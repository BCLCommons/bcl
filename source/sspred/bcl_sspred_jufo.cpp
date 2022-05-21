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
#include "sspred/bcl_sspred_jufo.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "sspred/bcl_sspred_jufo_ann.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

  //////////
  // data //
  //////////

    //! @brief returns total number of values used per residue
    //! @return total number of values used per residue
    size_t JUFO::GetNumberValuesPerResidue()
    {
      // initialize static size_t
      static const size_t s_number_values_per_residue( 27);

      // return
      return s_number_values_per_residue;
    }

    //! @brief returns average property vector to be used for undefined amino acids
    //! @return average property vector to be used for undefined amino acids
    linal::Vector< double> JUFO::CalculateAveragePropertyVector()
    {
      // initialize static average property array to be used for undefined amino acids
      const double average_7_property_array[ 7] =
        { 2.169, 0.146, 3.330, 0.372, 6.136, 0.000, 0.000}; //weighed averaged
      const linal::Vector< double> average_7_property_vector( 7, average_7_property_array);

      // initialize with zeroes
      linal::Vector< double> average_property_vector( GetNumberValuesPerResidue(), 0.0);

      // replace the 7 property part with the array
      // initialize static vector that holds these values
      average_property_vector.ReplaceElements( 0, average_7_property_vector);

      // return
      return average_7_property_vector;
    }

    //! @brief returns average property vector to be used for undefined amino acids
    //! @return average property vector to be used for undefined amino acids
    const linal::Vector< double> &JUFO::GetAveragePropertyVector()
    {
      // initialize static vector that holds these values
      static const linal::Vector< double> s_average_property_vector
      (
        CalculateAveragePropertyVector()
      );

      // return
      return s_average_property_vector;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    JUFO::JUFO() :
      m_Prediction( GetDefaultPredictionVector())
    {
    }

    //! @brief constructor from a linal::Vector3D
    //! @param VECTOR linal::Vector3D of probabilities
    JUFO::JUFO( const linal::Vector3D &VECTOR) :
      m_Prediction( VECTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new JUFO
    JUFO *JUFO::Clone() const
    {
      return new JUFO( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &JUFO::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &JUFO::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".jufo");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D JUFO::GetThreeStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> JUFO::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( m_Prediction, biol::GetEnvironmentTypes().e_Solution);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &JUFO::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      // initialize vector3d
      linal::Vector3D prediction;

      // temporary variable to store one state
      std::string one_state;

      // read one state
      ISTREAM >> one_state;

      // read COIL, HELIX and STRAND
      ISTREAM >> prediction( biol::GetSSTypes().COIL);
      ISTREAM >> prediction( biol::GetSSTypes().HELIX);
      ISTREAM >> prediction( biol::GetSSTypes().STRAND);

      // normalize the vector
      prediction.SetToSum();

      // set the predictions for this amino acid
      AMINO_ACID.SetSSPrediction( GetMethods().e_JUFO, JUFO( prediction));

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &JUFO::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      // call standard read function and return it
      return ReadStandardPredictionsForAASequence( ISTREAM, AA_SEQUENCE, GetMethods().e_JUFO);
    }

    //! @brief iterates over the sequence and calculates the jufo predictions for every residue in the sequence
    //! @param SEQUENCE AASequence for which JUFO will be calculated
    void JUFO::Calculate( biol::AASequence &SEQUENCE)
    {
      // make sure the blast profile exists for this sequence
      BCL_Assert
      (
        SEQUENCE.GetFirstAA()->GetBlastProfilePtr().IsDefined(),
        "Blast Profile is not available"
      );

      // initialize input, output vectors
      linal::Vector< double> input_data( 1053, 0.0);
      linal::Vector< double> output( 3);

      // for every residue in this protein model, calculate jufo
      for( size_t aa_ctr( 0); aa_ctr != SEQUENCE.GetSize(); ++aa_ctr)
      {
        size_t current_index( 0);
        // collect very far neighbor averaged information if there are residues earlier than 30 residues.
        if( aa_ctr > 30)
        {
          input_data.ReplaceElements( current_index, AverageData( SEQUENCE.SubSequence( 0, aa_ctr - 30)));
        }
        // else just use average values
        else
        {
          input_data.ReplaceElements( current_index, GetAveragePropertyVector());
        }
        current_index += GetNumberValuesPerResidue();

        // collect far neighbor grouped information
        for( size_t group_no( 0); group_no < 3; ++group_no)
        {
          // if there are enough residues to collect information
          if( aa_ctr > 15 + 5 * ( 3 - group_no - 1))
          {
            // find where the next group would begin
            size_t end( aa_ctr - 15 - 5 * ( 3 - group_no - 1));
            size_t begin( 0);
            // if there are enough residues to have a complete group of 5 residues, set the beginning to end - 5, otherwise begin from 0.
            if( end > 4)
            {  begin = end - 5;
            }
            input_data.ReplaceElements( current_index, AverageData( SEQUENCE.SubSequence( begin, end - begin)));
          }
          // else use average values
          else
          {
            input_data.ReplaceElements( current_index, GetAveragePropertyVector());
          }
          current_index += GetNumberValuesPerResidue();
        }

        // collect close neighbor information
        for( size_t num_neighbor( 15); num_neighbor > 0; --num_neighbor)
        {
          // if that neighbor exists
          if( aa_ctr >= num_neighbor)
          {
            input_data.ReplaceElements( current_index, GenerateData( *SEQUENCE.GetAA( aa_ctr - num_neighbor)));
          }
          // else use average values
          else
          {
            input_data.ReplaceElements( current_index, GetAveragePropertyVector());
          }
          current_index += GetNumberValuesPerResidue();
        }

        // collect information for this residue
        input_data.ReplaceElements( current_index, GenerateData( *SEQUENCE.GetAA( aa_ctr)));
        current_index += GetNumberValuesPerResidue();

        // collect close neighbor information
        for( size_t num_neighbor( 0); num_neighbor < 15; ++num_neighbor)
        {
          // if that neighbor exists
          if( SEQUENCE.GetSize() - aa_ctr - 1 > num_neighbor)
          {
            input_data.ReplaceElements( current_index, GenerateData( *SEQUENCE.GetAA( aa_ctr + num_neighbor + 1)));
          }
          // else use average values
          else
          {
            input_data.ReplaceElements( current_index, GetAveragePropertyVector());
          }
          current_index += GetNumberValuesPerResidue();
        }

        // collect far neighbor grouped information
        for( size_t group_no( 0); group_no < 3; ++group_no)
        {
          // if there are enough residues to collect information
          if( SEQUENCE.GetSize() - aa_ctr > 16 + 5 * group_no)
          {
            // find where this group would begin
            size_t begin( aa_ctr + 16 + 5 * ( group_no));
            size_t end( SEQUENCE.GetSize());
            // if there are enough residues to have a complete group of 5 residues, set the beginning to end - 5, otherwise begin from 0.
            if( begin < SEQUENCE.GetSize() - 4)
            {  end = begin + 5;
            }
            input_data.ReplaceElements( current_index, AverageData( SEQUENCE.SubSequence( begin, end - begin)));
          }
          // else use average values
          else
          {
            input_data.ReplaceElements( current_index, GetAveragePropertyVector());
          }
          current_index += GetNumberValuesPerResidue();
        }

        // collect very far neighbor averaged information
        if( SEQUENCE.GetSize() - aa_ctr > 30)
        {
          size_t begin( aa_ctr + 30);
          input_data.ReplaceElements( current_index, AverageData( SEQUENCE.SubSequence( begin, SEQUENCE.GetSize() - begin)));
        }
        // else use average values
        else
        {
          input_data.ReplaceElements( current_index, GetAveragePropertyVector());
        }
        current_index += GetNumberValuesPerResidue();

        // at this point, the input is ready so now let's call the network and get the result
        output = JUFOANN().F( input_data);
        double sum( 0.0);

        // normalize output
        for( size_t a( 0); a < 3; ++a)
        {
          output( a) += 0.1;
          output( a) /= 1.2;
          sum += output( a);
        }

        // divide by the sum
        output /= sum;

        // prepare the output vector
        linal::Vector< double> output_vector( biol::GetSSTypes().COIL.GetIndex() + 1, double( 0));
        output_vector.ReplaceElements( 0, output);

        // set the values
        util::ShPtr< biol::AABase> aa_with_jufo( SEQUENCE.GetAA( aa_ctr)->Clone());
        aa_with_jufo->SetSSPrediction
        (
          GetMethods().e_JUFO,
          JUFO( linal::Vector3D( output_vector.Begin()))
        );
        SEQUENCE.SetAA( aa_ctr, *aa_with_jufo);

      } // iterate over every residue

    }

    //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which JUFO will be calculated
    void JUFO::Calculate( assemble::ProteinModel &PROTEIN_MODEL)
    {
      // iterate over the Chains in the given model
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // calculate for the given sequence
        Calculate( *( *chain_itr)->GetSequence());
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &JUFO::Read( std::istream &ISTREAM)
    {
      // read members
      ISTREAM >> m_Prediction;

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &JUFO::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // end
      return OSTREAM;
    };

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generates the input data composed of 7(properties)+20(blast)=27 values for a single AA
    //! @param AMINO_ACID Amino acid of interest
    //! @return input data for the given amino acid
    linal::Vector< double> JUFO::GenerateData( const biol::AABase &AMINO_ACID)
    {
      // initialize a vector of size GetNumberValuesPerResidue()
      linal::Vector< double> numerical_representation( GetNumberValuesPerResidue());

      //Get 7 AA properties
      numerical_representation.ReplaceElements( 0, AMINO_ACID.GetType()->GetPropertiesForANN());
      numerical_representation.ReplaceElements
      (
        biol::AATypeData::s_NumberPropertyTypesForANN, AMINO_ACID.GetBlastProfile().GetProfile()
      );

      // return
      return numerical_representation;
    }

    //! @brief calculates the average values for input data for JUFO for a given sequence and returns it in a vector
    //! @param SEQUENCE AASequence of interest
    //! @return input data for given AASequence
    linal::Vector< double> JUFO::AverageData( const biol::AASequence &SEQUENCE)
    {
      linal::Vector< double> data( GetNumberValuesPerResidue(), 0.0);

      // sum up vectors of inputs
      for
      (
        biol::AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // sum up the data
        data += GenerateData( **aa_itr);
      }

      // now divide every value by the size of the sequence
      data /= double( SEQUENCE.GetSize());

      // for columns 8 to the end, convert to integer value
      for( size_t a( 7); a < 27; a++)
      {
        // convert to integer
        data( a) = int( data( a));
      }

      // return
      return data;
    }

  } // namespace sspred
} // namespace bcl
