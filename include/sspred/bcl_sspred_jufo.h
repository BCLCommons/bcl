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

#ifndef BCL_SSPRED_JUFO_H_
#define BCL_SSPRED_JUFO_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_method_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically
namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class JUFO
    //! @brief implements the JUFO secondary structure prediction algorithm
    //! @details stores and calculates JUFO secondary structure prediction using neural networks
    //!
    //! @see @link example_sspred_jufo.cpp @endlink
    //! @author karakam
    //! @date Jun 3, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API JUFO :
      public MethodInterface
    {

    private:

    //////////
    // data //
    //////////

      //! vector to store 3 state predictions
      linal::Vector3D m_Prediction;

      //! @brief returns total number of values used per residue
      //! @return total number of values used per residue
      static size_t GetNumberValuesPerResidue();

      //! @brief calculates average property vector to be used for undefined amino acids
      //! @return average property vector to be used for undefined amino acids
      static linal::Vector< double> CalculateAveragePropertyVector();

      //! @brief returns average property vector to be used for undefined amino acids
      //! @return average property vector to be used for undefined amino acids
      static const linal::Vector< double> &GetAveragePropertyVector();

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      JUFO();

      //! @brief constructor from a linal::Vector3D
      //! @param VECTOR linal::Vector3D of probabilities
      JUFO( const linal::Vector3D &VECTOR);

      //! @brief Clone function
      //! @return pointer to new JUFO
      JUFO *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get file extension associated with this Method
      //! @return file extension associated with this Method
      const std::string &GetFileExtension() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get three state, environment independent secondary structure prediction
      //! @return three state, environment independent secondary structure prediction
      linal::Vector3D GetThreeStatePrediction() const;

      //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
      //! @return three state secondary structure prediction
      linal::Matrix< double> GetNineStatePrediction() const;

      //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AMINO_ACID amino acid into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAA
      (
        std::istream &ISTREAM,
        biol::AABase &AMINO_ACID
      ) const;

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAASequence
      (
        std::istream &ISTREAM,
        biol::AASequence &AA_SEQUENCE
      ) const;

      //! @brief iterates over the sequence and calculates the jufo predictions for every residue in the sequence
      //! @param SEQUENCE AASequence for which JUFO will be calculated
      static void Calculate( biol::AASequence &SEQUENCE);

      //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
      //! @param PROTEIN_MODEL ProteinModel for which JUFO will be calculated
      static void Calculate( assemble::ProteinModel &PROTEIN_MODEL);

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief generates the input data composed of 7(properties)+20(blast)=27 values for a single AA
      //! @param AMINO_ACID Amino acid of interest
      //! @return input data for the given amino acid
      static linal::Vector< double> GenerateData( const biol::AABase &AMINO_ACID);

      //! @brief calculates the average values for input data for JUFO for a given sequence and returns it in a vector
      //! @param SEQUENCE AASequence of interest
      //! @return input data for given AASequence
      static linal::Vector< double> AverageData( const biol::AASequence &SEQUENCE);

    }; // class JUFO

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_JUFO_H_
