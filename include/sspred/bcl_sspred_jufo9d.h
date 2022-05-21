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

#ifndef BCL_SSPRED_JUFO9D_H_
#define BCL_SSPRED_JUFO9D_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "model/bcl_model.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_method_interface.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"

// external includes - sorted alphabetically
namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class JUFO9D
    //! @brief the JUFO9D secondary structure prediction algorithm
    //! @details secondary structure prediction class that implements Jufo9D which predicts 9-states.
    //!
    //! @see @link example_sspred_jufo9d.cpp @endlink
    //! @author karakam
    //! @date Jun 3, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API JUFO9D :
      public MethodInterface
    {

    private:

    //////////
    // data //
    //////////

      //! vector to store 3 state predictions
      linal::Matrix< double> m_Prediction;

      //! probability to pick a state at random
      static const double s_BaseProbability;

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
      JUFO9D();

      //! @brief constructor from a vector of 9 values
      //! @param VECTOR vector of values
      JUFO9D( const linal::Vector< double> &VECTOR);

      //! @brief constructor from a linal::Matrix
      //! @param MATRIX linal::Matrix of probabilities
      JUFO9D( const linal::Matrix< double> &MATRIX);

      //! @brief Clone function
      //! @return pointer to new JUFO
      JUFO9D *Clone() const;

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

      //! @brief probability that this is a TM-span of type given
      //! @return probability that this is a TM-span of type given
      double TMTypeProbability( const biol::SSType &SS_TYPE) const;

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

      //! @brief iterates over the sequences and calculates the jufo predictions for every residue in the sequence
      //! @param MODEL protein model of interest
      //! @param MULTIMER true if the sequence is multimeric
      static void Calculate( assemble::ProteinModel &MODEL, const bool MULTIMER = false);

      //! @brief iterates over the sequences and calculates the jufo predictions for every residue in the sequence
      //! @param MODEL protein model of interest
      //! @param MULTIMER true if the sequence is multimeric
      static void Calculate( assemble::ProteinModelWithCache &MODEL, const bool MULTIMER = false);

      //! @brief get the set of descriptors used for training and as input for Jufo9D ANNs
      //! @param MULTIMER true if the sequence is multimeric
      static util::Implementation< descriptor::Base< biol::AABase, float> > GetJufo9DANNDescriptors( const bool &MULTIMER);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator = for assigning given JUFO_RHS to this object
      //! @param JUFO_RHS JUFO to be assigned
      //! @return this object after assignment
      JUFO9D &operator =( const JUFO9D &JUFO_RHS);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief write two state TM-helix prediction from sequence to stream
      //! @param OSTREAM stream to write to
      //! @param AA_SEQUENCE sequence containing predictions
      //! @param SS_TYPE sstype for the 2 state prediction, must be helix or strand
      //! @return output stream
      static std::ostream &WriteTwoStateTMPredictions
      (
        std::ostream &OSTREAM,
        const biol::AASequence &AA_SEQUENCE,
        const biol::SSType &SS_TYPE
      );

      //! @brief write three state TM prediction from sequence to stream
      //! @param OSTREAM stream to write to
      //! @param AA_SEQUENCE sequence containing predictions
      //! @return output stream
      static std::ostream &WriteThreeStateTMPredictions( std::ostream &OSTREAM, const biol::AASequence &AA_SEQUENCE);

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class JUFO9D

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_JUFO9D_H_
