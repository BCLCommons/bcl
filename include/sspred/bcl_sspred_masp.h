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

#ifndef BCL_SSPRED_MASP_H_
#define BCL_SSPRED_MASP_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically

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
    //! @class MASP
    //! @brief stores prediction for MASP: Membrane-Association and secondary structure Predictor
    //! @details stores prediction for MASP, which predicts both TM-Strands for (TM-Beta proteins) and TM-Helices, as
    //!          well as soluble secondary structure
    //!
    //! @see @link example_sspred_masp.cpp @endlink
    //! @author mendenjl
    //! @date May 15, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MASP :
      public MethodInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Probability vectory for TM-prediction and secondary structure.  Only states MH, ME, SH, SE, SC are considered
      linal::Matrix< double> m_Prediction;

    public:

      //! single instance of that class
      static const util::ObjectInterface *s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MASP();

      //! @brief constructor from ss probabilities
      //! @param MEMBRANE_SS_PREDICTION the membrane secondary structure prediction for the residue
      //! @param SOLUBLE_SS_PREDICTION the soluble secondary structure prediction for the residue
      MASP
      (
        const linal::Vector3D &MEMBRANE_SS_PREDICTION,
        const linal::Vector3D &SOLUBLE_SS_PREDICTION
      );

      //! @brief Clone function
      //! @return pointer to new Masp
      MASP *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
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
      std::istream &ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const;

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const;

    ////////////////
    // Prediction //
    ////////////////

      //! @brief Predict whether a given protein is a beta barrel and add MASP-BB environment to the protein
      //!        The secondary structure and environment for all residues that are not predicted to be beta-barrels will
      //!        not be changed
      //! @param MODEL ProteinModel for which JUFO will be calculated
      //! @return true if the protein is predicted to be an integral, bacterial or mitochondrial outer-membrane, beta-barrel protein
      static bool PredictBetaBarrelEnvironmentSS( assemble::ProteinModel &MODEL);

      //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
      //!        The secondary structure and environment for all residues that are not predicted to be beta-barrels will
      //!        not be changed
      //! @param MODEL ProteinModel for which MASP-BB will be calculated
      //! @return true if the protein is predicted to be an integral, bacterial or mitochondrial outer-membrane, beta-barrel protein
      static bool PredictBetaBarrelEnvironmentSS( assemble::ProteinModelWithCache &MODEL);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MASP

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_MASP_H_
