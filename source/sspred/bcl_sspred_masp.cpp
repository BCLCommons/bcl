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
#include "sspred/bcl_sspred_masp.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "model/bcl_model.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "util/bcl_util_implementation.h"
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
    MASP::MASP() :
      m_Prediction( GetDefaultPredictionMatrix())
    {
    }

    //! @brief constructor from ss probabilities
    //! @param MEMBRANE_SS_PREDICTION the membrane secondary structure prediction for the residue
    //! @param SOLUBLE_SS_PREDICTION the soluble secondary structure prediction for the residue
    MASP::MASP
    (
      const linal::Vector3D &MEMBRANE_SS_PREDICTION,
      const linal::Vector3D &SOLUBLE_SS_PREDICTION
    ) :
      m_Prediction( 3, 3, 0.0)
    {
      // determine SS type values
      const size_t solution_index( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex());
      const size_t membrane_index( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex());
      m_Prediction( solution_index, 0) = SOLUBLE_SS_PREDICTION( 0);
      m_Prediction( solution_index, 1) = SOLUBLE_SS_PREDICTION( 1);
      m_Prediction( solution_index, 2) = SOLUBLE_SS_PREDICTION( 2);
      m_Prediction( membrane_index, 0) = MEMBRANE_SS_PREDICTION( 0);
      m_Prediction( membrane_index, 1) = MEMBRANE_SS_PREDICTION( 1);
      m_Prediction( membrane_index, 2) = MEMBRANE_SS_PREDICTION( 2);
      m_Prediction.AsVector().SetToSum( 1.0);
    }

    //! @brief Clone function
    //! @return pointer to new MASP
    MASP *MASP::Clone() const
    {
      return new MASP( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MASP::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &MASP::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".masp");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D MASP::GetThreeStatePrediction() const
    {
      return ConvertNineStateToThreeState( m_Prediction);
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> MASP::GetNineStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &MASP::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      // create a vector to hold the soluble-state ss prediction and the membrane-state ss prediction
      linal::Vector3D membrane_prediction, ss_prediction;
      std::string topology_code;
      ISTREAM >> topology_code;

      ISTREAM >> membrane_prediction( 0);
      ISTREAM >> membrane_prediction( 1);
      ISTREAM >> membrane_prediction( 2);
      ISTREAM >> ss_prediction( 0); // helix
      ISTREAM >> ss_prediction( 1); // strand
      ISTREAM >> ss_prediction( 2); // coil

      // call the reader function to read the actual ssprediction values for this amino acid
      AMINO_ACID.SetSSPrediction( GetMethods().e_MASP, MASP( membrane_prediction, ss_prediction));

      // clear any remaining numbers from the line
      std::getline( ISTREAM, topology_code);

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &MASP::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence
      (
        GetMethods().e_MASP, AA_SEQUENCE, **GetMethods().e_MASP
      );

      util::ChopHeader( ISTREAM);

      return ReadStandardPredictionsForAASequence( ISTREAM, AA_SEQUENCE, GetMethods().e_MASP);
    }

    //! @brief Predict whether a given protein is a beta barrel and add MASP-BB environment to the protein
    //!        The secondary structure and environment for all residues that are not predicted to be beta-barrels will
    //!        not be changed
    //! @param MODEL ProteinModel for which JUFO will be calculated
    //! @return true if the protein is predicted to be an integral, bacterial or mitochondrial outer-membrane, beta-barrel protein
    bool MASP::PredictBetaBarrelEnvironmentSS( assemble::ProteinModel &MODEL)
    {
      // create a protein-model-with-cache
      assemble::ProteinModelWithCache pmwc( MODEL, false);
      return PredictBetaBarrelEnvironmentSS( pmwc);
    }

    //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
    //!        The secondary structure and environment for all residues that are not predicted to be beta-barrels will
    //!        not be changed
    //! @param MODEL ProteinModel for which MASP-BB will be calculated
    //! @return true if the protein is predicted to be an integral, bacterial or mitochondrial outer-membrane, beta-barrel protein
    bool MASP::PredictBetaBarrelEnvironmentSS( assemble::ProteinModelWithCache &MODEL)
    {
      // make sure the blast profile exists for this sequence
      BCL_Assert
      (
        MODEL.GetIterator()->GetBlastProfilePtr().IsDefined(),
        "Blast Profile is not available"
      );

      // create the descriptor to generate the raw prediction values
      util::Implementation< descriptor::Base< biol::AABase, float> > aa_descriptor
      (
        "PredictionMean(storage=File(prefix=model," + model::Model::AddModelPath( "masp/bb/") + "))"
      );

      // set the object up
      aa_descriptor->SetObject( MODEL);

      // set the dimension (1 because we operate on elements of the sequence)
      aa_descriptor->SetDimension( 1);

      // create a descriptor iterator
      descriptor::Iterator< biol::AABase> itr( MODEL.GetIterator());

      // iterate through the protein model to get the per-residue beta-barrel prediction
      storage::Vector< linal::Vector< float> > predictions;
      predictions.AllocateMemory( MODEL.GetSize());
      for( ; itr.NotAtEnd(); ++itr)
      {
        predictions.PushBack( aa_descriptor->operator ()( itr));
        predictions.LastElement() += float( 0.125);
        predictions.LastElement() /= float( 1.25);
      }

      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MASP::Read( std::istream &ISTREAM)
    {
      // read members
      ISTREAM >> m_Prediction;

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MASP::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
