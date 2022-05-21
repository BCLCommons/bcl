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
#include "sspred/bcl_sspred_proftmb.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PROFTMB::PROFTMB() :
      m_Prediction( GetDefaultPredictionMatrix())
    {
    }

    //! @brief constructor from a linal::Matrix
    //! @param MATRIX linal::Matrix of probabilities
    PROFTMB::PROFTMB( const linal::Matrix< double> &MATRIX) :
      m_Prediction( MATRIX)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PROFTMB
    PROFTMB *PROFTMB::Clone() const
    {
      return new PROFTMB( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &PROFTMB::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &PROFTMB::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( "xxx.PROFTMB");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D PROFTMB::GetThreeStatePrediction() const
    {
      return ConvertNineStateToThreeState( m_Prediction);
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> PROFTMB::GetNineStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &PROFTMB::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      // initialize predictions vector
      linal::Matrix< double> predictions
      (
        biol::GetEnvironmentTypes().GetNumberReducedTypes(), biol::GetSSTypes().COIL.GetIndex() + 1, double( 0.0)
      );

      // initialize two doubles and a string
      double first, second;
      std::string tmp;

      // read the predictions for STRAND
      ISTREAM >> tmp >> tmp >> first >> second;
      predictions( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex(), biol::GetSSTypes().STRAND) = first + second;

      // read the predictions for COIL
      ISTREAM >> tmp >> tmp >> first >> second;
      predictions( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex(), biol::GetSSTypes().COIL) = first + second;

      // normalize predictions to 1
      predictions.AsVector().SetToSum( 1.0);

      // set the predictions for this amino acid
      AMINO_ACID.SetSSPrediction( GetMethods().e_PROFTMB, PROFTMB( predictions));

      // return
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &PROFTMB::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      // use read standard function for given sequence
      ReadStandardPredictionsForAASequence( ISTREAM, AA_SEQUENCE, GetMethods().e_PROFTMB);

      // end
      return ISTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PROFTMB::Read( std::istream &ISTREAM)
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
    std::ostream &PROFTMB::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // end
      return OSTREAM;

    };

  } // namespace sspred
} // namespace bcl
