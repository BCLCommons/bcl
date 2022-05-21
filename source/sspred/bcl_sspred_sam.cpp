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
#include "sspred/bcl_sspred_sam.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
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
    SAM::SAM() :
      m_Prediction( GetDefaultPredictionVector())
    {
    }

    //! @brief constructor from a linal::Vector3D
    //! @param VECTOR linal::Vector3D of probabilities
    SAM::SAM( const linal::Vector3D &VECTOR) :
      m_Prediction( VECTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SAM
    SAM *SAM::Clone() const
    {
      return new SAM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SAM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &SAM::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".rdb6");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D SAM::GetThreeStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> SAM::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( m_Prediction, biol::GetEnvironmentTypes().e_Solution);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &SAM::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      // initialize two temporary doubles since sam values come in pairs
      double first, second;

      // initialize vector for storing predictions
      linal::Vector3D prediction;

      // read and update STRAND values
      ISTREAM >> first >> second;
      prediction( biol::GetSSTypes().STRAND) = first + second;

      // read and update HELIX values
      ISTREAM >> first >> second;
      prediction( biol::GetSSTypes().HELIX) = first + second;

      // read and update COIL values
      ISTREAM >> first >> second;
      prediction( biol::GetSSTypes().COIL) = first + second;

      // normalize predictions to 1
      prediction.SetToSum();

      // set the predictions for this method
      AMINO_ACID.SetSSPrediction( GetMethods().e_SAM, SAM( prediction));

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &SAM::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      // initialize temporary string variable
      std::string line;

      // remove header lines
      util::ChopHeader( ISTREAM);

      // remove two empty lines
      std::getline( ISTREAM, line);
      std::getline( ISTREAM, line);

      // call standard read function and return it
      return ReadStandardPredictionsForAASequence( ISTREAM, AA_SEQUENCE, GetMethods().e_SAM);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SAM::Read( std::istream &ISTREAM)
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
    std::ostream &SAM::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // end
      return OSTREAM;

    };

  } // namespace sspred
} // namespace bcl
