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
#include "sspred/bcl_sspred_b2tmpred.h"

// includes from bcl - sorted alphabetically
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
    B2TMPRED::B2TMPRED() :
      m_Prediction( GetDefaultPredictionVector())
    {
    }

    //! @brief constructor from a linal::Vector3D
    //! @param VECTOR linal::Vector3D of probabilities
    B2TMPRED::B2TMPRED( const linal::Vector3D &VECTOR) :
      m_Prediction( VECTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new B2TMPRED
    B2TMPRED *B2TMPRED::Clone() const
    {
      return new B2TMPRED( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &B2TMPRED::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &B2TMPRED::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".hmm");

      // end
      return s_file_extension;
    }

    ////////////////
    // operations //
    ////////////////1.00

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D B2TMPRED::GetThreeStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> B2TMPRED::GetNineStatePrediction() const
    {
      // if tm strand was predicted
      if( m_Prediction( biol::GetSSTypes().STRAND) == 1)
      {
        return ConvertThreeStateToNineState( m_Prediction, biol::GetEnvironmentTypes().e_MembraneCore);
      }
      // otherwise
      else
      {
        return ConvertThreeStateToNineState( m_Prediction, biol::GetEnvironmentTypes().e_Solution);
      }
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &B2TMPRED::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      // issue warning
      BCL_MessageCrt( "B2TMPRED ReadPredictionsForAA should not have been called");

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &B2TMPRED::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      // static ShPtr to Method interface to be used
      static B2TMPRED s_default_strand_prediction( linal::Vector3D( 0.0, 1.0, 0.0));

      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence
      (
        GetMethods().e_B2TMPRED, AA_SEQUENCE, **GetMethods().e_B2TMPRED
      );

      // initialize necessary variables
      size_t start, end;
      storage::Vector< std::string> sub_strings;
      std::string line;

      // while reading lines
      while( std::getline( ISTREAM, line))
      {
        // if the line contains transmembrane
        if( line.find( "Transmembrane") != std::string::npos)
        {
          sub_strings = util::SplitString( util::SplitString( line, " ")( 3), ":");
          start = util::ConvertStringToNumericalValue< size_t>( sub_strings( 0));
          end   = util::ConvertStringToNumericalValue< size_t>( sub_strings( 1));

          // set the subsequence predictions
          MethodHandler::SetPredictionsForSubSequence( GetMethods().e_B2TMPRED, AA_SEQUENCE, s_default_strand_prediction, start - 1, end - 1);
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
    std::istream &B2TMPRED::Read( std::istream &ISTREAM)
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
    std::ostream &B2TMPRED::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // end
      return OSTREAM;

    };

  } // namespace sspred
} // namespace bcl
