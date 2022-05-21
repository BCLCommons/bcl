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
#include "sspred/bcl_sspred_conpred.h"

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
    CONPRED::CONPRED() :
      m_Prediction( GetDefaultPredictionVector())
    {
    }

    //! @brief constructor from a linal::Vector3D
    //! @param VECTOR linal::Vector3D of probabilities
    CONPRED::CONPRED( const linal::Vector3D &VECTOR) :
      m_Prediction( VECTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CONPRED
    CONPRED *CONPRED::Clone() const
    {
      return new CONPRED( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CONPRED::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &CONPRED::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( "_conpred.txt");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D CONPRED::GetThreeStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> CONPRED::GetNineStatePrediction() const
    {
      // if tm helix was predicted
      if( m_Prediction( biol::GetSSTypes().HELIX) == 1)
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
    std::istream &CONPRED::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      BCL_MessageCrt( "ReadPredictionsForAA function should not been called for CONPRED");

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &CONPRED::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      // static ShPtr to Method interface to be used
      static CONPRED s_default_helix_prediction( linal::Vector3D( 1.0, 0.0, 0.0));

      // skip header lines
      util::ChopHeader( ISTREAM);

      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence
      (
        GetMethods().e_CONPRED, AA_SEQUENCE, **GetMethods().e_CONPRED
      );

      // initialize necessary variables
      size_t start, end;
      std::string line;
      storage::Vector< std::string> sub_strings;

      // indicator if right location has been found
      bool right_location( false);

      // while reading lines
      while( std::getline( ISTREAM, line))
      {
        // update the boolean after an empty line is read after it is set
        if( right_location && util::TrimString( line).empty())
        {
          right_location = false;
        }
        if( right_location)
        {
          // split strings
          sub_strings = util::SplitString( line);

          // assert the the format is correct
          BCL_Assert
          (
            sub_strings.GetSize() >= 4,
            "not enough elements in conpred line: " + util::Format()( sub_strings.GetSize()) + " <= 4"
          );

          // update start and end
          start = util::ConvertStringToNumericalValue< size_t>( sub_strings( 1));
          end   = util::ConvertStringToNumericalValue< size_t>( sub_strings( 3));

          // save prediction
          MethodHandler::SetPredictionsForSubSequence
          (
            GetMethods().e_CONPRED, AA_SEQUENCE, s_default_helix_prediction, start - 1, end - 1
          );
        }
        // if begin can not be found
        if( line.find( "begin") != std::string::npos)
        {
          right_location = true;
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
    std::istream &CONPRED::Read( std::istream &ISTREAM)
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
    std::ostream &CONPRED::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // end
      return OSTREAM;

    };

  } // namespace sspred
} // namespace bcl
