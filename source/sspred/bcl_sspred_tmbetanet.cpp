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
#include "sspred/bcl_sspred_tmbetanet.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
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
    TMBETANET::TMBETANET() :
      m_Prediction( GetDefaultPredictionVector())
    {
    }

    //! @brief constructor from a linal::Vector3D
    //! @param VECTOR linal::Vector3D of probabilities
    TMBETANET::TMBETANET( const linal::Vector3D &VECTOR) :
      m_Prediction( VECTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new TMBETANET
    TMBETANET *TMBETANET::Clone() const
    {
      return new TMBETANET( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TMBETANET::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &TMBETANET::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".tmbetanet");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D TMBETANET::GetThreeStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> TMBETANET::GetNineStatePrediction() const
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
    std::istream &TMBETANET::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      // issue warning
      BCL_MessageCrt( "TMBETANET ReadPredictionsForAA should not have been called");

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &TMBETANET::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      // static ShPtr to Method interface to be used
      static const TMBETANET s_default_strand_prediction( linal::Vector3D( 0.0, 1.0, 0.0));

      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence
      (
        GetMethods().e_TMBETANET, AA_SEQUENCE, **GetMethods().e_TMBETANET
      );

      // initialize necessary variables
      std::string line;

      // while reading lines
      while( std::getline( ISTREAM, line))
      {
        if( line.empty())
        {
          // skip blank lines
          continue;
        }
        const storage::Vector< std::string> split_line( util::SplitString( line));

        if( split_line.GetSize() == size_t( 1) && split_line( 0) == "tm")
        {
          // tm header line
          continue;
        }

        if( split_line.GetSize() != size_t( 3))
        {
          BCL_MessageStd
          (
            "Skipping TMBetaNet prediction line: " + line + " because it did not have 3 entries "
          );
          continue;
        }

        const std::string &ss_type( split_line( 0));

        // check format
        if( ss_type == "tm")
        {
          // read start and end
          const int start( util::ConvertStringToNumericalValue< int>( split_line( 1)));
          const int end( util::ConvertStringToNumericalValue< int>( split_line( 2)));
          if( size_t( end) <= AA_SEQUENCE.GetSize())
          {
            MethodHandler::SetPredictionsForSubSequence
            (
              GetMethods().e_TMBETANET, AA_SEQUENCE, s_default_strand_prediction, start - 1, end - 1
            );
          }
          else
          {
            BCL_MessageCrt
            (
              "Warning, ignoring TMBetaNet prediction line out of range for this chain: " + line
            );
          }
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
    std::istream &TMBETANET::Read( std::istream &ISTREAM)
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
    std::ostream &TMBETANET::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
