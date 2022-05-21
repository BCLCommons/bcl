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
#include "sspred/bcl_sspred_boctopus.h"

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
    BOCTOPUS::BOCTOPUS() :
      m_Prediction( e_Outside)
    {
    }

    //! @brief constructor from the BOCTOPUS prediction
    //! @param PREDICTION prediction from BOCTOPUS file
    BOCTOPUS::BOCTOPUS( const char &PREDICTION) :
      m_Prediction( static_cast< Prediction>( PREDICTION))
    {
    }

    //! @brief Clone function
    //! @return pointer to new BOCTOPUS
    BOCTOPUS *BOCTOPUS::Clone() const
    {
      return new BOCTOPUS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BOCTOPUS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &BOCTOPUS::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".bocto_topo");

      // end
      return s_file_extension;

    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D BOCTOPUS::GetThreeStatePrediction() const
    {
      return m_Prediction == e_Membrane ? linal::Vector3D( 0, 1, 0) : linal::Vector3D( 0, 0, 1);
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> BOCTOPUS::GetNineStatePrediction() const
    {
      // if tm helix was predicted
      return ConvertThreeStateToNineState
             (
               GetThreeStatePrediction(),
               m_Prediction == e_Membrane
               ? biol::GetEnvironmentTypes().e_MembraneCore
               : biol::GetEnvironmentTypes().e_Solution
             );
    }

    //! @brief find the TMTypes with highest prediction and returns it
    //! @return TMTYpe with highest prediction
    biol::EnvironmentType BOCTOPUS::GetOneStateTMPrediction() const
    {
      return
        m_Prediction == e_Membrane
        ? biol::GetEnvironmentTypes().e_MembraneCore
          : m_Prediction == e_Inside
          ? biol::GetEnvironmentTypes().e_SolutionInside
          : biol::GetEnvironmentTypes().e_SolutionOutside;
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &BOCTOPUS::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      // issue warning
      BCL_MessageCrt( "BOCTOPUS ReadPredictionsForAA should not have been called");

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &BOCTOPUS::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // static ShPtr to Method interface to be used
      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence
      (
        GetMethods().e_BOCTOPUS, AA_SEQUENCE, **GetMethods().e_BOCTOPUS
      );

      // initialize necessary variables
      std::string line;
      bool read_flag( false);

      // create iterators on the sequence
      util::ShPtrVector< biol::AABase>::iterator aa_itr( AA_SEQUENCE.Begin());
      const util::ShPtrVector< biol::AABase>::iterator aa_itr_end( AA_SEQUENCE.End());

      // while reading lines
      while( std::getline( ISTREAM, line))
      {
        const std::string trimmed_line( util::TrimString( line));

        // if the line contains identifier to start reading
        if( util::EndsWith( trimmed_line, " predicted topology:") || util::StartsWith( line, ">"))
        {
          read_flag = true;
        }
        // line contains predictions
        else if( read_flag)
        {
          // iterate over the string
          for
          (
            std::string::const_iterator char_itr( trimmed_line.begin()),
              char_itr_end( trimmed_line.end());
            char_itr != char_itr_end && aa_itr != aa_itr_end; ++char_itr, ++aa_itr
          )
          {
            // get the prediction
            ( *aa_itr)->SetSSPrediction( GetMethods().e_BOCTOPUS, BOCTOPUS( *char_itr));
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
    std::istream &BOCTOPUS::Read( std::istream &ISTREAM)
    {
      // read members
      char prediction_char( ' ');
      ISTREAM >> prediction_char;
      m_Prediction = static_cast< Prediction>( prediction_char);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &BOCTOPUS::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << char( m_Prediction);

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
