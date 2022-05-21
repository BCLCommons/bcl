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
#include "sspred/bcl_sspred_kaksi.h"

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
    Kaksi::Kaksi() :
      m_Prediction( GetDefaultPredictionVector())
    {
    }

    //! @brief constructor from a linal::Vector3D
    //! @param VECTOR linal::Vector3D of probabilities
    Kaksi::Kaksi( const linal::Vector3D &VECTOR) :
      m_Prediction( VECTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Kaksi
    Kaksi *Kaksi::Clone() const
    {
      return new Kaksi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Kaksi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &Kaksi::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".kaksi");

      // end
      return s_file_extension;

    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D Kaksi::GetThreeStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> Kaksi::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( m_Prediction, biol::GetEnvironmentTypes().e_Solution);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Kaksi::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      // issue warning
      BCL_MessageCrt( "Kaksi ReadPredictionsForAA should not have been called");

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &Kaksi::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // static ShPtr to Method interface to be used
      static const Kaksi s_helix( linal::Vector3D( 1.0, 0.0, 0.0)), s_strand( linal::Vector3D( 0.0, 1.0, 0.0));
      static const Kaksi s_coil( linal::Vector3D( 0.0, 0.0, 1.0));
      static const Kaksi s_unknown( linal::Vector3D( util::GetUndefined< double>()));

      // initialize sequence predictions
      MethodHandler::InitializePredictionsForAASequence
      (
        GetMethods().e_Kaksi, AA_SEQUENCE, **GetMethods().e_Kaksi
      );

      // initialize necessary variables
      std::string line;

      // create iterators on the sequence
      util::ShPtrVector< biol::AABase>::iterator aa_itr( AA_SEQUENCE.Begin());
      const util::ShPtrVector< biol::AABase>::iterator aa_itr_end( AA_SEQUENCE.End());

      const char sequence_chain_id( AA_SEQUENCE.GetChainID());

      // track all chain ids seen in stride that are not == the sequence chain id
      std::string chain_ids_seen;

      // while reading lines
      bool in_right_chain( false);
      std::string kaksi_ss;
      kaksi_ss.reserve( AA_SEQUENCE.GetSize());
      while( std::getline( ISTREAM, line))
      {
        std::string trimmed_line( util::TrimString( line));

        if( trimmed_line.empty())
        {
          continue;
        }
        // check chain id
        if( trimmed_line[ 0] == '>')
        {
          if( trimmed_line.size() > 5)
          {
            if( trimmed_line[ 5] == sequence_chain_id)
            {
              in_right_chain = true;
            }
            else
            {
              chain_ids_seen += trimmed_line[ 5];
              // if we were previously in the right chain, then we are now on another chain, so stop looking for more
              // more residues
              if( in_right_chain)
              {
                break;
              }
              in_right_chain = false;
            }
          }
          continue;
        }
        if( !in_right_chain)
        {
          continue;
        }

        kaksi_ss += trimmed_line;
      }

      // determine offset from kaksi into aa vector
      if( kaksi_ss.size() < AA_SEQUENCE.GetSize() && !kaksi_ss.empty() && kaksi_ss[ 0] != 'X')
      {
        if
        (
          AA_SEQUENCE.GetFirstAA()->GetType() == biol::GetAATypes().XXX
          || AA_SEQUENCE.GetFirstAA()->GetType() == biol::GetAATypes().UNK
        )
        {
          kaksi_ss = 'X' + kaksi_ss;
        }
      }
      if( kaksi_ss.size() < AA_SEQUENCE.GetSize() && !kaksi_ss.empty() && kaksi_ss[ kaksi_ss.size() - 1] != 'X')
      {
        if
        (
          AA_SEQUENCE.GetLastAA()->GetType() == biol::GetAATypes().XXX
          || AA_SEQUENCE.GetLastAA()->GetType() == biol::GetAATypes().UNK
        )
        {
          kaksi_ss += 'X';
        }
      }
      if( kaksi_ss.size() > AA_SEQUENCE.GetSize() && kaksi_ss[ 0] != 'X' && kaksi_ss[ kaksi_ss.size() - 1] == 'X')
      {
        kaksi_ss.erase( kaksi_ss.size() - 1);
      }
      if( kaksi_ss.size() > AA_SEQUENCE.GetSize() && kaksi_ss[ 0] == 'X' && kaksi_ss[ kaksi_ss.size() - 1] != 'X')
      {
        kaksi_ss.erase( 0, 1);
      }

      if( kaksi_ss.size() != AA_SEQUENCE.GetSize())
      {
        BCL_MessageCrt
        (
          "Kaksi file contained different # of residues for chain " + util::Format()( sequence_chain_id)
          + ". " + util::Format()( kaksi_ss.size()) + " residues found; sequence "
          + AA_SEQUENCE.GetSequenceIdentification() + " had: "
          + util::Format()( AA_SEQUENCE.GetSize())
        );
      }
      // iterate over the string
      for
      (
        std::string::const_iterator char_itr( kaksi_ss.begin()), char_itr_end( kaksi_ss.end());
        char_itr != char_itr_end && aa_itr != aa_itr_end;
        ++char_itr, ++aa_itr
      )
      {
        // if this is a transmembrane helix
        switch( *char_itr)
        {
          case 'H': ( *aa_itr)->SetSSPrediction( GetMethods().e_Kaksi, s_helix); break;
          case 'b': ( *aa_itr)->SetSSPrediction( GetMethods().e_Kaksi, s_strand); break;
          case '.': ( *aa_itr)->SetSSPrediction( GetMethods().e_Kaksi, s_coil); break;
          case 'X':
          default: ( *aa_itr)->SetSSPrediction( GetMethods().e_Kaksi, s_unknown); break;
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
    std::istream &Kaksi::Read( std::istream &ISTREAM)
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
    std::ostream &Kaksi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
