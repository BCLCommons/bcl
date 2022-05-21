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
#include "sspred/bcl_sspred_talos.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> TALOS::s_Instance
    (
      GetObjectInstances().AddInstance( new TALOS())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    TALOS::TALOS() :
      m_Prediction( GetDefaultPredictionVector())
    {
    }

    //! @brief constructor from a linal::Vector3D
    //! @param VECTOR linal::Vector3D of probabilities
    TALOS::TALOS( const linal::Vector3D &VECTOR) :
      m_Prediction( VECTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new TALOS
    TALOS *TALOS::Clone() const
    {
      return new TALOS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TALOS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &TALOS::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( "SS.tab");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D TALOS::GetThreeStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> TALOS::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( m_Prediction, biol::GetEnvironmentTypes().e_Solution);
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &TALOS::ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const
    {
      // initialize vector3d
      linal::Vector3D prediction;

      // temporary variable to read in unused information
      std::string temp;

      // read in unused information
      ISTREAM >> temp >> temp;

      // read COIL, HELIX and STRAND (read in to temp so that nan's can be appropriately handled)
      ISTREAM >> temp;
      prediction( biol::GetSSTypes().HELIX) = util::ConvertStringToNumericalValue< double>( temp);
      ISTREAM >> temp;
      prediction( biol::GetSSTypes().STRAND) = util::ConvertStringToNumericalValue< double>( temp);
      ISTREAM >> temp;
      prediction( biol::GetSSTypes().COIL) = util::ConvertStringToNumericalValue< double>( temp);

      // read rest of line
      std::getline( ISTREAM, temp);

      // TALOS predictions can be (0.0, 0.0, 0.0) or nan
      if( !prediction.IsDefined() || ( prediction.X() == 0.0 && prediction.Y() == 0.0 && prediction.Z() == 0.0))
      {
        prediction( biol::GetSSTypes().HELIX) = 0.0;
        prediction( biol::GetSSTypes().STRAND) = 0.0;
        prediction( biol::GetSSTypes().COIL) = 1.0;
      }

      // normalize the vector
      prediction.SetToSum();

      // set the predictions for this amino acid
      AMINO_ACID.SetSSPrediction( GetMethods().e_TALOS, TALOS( prediction));

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &TALOS::ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const
    {
      // skip over header
      std::string this_line;
      while( !ISTREAM.eof())
      {
        std::getline( ISTREAM, this_line);
        const storage::Vector< std::string> split_line( util::SplitString( this_line, " "));
        const std::string first_entry( split_line.IsEmpty() ? "" : split_line.FirstElement());
        if( first_entry == "FORMAT")
        {
          break;
        }
      }
      std::getline( ISTREAM, this_line);

      // call standard read function and return it
      return ReadStandardPredictionsForAASequence( ISTREAM, AA_SEQUENCE, GetMethods().e_TALOS);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TALOS::Read( std::istream &ISTREAM)
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
    std::ostream &TALOS::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred

} // namespace bcl
