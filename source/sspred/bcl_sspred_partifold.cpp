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
#include "sspred/bcl_sspred_partifold.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PARTIFOLD::PARTIFOLD() :
      m_Prediction( GetDefaultPredictionVector())
    {
    }

    //! @brief constructor from a linal::Vector3D
    //! @param VECTOR linal::Vector3D of probabilities
    PARTIFOLD::PARTIFOLD( const linal::Vector3D &VECTOR) :
      m_Prediction( VECTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PARTIFOLD
    PARTIFOLD *PARTIFOLD::Clone() const
    {
      return new PARTIFOLD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PARTIFOLD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &PARTIFOLD::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".mfe");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D PARTIFOLD::GetThreeStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> PARTIFOLD::GetNineStatePrediction() const
    {
      // if tm strand was predicted
      if( m_Prediction( biol::GetSSTypes().STRAND) == 1.0)
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
    std::istream &PARTIFOLD::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      // initialize vector3d
      linal::Vector3D prediction;

      // initialize integer and string to store seqid and amino acid temporarily
      std::string one_state, contact_a, contact_b;

      // read in seqid, residue, and one_state
      ISTREAM >> one_state >> contact_a >> contact_b;

      // TM beta-barrel prediction states "M" for membrane and "C" for channel indicating side-chain orientation
      if( one_state == "M" || one_state == "C")
      {
        prediction( biol::GetSSTypes().STRAND) = 1.0;
      }
      else
      {
        prediction( biol::GetSSTypes().COIL) = 1.0;
      }

      // set the predictions for this amino acid
      AMINO_ACID.SetSSPrediction( GetMethods().e_PARTIFOLD, PARTIFOLD( prediction));

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &PARTIFOLD::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      // call standard read function and return it
      return ReadStandardPredictionsForAASequence( ISTREAM, AA_SEQUENCE, GetMethods().e_PARTIFOLD);

      // end
      return ISTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PARTIFOLD::Read( std::istream &ISTREAM)
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
    std::ostream &PARTIFOLD::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
