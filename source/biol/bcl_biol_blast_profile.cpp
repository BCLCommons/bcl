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
#include "biol/bcl_biol_blast_profile.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BlastProfile::BlastProfile() :
      m_Profile( AATypes::s_NumberStandardAATypes, util::GetUndefined< double>()),
      m_Probabilities( AATypes::s_NumberStandardAATypes, util::GetUndefined< double>()),
      m_MatchWeightRelativeToPseudocount( util::GetUndefined< double>()),
      m_Information( util::GetUndefined< double>())
    {
    }
    //! @brief constructor from a profile vector
    //! @param PROFILE blast profile vector
    BlastProfile::BlastProfile
    (
      const linal::Vector< double> &PROFILE,
      const double &MATCH_WEIGHT,
      const double &POSITION_INFORMATION
    ) :
      m_Profile( AATypes::s_NumberStandardAATypes, util::GetUndefined< double>()),
      m_Probabilities( AATypes::s_NumberStandardAATypes, util::GetUndefined< double>()),
      m_MatchWeightRelativeToPseudocount( MATCH_WEIGHT),
      m_Information( POSITION_INFORMATION)
    {
      SetProfile( PROFILE);
    }

    //! @brief virtual copy constructor
    BlastProfile *BlastProfile::Clone() const
    {
      return new BlastProfile( *this);
    }

  ///////////////////////
  // data access - Get //
  ///////////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BlastProfile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////////////
  // data access - Set //
  ///////////////////////

    //! @brief set BlastProfile to provided PROFILE
    //! @param PROFILE blast profile vector
    void BlastProfile::SetProfile( const linal::Vector< double> &PROFILE)
    {
      // assert the given PROFILE has correct size
      BCL_Assert
      (
        PROFILE.GetSize() == AATypes::s_NumberStandardAATypes,
        "The Blast profile vector passed should have a size of "
        + util::Format().W( 5)( size_t( AATypes::s_NumberStandardAATypes))
        + " but instead it is " + util::Format().W( 5)( PROFILE.GetSize()) + "!!"
      );

      // reset profile
      m_Profile = 0;

      // update profile
      m_Profile = PROFILE;

      // update probabilities
      SetProbabilities();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the conservation, sum of profile times probability
    //! @return the conservation
    double BlastProfile::CalculateConservation() const
    {
      // return the scalar product
      return linal::ScalarProduct( m_Profile, m_Probabilities);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads BlastProfile from ISTREAM
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BlastProfile::ReadProfile( std::istream &ISTREAM)
    {
      // iterate over profile
      for
      (
        double *ptr( m_Profile.Begin()), *const ptr_end( m_Profile.End());
        ptr != ptr_end;
        ++ptr
      )
      {
        // read value from stream
        ISTREAM >> *ptr;
      }

      // try reading in the probabilities, information per position, and gapless weight
      // some blast implementations (like deltablast) do not provide this information, so set probabilities up from
      // blast profile if the probabilities are not availability.
      std::string rest_of_line;
      std::getline( ISTREAM, rest_of_line);
      storage::Vector< double> remaining_numbers( util::SplitStringToNumerical< double>( rest_of_line));

      // four possibilities: blast probabilities given, information gain and match weight given, both, or neither
      if( ( remaining_numbers.GetSize() % 20) == size_t( 2))
      {
        // information and alignment weight were given, read them in
        m_Information = remaining_numbers( remaining_numbers.GetSize() - 2);
        m_MatchWeightRelativeToPseudocount = remaining_numbers( remaining_numbers.GetSize() - 1);
      }

      // update probabilities
      SetProbabilities();

      // return
      return ISTREAM;
    }

    //! @brief writes BlastProfile to OSTREAM
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &BlastProfile::WriteProfile( std::ostream &OSTREAM) const
    {
      // iterate over profile
      util::Format format;
      format.W( 7).FFP( 2);
      for
      (
        const double *ptr( m_Profile.Begin()), *const ptr_end( m_Profile.End());
        ptr != ptr_end;
        ++ptr
      )
      {
        // write value to stream
        OSTREAM << format( *ptr);
      }
      OSTREAM << format( m_Information) << format( m_MatchWeightRelativeToPseudocount);
      // return
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BlastProfile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Profile, ISTREAM);
      io::Serialize::Read( m_Probabilities, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &BlastProfile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // read members
      io::Serialize::Write( m_Profile, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Probabilities, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

    //! @brief set BlastProfile probabilities
    void BlastProfile::SetProbabilities()
    {
      // calculate probabilities by taking 10 ^ ( m_profile/10)
      m_Probabilities = double( 10) ^ ( m_Profile / double( 10));

      // normalize (typically by 20)
      m_Probabilities /= m_Probabilities.Sum();
    }

  } // namespace biol
} // namespace bcl
