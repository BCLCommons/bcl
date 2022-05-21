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

#ifndef BCL_BIOL_BLAST_PROFILE_H_
#define BCL_BIOL_BLAST_PROFILE_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_types.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BlastProfile
    //! @brief This class stores Blast information relevant to a single amino acid
    //! @details  This class stores the 20-value blast profile as the 20-value probability vector for a single amino
    //! acid as well as providing the Get/Set and Read/Write functionalities
    //!
    //! @see @link example_biol_blast_profile.cpp @endlink
    //! @author karakam
    //! @date 07/02/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BlastProfile :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! Position based scoring matrix
      linal::Vector< double> m_Profile;

      //! Position based scoring matrix converted to probabilities
      linal::Vector< double> m_Probabilities;

      //! Relative weight of gapless real matches to pseudocounts
      double m_MatchWeightRelativeToPseudocount;

      //! Information per position
      double m_Information;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BlastProfile();

      //! @brief constructor from a profile vector
      //! @param PROFILE blast profile vector
      BlastProfile
      (
        const linal::Vector< double> &PROFILE,
        const double &MATCH_WEIGHT = util::GetUndefined< double>(),
        const double &POSITION_INFORMATION = util::GetUndefined< double>()
      );

      //! @brief virtual copy constructor
      BlastProfile *Clone() const;

    ///////////////////////
    // data access - Get //
    ///////////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns Profile Vector for 20 standard AATypes
      //! @return Profile Vector for 20 standard AATypes
      const linal::Vector< double> &GetProfile() const
      {
        return m_Profile;
      }

      //! @brief returns the Profile value for the given AAType
      //! @return the Profile value for the given AAType
      const double &GetProfile( const AAType &TYPE) const
      {
        return m_Profile( TYPE.GetIndex());
      }

      //! @brief return BlastProfile as probabilities for 20 standard AATypes
      //! @return BlastProfile as probabilities for 20 standard AATypes
      const linal::Vector< double> &GetProbabilities() const
      {
        return m_Probabilities;
      }

      //! @brief return weight of alignment relative to pseudocount (0-1); higher values indicate more alignments at this position
      //! @return weight of alignment relative to pseudocount (0-1); higher values indicate more alignments at this position
      double GetAlignmentWeight() const
      {
        return m_MatchWeightRelativeToPseudocount;
      }

      //! @brief return information per position
      //! @return information per position
      double GetPositionInformation() const
      {
        return m_Information;
      }

    ///////////////////////
    // data access - Set //
    ///////////////////////

      //! @brief set BlastProfile to provided PROFILE
      //! @param PROFILE blast profile vector
      void SetProfile( const linal::Vector< double> &PROFILE);

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the conservation, sum of profile times probability
      //! @return the conservation
      double CalculateConservation() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief reads BlastProfile from ISTREAM
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadProfile( std::istream &ISTREAM);

      //! @brief writes BlastProfile to OSTREAM
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &WriteProfile( std::ostream &OSTREAM) const;

    protected:

      //! @brief reads blast
      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

      //! @brief set BlastProfile probabilities
      void SetProbabilities();

    }; //class BlastProfile

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_BLAST_PROFILE_H_
