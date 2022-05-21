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

#ifndef BCL_SCORE_PROTEIN_MODEL_GAP_H_
#define BCL_SCORE_PROTEIN_MODEL_GAP_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_loop.h"
#include "bcl_score_protein_model.h"
#include "biol/bcl_biol_aa_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelGap
    //! @brief Score protein models for structural gaps.
    //! @details This class scores protein models for structural gaps. This class is intended for usage during loop
    //! construction to enable fragment-based construction.
    //!
    //! @see @link example_score_protein_model_gap.cpp @endlink
    //! @author fischea
    //! @date Feb 27, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ProteinModelGap :
      public ProteinModel
    {

    //////////
    // data //
    //////////

    private:

      //! scheme of this score
      std::string m_Scheme;

      //! loop score for scoring gaps
      Loop m_Score;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from default values
      //! @param SCHEME scheme of the score
      ProteinModelGap( const std::string &SCHEME = GetDefaultScheme());

      //! @brief returns a pointer to a new ProteinModelGap
      //! @return pointer to a new ProteinModelGap
      ProteinModelGap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the scheme of this score
      //! @return the scheme of this score
      const std::string &GetScheme() const;

      //! @brief returns the default scheme of this score
      //! @return the default scheme of this score
      static const std::string &GetDefaultScheme();

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Score a protein model regarding gaps.
      //! @param PROTEIN_MODEL protein model for which to compute the gap score
      //! @return gap score of the protein model
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief Compute the gap score for the given sequence.
      //! @param SEQUENCE Sequence for which to compute the gap score.
      //! @return Gap score for the given sequence.
      double ComputeGapScore( const biol::AASequence &SEQUENCE) const;

      //! @brief Compute the gap score for the given endpoints of the gap.
      //! @param AA_N n-terminal endpoint of the gap
      //! @param AA_C c-terminal endpoint of the gap
      //! @return Gap score for the given sequence.
      double ComputeGapScore( const biol::AABase &AA_N, const biol::AABase &AA_C) const;

    }; // class ProteinModelGap

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_GAP_H_
