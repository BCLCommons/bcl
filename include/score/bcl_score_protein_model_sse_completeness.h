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

#ifndef BCL_SCORE_PROTEIN_MODEL_SSE_COMPLETENESS_H_
#define BCL_SCORE_PROTEIN_MODEL_SSE_COMPLETENESS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "biol/bcl_biol_aa_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelSSECompleteness
    //! @brief This class scores the completeness of a protein model
    //! @detail The completeness of the protein model is scored by computing the ratio of residues with defined
    //! coordinates and the total number of residues in the sequence.
    //!
    //! @see @link example_score_protein_model_sse_completeness.cpp @endlink
    //! @author mendenjl
    //! @date Jan 23, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ProteinModelSSECompleteness :
      public ProteinModel
    {

    //////////
    // data //
    //////////

      //! bool, whether to return # of aas missing rather than # of sses missing
      bool m_CountAAs;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_AAInstance;
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief whether to count missing aas or missing sses (default)
      ProteinModelSSECompleteness( const bool &COUNT_AAS = false) :
        m_CountAAs( COUNT_AAS)
      {
      }

      //! @brief returns a pointer to a new ProteinModelSSECompleteness
      //! @return pointer to a new ProteinModelSSECompleteness
      ProteinModelSSECompleteness *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the scheme of this score
      //! @return the scheme of this score
      const std::string &GetScheme() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief scores the completeness of the given protein model
      //! @detail the completeness is calculated by computing the ratio of residues with defined coordinates and
      //! the total number of residues in the sequence.
      //! @param PROTEIN_MODEL protein model for which to compute the completeness
      //! @return completeness of the given protein model with -1.0 being complete and 0.0 being empty
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    }; // class ProteinModelSSECompleteness

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_SSE_COMPLETENESS_H_
