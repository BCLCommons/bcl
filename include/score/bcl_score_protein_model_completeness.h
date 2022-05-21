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

#ifndef BCL_SCORE_PROTEIN_MODEL_COMPLETENESS_H_
#define BCL_SCORE_PROTEIN_MODEL_COMPLETENESS_H_

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
    //! @class ProteinModelCompleteness
    //! @brief This class scores the completeness of a protein model
    //! @detail The completeness of the protein model is scored by computing the ratio of residues with defined
    //! coordinates and the total number of residues in the sequence.
    //!
    //! @see @link example_score_protein_model_completeness.cpp @endlink
    //! @author fischea
    //! @date Dec 8, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ProteinModelCompleteness :
      public ProteinModel
    {

    //////////
    // data //
    //////////

    private:

      //! ignore terminal loops for the calculation of the completeness
      bool m_IgnoreTermLoops;

      //! scheme of this score
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from default values
      //! @param IGNORE_TERM_LOOPS ignore terminal loops for the calculation of the completeness
      //! @param SCHEME scheme of the score
      ProteinModelCompleteness( bool IGNORE_TERM_LOOPS = false, const std::string &SCHEME = GetDefaultScheme());

      //! @brief returns a pointer to a new ProteinModelCompleteness
      //! @return pointer to a new ProteinModelCompleteness
      ProteinModelCompleteness *Clone() const;

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

      //! @brief scores the completeness of the given protein model
      //! @detail the completeness is calculated by computing the ratio of residues with defined coordinates and
      //! the total number of residues in the sequence.
      //! @param PROTEIN_MODEL protein model for which to compute the completeness
      //! @return completeness of the given protein model with -1.0 being complete and 0.0 being empty
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief computes the number of residues in the given chain
      //! @param CHAIN chain for which to compute the number of residues
      //! @return number of residues with in the given chain
      size_t GetNumResidues( const assemble::Chain &CHAIN) const;

      //! @brief computes the number of residues in the given protein model
      //! @param MODEL protein model for which to compute the number of residues
      //! @return number of residues in the given protein model
      size_t GetNumResidues( const assemble::ProteinModel &MODEL) const;

      //! @brief computes the number of residues with defined backbone coordinates in a given SSE
      //! @param SSE SSE for which to compute the number of residues with defined backbone coordinates
      //! @return number of residues with defined backbone coordinates in the given SSE
      size_t GetNumResiduesDefined( const assemble::SSE &SSE) const;

      //! @brief computes the number of residues with defined backbone coordinates in a given chain
      //! @param CHAIN chain for which to compute the number of residues with defined backbone coordinates
      //! @return number of residues with defined backbone coordinates in the given chain
      size_t GetNumResiduesDefined( const assemble::Chain &CHAIN) const;

      //! @brief computes the number of residues with defined backbone coordinates in a given protein model
      //! @param MODEL protein model for which to compute the number of residues with defined backbone coordinates
      //! @return number of residues with defined backbone coordinates in the given protein model
      size_t GetNumResiduesDefined( const assemble::ProteinModel &MODEL) const;

    }; // class ProteinModelCompleteness

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_COMPLETENESS_H_
