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

#ifndef BCL_SCORE_PROTEIN_MODEL_SSE_CHIRALITY_H_
#define BCL_SCORE_PROTEIN_MODEL_SSE_CHIRALITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "scorestat/bcl_scorestat_protein_model_sse_triplet_chirality.h"
#include "util/bcl_util_sh_ptr.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelSSEChirality
    //! @brief score for proximity of imaginary linear loops between all neighborings sses to other sses in the model
    //!
    //! @see @link example_score_protein_model_sse_chirality.cpp @endlink
    //! @author mendenjl
    //! @date Mar 15, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelSSEChirality :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      double m_Pseudocount;                //!< Pseudocount to be applied to the raw counts

      bool   m_UseContactSpecificWeights;  //!< Use contact specific weights

      scorestat::ProteinModelSSETripletChirality m_HashMaker; //!< Object to make the hashes

      linal::Vector< double> m_Propensities; //!< Propensities for each possible type of interaction

      bool m_ScaleByAAs; //!< whether to scale by the number of aas in sses in the protein

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      //! @param CONSIDER_NON_ADJACENT whether to consider non-adjacent SSEs as well - this makes the score quadratically
      //!        slower to compute
      //! @param PSEUDOCOUNT pseudocount to use
      //! @param CONTACT_SPECIFIC true if log-odds should be computed specific for each contact type
      //! @param SCALE_BY_AAS true to scale by the number of AAs in the protein
      ProteinModelSSEChirality
      (
        const double &PSEUDOCOUNT = 20.0,
        const bool   &CONTACT_SPECIFIC = false,
        const bool   &SCALE_BY_AAS = true
      );

      //! virtual copy constructor
      ProteinModelSSEChirality *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the propensities vector
      const linal::Vector< double> &GetPropensities() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that scores the chain
      //! @param CHAIN the chain for which all neighbor scores are calculated
      //! @param MISSING_SSES sse pool of missing sses
      //! @return score
      double operator()
      (
        const assemble::Chain &CHAIN,
        const assemble::SSEPool &MISSING_SSES
      ) const;

      //! @brief operator that scores the Protein model
      //! @param PROTEIN_MODEL the protein model for which all neighbor scores are calculated
      //! @return score
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        std::ostream &OSTREAM
      ) const;

    protected:

      //! @brief helper function called by WriteDetailedSchemeAndValues and operator() so that the code remains in-sync
      //! @param CHAIN the chain of interest
      //! @param MISSING_SSES sse pool of missing sses
      //! @param OSTREAM the output stream to write the detailed scheme to for this chain
      //! @param DO_WRITE set to true to actually write to the output stream; otherwise, nothing will be written
      //! @return the final score
      double ScoreLoops
      (
        const assemble::Chain &CHAIN,
        const assemble::SSEPool &MISSING_SSES,
        std::ostream &OSTREAM,
        const bool &DO_WRITE
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @return ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class ProteinModelSSEChirality

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_PROTEIN_MODEL_SSE_CHIRALITY_H_
