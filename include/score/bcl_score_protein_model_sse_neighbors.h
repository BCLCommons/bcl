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

#ifndef BCL_SCORE_PROTEIN_MODEL_SSE_NEIGHBORS_H_
#define BCL_SCORE_PROTEIN_MODEL_SSE_NEIGHBORS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelSSENeighbors
    //! @brief score all pairs of neighboring sses within chains of a protein model
    //!
    //! @see @link example_score_protein_model_sse_neighbors.cpp @endlink
    //! @author woetzen, alexanns
    //! @date Nov 5, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelSSENeighbors :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! function to score neighboring pairs of sses
      util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> > m_ScorePair;

      //! normalize the final score by the total number of SSEs or sequences that have been considered
      bool m_Normalize;

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
      ProteinModelSSENeighbors();

      //! @brief construct from Pair function
      //! @param SP_PAIR_FUNCTION binary function to score a pair of sses
      //! @param NORMALIZE if true, final score will be normalized by the number of sses/sequences in the protein model
      ProteinModelSSENeighbors
      (
        const util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> > &SP_PAIR_FUNCTION,
        const bool NORMALIZE
      );

      //! virtual copy constructor
      ProteinModelSSENeighbors *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that scores the chain
      //! @param CHAIN the chain for which all neighbor scores are calculated
      //! @return score
      double operator()( const assemble::Chain &CHAIN) const;

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

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @param INDENT indentation
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ProteinModelSSENeighbors

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_PROTEIN_MODEL_SSE_NEIGHBORS_H_
