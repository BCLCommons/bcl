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

#ifndef BCL_SCORE_PROTEIN_MODEL_SSE_PAIRS_H_
#define BCL_SCORE_PROTEIN_MODEL_SSE_PAIRS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelSSEPairs
    //! @brief Iterate over all possible sse pairs within a protein model and call a binary function on SSEs for scoring
    //!
    //! @see @link example_score_protein_model_sse_pairs.cpp @endlink
    //! @author woetzen
    //! @date Aug 16, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelSSEPairs :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! function to score two pairs of sequences or sses
      util::Implementation< math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double> > m_SpScorePair;

      //! normalize the final score by the total number of SSEs or sequences that have been considered
      bool m_Normalize;

      //! score type
      ProteinModel::Type m_ScoreType;

      //! readable scheme
      std::string m_ReadableScheme;

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
      ProteinModelSSEPairs();

      //! @brief construct from Pair function
      //! @param SP_PAIR_FUNCTION binary function to score a pair of sses
      //! @param NORMALIZE if true, final score will be normalized by the number of sses/sequences in the protein model
      //! @param SCORE_TYPE score type
      //! @param READABLE_SCHEME scheme that is more human readable
      ProteinModelSSEPairs
      (
        const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double> > &SP_PAIR_FUNCTION,
        const bool NORMALIZE,
        const ProteinModel::Type &SCORE_TYPE = ProteinModel::e_Undefined,
        const std::string &READABLE_SCHEME = ""
      );

      //! virtual copy constructor
      ProteinModelSSEPairs *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get a more readable score scheme
      //! @return a more readable score scheme
      const std::string &GetReadableScheme() const
      {
        return m_ReadableScheme;
      }

      //! @brief get score type
      //! @return score type
      ProteinModel::Type GetType() const
      {
        return m_ScoreType;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

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

    }; // ProteinModelSSEPairs

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_SSE_PAIRS_H_
