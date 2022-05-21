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

#ifndef BCL_SCORE_PROTEIN_MODEL_SSE_H_
#define BCL_SCORE_PROTEIN_MODEL_SSE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelSSE
    //! @brief wrapper classes for scoring all SSEs in a given model
    //! @details This class collects all the SSEs from a given model and calculates the score for each individual SSE
    //! using the provided SSE scoring function stored in data member m_ScoreSingle. Depending on the initialization
    //! it can also normalized the sum of SSE scores by the count of SSEs
    //!
    //! @see @link example_score_protein_model_sse.cpp @endlink
    //! @author woetzen
    //! @date Aug 16, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelSSE :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! ShPtr to score to be used
      //! return first of the function is the actual score for the SSE, the second is the number of entities being scored
      util::Implementation< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > m_ScoreSingle;

      //! normalize by the number of the scored entities
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

      //! @brief default constructor
      ProteinModelSSE();

      //! @brief construct from a function
      //! @param SP_SINGLE_FUNCTION single function to be used
      //! @param NORMALIZE
      //! @param SCORE_TYPE score type
      //! @param READABLE_SCHEME scheme that is more human readable
      ProteinModelSSE
      (
        const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > &SP_SINGLE_FUNCTION,
        const bool NORMALIZE,
        const ProteinModel::Type &SCORE_TYPE = ProteinModel::e_Undefined,
        const std::string &READABLE_SCHEME = ""
      );

      //! @brief virtual copy constructor
      ProteinModelSSE *Clone() const;

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
      Type GetType() const
      {
        return m_ScoreType;
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an PROTEIN_MODEL and returning a t_ResultType object
      //! @param PROTEIN_MODEL Protein Model to be used to evaluate the function
      //! @return function value of the given argument
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

    }; // class ProteinModelSSE

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_PROTEIN_MODEL_SSE_H_
