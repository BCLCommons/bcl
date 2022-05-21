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

#ifndef BCL_SCORE_SSE_POOL_SSES_H_
#define BCL_SCORE_SSE_POOL_SSES_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPoolSSEs
    //! @brief scores all sses in the pool with the given sse score
    //! @details for each sse the given sse score is evaluted and evtl. normalized by the number of scored entities.
    //!          These scores are summed up and evtl. normalized by the number of sses that were scored.
    //!
    //! @see @link example_score_sse_pool_sses.cpp @endlink
    //! @author woetzen
    //! @date Jun 16, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPoolSSEs :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSEPool, biol::Membrane, double>
    {

    private:

    //////////
    // data //
    //////////

      //! @brief score for a single sse
      util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > m_SSEScore;

      //! @brief normalize by the number of scored entities in an sse
      bool m_NormalizeSSE;

      //! @brief normalize by the number of sses
      bool m_NormalizeNumberSSEs;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from sse score and normalization
      SSEPoolSSEs
      (
        const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > &SP_SCORE_SSE,
        const bool NORMALIZE_SSE,
        const bool NORMALIZE_NUMBER_SSES
      );

      //! @brief Clone function
      //! @return pointer to new SSEPoolSSEs
      SSEPoolSSEs *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief access to the sse score
      //! @return Function to score a single sse
      const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > &GetSSEScore() const
      {
        return m_SSEScore;
      }

      //! @brief set the score
      //! @param SP_SSE_FUNCTION score for a single sse
      void SetSSEScore
      (
        const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > &SP_SSE_FUNCTION
      )
      {
        m_SSEScore = SP_SSE_FUNCTION;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that scores the pool
      //! @param POOL the sse pool to be scored
      //! @param MEMBRANE membrane object
      //! @return the sum of the normalized scores
      double operator()( const assemble::SSEPool &POOL, const biol::Membrane &MEMBRANE) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param POOL Argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @param MEMBRANE membrane object
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::SSEPool &POOL,
        const biol::Membrane &MEMBRANE,
        std::ostream &OSTREAM
      ) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class SSEPoolSSEs

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_SSE_POOL_SSES_H_
