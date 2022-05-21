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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_SSE_TERM_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_SSE_TERM_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwiseSSETerm
    //! @brief Implementation of secondary structure term
    //! @details see Kazmier, K. Journal of Structural Biology
    //! "Algorithm for selection of optimized EPR distance restraints for de novo protein structure determination"
    //! Section 2.1.2
    //!
    //! @see @link example_score_data_set_pairwise_sse_term.cpp @endlink
    //! @author alexanns
    //! @date May 7, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseSSETerm :
      public math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double>
    {

    private:

    //////////
    // data //
    //////////

      //! pool defining the sses that could be present
      util::ShPtr< assemble::SSEPool> m_SSEPool;

      //! the scheme for this scoring function
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking optional scheme
      //! @param SCHEME the scheme for this scoring function
      explicit DataSetPairwiseSSETerm( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking optional scheme
      //! @param SSE_POOL pool to use as sse definitions
      //! @param SCHEME the scheme for this scoring function
      explicit DataSetPairwiseSSETerm
      (
        const util::ShPtr< assemble::SSEPool> &SSE_POOL, const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseSSETerm
      DataSetPairwiseSSETerm *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the score of a data set
      //! @param DATA data set to be scored
      //! @return the score of the current data set
      double operator()( const restraint::DataSetPairwise &DATA) const;

    //////////////////////
    // input and output //
    //////////////////////

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

    public:

      //! @brief for each sse gives the number of data points that fall within that sse
      //! @param SSES the sses which will be used to count how many data points fall within them
      //! @param DATA the data set which will be used to count how many data points fall within the SSES
      //! @return map which has for each sse the number of data points that fall within that sse
      static storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap>
      GetSSECounts
      (
        const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> &SSES,
        const restraint::DataSetPairwise &DATA
      );

      //! @brief calculates the average percentage of points positioned in each SSE up to the ideal value, Qprimeprime
      //! @param SSE_COUNT map with a count for each sse of the number of points in each sse
      //! @param Q_PRIME_PRIME see reference
      //! @param NUMBER_OF_POINTS the number of points in the data set i.e. 2 * (dataset size) ( this is l in reference)
      //! @return double which is the L component of the sse term score
      static double CalculateSSEScoreComponentL
      (
        const storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap> &SSE_COUNT,
        const size_t Q_PRIME_PRIME,
        const size_t NUMBER_OF_POINTS
      );

      //! @brief calculates S component of SSE score, derived from frachtion of SSEs containing exactly
      //! @param SSE_COUNT map with a count for each sse of the number of points in each sse
      //! @param Q_PRIME_PRIME see reference
      //! @param REMAINDER "R" in reference
      //! @param NUMBER_SSES "s" in reference
      //! @return double which is component S of the sse score
      static double CalculateSSEScoreComponentS
      (
        const storage::Map< util::SiPtr< const assemble::SSE>, size_t, assemble::SSELessThanNoOverlap> &SSE_COUNT,
        const size_t Q_PRIME,
        const size_t Q_PRIME_PRIME,
        const size_t REMAINDER,
        const size_t NUMBER_SSES
      );

    }; // class DataSetPairwiseSSETerm

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DATA_SET_PAIRWISE_SSE_TERM_H_ 
