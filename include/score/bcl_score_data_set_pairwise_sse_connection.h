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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_SSE_CONNECTION_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_SSE_CONNECTION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

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
    //! @class DataSetPairwiseSSEConnection
    //! @brief Implementation of element connection term
    //! @details see Kazmier, K. Journal of Structural Biology
    //! "Algorithm for selection of optimized EPR distance restraints for de novo protein structure determination"
    //! Section 2.1.3
    //!
    //! @see @link example_score_data_set_pairwise_sse_connection.cpp @endlink
    //! @author alexanns
    //! @date May 8, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseSSEConnection :
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
      explicit DataSetPairwiseSSEConnection( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking optional scheme
      //! @param SSE_POOL pool to use as sse definitions
      //! @param SCHEME the scheme for this scoring function
      explicit DataSetPairwiseSSEConnection
      (
        const util::ShPtr< assemble::SSEPool> &SSE_POOL, const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseSSEConnection
      DataSetPairwiseSSEConnection *Clone() const;

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

      //! @brief Gives information on how many times each SSE pair is connected by a data pair
      //! @param SSES the sses which will be used to count how many datapoints fall in each
      //! @param DATA the dataset which will be used to see how many of each data point are in each sse
      //! @return map of sse pairs and how many times they are connected
      static storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t>
      GetSSEPairCounts
      (
        const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> &SSES,
        const restraint::DataSetPairwise &DATA
      );

      //! @brief calculates the R component of the SSE connection score
      //! @param SSE_PAIRS_COUNTS map of sse pairs and how many times they are connected
      //! @param C_PRIME_PRIME maximum acceptable integer value for C, the ideal number of connections for each SSE pair
      //! @param NUMBER_RESTRAINTS this is "r" in the reference
      //! @return double which is the R component of the SSE connection score
      static double CalculateSSEConnectionScoreComponentR
      (
        const storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t> &SSE_PAIRS_COUNTS,
        const size_t C_PRIME_PRIME,
        const size_t NUMBER_RESTRAINTS
      );

      //! @brief sse connection score term C, fraction of SSE pairs that contain exactly the ideal number of restraints
      //! @param SSE_PAIRS_COUNTS map of sse pairs and how many times they are connected
      //! @param C_PRIME minimum acceptable integer value for C, the ideal number of connections for each SSE pair
      //! @param C_PRIME_PRIME maximum acceptable integer value for C, the ideal number of connections for each SSE pair
      //! @param NUM_SSE_PAIRS this is "p" in the reference
      //! @param REMAINDER this is "M" in the reference
      //! @return score for sse connection score component C
      static double CalculateSSEConnectionScoreComponentC
      (
        const storage::Map< storage::VectorND< 2, util::SiPtr< const assemble::SSE> >, size_t> &SSE_PAIRS_COUNTS,
        const size_t C_PRIME,
        const size_t C_PRIME_PRIME,
        const size_t NUM_SSE_PAIRS,
        const size_t REMAINDER
      );

    }; // class DataSetPairwiseSSEConnection

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DATA_SET_PAIRWISE_SSE_CONNECTION_H_ 
