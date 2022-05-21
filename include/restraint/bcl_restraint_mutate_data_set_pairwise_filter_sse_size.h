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

#ifndef BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_SSE_SIZE_H_
#define BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_SSE_SIZE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDataSetPairwiseFilterSSESize
    //! @brief Mutates a DataSetPairwise by removing data pairs with either point on an sse that is too small
    //! @details A pool of sses indicates the possible sses that are present so that it can be determined which sse
    //!          a given data point fall in. Minimum sse lengths are defined to determine if the data point falls within
    //!          an sse of desirable length or not. If there are overlapping sses in the pool, then a random
    //!          non-overlapping set is used.
    //!
    //! @see @link example_restraint_mutate_data_set_pairwise_filter_sse_size.cpp @endlink
    //! @author alexanns
    //! @date May 10, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDataSetPairwiseFilterSSESize :
      public math::MutateInterface< DataSetPairwise>
    {

    private:

    //////////
    // data //
    //////////

      //! the pool indicating the sses that exist
      util::ShPtr< assemble::SSEPool> m_SSEPool;

      //! the minimum lengths that each sse size should have for a data point to be kept if within that sse
      storage::Map< biol::SSType, size_t> m_MinSSELengths;

      //! the scheme of this mutate
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
      explicit MutateDataSetPairwiseFilterSSESize( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking optional scheme
      //! @param SSE_POOL pool to use as sse definitions
      //! @param MIN_SSE_LENGTHS minimum lengths each sse size should have for a data point to be kept if within an sse
      //! @param SCHEME the scheme for this scoring function
      MutateDataSetPairwiseFilterSSESize
      (
        const util::ShPtr< assemble::SSEPool> &SSE_POOL,
        const storage::Map< biol::SSType, size_t> &MIN_SSE_LENGTHS,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new MutateDataSetPairwiseFilterSSESize
      MutateDataSetPairwiseFilterSSESize *Clone() const;

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

      //! @brief virtual operator taking an DataSetPairwise and returning a mutated object
      //! @param DATA_SET DataSetPairwise of interest that will be mutated
      //! @return MutateResult that results from mutating to the argument DATA_SET
      math::MutateResult< DataSetPairwise> operator()( const DataSetPairwise &DATA_SET) const;

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

    private:

    }; // class MutateDataSetPairwiseFilterSSESize

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_FILTER_SSE_SIZE_H_ 
