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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SWAP_WITH_POOL_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SWAP_WITH_POOL_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_criteria_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSESwapWithPool
    //! @brief swaps the located SSE in ProteinModel with one from SSEPool
    //!
    //! @see @link example_fold_mutate_protein_model_sse_swap_with_pool.cpp @endlink
    //! @author karakam
    //! @date May 28, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MutateProteinModelSSESwapWithPool :
      public math::MutateInterface< assemble::ProteinModel>
    {
    private:

    //////////
    // data //
    //////////

      //! locates a sse in a protein model
      util::Implementation
      <
        find::LocatorCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        >
      > m_Locator;

      //! whether to bend the swapped SSE to match the phi/psi angles of the previous SSE at that location
      bool m_Bend;

      //! scheme
      std::string m_Scheme;

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
      MutateProteinModelSSESwapWithPool();

      //! @brief constructor from a Locator and a scheme
      //! @param LOCATOR Locator that locates an SSE from model to swap with pool
      //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSESwapWithPool
      (
        const find::LocatorCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        > &LOCATOR,
        const bool BEND = false,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSESwapWithPool>()
      );

      //! @brief constructor from a ShPtr to Locator and a scheme
      //! @param SP_LOCATOR ShPtr to locator that locates an SSE from model to swap with pool
      //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSESwapWithPool
      (
        const util::ShPtr
        <
          find::LocatorCriteriaInterface
          <
            util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
          >
        > &SP_LOCATOR,
        const bool BEND = false,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSESwapWithPool>()
      );

      //! @brief clone
      MutateProteinModelSSESwapWithPool *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
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

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateProteinModelSSESwapWithPool

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SWAP_WITH_POOL_H_
