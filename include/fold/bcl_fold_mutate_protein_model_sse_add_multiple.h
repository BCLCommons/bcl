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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_ADD_MULTIPLE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_ADD_MULTIPLE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_domain_interface.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSEAddMultiple
    //! @brief Mutate class to be used for adding several SSEs simultaneously
    //! @details Mutate class to be used for adding several SSEs simultaneously into
    //! a protein model using a SSEPool and PlacementDomainInterface
    //!
    //! @see @link example_fold_mutate_protein_model_sse_add_multiple.cpp @endlink
    //! @author weinerbe
    //! @date Feb 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelSSEAddMultiple :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! SSEPool
      util::ShPtr< assemble::SSEPool> m_SSEPool;

      //! Placement to be used to place the SSE in the model
      util::ShPtr< PlacementDomainInterface> m_Placement;

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
      MutateProteinModelSSEAddMultiple();

      //! @brief constructor from a SSEPool, a PlacementDomainInterface and a scheme
      //! @param SSE_POOL ShPtr to pool of SSEs to be used
      //! @param PLACEMENT placement for the SSEs in the pool into the model
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEAddMultiple
      (
        const util::ShPtr< assemble::SSEPool> &SSE_POOL,
        const PlacementDomainInterface &PLACEMENT,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEAddMultiple>()
      );

      //! @brief constructor from a ShPtr to a SSEPool, a ShPtr to PlacementDomainInterface and a scheme
      //! @param SSE_POOL ShPtr to pool of SSEs to be used
      //! @param SP_PLACEMENT ShPtr to placement for the SSEs in the pool into the model
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEAddMultiple
      (
        const util::ShPtr< assemble::SSEPool> &SSE_POOL,
        const util::ShPtr< PlacementDomainInterface> &SP_PLACEMENT,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEAddMultiple>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelSSEAddMultiple
      MutateProteinModelSSEAddMultiple *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a protein model and returning a mutate object of protein model type
      //! @param PROTEIN_MODEL protein model of interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

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

    }; // class MutateProteinModelSSEAddMultiple

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_ADD_MULTIPLE_H_ 
