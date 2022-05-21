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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SWAP_BODY_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SWAP_BODY_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"
#include "find/bcl_find_pick_criteria_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "restraint/bcl_restraint_body.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSESwapBody
    //! @brief swaps an SSE from one assemble::SSEGeometryInterface to another
    //! @details MutateProteinModelSSESwapBody swaps an SSE from one assemble::SSEGeometryInterface to another, where the
    //! coord::bodies are defined by a restraint::Body. The SSE does not necessarily need to start out in a
    //! restraint::Body but after the mutation the SSE will be in a assemble::SSEGeometryInterface defined by a
    //! restraint::Body.
    //!
    //! @see @link example_fold_mutate_protein_model_sse_swap_body.cpp @endlink
    //! @author alexanns
    //! @date February 2, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelSSESwapBody :
      public math::MutateInterface< assemble::ProteinModel>
    {
    private:

    //////////
    // data //
    //////////

      //! ShPtr to a restraint:Body which is the restraint which provides the bodies an sse could be swapped into
      util::ShPtr< restraint::Body> m_BodyRestraint; 

      //! ShPtr to a PickCriteriaInterface which determines how a assemble::SSEGeometryInterface is picked out of "m_BodyRestraint"
      //! which is then used as the assemble::SSEGeometryInterface an SSE is swapped into
      util::ShPtr
      <
        find::PickCriteriaInterface
        <
          util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface
        >
      > m_BodyPicker;

      //! ShPtr to a LocatorInterface determines which SSE in the protein model will be chosen to be swapped to a new
      //! body
      util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > m_Locator;

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
      MutateProteinModelSSESwapBody();

      //! @brief constructor taking each of the member variable types
      //! @param BODY_RESTRAINT ShPtr to a restraint::Body which will be "m_BodyRestraint"
      //! @param BODY_PICKER ShPtr to a PickerInterface which will be "m_BodyPicker"
      //! @param LOCATOR ShPtr to a LocatorInterface which will by "m_Locator"
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSESwapBody
      (
        const util::ShPtr< restraint::Body> &BODY_RESTRAINT,
        const util::ShPtr
        <
          find::PickCriteriaInterface
          <
            util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface
          >
        > &BODY_PICKER,
        const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &LOCATOR,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSESwapBody>()
      );

      //! @brief clone
      MutateProteinModelSSESwapBody *Clone() const;

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

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateProteinModelSSESwapBody

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_SWAP_BODY_H_
