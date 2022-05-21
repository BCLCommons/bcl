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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "fold/bcl_fold_mutate_protein_model_sse_swap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSESwap::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelSSESwap)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSESwap::MutateProteinModelSSESwap() :
      m_Locator(),
      m_Bend( false),
      m_Scheme( GetStaticClassName< MutateProteinModelSSESwap>())
    {
    }

    //! @brief constructor from a locator and a scheme
    //! @param LOCATOR locator to be used
    //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSESwap::MutateProteinModelSSESwap
    (
      const find::LocatorCriteriaInterface
      <
        util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
      > &LOCATOR,
      const bool BEND,
      const std::string &SCHEME
    ) :
      m_Locator( *LOCATOR.Clone()),
      m_Bend( BEND),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a ShPtr to a locator and a scheme
    //! @param SP_LOCATOR ShPtr to locator to be used
    //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSESwap::MutateProteinModelSSESwap
    (
      const util::ShPtr
      <
        find::LocatorCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        >
      > &SP_LOCATOR,
      const bool BEND,
      const std::string &SCHEME
    ) :
      m_Locator( *SP_LOCATOR),
      m_Bend( BEND),
      m_Scheme( SCHEME)
    {
    }

    //! @brief clone
    MutateProteinModelSSESwap *MutateProteinModelSSESwap::Clone() const
    {
      return new MutateProteinModelSSESwap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSESwap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSESwap::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // if it has no SSEs or only one SSE
      if
      (
        PROTEIN_MODEL.GetNumberSSE( biol::GetSSTypes().HELIX) < 2 &&
        PROTEIN_MODEL.GetNumberSSE( biol::GetSSTypes().STRAND) < 2
      )
      {
        // then return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // randomly locate a sse from the protein model, this sse will serve as criterion of finding an sse to swap with
      const util::SiPtr< const assemble::SSE> located_sse_a( assemble::LocatorSSERandom().Locate( PROTEIN_MODEL));

      // assert located_sse_a is defined
      BCL_Assert( located_sse_a.IsDefined(), "could not find sse a in model");

      // create temporary util::SiPtr< const SSE> of all the sses in protein model
      util::SiPtrVector< const assemble::SSE> sse_vector( PROTEIN_MODEL.GetSSEs());
      // convert to list
      util::SiPtrList< const assemble::SSE> sse_list( sse_vector.Begin(), sse_vector.End());

      // find the located sse in the list
      util::SiPtrList< const assemble::SSE>::iterator
        sse_list_itr( std::find( sse_list.Begin(), sse_list.End(), located_sse_a));

      // make sure it was found
      BCL_Assert( sse_list_itr != sse_list.End(), "The located SSE was not found in the model when searched");

      // remove it from the list
      sse_list.Remove( sse_list_itr);

      // locate sse that is going to be swapped with located_sse_a
      const util::SiPtr< const assemble::SSE> located_sse_b
      (
        m_Locator->Locate( sse_list, *located_sse_a)
      );

      // if located_sse_b is not defined (e.g. because located_sse_a was the only sse of its type in the model)
      if( !located_sse_b.IsDefined())
      {
        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // report selected sse from pool to be swapped
      BCL_MessageVrb( "selected sse A to be swapped: " + located_sse_a->GetIdentification());
      BCL_MessageVrb( "selected sse B to be swapped: " + located_sse_b->GetIdentification());

      // make copies of both SSEs
      util::ShPtr< assemble::SSE> new_sse_a( located_sse_a->Clone()), new_sse_b( located_sse_b->Clone());

      // if the SSEs need to be bent
      if( m_Bend)
      {
        // fit the sses
        new_sse_a->FitToSSE( *located_sse_b);
        new_sse_b->FitToSSE( *located_sse_a);
      }
      // don't bend the SSEs, keep their current phi/psi angles
      else
      {
        // transform sse_a to position of b
        math::TransformationMatrix3D transform_ab( math::Inverse( located_sse_a->GetOrientation()));
        transform_ab( located_sse_b->GetOrientation());

        // transform sseb to position of a
        math::TransformationMatrix3D transform_ba( math::Inverse( located_sse_b->GetOrientation()));
        transform_ba( located_sse_a->GetOrientation());

        new_sse_a->Transform( transform_ab);
        new_sse_b->Transform( transform_ba);
      }

      // make copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace the sse with the swapped hardcopies of sse a and b
      new_model->Replace( new_sse_a);
      new_model->Replace( new_sse_b);

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
