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
#include "fold/bcl_fold_mutate_protein_model_sse_swap_with_pool.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSESwapWithPool::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelSSESwapWithPool())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSESwapWithPool::MutateProteinModelSSESwapWithPool() :
      m_Locator(),
      m_Bend( false),
      m_Scheme( GetStaticClassName< MutateProteinModelSSESwapWithPool>())
    {
    }

    //! @brief constructor from a Locator and a scheme
    //! @param LOCATOR Locator that locates an SSE from model to swap with pool
    //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSESwapWithPool::MutateProteinModelSSESwapWithPool
    (
      const find::LocatorCriteriaInterface
      <
        util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
      > &LOCATOR,
      const bool BEND,
      const std::string &SCHEME
    ) :
      m_Locator( LOCATOR),
      m_Bend( BEND),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a ShPtr to Locator and a scheme
    //! @param SP_LOCATOR ShPtr to locator that locates an SSE from model to swap with pool
    //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSESwapWithPool::MutateProteinModelSSESwapWithPool
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
    MutateProteinModelSSESwapWithPool *MutateProteinModelSSESwapWithPool::Clone() const
    {
      return new MutateProteinModelSSESwapWithPool( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSESwapWithPool::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateProteinModelSSESwapWithPool::GetAlias() const
    {
      static const std::string s_alias( "MutateProteinModelSSESwapWithPool");
      return s_alias;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateProteinModelSSESwapWithPool::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Swaps two SSEs of the same type.");
      serializer.AddInitializer
      (
        "locator",
        "locates an SSE in a protein model",
        io::Serialization::GetAgent( &m_Locator)
      );
      serializer.AddInitializer
      (
        "bend",
        "bend SSEs after swapping",
        io::Serialization::GetAgent( &m_Bend)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSESwapWithPool::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // declare static random sse picker
      static const assemble::PickSSERandom random_sse_picker;

      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // if empty model
      if( PROTEIN_MODEL.GetSSEs().IsEmpty())
      {
        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // choose a random sse and make copy
      const util::SiPtr< const assemble::SSE> sse_from_model( assemble::LocatorSSERandom().Locate( PROTEIN_MODEL));

      // assert sse_from_model is defined
      BCL_Assert( sse_from_model.IsDefined(), "could not find sse in model");

      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );
      // if empty pool
      if( !sp_pool.IsDefined())
      {
        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // report selected sse from model to be swapped
      BCL_MessageVrb( "sse to be swapped with pool: " + sse_from_model->GetIdentification());

      // get non-overlapping SSEs from the pool
      util::SiPtrList< const assemble::SSE> non_overlapping_sses( sp_pool->GetNonOverlappingSSEs( PROTEIN_MODEL));

      // locate sse that is going to be swapped with located_sse_a
      const util::SiPtr< const assemble::SSE> sse_from_pool( m_Locator->Locate( non_overlapping_sses, *sse_from_model));

      // if the SSE was not defined
      if( !sse_from_pool.IsDefined())
      {
        BCL_MessageVrb( "Could not pick random non-overlapping sse from pool");
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // make a copy of the sse from the pool
      util::ShPtr< assemble::SSE> sse_to_be_inserted( sse_from_pool->Clone());

      // if the SSE should be bent
      if( m_Bend)
      {
        // fit the SSE
        sse_to_be_inserted->FitToSSE( *sse_from_model);
      }
      else
      // just place an ideal SSE at the location
      {
        // set the body of the sse_from_pool to located sse
        sse_to_be_inserted->Transform( sse_from_model->GetOrientation());
      }

      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // remove located sse
      new_model->Remove( *sse_from_model);

      // insert the sse from pool
      new_model->Insert( sse_to_be_inserted);

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
