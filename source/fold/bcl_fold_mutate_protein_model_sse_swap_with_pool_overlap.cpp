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
#include "fold/bcl_fold_mutate_protein_model_sse_swap_with_pool_overlap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "fold/bcl_fold_placement_sse_short_loop.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSESwapWithPoolOverlap::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance
      (
        new MutateProteinModelSSESwapWithPoolOverlap()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSESwapWithPoolOverlap::MutateProteinModelSSESwapWithPoolOverlap() :
      m_Collector(),
      m_SingleSwapsOnly(),
      m_Bend( false),
      m_Scheme( GetStaticClassName< MutateProteinModelSSESwapWithPoolOverlap>())
    {
    }

    //! @brief constructor from a Collector, single swap boolean and a scheme
    //! @param COLLECTOR collector to be used
    //! @param SINGLE_SWAPS_ONLY use single swaps only
    //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSESwapWithPoolOverlap::MutateProteinModelSSESwapWithPoolOverlap
    (
      const util::ShPtr
      <
        find::CollectorCriteriaInterface
        <
          util::SiPtrList< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        >
      > &COLLECTOR,
      const bool SINGLE_SWAPS_ONLY,
      const bool BEND,
      const std::string &SCHEME
    ) :
      m_Collector( *COLLECTOR),
      m_SingleSwapsOnly( SINGLE_SWAPS_ONLY),
      m_Bend( BEND),
      m_Scheme( SCHEME)
    {
    }

    //! @brief clone
    MutateProteinModelSSESwapWithPoolOverlap *MutateProteinModelSSESwapWithPoolOverlap::Clone() const
    {
      return new MutateProteinModelSSESwapWithPoolOverlap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSESwapWithPoolOverlap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateProteinModelSSESwapWithPoolOverlap::GetAlias() const
    {
      static const std::string s_alias( "MutateProteinModelSSESwapWithPoolOverlap");
      return s_alias;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateProteinModelSSESwapWithPoolOverlap::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Swaps two SSEs of the same type.");
      serializer.AddInitializer
      (
        "collector",
        "collect SSEs from a list",
        io::Serialization::GetAgent( &m_Collector)
      );
      serializer.AddInitializer
      (
        "bend",
        "bend SSEs after swapping",
        io::Serialization::GetAgent( &m_Bend)
      );
      serializer.AddInitializer
      (
        "single swaps",
        "disallow multiple swaps",
        io::Serialization::GetAgent( &m_SingleSwapsOnly)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSESwapWithPoolOverlap::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize static random picker
      static const assemble::PickSSERandom random_sse_picker;

      // initialize empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // if empty model
      if( PROTEIN_MODEL.GetSSEs().IsEmpty())
      {
        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );
      if( !sp_pool.IsDefined())
      {
        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // find sses in protein model that have at least one overlapping sse in the pool
      const util::SiPtrList< const assemble::SSE>
      eligible_sses
      (
        PROTEIN_MODEL.GetOverlappingSSEs( sp_pool->GetNonIdenticalSSEs( PROTEIN_MODEL))
      );

      // if there are no eligible sses
      if( eligible_sses.IsEmpty())
      {
        // warn user
        BCL_MessageVrb( "no sses in the protein model overlap with the pool");

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // chose one of these sses randomly
      const util::SiPtr< const assemble::SSE> located_sse( random_sse_picker.Pick( eligible_sses));

      // make sure the sse is found correctly
      BCL_Assert( located_sse.IsDefined(), "could not find sses in model using " + util::Format()( random_sse_picker));

      // report selected sse from pool to be swapped
      BCL_MessageVrb( "sse to be swapped with pool: " + located_sse->GetIdentification());

      // find overlapping sses from the protein model and collect the subset of the overlapping SSEs
      util::SiPtrList< const assemble::SSE> overlapping_sses
      (
        m_Collector->Collect( sp_pool->GetOverlappingSSEs( *located_sse), *located_sse)
      );

      // if none of the sses in overlapping_sses_vector fulfilled the criterion, return empty result
      if( overlapping_sses.IsEmpty())
      {
        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // make a copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // boolean to indicate if swapping with multiple overlapping was possible
      bool swap_multiple( false);

      // if there are more than one overlapping SSE
      if( !m_SingleSwapsOnly && overlapping_sses.GetSize() > 1)
      {
        // see if there are any pair of these overlapping SSEs
        // that do not overlap with each other and also do not overlap with anything else in the model

        // make a new list
        util::SiPtrList< const assemble::SSE> possible_sses;

        // iterate over the list
        for
        (
          util::SiPtrList< const assemble::SSE>::const_iterator
            pool_itr( overlapping_sses.Begin()), pool_itr_end( overlapping_sses.End());
          pool_itr != pool_itr_end; ++pool_itr
        )
        {
          // if this SSE overlaps with only one SSE, thus if it does not overlap with another SSE other than the located_sse
          if( PROTEIN_MODEL.GetOverlappingSSEs( **pool_itr).GetSize() == 1)
          {
            possible_sses.PushBack( *pool_itr);
          }
        }

        // if there are enough SSEs left
        if( possible_sses.GetSize() >= 2)
        {
          // make a list of possible pairs
          storage::Vector< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > possible_pairs;

          // iterate over the possible SSEs
          for
          (
            util::SiPtrList< const assemble::SSE>::const_iterator
              sse_itr_a( possible_sses.Begin()), sse_itr_end( possible_sses.End());
            sse_itr_a != sse_itr_end; ++sse_itr_a
          )
          {
             // initialize sse_itr_b
            util::SiPtrList< const assemble::SSE>::const_iterator sse_itr_b( sse_itr_a);
            ++sse_itr_b;

            // iterate while not reaching end
            while( sse_itr_b != sse_itr_end)
            {
              // if two sses do not overlap
              if( !biol::DoOverlap( **sse_itr_a, **sse_itr_b))
              {
                possible_pairs.PushBack
                (
                  storage::VectorND< 2, util::SiPtr< const assemble::SSE> >( **sse_itr_a, **sse_itr_b)
                );
                possible_pairs.PushBack
                (
                  storage::VectorND< 2, util::SiPtr< const assemble::SSE> >( **sse_itr_b, **sse_itr_a)
                );
              }
              ++sse_itr_b;

            } // sse_itr_b
          } // sse_itr_a

          // if such pairs were found
          if( !possible_pairs.IsEmpty())
          {
            // pick one randomly
            storage::VectorND< 2, util::SiPtr< const assemble::SSE> > sse_pair( possible_pairs.RemoveRandomElement());

            // make a copy of both SSEs
            util::ShPtr< assemble::SSE> sse_from_pool_a( sse_pair.First()->Clone());
            util::ShPtr< assemble::SSE> sse_from_pool_b( sse_pair.Second()->Clone());

            // if the SSE should be bent
            if( m_Bend)
            {
              // fit the sequence
              sse_from_pool_a->FitToSSE( *located_sse);
            }
            // just place an ideal SSE
            else
            {
              // put the first one in located SSEs place
              sse_from_pool_a->Transform( located_sse->GetOrientation());
            }

            // replace with located SSE
            new_model->ReplaceWithOverlapping( sse_from_pool_a);

            // now add the second one using short loop placement
            storage::Pair< math::TransformationMatrix3D, bool> transformation
            (
              PlacementSSEShortLoop().Place( *sse_from_pool_b, *new_model)
            );

            // if the placement was successful
            if( transformation.Second())
            {
              // apply the transformation
              sse_from_pool_b->Transform( transformation.First());

              // insert into model
              new_model->Insert( sse_from_pool_b);
            }

            // set the boolean
            swap_multiple = true;
          }
        }
      }

      // if swap multiple was not possible
      if( !swap_multiple)
      {

        // pick one randomly from overlapping_sses
        util::SiPtr< const assemble::SSE> overlapping_sse( random_sse_picker.Pick( overlapping_sses));

        // assert that one was picked correctly
        BCL_Assert( overlapping_sse.IsDefined(), "Could not pick random overlapping sse from pool");

        // make a copy of the sse from the pool
        util::ShPtr< assemble::SSE> sse_from_pool( overlapping_sse->Clone());

        // if the SSE should be bent
        if( m_Bend)
        {
          // fit the sequence
          sse_from_pool->FitToSSE( *located_sse);
        }
        // just place an ideal SSE
        else
        {
          // set the body of the sse_from_pool to located sse
          sse_from_pool->Transform( located_sse->GetOrientation());
        }

        // replace the located_sse with sse_from_pool
        new_model->ReplaceWithOverlapping( sse_from_pool);
      }

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
