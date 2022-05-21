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
#include "fold/bcl_fold_default_mutates.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_collector_sse_paired.h"
#include "assemble/bcl_assemble_collector_sse_size.h"
#include "assemble/bcl_assemble_locator_domain_random.h"
#include "assemble/bcl_assemble_locator_sse_furthest.h"
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_locator_sse_unpaired.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_pick_sse_short_loops.h"
#include "assemble/bcl_assemble_sse_compare_type.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "coord/bcl_coord_move_combine.h"
#include "coord/bcl_coord_move_rotate_defined.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_transform_random.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "find/bcl_find_collector_criteria_combined.h"
#include "find/bcl_find_collector_criteria_wrapper.h"
#include "find/bcl_find_locator.h"
#include "find/bcl_find_locator_criteria.h"
#include "find/bcl_find_locator_criteria_wrapper.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "fold/bcl_fold_mutate_domain_flip.h"
#include "fold/bcl_fold_mutate_domain_shuffle.h"
#include "fold/bcl_fold_mutate_protein_model_domain.h"
#include "fold/bcl_fold_mutate_protein_model_pair_strands.h"
#include "fold/bcl_fold_mutate_protein_model_sse.h"
#include "fold/bcl_fold_mutate_protein_model_sse_add.h"
#include "fold/bcl_fold_mutate_protein_model_sse_move.h"
#include "fold/bcl_fold_mutate_protein_model_sse_pair.h"
#include "fold/bcl_fold_mutate_protein_model_sse_pair_hinge.h"
#include "fold/bcl_fold_mutate_protein_model_sse_remove.h"
#include "fold/bcl_fold_mutate_protein_model_sse_resize.h"
#include "fold/bcl_fold_mutate_protein_model_sse_split.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap_with_pool.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap_with_pool_overlap.h"
#include "fold/bcl_fold_mutate_protein_model_strand_switch_sheet.h"
#include "fold/bcl_fold_mutate_sheet_cycle.h"
#include "fold/bcl_fold_mutate_sheet_divide.h"
#include "fold/bcl_fold_mutate_sheet_register_fix.h"
#include "fold/bcl_fold_mutate_sheet_register_shift.h"
#include "fold/bcl_fold_mutate_sheet_twist.h"
#include "fold/bcl_fold_mutate_sse_bend_ramachandran.h"
#include "fold/bcl_fold_mutate_sse_bend_random.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_placement_sse_next_to_sse.h"
#include "fold/bcl_fold_placement_sse_short_loop.h"
#include "fold/bcl_fold_setup.h"
#include "math/bcl_math_mutate_move_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DefaultMutates::s_Instance
    (
      GetObjectInstances().AddInstance( new DefaultMutates())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DefaultMutates::DefaultMutates()
    {
    }

    //! @brief Clone function
    //! @return pointer to new DefaultMutates
    DefaultMutates *DefaultMutates::Clone() const
    {
      return new DefaultMutates( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to single instance
    //! @return single instance of default mutates
    DefaultMutates &DefaultMutates::GetInstance()
    {
      static DefaultMutates s_single_instance;
      return s_single_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DefaultMutates::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the mutates and add them to Mutates enumerator
    void DefaultMutates::InitializeMutates()
    {
      // mutates already initialized
      if( e_AddSSENextToSSE.IsDefined())
      {
        return;
      }

      // create storage set of xyz axis
      storage::Set< coord::Axis> xyz_axes_set( coord::GetAxes().e_X, coord::GetAxes().e_Y, coord::GetAxes().e_Z);

    ///////////////////////////
    // collector and locator //
    ///////////////////////////

      // random sse locator
      const util::ShPtr< assemble::LocatorSSERandom> sp_locator_sse_random( new assemble::LocatorSSERandom());

      // helix collector
      const assemble::CollectorSSE collector_helix( biol::GetSSTypes().HELIX);

      // strand collector
      const assemble::CollectorSSE collector_strand( biol::GetSSTypes().STRAND);

      // random picker
      const assemble::PickSSERandom picker_sse_random;

      // random helix locator
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > sp_locator_helix_random
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          collector_helix, picker_sse_random
        )
      );

      // random strand locator
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > sp_locator_strand_random
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          collector_strand, picker_sse_random
        )
      );

      // pool picker
      const find::PickCriteriaWrapper
      <
        util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
      > picker_sse_pool( picker_sse_random);

      // random pool strand picker
      const find::PickCriteriaWrapper
      <
        util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
      >
      picker_strand_sse_pool_random( assemble::PickSSERandom( storage::Set< biol::SSType>( biol::GetSSTypes().STRAND)));

      //  short loop picker
      const size_t short_loop_length( 7);
      const assemble::PickSSEShortLoops picker_sse_pool_short_loop( short_loop_length);

      const find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE>
        locator_sse_random_crit_wrap( sp_locator_sse_random);

      // construct Set to stra
      storage::Set< contact::Type> pair_contact_types;
      pair_contact_types.Insert( contact::GetTypes().HELIX_HELIX);
      pair_contact_types.Insert( contact::GetTypes().HELIX_SHEET);
      pair_contact_types.Insert( contact::GetTypes().SHEET_HELIX);
      pair_contact_types.Insert( contact::GetTypes().HELIX_STRAND);
      pair_contact_types.Insert( contact::GetTypes().STRAND_HELIX);
      pair_contact_types.Insert( contact::GetTypes().UNDEFINED_HELIX_STRAND);
      pair_contact_types.Insert( contact::GetTypes().UNDEFINED_STRAND_HELIX);

      // create helix pair collector
      const util::ShPtr< assemble::CollectorSSEPaired> collector_sse_pair
      (
        new assemble::CollectorSSEPaired
        (
          assemble::GetSSEGeometryPackingPickers().e_BestInteractionWeight,
          pair_contact_types,
          12.0, // distance
          true  // orthogonal
        )
      );

      // create helix pair collector
      const util::ShPtr< assemble::CollectorSSEPaired> collector_helix_pair
      (
        new assemble::CollectorSSEPaired
        (
          assemble::GetSSEGeometryPackingPickers().e_BestInteractionWeight,
          storage::Set< contact::Type>( contact::GetTypes().HELIX_HELIX),
          12.0, // distance
          true  // orthogonal
        )
      );

      // create sheet collector
      const util::ShPtr< find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::ProteinModel> >
        collector_sheet( new assemble::CollectorSheet());

      // initialize first collector for use in locator 1
      // this collector collects all the sses that agree in type with the criterion sse
      const find::CollectorCriteriaCombined< assemble::SSE> collector_same_sstype
      (
        util::ShPtr< util::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, bool> >
        (
          new assemble::SSECompareType()
        )
      );

      // initialize random picker (this picker will be used with both collectors)
      const find::PickCriteriaWrapper
      <
        util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
      > picker( picker_sse_random);

      // create locator to locate sses of same type
      const find::LocatorCriteria
      <
        util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE,
        util::SiPtrList< const assemble::SSE>
      > locator_ss_type( util::CloneToShPtr( collector_same_sstype), util::CloneToShPtr( picker));

      // collector sse wrapped
      const assemble::CollectorSSE collector_sse;
      const find::CollectorCriteriaWrapper< util::SiPtrList< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>
      collector_sse_wrapped( collector_sse);

      // collector helix wrapped
      const find::CollectorCriteriaWrapper< util::SiPtrList< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>
      collector_helix_wrapped( assemble::CollectorSSE( biol::GetSSTypes().HELIX));

      // collector strand wrapped
      const find::CollectorCriteriaWrapper< util::SiPtrList< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>
      collector_strand_wrapped( assemble::CollectorSSE( biol::GetSSTypes().STRAND));

      // constructor locator for furthest SSE from center
      const assemble::LocatorSSEFurthest< assemble::DomainInterface> locator_sse_furthest( collector_sse_wrapped);

      // constructor locator for furthest helix from center
      const assemble::LocatorSSEFurthest< assemble::DomainInterface> locator_helix_furthest( collector_helix_wrapped);

      // constructor locator for furthest strand from center
      const assemble::LocatorSSEFurthest< assemble::DomainInterface> locator_strand_furthest( collector_strand_wrapped);

      // random helix locator criteria wrapped for domain
      const find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>
      locator_sse_random_wrapped
      (
        sp_locator_sse_random
      );

      // random helix locator criteria wrapped for domain
      const find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>
      locator_helix_random_wrapped
      (
        sp_locator_helix_random
      );

      // random strand locator criteria wrapped for domain
      const find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>
      locator_strand_random_wrapped
      (
        sp_locator_strand_random
      );

    ///////////////
    // placement //
    ///////////////

      // place an sse next to a located sse in the protein model with a distance deviation factor of 1.0
      const PlacementSSENextToSSE placement_next_to_sse( locator_sse_random_crit_wrap);

      // place an sse with short loop next to a located sse in the protein model with a distance deviation factor of 1.0
      const PlacementSSEShortLoop placement_short_loop( short_loop_length);

      // place a strand from loop next to an edge strand in a sheet
      const PlacementStrandNextToSheet placement_strand_next_to_sheet( 0.5);

    //////////////
    // sse adds //
    //////////////

      // adds a SSE from the pool next to an SSE in the model
      e_AddSSENextToSSE = GetMutates().AddMutate
      (
        MutateProteinModelSSEAdd
        (
          picker_sse_pool,
          placement_next_to_sse,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_sse_next_to_sse"
        )
      );

      // adds a SSE from pool next to an SSE in the model connected by a short loop
      e_AddSSEShortLoop = GetMutates().AddMutate
      (
        MutateProteinModelSSEAdd
        (
          picker_sse_pool_short_loop,
          placement_short_loop,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_sse_short_loop"
        )
      );

      // adds a strand from pool next to another strand/sheet in the model
      e_AddStrandNextToSheet = GetMutates().AddMutate
      (
        MutateProteinModelSSEAdd
        (
          picker_strand_sse_pool_random,
          placement_strand_next_to_sheet,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_strand_next_to_sheet"
        )
      );

    //////////////////
    // sse removals //
    //////////////////

      // remove a random SSE from model
      e_RemoveRandom = GetMutates().AddMutate
      (
        MutateProteinModelSSERemove
        (
          sp_locator_sse_random,
          MutateTree::GetMutateTypeName( MutateTree::e_Remove) + "_random"
        )
      );

      // remove an unpaired strand from model
      e_RemoveUnpairedStrand = GetMutates().AddMutate
      (
        MutateProteinModelSSERemove
        (
          util::CloneToShPtr( assemble::LocatorSSEUnpaired( assemble::PickSSEFurthestEuclidean(), contact::GetTypes().STRAND_STRAND, 20.0)),
          MutateTree::GetMutateTypeName( MutateTree::e_Remove) + "_unpaired_strand"
        )
      );

    ///////////////
    // sse swaps //
    ///////////////

      // swap a SSE of same type within the model
      e_SwapSSEs = GetMutates().AddMutate
      (
        MutateProteinModelSSESwap
        (
          locator_ss_type,
          DefaultFlags::GetFlagFitSwappedSSEs()->GetFlag(),
          MutateTree::GetMutateTypeName( MutateTree::e_Swap) + "_sses"
        )
      );

      // swap a SSE in model, with a non-overlapping SSE from the pool of same type
      e_SwapSSEWithPool = GetMutates().AddMutate
      (
        MutateProteinModelSSESwapWithPool
        (
          locator_ss_type,
          DefaultFlags::GetFlagFitSwappedSSEs()->GetFlag(),
          MutateTree::GetMutateTypeName( MutateTree::e_Swap) + "_sse_with_pool"
        )
      );

      // swap a SSE in model, with an overlapping SSE (or two if possible) from the pool of same type
      e_SwapSSEWithPoolOverlap = GetMutates().AddMutate
      (
        MutateProteinModelSSESwapWithPoolOverlap
        (
          util::CloneToShPtr( collector_same_sstype),
          false,
          DefaultFlags::GetFlagFitSwappedSSEs()->GetFlag(),
          MutateTree::GetMutateTypeName( MutateTree::e_Swap) + "_sse_with_pool_overlap"
        )
      );

    ////////////
    // resize //
    ////////////

      // resize move for SSEs cterm
      e_SSEResizeCTerm = GetMutates().AddMutate
      (
        MutateProteinModelSSEResize
        (
          sp_locator_sse_random,
          0.5, // extend prob
          math::Range< size_t>( 1, 3), // length range to extend/shrink
          biol::AASequenceFlexibility::e_CTerminal, // change both ends
          true, // recenter after resize
          assemble::SSEPool::GetCommandLineMinSSELengths(),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_resize_cterm"
        )
      );

      // resize move for SSEs nterm
      e_SSEResizeNTerm = GetMutates().AddMutate
      (
        MutateProteinModelSSEResize
        (
          sp_locator_sse_random,
          0.5, // extend prob
          math::Range< size_t>( 1, 3), // length range to extend/shrink
          biol::AASequenceFlexibility::e_NTerminal, // change both ends
          true, // recenter after resize
          assemble::SSEPool::GetCommandLineMinSSELengths(),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_resize_nterm"
        )
      );

    //////////////////////
    // single SSE moves //
    //////////////////////

      // move a SSE next to SSE
      e_SSEMoveNext = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_sse_random_wrapped,
          placement_next_to_sse,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_move_next"
        )
      );

      // move a SSE next to short loop SSE
      e_SSEMoveShortLoop = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_sse_random_wrapped,
          placement_short_loop,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_move_short_loop"
        )
      );

      // move the furthest SSE next to SSE
      e_SSEFurthestMoveNext = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_sse_furthest,
          placement_next_to_sse,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_furthest_move_next"
        )
      );

      // split coil length
      const math::Range< size_t> sse_split_coil_length_range( 1, 3);
      storage::Map< biol::SSType, size_t> sse_split_min_sizes;
      sse_split_min_sizes[ biol::GetSSTypes().HELIX] = 5;
      sse_split_min_sizes[ biol::GetSSTypes().STRAND] = 3;

      // ranges for the SSE locator by size
      storage::Map< biol::SSType, math::Range< size_t> > sse_size_ranges;
      sse_size_ranges[ biol::GetSSTypes().HELIX] = math::Range< size_t>( 15, util::GetUndefinedSize_t());
      sse_size_ranges[ biol::GetSSTypes().STRAND] = math::Range< size_t>( 9, util::GetUndefinedSize_t());

      // locator for large SSEs
      util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > sp_locator
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          util::ShPtr< assemble::CollectorSSESize>
          (
            new assemble::CollectorSSESize( sse_size_ranges)
          ),
          util::ShPtr< assemble::PickSSERandom>( new assemble::PickSSERandom())
        )
      );

      // iterate over ss prediction methods
      storage::Set< sspred::Method> ss_pred_methods( sspred::Methods::GetCommandLineMethods());
      for
      (
        storage::Set< sspred::Method>::const_iterator
          sspred_itr( ss_pred_methods.Begin()), sspred_itr_end( ss_pred_methods.End());
        sspred_itr != sspred_itr_end;
        ++sspred_itr
      )
      {
        // construct split moves for each method
        m_SSESplit[ *sspred_itr] = GetMutates().AddMutate
        (
          MutateProteinModelSSESplit
          (
            *sspred_itr,
            sp_locator,
            sse_split_min_sizes,
            sse_split_coil_length_range,
            MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_split_" + sspred_itr->GetName()
          )
        );
      }

      // bend a random SSE using a phi/psi angle pair picked from ramachandran distribution
      e_SSEBendRamachandran = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( MutateSSEBendRamachandran( math::Range< size_t>( 1, 2))),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_bend_ramachandran"
        )
      );

      // bend a random SSE with phi/psi changes within  +/- 5 degrees
      const math::Range< double> bend_small_range( -math::g_Pi / 36.0, math::g_Pi / 36.0);
      e_SSEBendRandomSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( MutateSSEBendRandom( bend_small_range, bend_small_range)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_bend_random_small"
        )
      );

      // bend a random SSE with phi/psi changes within +/- 20 degrees
      const math::Range< double> bend_large_range( -math::g_Pi/ 9.0, math::g_Pi/ 9.0);   // +- 20 degrees
      e_SSEBendRandomLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( MutateSSEBendRandom( bend_large_range, bend_large_range)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_bend_random_large"
        )
      );

      // init translation small and large along every axis
      const linal::Vector3D x_trans_s( 2.0, 0.0, 0.0);
      const linal::Vector3D y_trans_s( 0.0, 2.0, 0.0);
      const linal::Vector3D z_trans_s( 0.0, 0.0, 2.0);
      const linal::Vector3D x_trans_l( 6.0, 0.0, 0.0);
      const linal::Vector3D y_trans_l( 0.0, 6.0, 0.0);
      const linal::Vector3D z_trans_l( 0.0, 0.0, 6.0);

      // construct translation moves
      e_SSETranslateSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( 2.0, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_translate_small"
        )
      );
      e_SSETranslateXSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( x_trans_s, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_translate_x_small"
        )
      );
      e_SSETranslateYSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( y_trans_s, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_translate_y_small"
        )
      );
      e_SSETranslateZSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( z_trans_s, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_translate_z_small"
        )
      );
      e_SSETranslateLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( 2.0, 6.0, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_translate_large"
        )
      );
      e_SSETranslateXLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( x_trans_s, x_trans_l, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_translate_x_large"
        )
      );
      e_SSETranslateYLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( y_trans_s, y_trans_l, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_translate_y_large"
        )
      );
      e_SSETranslateZLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( z_trans_s, z_trans_l, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_translate_z_large"
        )
      );

      // initialize rotations
      const linal::Vector3D x_rot_s( math::g_Pi / 12.0, 0.0, 0.0);
      const linal::Vector3D y_rot_s( 0.0, math::g_Pi / 12.0, 0.0);
      const linal::Vector3D z_rot_s( 0.0, 0.0, math::g_Pi / 12.0);
      const linal::Vector3D x_rot_l( math::g_Pi / 4.0, 0.0, 0.0);
      const linal::Vector3D y_rot_l( 0.0, math::g_Pi / 4.0, 0.0);
      const linal::Vector3D z_rot_l( 0.0, 0.0, math::g_Pi / 4.0);

      // construct rotation moves
      e_SSERotateSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( math::g_Pi / 12.0, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_rotate_small"
        )
      );
      e_SSERotateXSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( x_rot_s, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_rotate_x_small"
        )
      );
      e_SSERotateYSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( y_rot_s, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_rotate_y_small"
        )
      );
      e_SSERotateZSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( z_rot_s, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_rotate_z_small"
        )
      );
      e_SSERotateLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( math::g_Pi / 12.0, math::g_Pi / 4.0, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_rotate_large"
        )
      );
      e_SSERotateXLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( x_rot_s, x_rot_l, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_rotate_x_large"
        )
      );
      e_SSERotateYLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( y_rot_s, y_rot_l, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_rotate_y_large"
        )
      );
      e_SSERotateZLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( z_rot_s, z_rot_l, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_rotate_z_large"
        )
      );

      // construct transformations
      e_SSETransformSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTransformRandom( 2.0, math::g_Pi / 12.0, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_transform_small"
        )
      );
      e_SSETransformLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTransformRandom( 2.0, 6.0, math::g_Pi / 12.0, math::g_Pi / 4.0, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_transform_large"
        )
      );

      // construct flips
      const coord::MoveRotateDefined x_flip( coord::MoveRotateDefined::GetFlipMove( coord::GetAxes().e_X));
      const coord::MoveRotateDefined y_flip( coord::MoveRotateDefined::GetFlipMove( coord::GetAxes().e_Y));
      const coord::MoveRotateDefined z_flip( coord::MoveRotateDefined::GetFlipMove( coord::GetAxes().e_Z));

      // add x, y and z flips
      e_SSEFlipX = GetMutates().AddMutate( MutateProteinModelSSE( sp_locator_sse_random, util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( x_flip), false)), MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_flip_x"));
      e_SSEFlipY = GetMutates().AddMutate( MutateProteinModelSSE( sp_locator_sse_random, util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( y_flip), false)), MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_flip_y"));
      e_SSEFlipZ = GetMutates().AddMutate( MutateProteinModelSSE( sp_locator_sse_random, util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( z_flip), false)), MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_flip_z"));

    ////////////////////////
    // single helix moves //
    ////////////////////////

      // move a helix next to SSE
      e_HelixMoveNext = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_helix_random_wrapped,
          placement_next_to_sse,
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_move_next"
        )
      );

      // move a helix next to short loop SSE
      e_HelixMoveShortLoop = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_helix_random_wrapped,
          placement_short_loop,
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_move_short_loop"
        )
      );

      // move the furthest helix next to SSE
      e_HelixFurthestMoveNext = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_helix_furthest,
          placement_next_to_sse,
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_furthest_move_next"
        )
      );

      // translation along xy axes 0 to 2 angstroms
      const linal::Vector3D xy_translation_small( 2.0, 2.0, 0.0);
      e_HelixTranslateXYSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( xy_translation_small, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_translate_xy_small"
        )
      );

      // translation along xy axes 2 to 4 angstroms
      const linal::Vector3D xy_translation_large( 4.0, 4.0, 0.0);
      e_HelixTranslateXYLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( xy_translation_small, xy_translation_large, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_translate_xy_large"
        )
      );

      // translation along z axis 0 to 2 angstroms
      const linal::Vector3D z_translation_small( 0.0, 0.0, 2.0);
      e_HelixTranslateZSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( z_translation_small, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_translate_z_small"
        )
      );

      // translation along z axis 2 to 4 angstroms
      const linal::Vector3D z_translation_large( 0.0, 0.0, 4.0);
      e_HelixTranslateZLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( z_translation_small, z_translation_large, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_translate_z_large"
        )
      );

      // rotation around xy axes 0 to 15 degrees
      const linal::Vector3D xy_rotation_small( math::g_Pi / 12.0, math::g_Pi / 12.0, 0.0);
      e_HelixRotateXYSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( xy_rotation_small, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_rotate_xy_small"
        )
      );

      // rotation around xy axes 15 to 45 degrees
      const linal::Vector3D xy_rotation_large( math::g_Pi / 4.0, math::g_Pi / 4.0, 0.0);
      e_HelixRotateXYLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( xy_rotation_small, xy_rotation_large, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_rotate_xy_large"
        )
      );

      // rotation around z axis 0 to 15 degrees
      const linal::Vector3D z_rotation_small( 0.0, 0.0, math::g_Pi / 12.0);
      e_HelixRotateZSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( z_rotation_small, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_rotate_z_small"
        )
      );

      // rotation around z axis 15 to 45 degrees
      const linal::Vector3D z_rotation_large( 0.0, 0.0, math::g_Pi / 4.0);
      e_HelixRotateZLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( z_rotation_small, z_rotation_large, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_rotate_z_large"
        )
      );

      // small transformation around xy axes
      e_HelixTransformXYSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTransformRandom( xy_translation_small, xy_rotation_small, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_transform_xy_small"
        )
      );

      // large transformation around xy axes
      e_HelixTransformXYLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTransformRandom( xy_translation_small, xy_translation_large, xy_rotation_small, xy_rotation_large, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_transform_xy_large"
        )
      );

      // small transformation around z axis
      e_HelixTransformZSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTransformRandom( z_translation_small, z_rotation_small, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_transform_z_small"
        )
      );

      // large transformation around z axis
      e_HelixTransformZLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_helix_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTransformRandom( z_translation_small, z_translation_large, z_rotation_small, z_rotation_large, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_transform_z_large"
        )
      );

      // construct xy flip which consists of a x flip and then a random rotation around Z axis between 0 and 360 degrees
      util::ShPtrList< coord::MoveInterface> moves_list;
      moves_list.PushBack( util::ShPtr< coord::MoveInterface>( x_flip.Clone()));
      moves_list.PushBack
      (
        util::ShPtr< coord::MoveInterface>( new coord::MoveRotateRandom( linal::Vector3D( 0.0, 0.0, 2.0 * math::g_Pi)))
      );
      coord::MoveCombine xy_flip( moves_list);

      // add xy and z flip
      e_HelixFlipXY = GetMutates().AddMutate( MutateProteinModelSSE( sp_locator_helix_random, util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( xy_flip), false)), MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_flip_xy"));
      e_HelixFlipZ = GetMutates().AddMutate( MutateProteinModelSSE( sp_locator_helix_random, util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( z_flip), false)), MutateTree::GetMutateTypeName( MutateTree::e_Helix) + "_flip_z"));

    /////////////////////////
    // single strand moves //
    /////////////////////////

      // move a strand next to a sheet
      e_StrandMoveNext = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_strand_random_wrapped,
          placement_next_to_sse, MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_move_next"
        )
      );

      // move the furthest strand next to a sheet
      e_StrandFurthestMoveNext = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_strand_furthest,
          placement_next_to_sse,
          MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_furthest_move_next"
        )
      );

      // move a strand next to a sheet
      e_StrandMoveSheet = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_strand_random_wrapped,
          placement_strand_next_to_sheet,
          MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_move_sheet"
        )
      );

      // move the furthest strand next to a sheet
      e_StrandFurthestMoveSheet = GetMutates().AddMutate
      (
        MutateProteinModelSSEMove
        (
          locator_strand_furthest,
          placement_strand_next_to_sheet,
          MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_furthest_move_sheet"
        )
      );

      // translation along z axis 0 to 2 angstroms
      e_StrandTranslateZSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_strand_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( z_translation_small, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_translate_z_small"
        )
      );

      // translation along z axis 2 to 4 angstroms
      e_StrandTranslateZLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_strand_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( z_translation_small, z_translation_large, true)), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_translate_z_large"
        )
      );

      // add flips
      e_StrandFlipX = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_strand_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( x_flip), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_flip_x"
        )
      );
      e_StrandFlipY = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_strand_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( y_flip), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_flip_y"
        )
      );
      e_StrandFlipZ = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_strand_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( z_flip), false)),
          MutateTree::GetMutateTypeName( MutateTree::e_Strand) + "_flip_z"
        )
      );

    ////////////////////
    // sse pair moves //
    ////////////////////

      // small translations along the hinge axis
      e_SSEPairTranslateNoHingeSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSEPairHinge
        (
          *collector_sse_pair,
          coord::MoveTranslateRandom( 2.0),
          false,
          MutateTree::GetMutateTypeName( MutateTree::e_SSEPair) + "_translate_no_hinge_small"
        )
      );

      // large translations along the hinge axis
      e_SSEPairTranslateNoHingeLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSEPairHinge
        (
          *collector_sse_pair,
          coord::MoveTranslateRandom( 2.0, 4.0),
          false,
          MutateTree::GetMutateTypeName( MutateTree::e_SSEPair) + "_translate_no_hinge_large"
        )
      );

      // translate sse pair small
      e_SSEPairTranslateSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSEPair
        (
          collector_sse_pair,
          math::Range< double>( 0.0, 1.0), // translation range
          math::Range< double>( 0.0, 0.0), // 0 to 15 degrees rotation
          MutateTree::GetMutateTypeName( MutateTree::e_SSEPair) + "_translate_small"
        )
      );

      // translate sse pair large
      e_SSEPairTranslateLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSEPair
        (
          collector_sse_pair,
          math::Range< double>( 1.0, 3.0), // translation range
          math::Range< double>( 0.0, 0.0), // 15 degree to 45 degree rotation
          MutateTree::GetMutateTypeName( MutateTree::e_SSEPair) + "_translate_large"
        )
      );

      // rotate sse pair small
      e_SSEPairRotateSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSEPair
        (
          collector_sse_pair,
          math::Range< double>( 0.0, 0.0), // translation range
          math::Range< double>( 0.0, math::g_Pi / 12.0), // 0 to 15 degrees rotation
          MutateTree::GetMutateTypeName( MutateTree::e_SSEPair) + "_rotate_small"
        )
      );

      // rotate sse pair large
      e_SSEPairRotateLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSEPair
        (
          collector_sse_pair,
          math::Range< double>( 0.0, 0.0), // translation range
          math::Range< double>( math::g_Pi / 18.0, math::g_Pi / 4.0), // 10 to 45 degree rotation
          MutateTree::GetMutateTypeName( MutateTree::e_SSEPair) + "_rotate_large"
        )
      );

      // translate sse pair small
      e_SSEPairTransformSmall = GetMutates().AddMutate
      (
        MutateProteinModelSSEPair
        (
          collector_sse_pair,
          math::Range< double>( 0.0, 1.0), // translation range
          math::Range< double>( 0.0, math::g_Pi / 18.0), // 0 to 10 degrees rotation
          MutateTree::GetMutateTypeName( MutateTree::e_SSEPair) + "_transform_small"
        )
      );

      // translate sse pair large
      e_SSEPairTransformLarge = GetMutates().AddMutate
      (
        MutateProteinModelSSEPair
        (
          collector_sse_pair,
          math::Range< double>( 1.0, 3.0), // translation range
          math::Range< double>( math::g_Pi / 18.0, math::g_Pi / 4.0), // 10 to 45 degrees rotation
          MutateTree::GetMutateTypeName( MutateTree::e_SSEPair) + "_transform_large"
        )
      );

    //////////////////////
    // helix pair moves //
    //////////////////////

      // small rotations around Z axis of the hinge
      e_HelixPairRotateZSmallNoHinge = GetMutates().AddMutate
      (
        MutateProteinModelSSEPairHinge
        (
          *collector_helix_pair,
          coord::MoveRotateRandom( z_rotation_small),
          false,
          MutateTree::GetMutateTypeName( MutateTree::e_HelixPair) + "_rotate_z_small_no_hinge"
        )
      );

      // small rotations around Z axis of the hinge, also rotation the hinge
      e_HelixPairRotateZSmallHinge = GetMutates().AddMutate
      (
        MutateProteinModelSSEPairHinge
        (
          *collector_helix_pair,
          coord::MoveRotateRandom( z_rotation_small),
          true,
          MutateTree::GetMutateTypeName( MutateTree::e_HelixPair) + "_rotate_z_small_hinge"
        )
      );

      // large rotations around Z axis of the hinge
      e_HelixPairRotateZLargeNoHinge = GetMutates().AddMutate
      (
        MutateProteinModelSSEPairHinge
        (
          *collector_helix_pair,
          coord::MoveRotateRandom( z_rotation_small, z_rotation_large),
          false,
          MutateTree::GetMutateTypeName( MutateTree::e_HelixPair) + "_rotate_z_large_no_hinge"
        )
      );

      // large rotations around Z axis of the hinge, also rotation the hinge
      e_HelixPairRotateZLargeHinge = GetMutates().AddMutate
      (
        MutateProteinModelSSEPairHinge
        (
          *collector_helix_pair,
          coord::MoveRotateRandom( z_rotation_small, z_rotation_large),
          true,
          MutateTree::GetMutateTypeName( MutateTree::e_HelixPair) + "_rotate_z_large_hinge"
        )
      );

    ////////////////////////
    // helix domain moves //
    ////////////////////////

      // initialize random helix domain locator of 2 to 4 helices
      assemble::LocatorDomainRandom helix_domain_locator( math::Range< size_t>( 2, 4), biol::GetSSTypes().HELIX);

      // shuffle helices in the domain
      e_HelixDomainShuffle = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          helix_domain_locator,
          MutateDomainShuffle( 2, DefaultFlags::GetFlagFitSwappedSSEs()->GetFlag()),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_shuffle"
        )
      );

      // translate the helix domain up to 2 angstroms
      e_HelixDomainTranslateSmall = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          helix_domain_locator,
          coord::MoveTranslateRandom( 2.0, true),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_translate_small"
        )
      );

      // translate the helix domain up to 6 angstroms
      e_HelixDomainTranslateLarge = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          helix_domain_locator,
          coord::MoveTranslateRandom( 2.0, 6.0, true),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_translate_large"
        )
      );

      // rotate the helix domain 0 to 15 degrees
      e_HelixDomainRotateSmall = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          helix_domain_locator,
          coord::MoveRotateRandom( math::g_Pi / 12.0, true),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_rotate_small"
        )
      );

      // rotate the helix domain 15 to 45 degrees
      e_HelixDomainRotateLarge = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          helix_domain_locator,
          coord::MoveRotateRandom( math::g_Pi / 12.0, math::g_Pi / 4.0, true),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_rotate_large"
        )
      );

      // transform the helix domain with 0 to 2 angstrom translation and 0 to 15 degree rotation
      e_HelixDomainTransformSmall = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          helix_domain_locator,
          coord::MoveTransformRandom( 2.0, math::g_Pi / 6.0),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_transform_small"
        )
      );

      // transform the helix domain with 2 to 6 angstrom translation and 15 to 45 degree rotation
      e_HelixDomainTransformLarge = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          helix_domain_locator,
          coord::MoveTransformRandom( 2.0, 6.0, math::g_Pi / 12.0, math::g_Pi / 4.0),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_transform_large"
        )
      );

      // flip helix domain externally around a randomly determined axis
      e_HelixDomainFlipExt = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateDomainFlip( xyz_axes_set, false),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_flip_ext"
        )
      );

      // flip each helix in the domain internally around an axis
      e_HelixDomainFlipInt = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateDomainFlip( xyz_axes_set, true, true),
          MutateTree::GetMutateTypeName( MutateTree::e_HelixDomain) + "_flip_int"
        )
      );

    /////////////////
    // sheet moves //
    /////////////////

      // rotate the sheet up to 15 degrees
      e_SheetRotateSmall = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          coord::MoveRotateRandom( math::g_Pi / 12.0, true),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_rotate_small"
        )
      );

      // rotate the sheet between 15 and 45 degrees
      e_SheetRotateLarge = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          coord::MoveRotateRandom( math::g_Pi / 12.0, math::g_Pi / 4.0, true),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_rotate_large"
        )
      );

      // translate sheet upto 2 angstroms
      e_SheetTranslateSmall = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          coord::MoveTranslateRandom( 2.0, true),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_translate_small"
        )
      );

      // translate sheet between 2 and 4 angstroms
      e_SheetTranslateLarge = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          coord::MoveTranslateRandom( 2.0, 4.0, true),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_translate_large"
        )
      );

      // transform sheet with upto 2 angstrom translation and 15 degrees rotation
      e_SheetTransformSmall = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          coord::MoveTransformRandom( 2.0, math::g_Pi / 12.0),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_transform_small"
        )
      );

      // transform sheet with upto 2 angstrom translation and 15 degrees rotation
      e_SheetTransformLarge = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          coord::MoveTransformRandom( 2.0, 4.0, math::g_Pi / 12.0, math::g_Pi / 4.0),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_transform_large"
        )
      );

      // pair strands
      e_SheetPairStrands = GetMutates().AddMutate
      (
        MutateProteinModelPairStrands( MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_pair_strands")
      );

      // move strand from one sheet to another
      e_SheetSwitchStrand = GetMutates().AddMutate
      (
        MutateProteinModelStrandSwitchSheet
        (
          placement_strand_next_to_sheet, MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_switch_strand"
        )
      );

      e_SheetFlipExt = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateDomainFlip( xyz_axes_set, false),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_flip_ext"
        )
      );

      e_SheetFlipInt = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateDomainFlip( xyz_axes_set, true, true),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_flip_int"
        )
      );

      e_Sheet_flip_int_sub = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateDomainFlip( xyz_axes_set, true, false, false),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_flip_int_sub"
        )
      );

      e_Sheet_flip_int_sub_diff = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateDomainFlip( xyz_axes_set, true, false, true),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_flip_int_sub_diff"
        )
      );

      // divide a sheet into two pieces
      e_Sheet_divide = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetDivide( 4, 2, false),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_divide"
        )
      );

      // divide a sheet into two pieces as a sandwich
      e_Sheet_divide_sandwich = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetDivide( 4, 2, true),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_divide_sandwich"
        )
      );

      // adjusts twist angles each pair strand in the selected sheet between 0 and 5 degrees
      e_Sheet_twist_small = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetTwist( math::Range< double>( 0.0, math::Angle::Radian( 2.0))),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_twist_small"
        )
      );

      // adjusts twist angles each pair strand in the selected sheet between 0 and 5 degrees
      e_Sheet_twist_large = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetTwist( math::Range< double>( 0.0, math::Angle::Radian( 10.0))),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_twist_large"
        )
      );

      // does 1 to 2 swaps within the sheet
      e_Sheet_shuffle = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateDomainShuffle( 2, DefaultFlags::GetFlagFitSwappedSSEs()->GetFlag()),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_shuffle"
        )
      );

      // build nr_rotations
      const math::Range< size_t> nr_rotations( 1, 3);
      const math::Range< size_t> subset_size( 2, 4);

      // cycles all the strands in the sheet by shifting 1 to 3 positions
      e_Sheet_cycle = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetCycle( true, false, false, nr_rotations),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_cycle"
        )
      );

      // cycles all the strands in the sheet shifting 1 to 3 positions while keeping pairwise orientations intact
      e_Sheet_cycle_intact = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetCycle( true, true, false, nr_rotations),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_cycle_intact"
        )
      );

      // cycles the positions, shifting 1 to 3 positions, of a subset (2 to 4) strands in a given sheet
      e_Sheet_cycle_subset = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetCycle( true, false, true, nr_rotations, subset_size),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_cycle_subset"
        )
      );

      // cycles the positions, shifting 1 to 3 positions, of a subset (2 to 4) strands in a given sheet
      // while keeping pairwise orientations intact
      e_Sheet_cycle_subset_intact = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetCycle( true, true, true, nr_rotations, subset_size),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_cycle_subset_intact"
        )
      );

      // fixes the hyrogen bonding pattern of strands in a sheet
      e_Sheet_register_fix = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetRegisterFix(),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_register_fix"
        )
      );

      // shifts the register of strands in a sheet by translations
      e_Sheet_register_shift = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetRegisterShift( storage::VectorND< 2, double>( 1.0, 0.0)),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_register_shift"
        )
      );

      // shifts the register of strands in a sheet by flips, half the translation and a flip around z axis
      e_Sheet_register_shift_flip = GetMutates().AddMutate
      (
        MutateProteinModelDomain
        (
          *collector_sheet,
          MutateSheetRegisterShift( storage::VectorND< 2, double>( 0.0, 1.0)),
          MutateTree::GetMutateTypeName( MutateTree::e_Sheet) + "_register_shift_flip"
        )
      );
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void DefaultMutates::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      // // get the pool using the empty model
      util::ShPtr< assemble::SSEPool> sp_pool
      (
        GetSetup().GetEmptyModel()->GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );
      BCL_Assert( sp_pool.IsDefined(), "pool is not initialized!");

      // get the average number of helices and strands, and calculate the ratios
      const storage::Pair< double, double> avg_sse_counts( sp_pool->CalculateAverageHelixStrandCounts());
      const double avg_helix_count( avg_sse_counts.First());
      const double avg_strand_count( avg_sse_counts.Second());
      const double avg_sse_count( avg_helix_count + avg_strand_count);
      const double avg_helix_ratio( avg_helix_count / avg_sse_count);
      const double avg_strand_ratio( avg_strand_count / avg_sse_count);

      BCL_MessageStd
      (
        "average helix/strand/sse counts from pool: " + util::Format()( avg_helix_count) + " / " +
         util::Format()( avg_strand_count) + " / " + util::Format()( avg_sse_count)
      );

      // initialize variables that are based on the counts
      const bool has_helix( avg_helix_count > 0);
      const bool has_strand( avg_strand_count > 0);
      const bool has_helixdomain( avg_helix_count >= 1.5);
      const bool has_helixpair( avg_helix_count >= 1.5);
      const bool has_sheet( avg_strand_count >= 1.5);

      // first choose probabilities for add, remove, swap, move
      MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_Add,         0.075);
      MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_Remove,      0.025);
      MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_Swap,        0.10);
      MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_SSE,         0.10);
      MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_Helix,       0.20 * avg_helix_ratio);
      MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_Strand,      0.20 * avg_strand_ratio);
      const double ssepair_weight(                                     0.20);
      MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_SSEPair,     ssepair_weight);
      const double domain_weight(                                      0.30);

      // if has a helix pair
      if( has_helixpair)
      {
        // update helix pair weight
        MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_HelixPair, ssepair_weight * avg_helix_ratio);
        MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_SSEPair,   ssepair_weight * avg_strand_ratio);
      }

      // update domain weights
      if( has_helixdomain)
      {
        MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_HelixDomain, domain_weight * avg_helix_ratio);
      }

      if( has_sheet)
      {
        MUTATE_TREE.SetDefaultMutateTypeProbability( MutateTree::e_Sheet,       domain_weight * avg_strand_ratio);
      }

      // adds
      MUTATE_TREE.SetDefaultMutateProbability( e_AddSSENextToSSE     , 0.75);
      MUTATE_TREE.SetDefaultMutateProbability( e_AddSSEShortLoop     , 0.25);
      MUTATE_TREE.SetDefaultMutateProbability( e_AddStrandNextToSheet, 1.00 * avg_strand_ratio);
      // remove
      MUTATE_TREE.SetDefaultMutateProbability( e_RemoveRandom        , 0.5);
      MUTATE_TREE.SetDefaultMutateProbability( e_RemoveUnpairedStrand, 0.5);

      // swap
      MUTATE_TREE.SetDefaultMutateProbability( e_SwapSSEs       , 0.8);
      MUTATE_TREE.SetDefaultMutateProbability( e_SwapSSEWithPool, 0.2);
      // if the pool is overlapping
      if( sp_pool->IsOverlapping())
      {
        MUTATE_TREE.SetDefaultMutateProbability( e_SwapSSEWithPoolOverlap, 0.2);
      }

      // single SSE moves
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEMoveNext,         3.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEMoveShortLoop,    3.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEFurthestMoveNext, 3.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEBendRamachandran, 1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEBendRandomSmall,  1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEBendRandomLarge,  1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETranslateSmall,   3.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETranslateXSmall,  1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETranslateYSmall,  1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETranslateZSmall,  1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETranslateLarge,   3.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETranslateXLarge,  1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETranslateYLarge,  1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETranslateZLarge,  1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSERotateSmall,      3.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSERotateXSmall,     1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSERotateYSmall,     1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSERotateZSmall,     1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSERotateLarge,      3.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSERotateXLarge,     1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSERotateYLarge,     1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSERotateZLarge,     1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETransformSmall,   6.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSETransformLarge,   6.0);

      // if secondary structure prediction were provided and sse resizing is not disabled
      if
      (
        sspred::Methods::GetFlagReadSSPredictions()->GetFlag() &&
        DefaultFlags::GetFlagEnableSSEResize()->GetFlag()
      )
      {
        // average number of non overlapping sses
        const double number_non_overlapping_long_sses( sp_pool->CountLongNonOverlappingSSEs( 30, 20));
        const double split_probability( ( 2.0 + 2 * number_non_overlapping_long_sses) / m_SSESplit.GetSize());

        // resize move
        MUTATE_TREE.SetDefaultMutateProbability( e_SSEResizeCTerm, 1.5);
        MUTATE_TREE.SetDefaultMutateProbability( e_SSEResizeNTerm, 1.5);

        // iterate over ss prediction methods
        for
        (
          storage::Map< sspred::Method, Mutate>::const_iterator
            sspred_itr( m_SSESplit.Begin()), sspred_itr_end( m_SSESplit.End());
          sspred_itr != sspred_itr_end;
          ++sspred_itr
        )
        {
          MUTATE_TREE.SetDefaultMutateProbability( sspred_itr->second, split_probability);
        }
      }

      // single helix moves
      if( has_helix)
      {
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixMoveNext,             2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixMoveShortLoop,        2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixFurthestMoveNext,     2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixTranslateXYSmall,     1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixTranslateXYLarge,     2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixTranslateZSmall,      1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixTranslateZLarge,      2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixRotateXYSmall,        1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixRotateXYLarge,        2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixRotateZSmall,         1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixRotateZLarge,         2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixTransformXYSmall,     1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixTransformXYLarge,     2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixTransformZSmall,      1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixTransformZLarge,      2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixFlipXY,               2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixFlipZ,                2.0);
      }

      // single strand moves
      if( has_strand)
      {
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandMoveNext,          2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandFurthestMoveNext,  2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandMoveSheet,         2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandFurthestMoveSheet, 2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandTranslateZSmall,   1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandTranslateZLarge,   2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandFlipX,             2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandFlipY,             2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_StrandFlipZ,             2.0);
      }

      // sse pair moves
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEPairTranslateNoHingeSmall, 1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEPairTranslateNoHingeLarge, 2.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEPairTranslateSmall,        1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEPairTranslateLarge,        2.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEPairRotateSmall,           1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEPairRotateLarge,           2.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEPairTransformSmall,        1.0);
      MUTATE_TREE.SetDefaultMutateProbability( e_SSEPairTransformLarge,        2.0);

      // helix pair moves
      if( has_helixpair)
      {
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixPairRotateZSmallNoHinge, 1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixPairRotateZSmallHinge,   1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixPairRotateZLargeNoHinge, 1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixPairRotateZLargeHinge,   1.0);
      }

      // helix domain moves
      if( has_helixdomain)
      {
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainShuffle,        8.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainTranslateSmall, 1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainTranslateLarge, 2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainRotateSmall,    1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainRotateLarge,    2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainTransformSmall, 1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainTransformLarge, 2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainFlipExt,        2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_HelixDomainFlipInt,        2.0);
      }

      // sheet moves
      if( has_sheet)
      {
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetRotateSmall,           1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetRotateLarge,           2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetTranslateSmall,        1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetTranslateLarge,        2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetTransformSmall,        1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetTransformLarge,        2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetPairStrands,           2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetSwitchStrand,          4.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetFlipExt,               2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_SheetFlipInt,               2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_flip_int_sub,         2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_flip_int_sub_diff,    2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_divide,               2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_divide_sandwich,      2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_twist_small,          2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_twist_large,          2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_shuffle,              8.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_cycle,                4.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_cycle_intact,         4.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_cycle_subset,         4.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_cycle_subset_intact,  4.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_register_fix,         1.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_register_shift,       2.0);
        MUTATE_TREE.SetDefaultMutateProbability( e_Sheet_register_shift_flip,  2.0);
      }

    }

    //! @brief initializes the mutate tree for the given protein model
    //! @param TREE mutate tree that shall be optimized
    //! @param MODEL protein model for which the mutate tree shall be optimized
    void DefaultMutates::InitializeMutateTree( MutateTree &TREE, const assemble::ProteinModel &MODEL) const
    {
      // get the pool of the protein model
      util::ShPtr< assemble::SSEPool> sp_pool
      (
        MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );

      // get the average number of helices and strands, and calculate the ratios
      const storage::Pair< double, double> avg_sse_counts
      (
        sp_pool.IsDefined() ? sp_pool->CalculateAverageHelixStrandCounts() : storage::Pair< double, double>( 5, 5)
      );
      const double avg_helix_count( avg_sse_counts.First());
      const double avg_strand_count( avg_sse_counts.Second());
      const double avg_sse_count( avg_helix_count + avg_strand_count);
      const double avg_helix_ratio( avg_helix_count / avg_sse_count);
      const double avg_strand_ratio( avg_strand_count / avg_sse_count);

      BCL_MessageStd
      (
        "average helix/strand/sse counts from pool: " + util::Format()( avg_helix_count) + " / " +
         util::Format()( avg_strand_count) + " / " + util::Format()( avg_sse_count)
      );

      // initialize variables that are based on the counts
      const bool has_helix( avg_helix_count > 0);
      const bool has_strand( avg_strand_count > 0);
      const bool has_helixdomain( avg_helix_count >= 1.5);
      const bool has_helixpair( avg_helix_count >= 1.5);
      const bool has_sheet( avg_strand_count >= 1.5);

      // first choose probabilities for add, remove, swap, move
      TREE.SetDefaultMutateTypeProbability( MutateTree::e_Add,         0.075);
      TREE.SetDefaultMutateTypeProbability( MutateTree::e_Remove,      0.025);
      TREE.SetDefaultMutateTypeProbability( MutateTree::e_Swap,        0.10);
      TREE.SetDefaultMutateTypeProbability( MutateTree::e_SSE,         0.10);
      TREE.SetDefaultMutateTypeProbability( MutateTree::e_Helix,       0.20 * avg_helix_ratio);
      TREE.SetDefaultMutateTypeProbability( MutateTree::e_Strand,      0.20 * avg_strand_ratio);
      const double ssepair_weight(                                     0.20);
      TREE.SetDefaultMutateTypeProbability( MutateTree::e_SSEPair,     ssepair_weight);
      const double domain_weight(                                      0.30);

      // if has a helix pair
      if( has_helixpair)
      {
        // update helix pair weight
        TREE.SetDefaultMutateTypeProbability( MutateTree::e_HelixPair, ssepair_weight * avg_helix_ratio);
        TREE.SetDefaultMutateTypeProbability( MutateTree::e_SSEPair,   ssepair_weight * avg_strand_ratio);
      }

      // update domain weights
      if( has_helixdomain)
      {
        TREE.SetDefaultMutateTypeProbability( MutateTree::e_HelixDomain, domain_weight * avg_helix_ratio);
      }
      if( has_sheet)
      {
        TREE.SetDefaultMutateTypeProbability( MutateTree::e_Sheet,       domain_weight * avg_strand_ratio);
      }

      // adds
      TREE.SetDefaultMutateProbability( e_AddSSENextToSSE     , 0.75);
      TREE.SetDefaultMutateProbability( e_AddSSEShortLoop     , 0.25);
      TREE.SetDefaultMutateProbability( e_AddStrandNextToSheet, 1.00 * avg_strand_ratio);
      // remove
      TREE.SetDefaultMutateProbability( e_RemoveRandom        , 0.5);
      TREE.SetDefaultMutateProbability( e_RemoveUnpairedStrand, 0.5);

      // swap
      TREE.SetDefaultMutateProbability( e_SwapSSEs       , 0.8);
      TREE.SetDefaultMutateProbability( e_SwapSSEWithPool, 0.2);
      // if the pool is overlapping
      if( !sp_pool.IsDefined() || sp_pool->IsOverlapping())
      {
        TREE.SetDefaultMutateProbability( e_SwapSSEWithPoolOverlap, 0.2);
      }

      // single SSE moves
      TREE.SetDefaultMutateProbability( e_SSEMoveNext,         3.0);
      TREE.SetDefaultMutateProbability( e_SSEMoveShortLoop,    3.0);
      TREE.SetDefaultMutateProbability( e_SSEFurthestMoveNext, 3.0);
      TREE.SetDefaultMutateProbability( e_SSEBendRamachandran, 1.0);
      TREE.SetDefaultMutateProbability( e_SSEBendRandomSmall,  1.0);
      TREE.SetDefaultMutateProbability( e_SSEBendRandomLarge,  1.0);
      TREE.SetDefaultMutateProbability( e_SSETranslateSmall,   3.0);
      TREE.SetDefaultMutateProbability( e_SSETranslateXSmall,  1.0);
      TREE.SetDefaultMutateProbability( e_SSETranslateYSmall,  1.0);
      TREE.SetDefaultMutateProbability( e_SSETranslateZSmall,  1.0);
      TREE.SetDefaultMutateProbability( e_SSETranslateLarge,   3.0);
      TREE.SetDefaultMutateProbability( e_SSETranslateXLarge,  1.0);
      TREE.SetDefaultMutateProbability( e_SSETranslateYLarge,  1.0);
      TREE.SetDefaultMutateProbability( e_SSETranslateZLarge,  1.0);
      TREE.SetDefaultMutateProbability( e_SSERotateSmall,      3.0);
      TREE.SetDefaultMutateProbability( e_SSERotateXSmall,     1.0);
      TREE.SetDefaultMutateProbability( e_SSERotateYSmall,     1.0);
      TREE.SetDefaultMutateProbability( e_SSERotateZSmall,     1.0);
      TREE.SetDefaultMutateProbability( e_SSERotateLarge,      3.0);
      TREE.SetDefaultMutateProbability( e_SSERotateXLarge,     1.0);
      TREE.SetDefaultMutateProbability( e_SSERotateYLarge,     1.0);
      TREE.SetDefaultMutateProbability( e_SSERotateZLarge,     1.0);
      TREE.SetDefaultMutateProbability( e_SSETransformSmall,   6.0);
      TREE.SetDefaultMutateProbability( e_SSETransformLarge,   6.0);

      // if secondary structure prediction were provided and sse resizing is not disabled
      if
      (
        sspred::Methods::GetFlagReadSSPredictions()->GetFlag() &&
        DefaultFlags::GetFlagEnableSSEResize()->GetFlag()
      )
      {
        // average number of non overlapping sses
        const double number_non_overlapping_long_sses
        (
          sp_pool.IsDefined() ? sp_pool->CountLongNonOverlappingSSEs( 30, 20) : 2.0
        );
        const double split_probability( ( 2.0 + 2 * number_non_overlapping_long_sses) / m_SSESplit.GetSize());

        // resize move
        TREE.SetDefaultMutateProbability( e_SSEResizeCTerm, 1.5);
        TREE.SetDefaultMutateProbability( e_SSEResizeNTerm, 1.5);

        // iterate over ss prediction methods
        for
        (
          storage::Map< sspred::Method, Mutate>::const_iterator
            sspred_itr( m_SSESplit.Begin()), sspred_itr_end( m_SSESplit.End());
          sspred_itr != sspred_itr_end;
          ++sspred_itr
        )
        {
          TREE.SetDefaultMutateProbability( sspred_itr->second, split_probability);
        }
      }

      // single helix moves
      if( has_helix)
      {
        TREE.SetDefaultMutateProbability( e_HelixMoveNext,             2.0);
        TREE.SetDefaultMutateProbability( e_HelixMoveShortLoop,        2.0);
        TREE.SetDefaultMutateProbability( e_HelixFurthestMoveNext,     2.0);
        TREE.SetDefaultMutateProbability( e_HelixTranslateXYSmall,     1.0);
        TREE.SetDefaultMutateProbability( e_HelixTranslateXYLarge,     2.0);
        TREE.SetDefaultMutateProbability( e_HelixTranslateZSmall,      1.0);
        TREE.SetDefaultMutateProbability( e_HelixTranslateZLarge,      2.0);
        TREE.SetDefaultMutateProbability( e_HelixRotateXYSmall,        1.0);
        TREE.SetDefaultMutateProbability( e_HelixRotateXYLarge,        2.0);
        TREE.SetDefaultMutateProbability( e_HelixRotateZSmall,         1.0);
        TREE.SetDefaultMutateProbability( e_HelixRotateZLarge,         2.0);
        TREE.SetDefaultMutateProbability( e_HelixTransformXYSmall,     1.0);
        TREE.SetDefaultMutateProbability( e_HelixTransformXYLarge,     2.0);
        TREE.SetDefaultMutateProbability( e_HelixTransformZSmall,      1.0);
        TREE.SetDefaultMutateProbability( e_HelixTransformZLarge,      2.0);
        TREE.SetDefaultMutateProbability( e_HelixFlipXY,               2.0);
        TREE.SetDefaultMutateProbability( e_HelixFlipZ,                2.0);
      }

      // single strand moves
      if( has_strand)
      {
        TREE.SetDefaultMutateProbability( e_StrandMoveNext,          2.0);
        TREE.SetDefaultMutateProbability( e_StrandFurthestMoveNext,  2.0);
        TREE.SetDefaultMutateProbability( e_StrandMoveSheet,         2.0);
        TREE.SetDefaultMutateProbability( e_StrandFurthestMoveSheet, 2.0);
        TREE.SetDefaultMutateProbability( e_StrandTranslateZSmall,   1.0);
        TREE.SetDefaultMutateProbability( e_StrandTranslateZLarge,   2.0);
        TREE.SetDefaultMutateProbability( e_StrandFlipX,             2.0);
        TREE.SetDefaultMutateProbability( e_StrandFlipY,             2.0);
        TREE.SetDefaultMutateProbability( e_StrandFlipZ,             2.0);
      }

      // sse pair moves
      TREE.SetDefaultMutateProbability( e_SSEPairTranslateNoHingeSmall, 1.0);
      TREE.SetDefaultMutateProbability( e_SSEPairTranslateNoHingeLarge, 2.0);
      TREE.SetDefaultMutateProbability( e_SSEPairTranslateSmall,        1.0);
      TREE.SetDefaultMutateProbability( e_SSEPairTranslateLarge,        2.0);
      TREE.SetDefaultMutateProbability( e_SSEPairRotateSmall,           1.0);
      TREE.SetDefaultMutateProbability( e_SSEPairRotateLarge,           2.0);
      TREE.SetDefaultMutateProbability( e_SSEPairTransformSmall,        1.0);
      TREE.SetDefaultMutateProbability( e_SSEPairTransformLarge,        2.0);

      // helix pair moves
      if( has_helixpair)
      {
        TREE.SetDefaultMutateProbability( e_HelixPairRotateZSmallNoHinge, 1.0);
        TREE.SetDefaultMutateProbability( e_HelixPairRotateZSmallHinge,   1.0);
        TREE.SetDefaultMutateProbability( e_HelixPairRotateZLargeNoHinge, 1.0);
        TREE.SetDefaultMutateProbability( e_HelixPairRotateZLargeHinge,   1.0);
      }

      // helix domain moves
      if( has_helixdomain)
      {
        TREE.SetDefaultMutateProbability( e_HelixDomainShuffle,        8.0);
        TREE.SetDefaultMutateProbability( e_HelixDomainTranslateSmall, 1.0);
        TREE.SetDefaultMutateProbability( e_HelixDomainTranslateLarge, 2.0);
        TREE.SetDefaultMutateProbability( e_HelixDomainRotateSmall,    1.0);
        TREE.SetDefaultMutateProbability( e_HelixDomainRotateLarge,    2.0);
        TREE.SetDefaultMutateProbability( e_HelixDomainTransformSmall, 1.0);
        TREE.SetDefaultMutateProbability( e_HelixDomainTransformLarge, 2.0);
        TREE.SetDefaultMutateProbability( e_HelixDomainFlipExt,        2.0);
        TREE.SetDefaultMutateProbability( e_HelixDomainFlipInt,        2.0);
      }

      // sheet moves
      if( has_sheet)
      {
        TREE.SetDefaultMutateProbability( e_SheetRotateSmall,           1.0);
        TREE.SetDefaultMutateProbability( e_SheetRotateLarge,           2.0);
        TREE.SetDefaultMutateProbability( e_SheetTranslateSmall,        1.0);
        TREE.SetDefaultMutateProbability( e_SheetTranslateLarge,        2.0);
        TREE.SetDefaultMutateProbability( e_SheetTransformSmall,        1.0);
        TREE.SetDefaultMutateProbability( e_SheetTransformLarge,        2.0);
        TREE.SetDefaultMutateProbability( e_SheetPairStrands,           2.0);
        TREE.SetDefaultMutateProbability( e_SheetSwitchStrand,          4.0);
        TREE.SetDefaultMutateProbability( e_SheetFlipExt,               2.0);
        TREE.SetDefaultMutateProbability( e_SheetFlipInt,               2.0);
        TREE.SetDefaultMutateProbability( e_Sheet_flip_int_sub,         2.0);
        TREE.SetDefaultMutateProbability( e_Sheet_flip_int_sub_diff,    2.0);
        TREE.SetDefaultMutateProbability( e_Sheet_divide,               2.0);
        TREE.SetDefaultMutateProbability( e_Sheet_divide_sandwich,      2.0);
        TREE.SetDefaultMutateProbability( e_Sheet_twist_small,          2.0);
        TREE.SetDefaultMutateProbability( e_Sheet_twist_large,          2.0);
        TREE.SetDefaultMutateProbability( e_Sheet_shuffle,              8.0);
        TREE.SetDefaultMutateProbability( e_Sheet_cycle,                4.0);
        TREE.SetDefaultMutateProbability( e_Sheet_cycle_intact,         4.0);
        TREE.SetDefaultMutateProbability( e_Sheet_cycle_subset,         4.0);
        TREE.SetDefaultMutateProbability( e_Sheet_cycle_subset_intact,  4.0);
        TREE.SetDefaultMutateProbability( e_Sheet_register_fix,         1.0);
        TREE.SetDefaultMutateProbability( e_Sheet_register_shift,       2.0);
        TREE.SetDefaultMutateProbability( e_Sheet_register_shift_flip,  2.0);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DefaultMutates::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DefaultMutates::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
