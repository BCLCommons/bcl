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
#include "fold/bcl_fold_mutate_protein_model_sse_add.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "biol/bcl_biol_membrane.h"
#include "fold/bcl_fold_mutate_protein_model_sse_remove.h"
#include "fold/bcl_fold_mutate_sse_bend_ramachandran.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEAdd::MutateProteinModelSSEAdd() :
      m_SSEPoolPicker(),
      m_Placement(),
      m_Scheme( GetStaticClassName< MutateProteinModelSSEAdd>())
    {
    }

    //! @brief constructor from a PickCriteriaInterface, PlacementInterface and a scheme
    //! @param PICKER picker to be used for picking sses from the pool
    //! @param PLACEMENT placement to place the picked sse in the model
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEAdd::MutateProteinModelSSEAdd
    (
      const find::PickCriteriaInterface
      <
        util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
      > &PICKER,
      const PlacementInterface< assemble::SSE, assemble::ProteinModel> &PLACEMENT,
      const std::string &SCHEME
    ) :
      m_SSEPoolPicker( PICKER.Clone()),
      m_Placement( PLACEMENT),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a PickCriteriaInterface, PlacementInterface amd a scheme
    //! @param SP_PICKER ShPtr to picker to be used for picking sses from the pool
    //! @param SP_PLACEMENT ShPtr to placement to place the picked sse in the model
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEAdd::MutateProteinModelSSEAdd
    (
      const util::ShPtr
      <
        find::PickCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        >
      > &SP_PICKER,
      const util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > &SP_PLACEMENT,
      const std::string &SCHEME
    ) :
      m_SSEPoolPicker( *SP_PICKER),
      m_Placement( *SP_PLACEMENT),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSEAdd
    MutateProteinModelSSEAdd *MutateProteinModelSSEAdd::Clone() const
    {
      return new MutateProteinModelSSEAdd( *this);
    };

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEAdd::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelSSEAdd)
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelSSEAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEAdd::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

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

      // choose a random sse from the SSEs that do not overlap with SSEs in PROTEIN_MODEL
      util::SiPtrList< const assemble::SSE> non_overlapping_sses( sp_pool->GetNonOverlappingSSEs( PROTEIN_MODEL));

      // if there are no overlapping sses left display messsage and return
      if( non_overlapping_sses.IsEmpty())
      {
        // warn user
        BCL_MessageVrb( "skipping adding, no non-overlapping sses left in pool");
//        static fold::MutateProteinModelSSERemove s_remover
//        (
//          MutateProteinModelSSERemove
//          (
//            util::ShPtr< assemble::LocatorSSERandom>( new assemble::LocatorSSERandom()),
//            "remove_random"
//          )
//        );
//        math::MutateResult< assemble::ProteinModel> result( s_remover( PROTEIN_MODEL));
//        if( result.GetArgument().IsDefined())
//        {
//          return operator()( *result.GetArgument());
//        }

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // called the picker to pick an sse from the pool to be placed in the protein model
      const util::SiPtr< const assemble::SSE> sse_from_pool( m_SSEPoolPicker->Pick( non_overlapping_sses, PROTEIN_MODEL));

      // if pick could not find one just return
      if( !sse_from_pool.IsDefined())
      {
        // output the picker that failed
        BCL_MessageVrb( "undefined sse picked by " + util::Format()( *m_SSEPoolPicker));

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // report selected sse from pool to be added
      BCL_MessageVrb( "sse from pool to be added " + sse_from_pool->GetIdentification());

      // initialize transformation
      storage::Pair< math::TransformationMatrix3D, bool> transformation;

      // if protein model is empty insert this sse at the center
      if( PROTEIN_MODEL.GetSSEs().IsEmpty())
      {
        transformation = m_Placement->Place( *sse_from_pool);
      }
      // if PROTEIN MODEL is not empty
      else
      {
        // get the placement
        transformation = m_Placement->Place( *sse_from_pool, PROTEIN_MODEL);
      }

      // if the placement was not successful in finding a good placement
      if( !transformation.Second())
      {
        BCL_MessageStd( "The provided transformation is not successful, skipping the add!");

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // make a copy of the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // get membrane for current protein model
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // create a copy of the selected sse assuming the sse from the pool is already
      // in the ideal conformation and in origin
      util::ShPtr< assemble::SSE> this_sse( sse_from_pool->Clone());

      // transform with the assemble transformation matrix3D
      this_sse->Transform( transformation.First());

      static MutateSSEBendRamachandran s_bender( math::Range< size_t>( 1000, 1000), biol::AASequenceFlexibility::e_Bidirectional);
      this_sse = s_bender( *this_sse, sp_membrane).GetArgument();

      // print sses before the add and after the add
      BCL_MessageDbg( "#SSEs in model before: " + util::Format()( new_model->GetNumberSSEs()));

      // add to protein model
      new_model->Insert( this_sse);

      // print sses before the add and after the add
      BCL_MessageDbg( "#SSEs in model after: " + util::Format()( new_model->GetNumberSSEs()));

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
