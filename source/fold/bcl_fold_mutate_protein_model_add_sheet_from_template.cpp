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
#include "fold/bcl_fold_mutate_protein_model_add_sheet_from_template.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sses_random.h"
#include "assemble/bcl_assemble_sheet_template_handler.h"
#include "assemble/bcl_assemble_sse_pool.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelAddSheetFromTemplate::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelAddSheetFromTemplate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from placer
    //! @param PLACER object used to place the sheet in the protein model
    MutateProteinModelAddSheetFromTemplate::MutateProteinModelAddSheetFromTemplate
    (
      const PlacementInterface< assemble::Domain, assemble::ProteinModel> &PLACER
    ) :
      m_Placer( util::CloneToShPtr( PLACER))
    {
    }

    //! @brief returns a pointer to a new MutateProteinModelAddSheetFromTemplate
    //! @return pointer to a new MutateProteinModelAddSheetFromTemplate
    MutateProteinModelAddSheetFromTemplate *MutateProteinModelAddSheetFromTemplate::Clone() const
    {
      return new MutateProteinModelAddSheetFromTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateProteinModelAddSheetFromTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateProteinModelAddSheetFromTemplate::GetScheme() const
    {
      static const std::string s_scheme( "add_sheet_from_template");
      return s_scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief adds a sheet to the given protein model and returns the mutated model
    //! @detail a random number of strands is selected from the SSE pool of the given model and fitted to a sheet
    //! template in the library before being added to the given protein model
    //! @param MODEL protein model to add the sheet to
    //! @return the mutated protein model with additional information regarding the mutation
    math::MutateResult< assemble::ProteinModel> MutateProteinModelAddSheetFromTemplate::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // undefined protein model which will be returned if application of the mutate fails
      util::ShPtr< assemble::ProteinModel> sp_result_model;

      // get the strands from the SSE pool which are not overlapping with SSEs in the protein model
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );
      util::SiPtrList< const assemble::SSE> non_overlapping_strands
      (
        sp_pool->GetNonOverlappingSSEs( MODEL, biol::GetSSTypes().STRAND)
      );

      //! a sheet can only be formed if there are at least two strands available for insertion
      const size_t number_strands( non_overlapping_strands.GetSize());
      if( number_strands < 2)
      {
        return math::MutateResult< assemble::ProteinModel>( sp_result_model, *this);
      }

      // randomly choose how many strands from the pool to form the sheet from and pick them from the list
      const size_t num_strands_to_pick( random::GetGlobalRandom().Random< size_t>( 2, number_strands));
      const storage::Set< biol::SSType> sse_types_to_pick( biol::GetSSTypes().STRAND);
      const util::SiPtrList< const assemble::SSE> picked_strands
      (
        assemble::PickSSEsRandom( sse_types_to_pick, num_strands_to_pick).Pick( non_overlapping_strands)
      );

      // pick a random template from the library that has a proper number of strands of the right sequence length
      util::SiPtrVector< const assemble::SSE> targets;
      for
      (
        util::SiPtrList< const assemble::SSE>::const_iterator it( picked_strands.Begin()), it_end( picked_strands.End());
        it != it_end;
        ++it
      )
      {
        targets.PushBack( *it);
      }
      const assemble::FoldTemplate &fold_template( assemble::SheetTemplateHandler::GetRandomTemplate( targets));

      // application of the mutate failed if no adequate template could be found
      if( fold_template.GetGeometries().IsEmpty())
      {
        return math::MutateResult< assemble::ProteinModel>( sp_result_model, *this);
      }

      // fit the selected strands to the template and compute where to place the sheet in the model
      util::ShPtr< assemble::Domain> sp_new_sheet( fold_template.FitSSEs( targets).Clone());
      const storage::Pair< math::TransformationMatrix3D, bool> placement( m_Placer->Place( *sp_new_sheet, MODEL));
      const math::TransformationMatrix3D &trans_matrix( placement.First());
      sp_new_sheet->Transform( trans_matrix);
      const util::SiPtrVector< const assemble::SSE> sses_to_insert( sp_new_sheet->GetSSEs());

      // add the sheet to the model
      sp_result_model = util::ShPtr< assemble::ProteinModel>( MODEL.Clone());
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator it( sses_to_insert.Begin()), it_end( sses_to_insert.End());
        it != it_end;
        ++it
       )
      {
        sp_result_model->Insert( util::CloneToShPtr( **it));
      }

      return math::MutateResult< assemble::ProteinModel>( sp_result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read object from input stream
    //! @param ISTREAM input stream to read object from
    //! @return input stream which was read from
    std::istream &MutateProteinModelAddSheetFromTemplate::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write object into  output stream
    //! @param OSTREAM output stream to write object into
    //! @param INDENT number of indentations to separate members
    //! @return output stream object was written into
    std::ostream &MutateProteinModelAddSheetFromTemplate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
