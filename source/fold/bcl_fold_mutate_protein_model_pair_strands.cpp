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
#include "fold/bcl_fold_mutate_protein_model_pair_strands.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_placement_strand_next_to_sheet.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelPairStrands::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelPairStrands())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a scheme
    MutateProteinModelPairStrands::MutateProteinModelPairStrands() :
      m_Scheme( GetStaticClassName< MutateProteinModelPairStrands>())
    {
    }

    //! @brief constructor from a scheme
    //! @param SCHEME scheme of this mutate
    MutateProteinModelPairStrands::MutateProteinModelPairStrands( const std::string &SCHEME) :
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelPairStrands
    MutateProteinModelPairStrands *MutateProteinModelPairStrands::Clone() const
    {
      return new MutateProteinModelPairStrands( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelPairStrands::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelPairStrands::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize static failure return value
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // if no or only one strand in the model
      if( PROTEIN_MODEL.GetNumberSSE( biol::GetSSTypes().STRAND) <= 1)
      {
        BCL_MessageVrb( "One or less strands in this model yet, skipping pairing strands");
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }
      // collect all the Sheets in the model
      util::ShPtrVector< assemble::Domain> sheets_unfiltered( assemble::CollectorSheet().Collect( PROTEIN_MODEL));
      util::ShPtrVector< assemble::Domain> sheets;

      // iterate over the sheet
      for
      (
        util::ShPtrVector< assemble::Domain>::iterator
          sheet_itr( sheets_unfiltered.Begin()), sheet_itr_end( sheets_unfiltered.End());
        sheet_itr != sheet_itr_end; ++sheet_itr
      )
      {
        // if the sheet is of type beta-barrel then remove it
        if( ( *sheet_itr)->GetTopology()->GetType() == assemble::Topology::e_Sheet)
        {
          sheets.PushBack( *sheet_itr);
        }
      }

      // if there is only one sheet ( which will happen with one sheet with at least two strands in it
      if( sheets.GetSize() == 1)
      {
        BCL_MessageVrb( "Only one sheet found in this model yet, skipping pairing strands");
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // initialize vector to store individual strands that are not in any sheet
      util::ShPtrVector< assemble::Domain> sheets_with_single_strands;

      // iterate over the sheets
      for
      (
        util::ShPtrVector< assemble::Domain>::iterator sheet_itr( sheets.Begin()), sheet_itr_end( sheets.End());
        sheet_itr != sheet_itr_end; ++sheet_itr
      )
      {
        // if there is only one strand in the Sheet
        if( ( *sheet_itr)->GetNumberSSEs() == 1)
        {
          sheets_with_single_strands.PushBack( *sheet_itr);
        }
      }

      // if there are no unpaired strands
      if( sheets_with_single_strands.IsEmpty())
      {
        BCL_MessageVrb( "No unpaired strands in this model yet, skipping pairing strands");
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // now we know that there are at least one unpaired strands
      // so pick one randomly
      util::ShPtrVector< assemble::Domain>::iterator random_itr
      (
        random::GetGlobalRandom().Iterator
        (
          sheets_with_single_strands.Begin(),
          sheets_with_single_strands.End(),
          sheets_with_single_strands.GetSize()
        )
      );

      // make sure the random iterator returned correctly
      BCL_Assert
      (
        random_itr != sheets_with_single_strands.End(),
        "Error occurred in picking a random sheet with a single unpaired strand"
      );

      // make a copy of this SSE
      util::ShPtr< assemble::SSE> new_sse( ( *random_itr)->GetSSEs().FirstElement()->Clone());

      // idealize the new SSE
      new_sse->SetToIdealConformationAtOrigin();

      // now we have to pick the sheet to move this pair next to
      // first create the iterator to pick the sheet
      util::ShPtrVector< assemble::Domain>::const_iterator random_sheet_itr;

      // if there are more unpaired strands
      if( sheets_with_single_strands.GetSize() >= 2)
      {
        // remove the picked strand
        sheets_with_single_strands.Remove( random_itr);

        // pick another sheet with single strand randomly
        random_sheet_itr =
          random::GetGlobalRandom().Iterator
          (
            sheets_with_single_strands.Begin(),
            sheets_with_single_strands.End(),
            sheets_with_single_strands.GetSize()
          );
        // make sure the random iterator returned correctly
        BCL_Assert
        (
          random_sheet_itr != sheets_with_single_strands.End(),
          "Error occurred in picking a random sheet to pair the strand with"
        );

      }
      // otherwise, if there was only one unpaired strand then we have to pick one of the sheets with multiple strands
      else
      {
        // remove the picked single strand from sheets list
        sheets.Remove( std::find( sheets.Begin(), sheets.End(), *random_itr));

        // now get a random sheet to place next to
        random_sheet_itr =
          random::GetGlobalRandom().Iterator
          (
            sheets.Begin(),
            sheets.End(),
            sheets.GetSize()
          );
        // make sure the random iterator returned correctly
        BCL_Assert
        (
          random_sheet_itr != sheets.End(),
          "Error occurred in picking a random sheet to pair the strand with"
        );
      }

      // calculate the placement
      storage::Pair< math::TransformationMatrix3D, bool> transformation
      (
        PlacementStrandNextToSheet().Place( *new_sse, **random_sheet_itr)
      );

      // make sure the placement was calculated correctly
      BCL_Assert
      (
        transformation.Second(),
        "The placement failed for placing unpaired strand next to sheet"
      );

      // apply the transformation to new sse
      new_sse->Transform( transformation.First());

      // make a new model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace the one in model with the new SSE
      new_model->Replace( new_sse);

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelPairStrands::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelPairStrands::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
