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
#include "fold/bcl_fold_mutate_protein_model_strand_switch_sheet.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelStrandSwitchSheet::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelStrandSwitchSheet())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelStrandSwitchSheet::MutateProteinModelStrandSwitchSheet() :
      m_Placement(),
      m_Scheme( GetStaticClassName< MutateProteinModelStrandSwitchSheet>())
    {
    }

    //! @brief constructor from a placement and a scheme
    //! @param PLACEMENT reference to the PlacementStrandNextToSheet to be used
    //! @param SCHEME Scheme to be used
    MutateProteinModelStrandSwitchSheet::MutateProteinModelStrandSwitchSheet
    (
      const PlacementStrandNextToSheet &PLACEMENT,
      const std::string &SCHEME
    ) :
      m_Placement( PLACEMENT),
      m_Scheme( SCHEME)
    {

    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelStrandSwitchSheet
    MutateProteinModelStrandSwitchSheet *MutateProteinModelStrandSwitchSheet::Clone() const
    {
      return new MutateProteinModelStrandSwitchSheet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelStrandSwitchSheet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an ProteinModel and returning a mutated ProteinModel
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelStrandSwitchSheet::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // make an empty return pointer
      static util::ShPtr< assemble::ProteinModel> s_empty_result;

      // initialize the static sheet collector
      static const assemble::CollectorSheet s_sheet_collector;

      // collect all the sheets in the model
      const util::ShPtrVector< assemble::Domain> sheet_vector( s_sheet_collector.Collect( PROTEIN_MODEL));

      // if there are not at least two sheets
      if( sheet_vector.GetSize() < 2)
      {
        // then this move can't be applied
        return math::MutateResult< assemble::ProteinModel>( s_empty_result, *this);
      }

      // create a SiPtrVector of sheets
      util::SiPtrVector< const assemble::Domain> possible_sheets( sheet_vector);

      // pick two randomly
      const util::SiPtr< const assemble::Domain> receiver_sheet( possible_sheets.RemoveRandomElement());
      const util::SiPtr< const assemble::Domain> donator_sheet( possible_sheets.RemoveRandomElement());

      // now pick one of the edge from the donator sheet
      util::SiPtr< const assemble::SSEGeometryInterface> sse_to_move
      (
        random::GetGlobalRandom().Boolean() ?
          donator_sheet->GetTopology()->GetElements().FirstElement() :
          donator_sheet->GetTopology()->GetElements().LastElement()
      );

      // make a copy of the SSE
      util::ShPtr< assemble::SSE> new_sse( sse_to_move->Clone());
      // check the dynamic cast
      BCL_Assert( new_sse.IsDefined(), "The dynamic cast has failed from SSEGeometryInterface to SSE")

      // get the placement
      storage::Pair< math::TransformationMatrix3D, bool> placement
      (
        m_Placement.Place( *new_sse, *receiver_sheet)
      );

      // if the placement failed
      if( !placement.Second())
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_result, *this);
      }

      // build up the transformation that moves it to origin and then moves to the placement determined
      math::TransformationMatrix3D transformation( math::Inverse( new_sse->GetOrientation()));
      transformation( placement.First());

      // apply the transformation
      new_sse->Transform( transformation);

      // make a copy of the model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace the sse
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
    std::istream &MutateProteinModelStrandSwitchSheet::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Placement, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelStrandSwitchSheet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Placement, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
