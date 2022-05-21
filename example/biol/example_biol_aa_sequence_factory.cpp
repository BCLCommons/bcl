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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "biol/bcl_biol_aa_sequence_factory.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_phi_psi.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_sequence_factory.cpp
  //!
  //! @author weinerbe
  //! @date Jan 25, 2011
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAASequenceFactory :
    public ExampleInterface
  {
  public:

    ExampleBiolAASequenceFactory *Clone() const
    {
      return new ExampleBiolAASequenceFactory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // get helix
      util::ShPtr< assemble::SSE> sp_helix
      (
        Proteins::GetSSE( pdb_filename, 'A', 23, 34, biol::GetAAClasses().e_AABackBone)
      );
      const biol::AASequencePhiPsi phi_psi_helix( *sp_helix);

      // get 2 strands
      util::ShPtr< assemble::SSE> sp_strand_large
      (
        Proteins::GetSSE( pdb_filename, 'A', 10, 17, biol::GetAAClasses().e_AABackBone)
      );
      util::ShPtr< assemble::SSE> sp_strand_small
      (
        Proteins::GetSSE( pdb_filename, 'A', 40, 45, biol::GetAAClasses().e_AABackBone)
      );
      const biol::AASequencePhiPsi phi_psi_strand_small( *sp_strand_small);
      const biol::AASequencePhiPsi phi_psi_strand_large( *sp_strand_large);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      biol::AASequenceFactory def_construct;

    ////////////////
    // operations //
    ////////////////

      // split helix
      util::ShPtr< assemble::SSE> sp_helix_left( new assemble::SSE( sp_helix->SubSequence( 0, 3), sp_helix->GetType()));
      util::ShPtr< assemble::SSE> sp_helix_right( new assemble::SSE( sp_helix->SubSequence( 3, 3), sp_helix->GetType()));

      // transformation for prepending
      const math::TransformationMatrix3D trans_prepend
      (
        biol::AASequenceFactory::TransformationPrepend
        (
          *sp_helix_left->GetLastAA(),
          *sp_helix_right,
          sp_helix_right->GetFirstAA()->CalculatePhi( sp_helix_left->GetLastAA()->GetAtom( biol::GetAtomTypes().C))
        )
      );

      BCL_MessageVrb
      (
        "transformation for prepending left to right helix: " + sp_helix_left->GetIdentification() + '\t' +
        sp_helix_right->GetIdentification() + '\n' + util::Format()( trans_prepend)
      );

      // transformation for appending
      const math::TransformationMatrix3D trans_append
      (
        biol::AASequenceFactory::TransformationAppend
        (
          *sp_helix_left,
          *sp_helix_right->GetFirstAA(),
          sp_helix_right->GetFirstAA()->CalculatePhi( sp_helix_left->GetLastAA()->GetAtom( biol::GetAtomTypes().C))
        )
      );

      BCL_MessageVrb
      (
        "transformation for appending right to left helix: " + sp_helix_right->GetIdentification() + '\t' +
        sp_helix_left->GetIdentification() + '\n' + util::Format()( trans_append)
      );

      // fit the helix to itself after idealizing
      sp_helix->SetToIdealConformationAtOrigin();
      biol::AASequenceFactory::FitSequence( *sp_helix, phi_psi_helix, biol::GetSSTypes().HELIX);
      sp_helix->SetGeometry();

      // fit the large strand to the small strand phi/psis
      sp_strand_large->SetToIdealConformationAtOrigin();
      biol::AASequenceFactory::FitSequence( *sp_strand_large, phi_psi_strand_small, biol::GetSSTypes().STRAND);
      sp_strand_large->SetGeometry();

      // fit the small strand to the large strand phi/psis
      sp_strand_small->SetToIdealConformationAtOrigin();
      biol::AASequenceFactory::FitSequence( *sp_strand_small, phi_psi_strand_large, biol::GetSSTypes().STRAND);
      sp_strand_small->SetGeometry();

      // write out the model
      assemble::Chain test_chain( *( protein_model.GetEmptyChains().FirstElement()));
      test_chain.Insert( sp_helix);
      test_chain.Insert( sp_strand_large);
      test_chain.Insert( sp_strand_small);
      Proteins::WriteChainToPDB( test_chain, AddExampleOutputPathToFilename( def_construct, "fit_sequences.pdb"));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAASequenceFactory

  const ExampleClass::EnumType ExampleBiolAASequenceFactory::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAASequenceFactory())
  );

} // namespace bcl
