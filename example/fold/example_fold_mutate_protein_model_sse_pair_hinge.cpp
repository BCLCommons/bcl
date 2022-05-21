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
#include "fold/bcl_fold_mutate_protein_model_sse_pair_hinge.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sse_paired.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "coord/bcl_coord_move_rotate_defined.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_pair_hinge.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSEPairHinge :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSEPairHinge *Clone() const
    { return new ExampleFoldMutateProteinModelSSEPairHinge( *this);}

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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      //build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create sse pair collector
      util::ShPtr
      <
        find::CollectorInterface
        <
          storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface
        >
      >
      sp_sse_pair_collect
      (
        new assemble::CollectorSSEPaired
        (
          assemble::GetSSEGeometryPackingPickers().e_BestInteractionWeight,
          storage::Set< contact::Type>( contact::GetTypes().Begin(), contact::GetTypes().End()),
          12.0, // distance
          true  // orthogonal
        )
      );

      // move rotate defined
      util::ShPtr< coord::MoveInterface> rotate_defined
      (
        new coord::MoveRotateDefined( math::g_Pi / 5, coord::GetAxes().e_Z)
      );

      // create mutate object "mutate_sse_pair_with_hinge" from max translation and rotation and CollectorSSEPaired
      fold::MutateProteinModelSSEPairHinge mutate_sse_pair_with_hinge
      (
        sp_sse_pair_collect,
        rotate_defined,
        true // boolean to move hinge too
      );

      // create mutate object "mutate_sse_pair_without_hinge" from max translation and rotation and CollectorSSEPaired
      fold::MutateProteinModelSSEPairHinge mutate_sse_pair_without_hinge
      (
        sp_sse_pair_collect,
        rotate_defined,
        false // boolean to move hinge too
      );

      // mutate the protein model (random sse will be moved
      assemble::ProteinModel mutated_model_a( *mutate_sse_pair_with_hinge( protein_model).GetArgument());

      // mutate the protein model (random sse will be moved
      assemble::ProteinModel mutated_model_b( *mutate_sse_pair_without_hinge( protein_model).GetArgument());

      const double rmsd_start( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, protein_model, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      const double rmsd_mutate_a( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model_a, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      const double rmsd_mutate_b( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model_b, biol::GetAtomTypes().GetBackBoneAtomTypes()));

      BCL_MessageStd( "rmsd before mutating sses " + util::Format()( rmsd_start));
      BCL_MessageStd( "rmsd after mutating sses with hinge " + util::Format()( rmsd_mutate_a));
      BCL_MessageStd( "rmsd after mutating sses without hinge" + util::Format()( rmsd_mutate_b));

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write mutate_sse_pair_with_hinge.pdb");
      std::string out_filename
      (
        AddExampleOutputPathToFilename( mutate_sse_pair_with_hinge, "mutate_sse_pair_with_hinge.pdb")
      );
      Proteins::WriteModelToPDB( mutated_model_a, out_filename);

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write mutated_domain_without_hinge.pdb");
      out_filename = AddExampleOutputPathToFilename( mutate_sse_pair_with_hinge, "mutate_sse_pair_without_hinge.pdb");
      Proteins::WriteModelToPDB( mutated_model_b, out_filename);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSEPairHinge

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSEPairHinge::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSEPairHinge())
  );

} // namespace bcl
