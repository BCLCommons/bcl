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
#include "fold/bcl_fold_mutate_protein_model_sse_pair_align_and_pull.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sse_paired.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_pair_align_and_pull.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSEPairAlignAndPull :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSEPairAlignAndPull *Clone() const
    { return new ExampleFoldMutateProteinModelSSEPairAlignAndPull( *this);}

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
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

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

      // create mutate object "mutate_sse_pair" from max translation and rotation and CollectorSSEPaired
      fold::MutateProteinModelSSEPairAlignAndPull mutate_sse_pair
      (
        sp_sse_pair_collect,
        false
      );

      // mutate the protein model (random sse will be moved
      assemble::ProteinModel mutated_model( *mutate_sse_pair( protein_model).GetArgument());

      const double rmsd_start( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, protein_model, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      const double rmsd_mutate( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      BCL_MessageStd( "rmsd before mutating sses " + util::Format()( rmsd_start));
      BCL_MessageStd( "rmsd after mutating sses " + util::Format()( rmsd_mutate));

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write mutated_sse_pair.pdb");
      Proteins::WriteModelToPDB( mutated_model, AddExampleOutputPathToFilename( mutate_sse_pair, "mutated_sse_pair.pdb"));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSEPairAlignAndPull

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSEPairAlignAndPull::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSEPairAlignAndPull())
  );

} // namespace bcl
