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
#include "fold/bcl_fold_mutate_protein_model_sse_remove.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_remove.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSERemove :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSERemove *Clone() const
    { return new ExampleFoldMutateProteinModelSSERemove( *this);}

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
      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create mutate object that removes sse A 1 7
      fold::MutateProteinModelSSERemove mutate( util::CloneToShPtr( assemble::LocatorSSE( 'A', 1, 7)));
      // create mutate object that removes random sse
      fold::MutateProteinModelSSERemove mutate_random( util::CloneToShPtr( assemble::LocatorSSERandom()));
//      // create mutate object that removes furthest SSE from center of protein model
//      assemble::MutateProteinModelSSERemove mutate_furthest( ( assemble::LocatorSSEFurthest< assemble::ProteinModel>()));

      // mutate the protein model A 1 7 sse will be removed
      assemble::ProteinModel mutated_model( *mutate( protein_model).GetArgument());
      // create ProteinModel "mutated_model_random" removing a random sse
      assemble::ProteinModel mutated_model_random( *mutate_random( protein_model).GetArgument());
//      // create ProteinModel "mutated_model_furthest" removing the sse furthest apart from center
//      assemble::ProteinModel mutated_model_furthest( mutate_furthest( protein_model));

      const double rmsd_start( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, protein_model, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      const double rmsd_mutate( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      const double rmsd_random( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model_random, biol::GetAtomTypes().GetBackBoneAtomTypes()));
//      const double rmsd_furthest( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model_furthest));

      const size_t number_aas_start( protein_model.GetNumberAAs());
      const size_t number_aas_mutate( mutated_model.GetNumberAAs());
      const size_t number_aas_random( mutated_model_random.GetNumberAAs());
//      const size_t number_aas_furthest( mutated_model_furthest.GetNumberAAs());

      BCL_MessageStd( "rmsd before removal " + util::Format()( rmsd_start));
      BCL_MessageStd( "number aas before removal " + util::Format()( number_aas_start));
      BCL_MessageStd( "rmsd after removal of sse 1-7 " + util::Format()( rmsd_mutate));
      BCL_MessageStd( "number aas after removal of sse 1-7 " + util::Format()( number_aas_mutate));
      BCL_MessageStd( "rmsd after removal of random sse " + util::Format()( rmsd_random));
      BCL_MessageStd( "number aas after removal of random sse " + util::Format()( number_aas_random));
//      BCL_MessageStd( "rmsd after removal of furthest sse " + util::Format()( rmsd_furthest));
//      BCL_MessageStd( "number aas after removal of furthest sse " + util::Format()( number_aas_furthest));

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write mutated_remove_first.pdb");
      std::string out_filename( AddExampleOutputPathToFilename( mutate, "mutated_remove_first.pdb"));
      Proteins::WriteModelToPDB( mutated_model, out_filename);
      // write mutated protein model to an example pdb
      BCL_MessageStd( "write mutated_remove_random.pdb");
      out_filename = AddExampleOutputPathToFilename( mutate, "mutated_remove_random.pdb");
      Proteins::WriteModelToPDB( mutated_model_random, out_filename);
//      // write mutated protein model to an example pdb
//      BCL_MessageStd( "write mutated_remove_furthest.pdb");
//      out_filename = AddExampleOutputPathToFilename( mutate, "mutated_remove_furthest.pdb");
//      Proteins::WriteModelToPDB( mutated_model_furthest, out_filename);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSERemove

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSERemove::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSERemove())
  );

} // namespace bcl
