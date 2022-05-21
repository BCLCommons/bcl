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
#include "fold/bcl_fold_mutate_protein_model_sse_pair_clash.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sse_paired.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_pair_clash.cpp
  //!
  //! @author mendenjl
  //! @date Jan 08, 2018
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSEPairClash :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSEPairClash *Clone() const
    { return new ExampleFoldMutateProteinModelSSEPairClash( *this);}

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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1X91A.ideal.pdb"));
      //build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // create mutate object "mutate_sse_pair" from max translation and rotation and CollectorSSEPaired
      fold::MutateProteinModelSSEPairClash mutate_sse_clash_remover;

      // mutate the protein model (random sse will be moved
      assemble::ProteinModel mutated_model( *mutate_sse_clash_remover( protein_model).GetArgument());

      score::AAPairHiResClash clash_score;
      BCL_ExampleCheckWithinAbsTolerance( clash_score( protein_model), 0.390367, 1e-4);
      BCL_ExampleCheckWithinAbsTolerance( clash_score( protein_model = *mutate_sse_clash_remover( protein_model).GetArgument()), 0.0, 1e-4);

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write mutated_sse_clash_remover.pdb");
      Proteins::WriteModelToPDB( mutated_model, AddExampleOutputPathToFilename( mutate_sse_clash_remover, "mutated_sse_clash_remover.pdb"));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSEPairClash

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSEPairClash::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSEPairClash())
  );

} // namespace bcl
