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
#include "fold/bcl_fold_mutate_protein_model_replicate_conformation.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_protein_model_conformations.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_replicate_conformation.cpp
  //!
  //! @author alexanns
  //! @date May 7, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelReplicateConformation :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelReplicateConformation *Clone() const
    {
      return new ExampleFoldMutateProteinModelReplicateConformation( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor taking parameters
      BCL_MessageStd( "testing parameter constructor");
      // number of times to copy the model
      const int num_replications( 5);

      // create object for replicating the models in the ensemble
      const fold::MutateProteinModelReplicateConformation mutate
      (
        util::ShPtr< assemble::CollectorProteinModelConformations>
        (
          new assemble::CollectorProteinModelConformations( true)
        ),
        num_replications,
        "replicate_scheme"
      );

      // construct min_sse_sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 5;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 9;

      // initialize pdb filename
      const std::string pdb_file( AddExampleInputPathToFilename( e_Biology, "2yv8_ideal.pdb"));

      // get the protein model
      const assemble::ProteinModel model( Proteins::GetModel( pdb_file, biol::GetAAClasses().e_AABackBone, min_sse_sizes));

      // create ensemble
      const assemble::ProteinEnsemble ensemble( model);

      // copy the models
      math::MutateResult< assemble::ProteinModel> mutate_result( mutate( ensemble));

      // make sure the ensemble has the correct size
      BCL_ExampleCheck( mutate_result.GetArgument()->GetConformationalEnsemble().GetSize(), 5);

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // operator ()
      {
        BCL_MessageStd( "testing operator");
        // number of times to copy the model
        const int num_replications( 4);

        // create object for replicating the models in the ensemble
        const fold::MutateProteinModelReplicateConformation mutate
        (
          util::ShPtr< assemble::CollectorProteinModelConformations>
          (
            new assemble::CollectorProteinModelConformations( true)
          ),
          num_replications,
          "replicate_scheme"
        );

        // construct min_sse_sizes
        storage::Map< biol::SSType, size_t> min_sse_sizes;
        min_sse_sizes[ biol::GetSSTypes().STRAND] = 5;
        min_sse_sizes[ biol::GetSSTypes().HELIX] = 9;

        // initialize pdb filename
        const std::string pdb_file( AddExampleInputPathToFilename( e_Biology, "2yv8_ideal.pdb"));

        // get the protein model
        const assemble::ProteinModel model( Proteins::GetModel( pdb_file, biol::GetAAClasses().e_AABackBone, min_sse_sizes));

        // create ensemble
         assemble::ProteinEnsemble ensemble( model);

        // seed the ensemble
        ensemble.InsertElement( util::ShPtr< assemble::ProteinModel>( ensemble.HardCopy()));
        ensemble.InsertElement( util::ShPtr< assemble::ProteinModel>( ensemble.HardCopy()));

        // copy the models
        math::MutateResult< assemble::ProteinModel> mutate_result( mutate( ensemble));

        // make sure the ensemble has the correct size
        BCL_ExampleCheck( mutate_result.GetArgument()->GetConformationalEnsemble().GetSize(), 14);
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "Testing read write")
      WriteBCLObject( mutate);
      fold::MutateProteinModelReplicateConformation mutate_read;
      ReadBCLObject( mutate_read);
      {
        BCL_MessageStd( "read in object is \n" + util::Format()( mutate_read));
        // copy the models
        math::MutateResult< assemble::ProteinModel> mutate_result( mutate_read( ensemble));

        // make sure the ensemble has the correct size
        BCL_ExampleCheck( mutate_result.GetArgument()->GetConformationalEnsemble().GetSize(), 5);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelReplicateConformation

  const ExampleClass::EnumType ExampleFoldMutateProteinModelReplicateConformation::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelReplicateConformation())
  );
  
} // namespace bcl
