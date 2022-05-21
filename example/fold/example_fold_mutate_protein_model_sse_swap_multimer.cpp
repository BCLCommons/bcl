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
#include "fold/bcl_fold_mutate_protein_model_sse_swap_multimer.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_swap_multimer.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSESwapMultimer :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSESwapMultimer *Clone() const
    {
      return new ExampleFoldMutateProteinModelSSESwapMultimer( *this);
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
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // add the multiplier information to the protein model data
      protein_model.Transform( math::Inverse( protein_model.GetOrientation()));
      protein_model.Translate( linal::Vector3D( 16.0, 0.0, 0.0));
      util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        new assemble::ProteinModelMultiplier( linal::Vector3D( 0.0, 0.0, 1.0), 5, protein_model)
      );
      util::ShPtr< assemble::ProteinModelData> model_data( protein_model.GetProteinModelData());
      model_data->Insert( assemble::ProteinModelData::e_Multiplier, sp_multiplier);
      protein_model.SetProteinModelData( model_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      fold::MutateProteinModelSSESwapMultimer def_construct;
      BCL_ExampleIndirectCheck
      (
        def_construct.GetScheme(),
        "bcl::fold::MutateProteinModelSSESwapMultimer",
        "default constructor"
      );

      // test constructor from scheme
      fold::MutateProteinModelSSESwapMultimer multimer_construct( true, "example_swap_multimer");

      // test Clone
      util::ShPtr< fold::MutateProteinModelSSESwapMultimer> clone_construct( multimer_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< fold::MutateProteinModelSSESwapMultimer>(),
        clone_construct->GetClassIdentifier()
      );

      // check GetScheme
      BCL_ExampleCheck( multimer_construct.GetScheme(), "example_swap_multimer");

    ///////////////
    // operators //
    ///////////////

      // test () operator
      util::ShPtr< assemble::ProteinModel> sp_mutated_model( multimer_construct( protein_model).GetArgument());
      Proteins::WriteModelToPDB
      (
        *sp_mutated_model,
        AddExampleOutputPathToFilename( multimer_construct, "mutate_sse_swap_multimer.pdb")
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( multimer_construct);

      // read the object back in
      fold::MutateProteinModelSSESwapMultimer mutate_read;
      ReadBCLObject( mutate_read);

      BCL_ExampleIndirectCheck
      ( multimer_construct.GetScheme(), mutate_read.GetScheme(), "read and write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSESwapMultimer

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSESwapMultimer::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSESwapMultimer())
  );
  
} // namespace bcl
