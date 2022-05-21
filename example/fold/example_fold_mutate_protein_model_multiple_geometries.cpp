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
#include "fold/bcl_fold_mutate_protein_model_multiple_geometries.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_multiple_geometries.cpp
  //! @details checks mutate and writes out example PDB
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelMultipleGeometries :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelMultipleGeometries *Clone() const
    {
      return new ExampleFoldMutateProteinModelMultipleGeometries( *this);
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

      // create the SSE pool and insert into the model
      util::ShPtr< assemble::SSEPool> sp_pool( new assemble::SSEPool( protein_model.GetSSEs()));
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData());
      sp_model_data->Insert( assemble::ProteinModelData::e_Pool, sp_pool);
      protein_model.SetProteinModelData( sp_model_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "testing default constructor");
      fold::MutateProteinModelMultipleGeometries def_construct;

      // test constructor from a SSEPool and a fold template
      BCL_MessageStd( "testing constructor from fold template");
      const assemble::FoldTemplate fold_template
      (
        protein_model,
        util::ShPtr< assemble::CollectorTopologyInterface>( new assemble::CollectorTopologyCombined()),
        "1UBI",
        false
      );
      const fold::MutateProteinModelMultipleGeometries template_construct( fold_template);

      // test copy constructor
      BCL_MessageStd( "testing copy constructor");
      const fold::MutateProteinModelMultipleGeometries copy_construct( template_construct);

      // test clone function
      BCL_MessageStd( "testing clone function");
      const util::ShPtr< fold::MutateProteinModelMultipleGeometries> clone_construct( template_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::fold::MutateProteinModelMultipleGeometries");
      BCL_Example_Check
      (
        GetStaticClassName< fold::MutateProteinModelMultipleGeometries>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< fold::MutateProteinModelMultipleGeometries>() +
        " but should give " + correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< fold::MutateProteinModelMultipleGeometries>() == clone_construct->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_construct->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // check GetScheme
      BCL_Example_Check
      (
        def_construct.GetScheme() == correct_static_class_name,
        "GetScheme gives " + def_construct.GetScheme() + " but should give " + correct_static_class_name
      );

    ///////////////
    // operators //
    ///////////////

      // check operator ()
      BCL_MessageStd( "testing () operator");

      // test w/ same size template
      assemble::ProteinModel empty_model( protein_model.GetEmptyChains());
      empty_model.SetProteinModelData( protein_model.GetProteinModelData());
      const assemble::ProteinModel mutated_model( *template_construct( empty_model).GetArgument());
      BCL_ExampleIndirectCheck( mutated_model.GetSSEs().GetSize(), 5, "() operator with same size template");
      Proteins::WriteModelToPDB
      (
        mutated_model, AddExampleOutputPathToFilename( def_construct, "1ubi_mutate_multiple_geometries.pdb")
      );

      // test w/ large template
      const fold::MutateProteinModelMultipleGeometries large_mutate( storage::VectorND< 3, double>( 0.0, 0.0, 1.0));
      const assemble::ProteinModel mutated_model_large( *large_mutate( empty_model).GetArgument());
      BCL_ExampleIndirectCheck( mutated_model_large.GetSSEs().GetSize(), 5, "() operator with larger template");

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( def_construct);

      // read the object back in
      fold::MutateProteinModelMultipleGeometries mutate_read;
      ReadBCLObject( mutate_read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelMultipleGeometries

  const ExampleClass::EnumType ExampleFoldMutateProteinModelMultipleGeometries::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelMultipleGeometries())
  );

} // namespace bcl
