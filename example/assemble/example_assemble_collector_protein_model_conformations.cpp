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
#include "assemble/bcl_assemble_collector_protein_model_conformations.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_ensemble.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_protein_model_conformations.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date May 7, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorProteinModelConformations :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorProteinModelConformations *Clone() const
    {
      return new ExampleAssembleCollectorProteinModelConformations( *this);
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

      // default constructor
      assemble::CollectorProteinModelConformations default_constr;

      // constructor taking parameters
      assemble::CollectorProteinModelConformations param_constr( true);

      // clone constructor
      util::ShPtr< assemble::CollectorProteinModelConformations> clone_constr( param_constr.Clone());

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // Collect
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

      // test using default constructor
      {
        const util::SiPtrList< const assemble::ProteinModel> collected( default_constr.Collect( ensemble));

        BCL_ExampleCheck( collected.GetSize(), 2);
      }

      // test using parameter constructor
      {
        const util::SiPtrList< const assemble::ProteinModel> collected( param_constr.Collect( ensemble));

        BCL_ExampleCheck( collected.GetSize(), 3);
      }

      // test using clone constructor
      {
        const util::SiPtrList< const assemble::ProteinModel> collected( clone_constr->Collect( ensemble));

        BCL_ExampleCheck( collected.GetSize(), 3);
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "Testing read write")
      WriteBCLObject( param_constr);
      assemble::CollectorProteinModelConformations collector_read;
      ReadBCLObject( collector_read);
      {
        const util::SiPtrList< const assemble::ProteinModel> collected( collector_read.Collect( ensemble));

        BCL_ExampleCheck( collected.GetSize(), 3);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorProteinModelConformations

  const ExampleClass::EnumType ExampleAssembleCollectorProteinModelConformations::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorProteinModelConformations())
  );
  
} // namespace bcl
