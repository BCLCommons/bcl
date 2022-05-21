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
#include "fold/bcl_fold_mutate_protein_model_grow_sse.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_phi_psi_generator_ramachandran.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_grow_sse.cpp
  //!
  //! @author alexanns
  //! @date January 15, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelGrowSSE :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelGrowSSE *Clone() const
    {
      return new ExampleFoldMutateProteinModelGrowSSE( *this);
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
      // create protein model
      const assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")));

      const util::ShPtr
      <
        find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface>
      > loop_locator( new assemble::LocatorSSE( 'A', 18, 22));

      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> >
      > phi_psi_gen( new fold::PhiPsiGeneratorRamachandran( biol::Ramachandran::GetDefaultHistogramFilename()));

      const util::ShPtr
      <
        find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel>
      > anchor_aa_locator( new assemble::LocatorAA( 'A', 17));

      const util::ShPtr
      <
        find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface>
      > anchor_sse_locator( new assemble::LocatorSSE( 'A', 10, 17));

      const util::ShPtr
      <
        find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel>
      > target_residue_locator( new assemble::LocatorAA( 'A', 23));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const fold::MutateProteinModelGrowSSE def_constr;

      // this suffix is used for all the correctly mutated files for this system architecture
      const std::string correct_suffix( ".correct");

      // an input stream, used to load the correct file
      io::IFStream input;

      // constructor taking parameters
      const fold::MutateProteinModelGrowSSE param_constr
      (
        loop_locator, phi_psi_gen, anchor_aa_locator, biol::AASequenceFlexibility::e_CTerminal
      );
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( param_constr, "mutate_protein_model_grow_sse_a.pdb")
        );
        const util::ShPtr< assemble::ProteinModel> new_model( param_constr( model).GetArgument());

        // write the model out to a string stream
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);

        BCL_ExampleCheck( io::File::FilesMatch( out_pdb_name, out_pdb_name + correct_suffix), true);
      }

      // clone constructor
      const util::ShPtr< fold::MutateProteinModelGrowSSE> clone_constr( param_constr.Clone());
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( param_constr, "mutate_protein_model_grow_sse_b.pdb")
        );
        const util::ShPtr< assemble::ProteinModel> new_model( clone_constr->operator()( model).GetArgument());

        // write the model out to a string stream
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( out_pdb_name, out_pdb_name + correct_suffix), true);
      }

    ///////////////
    // operators //
    ///////////////

      // operator()() for mutating
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( param_constr, "mutate_protein_model_grow_sse_c.pdb")
        );
        const util::ShPtr< assemble::ProteinModel> new_model( param_constr( model).GetArgument());

        // write the model out to a string stream
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( out_pdb_name, out_pdb_name + correct_suffix), true);
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( param_constr);

      // read the object back in
      fold::MutateProteinModelGrowSSE mutate_read;
      ReadBCLObject( mutate_read);

      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( param_constr, "mutate_protein_model_grow_sse_d.pdb")
        );
        const util::ShPtr< assemble::ProteinModel> new_model( param_constr( model).GetArgument());

        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        // close the input stream
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( out_pdb_name, out_pdb_name + correct_suffix), true);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelGrowSSE

  const ExampleClass::EnumType ExampleFoldMutateProteinModelGrowSSE::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelGrowSSE())
  );

} // namespace bcl
