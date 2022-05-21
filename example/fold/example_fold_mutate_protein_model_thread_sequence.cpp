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
#include "fold/bcl_fold_mutate_protein_model_thread_sequence.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "align/bcl_align_handler_pir.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    // explicit instantiation of AlignmentNode< biol::AABase>
    template class AlignmentNode< biol::AABase>;
  } // namespace align

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_thread_sequence.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Sep 18, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelThreadSequence :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelThreadSequence *Clone() const
    {
      return new ExampleFoldMutateProteinModelThreadSequence( *this);
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
      // read alignment
      const align::HandlerPIR< biol::AABase> pir_reader;
      const std::string alignment_file( AddExampleInputPathToFilename( e_Biology, "2lzm_1ubiA.pir"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, alignment_file);
      util::ShPtr< align::AlignmentNode< biol::AABase> > alignment( pir_reader.ReadAlignment( read, biol::AASequence()));
      io::File::CloseClearFStream( read);

      // create map of chain and alignment
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > chain_align;
      chain_align.Insert( std::pair< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >( 'A', alignment));

      // get the coordinates
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb"));
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageDbg( "test default constructor");
      fold::MutateProteinModelThreadSequence def_constr;
      const std::string correct_pdb
      (
        AddExampleOutputPathToFilename( def_constr, "MutateProteinModelThreadSequence_correct.pdb")
      );

      // constructor taking parameters
      BCL_MessageDbg( "test constructor taking parameters");
      fold::MutateProteinModelThreadSequence param_constr( chain_align);
      {
        math::MutateResult< assemble::ProteinModel> mutate_result( param_constr( protein_model));

        BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);
        const std::string outpdb_filename
        (
          AddExampleOutputPathToFilename( param_constr, "MutateProteinModelThreadSequence_a.pdb")
        );
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, outpdb_filename);
        pdb::Factory().WriteModelToPDB( *mutate_result.GetArgument(), write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( outpdb_filename, correct_pdb), true);
      }

      // clone constructor
      BCL_MessageDbg( "test clone");
      util::ShPtr< fold::MutateProteinModelThreadSequence> clone_constr( param_constr.Clone());
      {
        math::MutateResult< assemble::ProteinModel> mutate_result( clone_constr->operator()( protein_model));

        BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);
        const std::string outpdb_filename
        (
          AddExampleOutputPathToFilename( param_constr, "MutateProteinModelThreadSequence_b.pdb")
        );
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, outpdb_filename);
        pdb::Factory().WriteModelToPDB( *mutate_result.GetArgument(), write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( outpdb_filename, correct_pdb), true);
      }

    /////////////////
    // data access //
    /////////////////

      // get scheme
      BCL_MessageDbg( "test scheme");
      fold::MutateProteinModelThreadSequence scheme_constr( chain_align, "new_scheme");
      BCL_ExampleCheck( scheme_constr.GetScheme(), "new_scheme");

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // operator
      BCL_MessageDbg( "test operator");
      {
        math::MutateResult< assemble::ProteinModel> mutate_result( scheme_constr( protein_model));

        BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);
        const std::string outpdb_filename
        (
          AddExampleOutputPathToFilename( param_constr, "MutateProteinModelThreadSequence_c.pdb")
        );
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, outpdb_filename);
        pdb::Factory().WriteModelToPDB( *mutate_result.GetArgument(), write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck( io::File::FilesMatch( outpdb_filename, correct_pdb), true);
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageDbg( "test write");
      WriteBCLObject( *clone_constr);
      // read is not implemented for alignments

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelThreadSequence

  const ExampleClass::EnumType ExampleFoldMutateProteinModelThreadSequence::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelThreadSequence())
  );

} // namespace bcl
