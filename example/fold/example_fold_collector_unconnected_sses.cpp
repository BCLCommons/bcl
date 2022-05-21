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
#include "fold/bcl_fold_collector_unconnected_sses.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_collector_unconnected_sses.cpp
  //! @brief this example tests the implementation of fold::CollectorUnconnectedSSEs which collects unconnected SSEs
  //! from a given domain.
  //!
  //! @author fischea
  //! @date Dec 17, 2014
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleFoldCollectorUnconnectedSSEs :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief clone function
    //! @return pointer to a new ExampleFoldCollectorUnconnectedSSEs
    ExampleFoldCollectorUnconnectedSSEs *Clone() const
    {
      return new ExampleFoldCollectorUnconnectedSSEs( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! @detail this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test construction methods
      const fold::CollectorUnconnectedSSE collector_default;
      const fold::CollectorUnconnectedSSE collector_1
      (
        biol::AASequenceFlexibility::e_NTerminal,
        true,
        storage::Set< biol::SSType>( biol::GetSSTypes().COIL),
        false
      );

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck
      (
        collector_default.GetClassIdentifier(), ( GetStaticClassName< fold::CollectorUnconnectedSSE>())
      );

    ////////////////
    // operations //
    ////////////////

      // test using a pdb files with different numbers of disconnected SSEs
      const std::string test_pdb_no_loops( AddExampleInputPathToFilename( e_Biology, "2AP3_no_loops.pdb"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, test_pdb_no_loops);
      const pdb::Handler pdb_no_loops( read);
      io::File::CloseClearFStream( read);
      const assemble::ProteinModel model_no_loops( pdb::Factory().ProteinModelFromPDB( pdb_no_loops));
      const std::string test_pdb_one_loops( AddExampleInputPathToFilename( e_Biology, "2AP3_one_loop.pdb"));
      BCL_ExampleMustOpenInputFile( read, test_pdb_one_loops);
      const pdb::Handler pdb_one_loop( read);
      io::File::CloseClearFStream( read);
      const assemble::ProteinModel model_one_loop( pdb::Factory().ProteinModelFromPDB( pdb_one_loop));

      // preparation of the test files
      const util::SiPtrVector< const assemble::SSE> loops( model_no_loops.GetSSEs( biol::GetSSTypes().COIL));
      const util::SiPtrList< const assemble::SSE> all_loops( loops.Begin(), loops.End());
      const util::SiPtrList< const assemble::SSE> all_loops_2( all_loops.Begin(), --all_loops.End());

      // test collection method
      // TODO needs better testing, direct list comparison impossible due to NaNs in loop coordinates
      const util::SiPtrList< const assemble::SSE> result_1( collector_1.Collect( model_no_loops));
      BCL_ExampleCheck( result_1.GetSize(), all_loops.GetSize());
      const util::SiPtrList< const assemble::SSE> result_2( collector_1.Collect( model_one_loop));
      BCL_ExampleCheck( result_2.GetSize(), all_loops_2.GetSize());

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( collector_1);

      // read object
      fold::CollectorUnconnectedSSE collector_read;
      ReadBCLObject( collector_read);

      // compare behavior of the read in object to the written out object
      BCL_ExampleCheck
      (
        ( collector_1.Collect( model_no_loops)).GetSize(), ( collector_read.Collect( model_no_loops)).GetSize()
      );

      return 0;
    }

  }; // class ExampleFoldCollectorUnconnectedSSEs

  //! single instance of this class
  const ExampleClass::EnumType ExampleFoldCollectorUnconnectedSSEs::s_Instance
  (
     GetExamples().AddEnum( ExampleFoldCollectorUnconnectedSSEs())
  );

} // namespace bcl
