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
#include "assemble/bcl_assemble_sse_pool.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_line.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_pool.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEPool :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEPool *Clone() const
    {
      return new ExampleAssembleSSEPool( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));

      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create the sse pool
      BCL_MessageStd( "Creating the SSEPool");
      assemble::SSEPool sse_pool( protein_model.GetSSEs());

      // print the sse pool
      BCL_MessageStd( "Printing the SSEPool");
      sse_pool.WriteSSEPool( util::GetLogger());

    /////////////////
    // data access //
    /////////////////

      // write pool to file
      io::OFStream write;
      BCL_MessageStd( "write 1IE9.pool");
      std::string out_filename( AddExampleOutputPathToFilename( sse_pool, "1IE9.pool"));
      BCL_ExampleMustOpenOutputFile( write, out_filename);
      sse_pool.WriteSSEPool( write);
      io::File::CloseClearFStream( write);

      // read from the pool
      BCL_MessageStd( "read 1IE9.pool");
      assemble::SSEPool sse_pool_from_file;
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, out_filename);
      sse_pool_from_file.ReadSSEPool( read, protein_model, 3, 3);
      io::File::CloseClearFStream( read);

      BCL_MessageStd( "test write and read ");
      // write pool to bcl file
      WriteBCLObject( sse_pool);
      // read from the pool
      assemble::SSEPool sse_pool_from_bcl_file;
      ReadBCLObject( sse_pool_from_bcl_file);

      BCL_MessageStd( "iterating over both pools to ensure read functions correctly");
      BCL_ExampleIndirectAssert( sse_pool.GetSize(), sse_pool_from_file.GetSize(), "I/O");
      BCL_ExampleIndirectAssert( sse_pool.GetSize(), sse_pool_from_bcl_file.GetSize(), "I/O");
      // iterate over both pools at the same time
      for
      (
        assemble::SSEPool::const_iterator
          sse_itr_orig( sse_pool.Begin()), sse_itr_orig_end( sse_pool.End()),
          sse_itr_read( sse_pool_from_file.Begin()),
          sse_itr_bcl_read( sse_pool_from_bcl_file.Begin());
        sse_itr_orig != sse_itr_orig_end;
        ++sse_itr_orig, ++sse_itr_read, ++sse_itr_bcl_read
      )
      {
        //check that reading from normal input file was correct
        BCL_ExampleIndirectCheck( **sse_itr_orig, **sse_itr_read, "I/O");
        // check that reading from a bcl stype input file was correct
        BCL_ExampleIndirectCheck( **sse_itr_orig, **sse_itr_bcl_read, "I/O");
      }

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test ReadSSEPoolInformation function
      BCL_MessageStd( "test ReadSSEPoolInformation function");
      const storage::Vector< pdb::Line> pool_info
      (
        sse_pool.ReadSSEPoolInformation( out_filename)
      );
      // check "pool_info" has the correct number of pdb::Lines in it
      BCL_ExampleIndirectCheck( pool_info.GetSize(), 12, "sse_pool.ReadSSEPoolInformation( out_filename)");

      // test PeakChains function
      BCL_MessageStd( "test GetChainsRepresented function");
      std::string filename( AddExampleInputPathToFilename( e_Biology, "test.pool"));
      const storage::Set< char> chain_ids( sse_pool.GetChainsRepresented( filename));
      BCL_MessageStd( "The chains are : \n" + util::Format()( chain_ids));
      // check "chain_ids" has the correct number of chainids in it
      const size_t number_ids( chain_ids.GetSize());
      const size_t correct_number_ids( 2);
      BCL_Example_Check
      (
        number_ids == number_ids,
        "There should have been " + util::Format()( correct_number_ids) + " chain ids read in from "
        + filename + " but instead only " + util::Format()( number_ids) + " were read in."
      );
      // check that "chain_ids" has the correct chains in it
      BCL_Example_Check
      (
        *chain_ids.Begin() == 'A',
        "The first chain should be " + util::Format()( 'A') + " but instead is " + util::Format()( *chain_ids.Begin())
      );
      BCL_Example_Check
      (
        *++chain_ids.Begin() == 'B',
        "The second chain should be " + util::Format()( 'B') + " but instead is " + util::Format()( *++chain_ids.Begin())
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEPool

  const ExampleClass::EnumType ExampleAssembleSSEPool::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEPool())
  );

} // namespace bcl

