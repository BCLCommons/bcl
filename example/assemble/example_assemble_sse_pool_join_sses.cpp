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
#include "assemble/bcl_assemble_sse_pool_join_sses.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_pool_join_sses.cpp
  //! @brief this example demonstrates how a pool is mutated by joining adjacent SSEs of the same sstype
  //!
  //! @author woetzen
  //! @date Jun 25, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEPoolJoinSSEs :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEPoolJoinSSEs *Clone() const
    {
      return new ExampleAssembleSSEPoolJoinSSEs( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));

      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AA, ssetype_min_size)
      );

      // create the sse pool
      BCL_MessageStd( "Creating the SSEPool");
      assemble::SSEPool sse_pool( protein_model.GetSSEs(), false);

      // mutate and pick
      util::ShPtr< assemble::LocatorSSERandom> sp_locate_sse( new assemble::LocatorSSERandom());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from pick and mutate
      assemble::SSEPoolJoinSSEs ssepool_join( sp_locate_sse);

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( ssepool_join.GetClassIdentifier(), GetStaticClassName< assemble::SSEPoolJoinSSEs>());

    ///////////////
    // operators //
    ///////////////

      // mutate the pool
      util::ShPtr< assemble::SSEPool> sp_pool( ssepool_join( sse_pool).GetArgument());

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( ssepool_join, "pool_joined.bcl"));
      sp_pool->WriteSSEPool( write);
      io::File::CloseClearFStream( write);

      // write bcl object
      WriteBCLObject( ssepool_join);
      assemble::SSEPoolJoinSSEs ssepool_join_read
      (
        ( util::ShPtr< assemble::LocatorSSERandom>())
      );
      ReadBCLObject( ssepool_join_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEPoolJoinSSEs

  const ExampleClass::EnumType ExampleAssembleSSEPoolJoinSSEs::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEPoolJoinSSEs())
  );

} // namespace bcl
