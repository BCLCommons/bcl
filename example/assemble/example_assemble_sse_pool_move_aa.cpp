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
#include "assemble/bcl_assemble_sse_pool_move_aa.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_pool_move_aa.cpp
  //! @brief example that demonstrates how to use the SSEPoolMoveAA
  //!
  //! @author woetzen
  //! @date Jun 18, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEPoolMoveAA :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEPoolMoveAA *Clone() const
    {
      return new ExampleAssembleSSEPoolMoveAA( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from pick and mutate
      assemble::SSEPoolMoveAA ssepool_move_aa( math::Range< size_t>( 5, 5));

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( ssepool_move_aa.GetClassIdentifier(), GetStaticClassName< assemble::SSEPoolMoveAA>());

      // scheme
      BCL_ExampleCheck( ssepool_move_aa.GetScheme(), assemble::SSEPoolMoveAA::GetDefaultScheme());

    ///////////////
    // operators //
    ///////////////

      // mutate the pool
      util::ShPtr< assemble::SSEPool> sp_pool( ssepool_move_aa( sse_pool).GetArgument());

    //////////////////////
    // input and output //
    //////////////////////

      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( ssepool_move_aa, "pool_move_aa.bcl"));
      sp_pool->WriteSSEPool( write);
      io::File::CloseClearFStream( write);

      // write bcl object
      WriteBCLObject( ssepool_move_aa);
      assemble::SSEPoolMoveAA ssepool_move_aa_read( math::Range< size_t>( 0, 1));
      ReadBCLObject( ssepool_move_aa_read);

      while( ( sp_pool = ssepool_move_aa( *sp_pool).GetArgument()).IsDefined() && sp_pool->GetSize() > 18)
      {
        ;
        BCL_MessageDbg( "size of pool: " + util::Format()( sp_pool->GetSize()));
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEPoolMoveAA

  const ExampleClass::EnumType ExampleAssembleSSEPoolMoveAA::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEPoolMoveAA())
  );

} // namespace bcl
