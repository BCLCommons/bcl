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
#include "assemble/bcl_assemble_sse_pool_split_sse.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_pool_split_sse.cpp
  //! @brief example that demonstrates how to use the SSEPoolSplitSSE
  //!
  //! @author woetzen
  //! @date Jun 18, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEPoolSplitSSE :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEPoolSplitSSE *Clone() const
    {
      return new ExampleAssembleSSEPoolSplitSSE( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AA, ssetype_min_size)
      );

      // create the sse pool
      BCL_MessageStd( "Creating the SSEPool");
      assemble::SSEPool sse_pool( protein_model.GetSSEs(), false);

      // pick
      util::ShPtr< assemble::LocatorSSERandom> sp_locate_sse( new assemble::LocatorSSERandom());

      // method
      sspred::Method method( sspred::GetMethods().e_JUFO);
      sspred::MethodHandler::ReadPredictionsForProteinModel( storage::Set< sspred::Method>( method), protein_model, "1ubi", AddExampleInputPathToFilename( e_Biology, ""));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from pick and mutate
      assemble::SSEPoolSplitSSE ssepool_split
      (
        method,
        sp_locate_sse,
        util::ShPtr< math::MutateInterface< assemble::SSE> >()
      );

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( ssepool_split.GetClassIdentifier(), GetStaticClassName< assemble::SSEPoolSplitSSE>());

      // scheme
      BCL_ExampleCheck( ssepool_split.GetScheme(), assemble::SSEPoolSplitSSE::GetDefaultScheme());

    ///////////////
    // operators //
    ///////////////

      // mutate the pool
      util::ShPtr< assemble::SSEPool> sp_pool( ssepool_split( sse_pool).GetArgument());

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( ssepool_split, "pool_split.bcl"));
      sp_pool->WriteSSEPool( write);
      io::File::CloseClearFStream( write);

      // write bcl object
      WriteBCLObject( ssepool_split);
      assemble::SSEPoolSplitSSE ssepool_split_read
      (
        sspred::GetMethods().e_Undefined,
        ( util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> >()),
        ( util::ShPtr< math::MutateInterface< assemble::SSE> >())
      );
      ReadBCLObject( ssepool_split_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEPoolSplitSSE

  const ExampleClass::EnumType ExampleAssembleSSEPoolSplitSSE::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEPoolSplitSSE())
  );

} // namespace bcl
