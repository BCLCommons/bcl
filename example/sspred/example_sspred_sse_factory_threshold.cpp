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
#include "sspred/bcl_sspred_sse_factory_threshold.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "io/bcl_io_file.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sspred_sse_factory_threshold.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredSseFactoryThreshold :
    public ExampleInterface
  {
  public:

    ExampleSspredSseFactoryThreshold *Clone() const
    {
      return new ExampleSspredSseFactoryThreshold( *this);
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
      // create AAsequence sequence
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // initialize ssmethod
      const sspred::Method method( sspred::GetMethods().e_JUFO);

      BCL_MessageStd( "get sequence from the pdb");
      biol::AASequence sequence( *Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AACaCb).GetSequences()( 0));

      // read ss predictions
      sspred::MethodHandler::ReadPredictionsForAASequence
      (
        storage::Set< sspred::Method>( method), sequence, "1ubi", AddExampleInputPathToFilename( e_Biology, "")
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a SSEFactoryConformation behind a pointer
      BCL_MessageStd( "Creating a SSEFactorySSPrediction using constructor from sspred::methods");
      sspred::SSEFactoryThreshold sse_factory( method);

      // threshold map for predictions
      storage::Map< biol::SSType, double> threshold_map;
      threshold_map[ biol::GetSSTypes().HELIX] = 0.5;
      threshold_map[ biol::GetSSTypes().STRAND] = 0.5;
      threshold_map[ biol::GetSSTypes().COIL] = 0.5;

      // factory with all constructor arguments
      sspred::SSEFactoryThreshold sse_factory_jufo( method, threshold_map, false);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // call the operator and construct the sse pool
      BCL_MessageStd( "Calling the operator to get a SSEPool from SSEFactory with " + method.GetName());
      assemble::SSEPool sse_pool( sse_factory( sequence));
      assemble::SSEPool sse_pool_jufo( sse_factory_jufo( sequence));

      // output number of sses found
      BCL_MessageStd( "number of SSEs found: " + util::Format()( sse_pool.GetSize()));
      BCL_MessageStd( "number of SSEs with " + method.GetName() + " found: " + util::Format()( sse_pool_jufo.GetSize()));

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( sse_factory_jufo, "pool_factory_threshold_" + method.GetName() + ".bcl"));
      sse_pool_jufo.WriteSSEPool( write);
      io::File::CloseClearFStream( write);

      // write bcl object
      WriteBCLObject( sse_factory_jufo);
      sspred::SSEFactoryThreshold factory_sspred_read;
      ReadBCLObject( factory_sspred_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredSEFactoryThreshold

  const ExampleClass::EnumType ExampleSspredSseFactoryThreshold::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredSseFactoryThreshold())
  );

} // namespace bcl

