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
#include "assemble/bcl_assemble_sse_factory_mc.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "assemble/bcl_assemble_sse_pool_agreement.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_factory_mc.cpp
  //! @brief this example shows how the ssefactorymc generates sse pools
  //!
  //! @author woetzen
  //! @date Jun 19, 2011
  //! @remarks status empty
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEFactoryMC :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEFactoryMC *Clone() const
    {
      return new ExampleAssembleSSEFactoryMC( *this);
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

    /////////////////
    // preparation //
    /////////////////

      const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_file_name, biol::GetAAClasses().e_AA)
      );

      const storage::Set< sspred::Method> methods( sspred::GetMethods().e_JUFO, sspred::GetMethods().e_PSIPRED);
      sspred::MethodHandler::ReadPredictionsForProteinModel( methods, model, "1ubi", AddExampleInputPathToFilename( e_Biology, ""));

      const assemble::SSEPool native_pool( model.GetSSEs(), false);

      const biol::Membrane membrane;

      // iterate over methods
      for
      (
        storage::Set< sspred::Method>::const_iterator method_itr( methods.Begin()), method_itr_end( methods.End());
        method_itr != method_itr_end;
        ++method_itr
      )
      {
        const double confidence_threshold( 0.5);
        const sspred::Method method( *method_itr);

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        // default constructor
        assemble::SSEFactoryMC factory_mc( method, confidence_threshold);

      /////////////////
      // data access //
      /////////////////

        // set the method
        factory_mc.SetMethod( method);

        // set the scoring function
        factory_mc.SetScoringFunction( factory_mc.DefaultScoringFunction());

        // set the mutate
        factory_mc.SetMutate( factory_mc.DefaultMutate());

      ///////////////
      // operators //
      ///////////////

        BCL_MessageStd
        (
          "score of native, Q3 and agreement score between native and native: " +
          method.GetName() + "\t" +
          util::Format()( factory_mc.DefaultScoringFunction()->operator ()( native_pool, membrane)) + "\t" +
          util::Format()( assemble::SSEPoolAgreement().Q3Score( native_pool, native_pool)) + "\t" +
          util::Format()( assemble::SSEPoolAgreement()( native_pool, native_pool))
        );

        // create pool
        factory_mc.SetMaxNumberIterations( 1000);
        assemble::SSEPool new_pool( factory_mc( *model.GetChain( 'A')->GetSequence()));
        BCL_MessageStd
        (
          "score of new, Q3 and agreement score between native and new pool: " +
          method.GetName() + "\t" +
          util::Format()( factory_mc.DefaultScoringFunction()->operator ()( new_pool, membrane)) + "\t" +
          util::Format()( assemble::SSEPoolAgreement().Q3Score( native_pool, new_pool)) + "\t" +
          util::Format()( assemble::SSEPoolAgreement()( native_pool, new_pool))
        );

        // decrease iterations and increase number optimizations
        factory_mc.SetMaxNumberIterations( 750);
        factory_mc.SetNumberOptimizations( 2);
        assemble::SSEPool min_pool( factory_mc( *model.GetChain( 'A')->GetSequence()));
        BCL_MessageStd
        (
          "score of min, Q3 and agreement score between native and min pool: " +
          method.GetName() + "\t" +
          util::Format()( factory_mc.DefaultScoringFunction()->operator ()( min_pool, membrane)) + "\t" +
          util::Format()( assemble::SSEPoolAgreement().Q3Score( native_pool, min_pool)) + "\t" +
          util::Format()( assemble::SSEPoolAgreement()( native_pool, min_pool))
        );

      ////////////////
      // operations //
      ////////////////

      //////////////////////
      // input and output //
      //////////////////////

        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( new_pool, "pool_factory_mc_" + method.GetName() + ".bcl"));
        new_pool.WriteSSEPool( write);
        io::File::CloseClearFStream( write);

        storage::Map< biol::SSType, size_t> min_sse_size;
        min_sse_size[ biol::GetSSTypes().HELIX] = 5;
        min_sse_size[ biol::GetSSTypes().STRAND] = 3;
        while( new_pool.Join( min_sse_size));

        BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( new_pool, "pool_factory_mc_joined_" + method.GetName() + ".bcl"));
        new_pool.WriteSSEPool( write);
        io::File::CloseClearFStream( write);

      //////////////////////
      // helper functions //
      //////////////////////

      } // iterate over methods

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEFactoryMC

  const ExampleClass::EnumType ExampleAssembleSSEFactoryMC::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEFactoryMC())
  );

} // namespace bcl
