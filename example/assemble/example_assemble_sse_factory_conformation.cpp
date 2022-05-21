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
#include "assemble/bcl_assemble_sse_factory_conformation.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_factory_conformation.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEFactoryConformation :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEFactoryConformation *Clone() const
    {
      return new ExampleAssembleSSEFactoryConformation( *this);
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
      util::ShPtr< biol::AASequence> sp_sequence
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone).GetChain( 'A')->GetSequence()
      );

      assemble::Chain chain( assemble::ConstructChainWithSSEsFromConformation( sp_sequence));

      BCL_MessageStd( util::Format()( chain.GetNumberSSEs()) + " sses found.");
      BCL_Example_Check
      (
        chain.GetNumberSSE( biol::GetSSTypes().HELIX) == 1 &&
        chain.GetNumberSSE( biol::GetSSTypes().STRAND) == 4,
        "Incorrect number of helices or strands calculated"
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a SSEFactoryConformation behind a pointer
      BCL_MessageStd( "Creating a SSEFactoryConformation");
      util::ShPtr< assemble::SSEFactoryInterface> sse_factory( new assemble::SSEFactoryConformation());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // call the operator and construct the sse pool
      BCL_MessageStd( "Calling the operator to get a SSEPool");
      assemble::SSEPool sse_pool( sse_factory->operator()( *sp_sequence));

      // iterate over the sses in the pool
      BCL_MessageStd( "iterating over the sses in the SSEPool");
      for
      (
        assemble::SSEPool::const_iterator
          sse_itr( sse_pool.Begin()),
          sse_itr_end( sse_pool.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // print out the sses in the pool
        BCL_MessageStd
        (
          util::Format()( ( *sse_itr)->GetType().GetName()) + "\t" +
          util::Format().W( 3)( ( *sse_itr)->GetFirstAA()->GetSeqID()) + " to " +
          util::Format().W( 3)( ( *sse_itr)->GetLastAA()->GetSeqID())
        );
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEFactoryConformation

  const ExampleClass::EnumType ExampleAssembleSSEFactoryConformation::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEFactoryConformation())
  );

} // namespace bcl

